"""NIST SRM 1950 / SRM1950-DB adapter (metabolite arm — structure oracle).

SRM1950-DB (Mandal et al. 2025, Anal. Chem., DOI 10.1021/acs.analchem.4c05018) is the certified
NIST SRM 1950 human-plasma reference set — 1,058 metabolites at srm1950-data.wishartlab.com. The
CSV delivery (metabolites.csv) ships ``HMDB_ID``, ``NAME``, ``SMILES`` and ``INCHIKEY`` columns,
but at acquisition the **INCHIKEY column is empty** while SMILES is ~95% populated.

Design mirrors ``necs_metabolon`` (small enough to load in full) with ONE acquisition-driven
difference: the independent structure-oracle InChIKey is **derived from the certified SMILES** via
RDKit (deterministic, standard cheminformatics; zero shared infra with BioMapper's resolver) when
the delivery's INCHIKEY column is empty. An explicit INCHIKEY (should a future delivery populate
it) is preferred verbatim. Rows whose SMILES fails to parse (or is absent) yield no gold structure
and are retained as coverage-only — excluded from the accuracy denominator by the scorer.

Network is isolated behind ``fetch_supplement`` so the transform is fully unit-testable offline.
"""

from __future__ import annotations

import hashlib
import io
import re
from dataclasses import dataclass
from typing import Any

import pandas as pd

from biomapper.benchmarks.config import SRM1950, DatasetConfig

HAS_STRUCTURE_COL = "has_gold_structure"

_ACCESSION_DIGITS = re.compile(r"(\d+)")


class RowIndexGoldColumnError(ValueError):
    """A gold accession column whose values are really the row number.

    The delivery's identifier column ran in file order, one value per row, with the numeric parts
    forming the sequence one to n against chemically unrelated names. It has an accession's exact
    format, so it reads as gold to anything that greps for one, and the identifier-based coverage
    figure derived from it was an artefact of the synthetic column rather than a resolver result.
    The run refuses it at acquisition rather than scoring against it.
    """


# Canonical held-out column -> candidate raw headers (case-insensitive, exact after strip).
QUERY_CANDIDATES: tuple[str, ...] = ("NAME", "Name", "metabolite_name", "Metabolite")
SMILES_CANDIDATES: tuple[str, ...] = ("SMILES", "Smiles", "Canonical SMILES")
INCHIKEY_CANDIDATES: tuple[str, ...] = ("INCHIKEY", "InChIKey", "INCHI_KEY", "InChI Key")
HMDB_CANDIDATES: tuple[str, ...] = ("HMDB_ID", "HMDB", "HMDBID", "HMDB ID")


def sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def fetch_supplement(url: str, *, timeout: float = 60.0) -> bytes:
    """Fetch the SRM1950-DB metabolites.csv bytes (network). Isolated so tests never hit it."""
    import requests

    resp = requests.get(url, timeout=timeout)
    resp.raise_for_status()
    return resp.content


def parse_csv(raw: bytes) -> pd.DataFrame:
    """Parse the delivery bytes into a raw DataFrame (all cells as strings, blanks preserved)."""
    text = raw.decode("utf-8-sig")
    return pd.read_csv(io.StringIO(text), dtype=str).fillna("")


def _norm(value: Any) -> str:
    return "" if value is None else str(value).strip()


def _resolve_column(raw_df: pd.DataFrame, candidates: tuple[str, ...]) -> str | None:
    lookup = {str(c).strip().lower(): c for c in raw_df.columns}
    for cand in candidates:
        hit = lookup.get(cand.strip().lower())
        if hit is not None:
            return hit
    return None


def inchikey_from_smiles(smiles: Any) -> str:
    """Standard InChIKey from a SMILES via RDKit, or "" when absent/unparseable.

    The certified structure is the dataset's own SMILES; RDKit's conversion is deterministic and
    shares no infrastructure with BioMapper's resolver, preserving oracle independence.
    """
    s = _norm(smiles)
    if not s:
        return ""
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")  # type: ignore[attr-defined]
    mol = Chem.MolFromSmiles(s)
    if mol is None:
        return ""
    return Chem.MolToInchiKey(mol)


def is_row_index_column(values: Any) -> bool:
    """True when a column of accession-shaped strings is really the row number in disguise.

    The rule is deliberately narrow, because the cost of a false positive is a refused run. All
    four must hold over the FULL column: every value populated, numeric parts unique, monotonically
    increasing, and exactly the consecutive sequence starting at one. A genuine accession set is
    unique but not consecutive; a sorted genuine set has gaps; a filtered subset of a corrupt
    column also has gaps, which is why this must be evaluated on the raw delivery and never on a
    downstream slice.

    A single row satisfies "the consecutive sequence starting at one" trivially and is not
    evidence, so the guard needs more than one row before it will fire.
    """
    parsed: list[int] = []
    for value in list(values):
        text = _norm(value)
        if not text:
            return False
        match = _ACCESSION_DIGITS.search(text)
        if match is None:
            return False
        parsed.append(int(match.group(1)))
    if len(parsed) < 2:
        return False
    return parsed == list(range(1, len(parsed) + 1))


def _refuse_row_index_gold_columns(raw_df: pd.DataFrame) -> None:
    """Fail the run loudly when the delivery ships an accession column that is a row index.

    Complements (and does not replace) the generic uniqueness-and-monotonicity quarantine that
    handles *unknown* columns: this one refuses a known-bad column outright at acquisition, before
    any figure can be computed from it.
    """
    hmdb_raw = _resolve_column(raw_df, HMDB_CANDIDATES)
    if hmdb_raw is None:
        return
    if is_row_index_column(raw_df[hmdb_raw].tolist()):
        raise RowIndexGoldColumnError(
            f"{hmdb_raw!r} in the delivery is a row index wearing an accession's format: its "
            f"numeric parts are unique, monotonic, and exactly the consecutive sequence starting "
            f"at one over all rows. Refusing to build the input rather than scoring coverage "
            f"against a synthetic gold column."
        )


def build_input_df(raw_df: pd.DataFrame, config: DatasetConfig = SRM1950) -> pd.DataFrame:
    """Build the mapper-ready input_df: name query + held-out gold columns + structure flag.

    Gold InChIKey = the delivery's INCHIKEY when present, else derived from the certified SMILES.
    The gold columns ride along untouched into the mapper output (``provided_id_columns=[]``) and
    are consumed only by the scorer.
    """
    query_raw = _resolve_column(raw_df, QUERY_CANDIDATES)
    if query_raw is None:
        raise KeyError(
            f"SRM1950 delivery is missing a recognizable NAME column; tried {QUERY_CANDIDATES!r} "
            f"against {list(raw_df.columns)!r}"
        )
    smiles_raw = _resolve_column(raw_df, SMILES_CANDIDATES)
    inchikey_raw = _resolve_column(raw_df, INCHIKEY_CANDIDATES)
    _refuse_row_index_gold_columns(raw_df)

    out = pd.DataFrame()
    out[config.name_column] = raw_df[query_raw].map(_norm)
    smiles = raw_df[smiles_raw].map(_norm) if smiles_raw is not None else pd.Series([""] * len(raw_df))
    explicit_ik = raw_df[inchikey_raw].map(_norm) if inchikey_raw is not None else pd.Series([""] * len(raw_df))
    assert config.gold_smiles_column is not None  # SRM1950 config carries a gold SMILES column
    out[config.gold_smiles_column] = smiles.values
    # Prefer an explicit delivery InChIKey; otherwise derive from the certified SMILES.
    out[config.gold_inchikey_column] = [
        ik if ik else inchikey_from_smiles(sm) for ik, sm in zip(explicit_ik.values, smiles.values)
    ]
    # The delivery's identifier column is NOT emitted. See ``RowIndexGoldColumnError``: it was a row
    # index in accession clothing, and a quarantined-but-present gold column is a trap for the next
    # person who greps for a gold identifier. The structure oracle is the certified SMILES-derived
    # InChIKey, which never read this column, so accuracy is unaffected by the drop.
    out[HAS_STRUCTURE_COL] = out[config.gold_inchikey_column].map(lambda s: bool(_norm(s)))
    return out


def build_card(
    raw_df: pd.DataFrame,
    source_sha: str,
    config: DatasetConfig = SRM1950,
) -> dict[str, Any]:
    """Build the dataset_card: N, input_type, per-column coverage, oracle provenance, SHA, license."""
    input_df = build_input_df(raw_df, config)
    n = len(input_df)
    coverage: dict[str, dict[str, Any]] = {}
    for namespace, column in config.gold_coverage_columns:
        present = int((input_df.get(column, pd.Series([""] * n)).map(_norm) != "").sum())
        coverage[namespace] = {"n": present, "fraction": (present / n) if n else 0.0}
    return {
        "dataset": config.key,
        "arm": config.arm,
        "entity_type": config.entity_type,
        "input_type": config.input_type,
        "target_vocabs": list(config.target_vocabs),
        "n_rows": n,
        "coverage": coverage,
        "structure_oracle_column": config.gold_inchikey_column,
        # Load-bearing provenance: the delivery's InChIKey column is empty, so the oracle InChIKey
        # is derived from the certified SMILES (recorded so a reviewer knows the oracle's origin).
        "structure_oracle_source": "derived_from_certified_smiles",
        "source_doi": config.source_doi,
        "source_url": config.source_url,
        "source_sha256": source_sha,
        "license": config.license,
    }


@dataclass(frozen=True)
class SRM1950Bundle:
    input_df: pd.DataFrame
    card: dict[str, Any]


def load_srm1950(source: bytes | str | pd.DataFrame, config: DatasetConfig = SRM1950) -> SRM1950Bundle:
    """Load SRM1950 from raw CSV bytes (SHA pinned), a URL (fetched), or a DataFrame (tests).

    When ``source`` is a DataFrame the card's ``source_sha256`` is computed over its canonical CSV
    bytes so the pin is deterministic for tests.
    """
    if isinstance(source, pd.DataFrame):
        raw_df = source
        raw_bytes = raw_df.to_csv(index=False).encode("utf-8")
    elif isinstance(source, bytes):
        raw_bytes = source
        raw_df = parse_csv(raw_bytes)
    elif isinstance(source, str):
        raw_bytes = fetch_supplement(source)
        raw_df = parse_csv(raw_bytes)
    else:
        raise TypeError(f"unsupported source type {type(source)!r}")

    sha = sha256_bytes(raw_bytes)
    return SRM1950Bundle(input_df=build_input_df(raw_df, config), card=build_card(raw_df, sha, config))
