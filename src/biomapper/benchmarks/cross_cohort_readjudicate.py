"""Re-adjudicate a cross-cohort certificate verdict against an OUTSIDE source.

A ``refuted`` verdict from :mod:`biomapper.benchmarks.cross_cohort_certify` is not a finding. Two
independent reasons:

1. The comparison is an InChIKey **first block**, which is the connectivity layer. It is neither
   tautomer- nor charge-invariant, so a tautomer, a charge state or a salt form can differ in the
   first block while being the same compound for harmonization purposes. The first-block rule
   therefore over-flags tautomers as refuted. It also fails the other way, silently accepting a
   stereoisomer error whenever one side has no stereo layer, which no amount of re-adjudication
   here can recover.
2. The NECS curated gold carries roughly 5% InChIKey errors. A disagreement with it is evidence
   that one of the two sides is wrong, not evidence about which.

So each non-certified case is checked against PubChem's **name** index, which is a different lookup
than either side used: the NECS side came from the MOESM5 curated key, the cohort side from the
cohort's own vendor identifier. Formula and monoisotopic mass come back with it, which is what
separates a tautomer or charge artifact from a genuinely different molecule.

Outcomes, and what each licenses:

``tautomer_or_charge_or_salt_artifact``
    Both sides resolve to the same molecular formula, or to formulas differing only in hydrogen
    count, and their masses agree within tolerance. The blocks differ because the first block is not
    invariant to that. Not a mapping error. Do not quote it as one.
``necs_gold_suspect``
    The outside source agrees with the cohort side and disagrees with the NECS gold. This is the
    documented ~5% gold defect showing up, not a BioMapper error.
``cohort_id_suspect``
    The outside source agrees with the NECS gold and disagrees with the cohort's vendor identifier.
    A defect in the cohort's published annotation.
``both_sides_disagree_with_outside``
    Neither side matches the outside source. Something upstream of both is wrong; escalate by hand.
``genuine_structural_disagreement``
    Formulas genuinely differ and the outside source corroborates one side. The link is wrong.
``outside_source_unresolved``
    PubChem could not resolve the name, so nothing is concluded. Refusal stands and is reported as
    unadjudicated, never quietly folded into either direction.

Every outcome is emitted per case with both blocks, the outside block, and the formulas, so a
reviewer can check the call rather than take it.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import quote

import pandas as pd

from biomapper.benchmarks.adapters.metlinkr import force_ipv4
from biomapper.benchmarks.scorers.independent_inchikey import _PUG_REST, _first_block

# Monoisotopic masses agree to well under this when two records describe one compound; a real
# difference in composition is far larger. Deliberately loose, because the question here is "same
# molecule or not", not "same to the millidalton".
MASS_TOLERANCE_DA = 0.02

# Monoisotopic mass of hydrogen. Two records that differ by a protonation state differ in mass by
# very close to a multiple of this, so the mass check has to EXPECT that difference rather than
# demand equality. Demanding equality would reject the exact case this module exists to detect.
HYDROGEN_MASS_DA = 1.007825

# How many hydrogens a tautomer, protonation state or salt form may differ by before the two records
# stop being one compound. Keto-enol and zwitterion cases are 1 to 2; a larger gap is a reduction or
# a different molecule, and excusing it would let a real mapping error through as an artifact.
MAX_HYDROGEN_DELTA = 2

PROPERTIES = "InChIKey,MolecularFormula,MonoisotopicMass"


@dataclass(frozen=True)
class OutsideRecord:
    """What an outside source says about a name. ``status`` separates a miss from a failure."""

    block: str | None
    formula: str | None
    mass: float | None
    status: str  # success | clean_miss | lookup_failed

    @property
    def resolved(self) -> bool:
        return self.status == "success" and self.block is not None


def _parse_formula(formula: str | None) -> dict[str, int] | None:
    """Element counts from a molecular formula, or ``None`` if it does not parse.

    Charge suffixes (``+``, ``-``, ``2-``) are dropped: a charge state is exactly the kind of
    difference this module exists to classify as an artifact rather than a different compound.
    """
    if not formula:
        return None
    import re

    cleaned = re.sub(r"[+-]\d*$|\d*[+-]$", "", str(formula).strip())
    if not cleaned or not re.fullmatch(r"(?:[A-Z][a-z]?\d*)+", cleaned):
        return None
    counts: dict[str, int] = {}
    for element, digits in re.findall(r"([A-Z][a-z]?)(\d*)", cleaned):
        counts[element] = counts.get(element, 0) + (int(digits) if digits else 1)
    return counts


def formulas_differ_only_in_hydrogen(left: str | None, right: str | None) -> bool:
    """True when two formulas match on every element except hydrogen.

    That is the signature of a tautomer, a protonation state or a charge variant, which the
    first-block comparison reports as refuted even though the compounds are the same for
    harmonization. Returns False when either formula does not parse, so an unparseable formula never
    buys a free pass.
    """
    a, b = _parse_formula(left), _parse_formula(right)
    if a is None or b is None:
        return False
    return {k: v for k, v in a.items() if k != "H"} == {k: v for k, v in b.items() if k != "H"}


def compositions_equivalent(left: OutsideRecord, right: OutsideRecord) -> bool:
    """True when two records describe one compound in a different protonation or tautomer form.

    The heavy-atom composition must match exactly. Hydrogens may differ by at most
    :data:`MAX_HYDROGEN_DELTA`, and when both masses are known the observed mass gap must match the
    gap that hydrogen count implies. Checking the mass against the EXPECTED difference rather than
    against zero is the point: a protonated and a neutral record of one compound differ by very
    close to one hydrogen mass, so demanding equal masses would reject every case this is meant
    to catch,
    while ignoring mass entirely would let two genuinely different compounds that happen to share a
    heavy-atom skeleton pass as an artifact.
    """
    a, b = _parse_formula(left.formula), _parse_formula(right.formula)
    if a is None or b is None:
        return False
    if {k: v for k, v in a.items() if k != "H"} != {k: v for k, v in b.items() if k != "H"}:
        return False
    hydrogen_delta = abs(a.get("H", 0) - b.get("H", 0))
    if hydrogen_delta > MAX_HYDROGEN_DELTA:
        return False
    if left.mass is None or right.mass is None:
        return True  # formula-only evidence; stated as such by the caller's rationale
    expected_gap = hydrogen_delta * HYDROGEN_MASS_DA
    return abs(abs(left.mass - right.mass) - expected_gap) <= MASS_TOLERANCE_DA


class OutsideResolver:
    """PubChem name-index lookups returning block, formula and mass together.

    A different lookup than either side of the certificate used, which is the whole point: the NECS
    side came from the curated supplement key and the cohort side from a vendor identifier, so the
    name index is independent of both. Cached, IPv4-forced (the desktop IPv6 route to some CDNs is
    broken), and fail-soft.
    """

    def __init__(self, *, timeout: float = 20.0, session: Any | None = None) -> None:  # noqa: ANN401
        import requests

        self._timeout = timeout
        self._session = session or requests.Session()
        self._cache: dict[str, OutsideRecord] = {}

    def by_name(self, name: str) -> OutsideRecord:
        return self._resolve(f"name:{name}", f"compound/name/{quote(name, safe='')}")

    def by_cid(self, cid: str) -> OutsideRecord:
        return self._resolve(f"cid:{cid}", f"compound/cid/{quote(cid, safe='')}")

    def _resolve(self, cache_key: str, path: str) -> OutsideRecord:
        cached = self._cache.get(cache_key)
        # A transient failure is never cached as terminal: a resumed run must retry it rather than
        # bank a refusal the service would not have given a minute later.
        if cached is not None and cached.status != "lookup_failed":
            return cached
        url = f"{_PUG_REST}/{path}/property/{PROPERTIES}/JSON"
        try:
            with force_ipv4():
                response = self._session.get(url, timeout=self._timeout)
        except Exception:  # noqa: BLE001 - an outside lookup must never abort adjudication
            record = OutsideRecord(None, None, None, "lookup_failed")
            self._cache[cache_key] = record
            return record
        if response.status_code == 404:
            record = OutsideRecord(None, None, None, "clean_miss")
        elif response.status_code != 200:
            record = OutsideRecord(None, None, None, "lookup_failed")
        else:
            record = self._parse(response.text)
        self._cache[cache_key] = record
        return record

    @staticmethod
    def _parse(body: str) -> OutsideRecord:
        try:
            properties = json.loads(body)["PropertyTable"]["Properties"]
        except Exception:  # noqa: BLE001
            return OutsideRecord(None, None, None, "lookup_failed")
        if not properties:
            return OutsideRecord(None, None, None, "clean_miss")
        # A name can hit several CIDs. The first is PubChem's own ranking; taking it is a choice, so
        # an ambiguous name is recorded as a clean miss rather than silently adjudicated off rank 1.
        if len(properties) > 1:
            blocks = {_first_block(str(p.get("InChIKey", ""))) for p in properties}
            blocks.discard(None)
            if len(blocks) > 1:
                return OutsideRecord(None, None, None, "clean_miss")
        entry = properties[0]
        mass = entry.get("MonoisotopicMass")
        try:
            mass_value = float(mass) if mass is not None else None
        except (TypeError, ValueError):
            mass_value = None
        return OutsideRecord(
            block=_first_block(str(entry.get("InChIKey", ""))) or None,
            formula=str(entry.get("MolecularFormula", "")) or None,
            mass=mass_value,
            status="success",
        )


def classify(
    necs_block: str,
    cohort_block: str,
    outside_necs: OutsideRecord,
    outside_cohort: OutsideRecord,
) -> tuple[str, str]:
    """Return ``(outcome, rationale)`` for one non-certified case.

    The order of the checks is the substance of this function. The outside source is asked FIRST
    whether it backs each side, because "the two recorded keys disagree but the outside source gives
    the same structure for both names" is a defect in one of the recorded keys, not a tautomer
    artifact. Only once the outside source has corroborated BOTH sides at different structures does
    the composition test get to excuse the difference. Running the artifact test first would
    mis-label every gold defect whose true structure happens to share a formula with the cohort's.
    """
    if not outside_necs.resolved and not outside_cohort.resolved:
        return (
            "outside_source_unresolved",
            "PubChem's name index resolved neither side, so nothing is concluded and the refusal "
            "stands as unadjudicated",
        )

    necs_agrees = outside_necs.resolved and outside_necs.block == necs_block
    cohort_agrees = outside_cohort.resolved and outside_cohort.block == cohort_block

    if necs_agrees and cohort_agrees:
        if necs_block == cohort_block:
            return (
                "outside_source_unresolved",
                "the outside source backs both sides at the SAME structure, so this case was not a "
                "structural disagreement and its refusal came from somewhere else",
            )
        if compositions_equivalent(outside_necs, outside_cohort):
            return (
                "tautomer_or_charge_or_salt_artifact",
                f"the outside source corroborates both sides, and they share a heavy-atom "
                f"composition ({outside_necs.formula} vs {outside_cohort.formula}) with a mass gap "
                f"consistent with a hydrogen count difference. The InChIKey first blocks differ "
                f"because the first block is not invariant to tautomer, protonation or salt form",
            )
        return (
            "genuine_structural_disagreement",
            f"the outside source corroborates BOTH sides at different compositions "
            f"({outside_necs.formula} vs {outside_cohort.formula}), so the two cohorts genuinely "
            "measured different molecules and the link is wrong",
        )

    if cohort_agrees and outside_necs.resolved and not necs_agrees:
        return (
            "necs_gold_suspect",
            f"the outside source gives {outside_necs.block} for the NECS name while the curated "
            f"gold says {necs_block or '(none)'}; the cohort side is corroborated. This is the "
            "documented NECS gold defect showing up, not a BioMapper error",
        )
    if necs_agrees and outside_cohort.resolved and not cohort_agrees:
        return (
            "cohort_id_suspect",
            f"the outside source gives {outside_cohort.block} for the cohort name while its vendor "
            f"identifier resolves to {cohort_block or '(none)'}; the NECS side is corroborated",
        )
    if outside_necs.resolved and outside_cohort.resolved:
        return (
            "both_sides_disagree_with_outside",
            "neither side matches the outside source, so the disagreement is upstream of both and "
            "needs a hand check",
        )
    return (
        "outside_source_unresolved",
        "only one side resolved against the outside source, which is not enough to place the "
        "defect",
    )


def readjudicate(cases: pd.DataFrame, resolver: OutsideResolver) -> pd.DataFrame:
    """Re-adjudicate every non-certified case. Returns the cases with outcome columns appended."""
    open_cases = cases[cases["verdict"] != "certified"].copy()
    rows: list[dict[str, Any]] = []
    for _, case in open_cases.iterrows():
        # Series.to_dict() is keyed Hashable; the case table's columns are all strings, so re-key
        # explicitly rather than assuming it.
        case_fields: dict[str, Any] = {str(k): v for k, v in case.to_dict().items()}
        refusal_class = str(case.get("refusal_class", ""))
        # A refusal from a panel that carries no structure-resolvable identifier is a property of
        # the source, not a case to chase. Recorded as such rather than sent to PubChem, which
        # would invent an adjudication the certificate never made.
        if refusal_class in (
            "cohort_panel_has_no_structure_resolvable_id",
            "no_independent_structure_either_side",
        ):
            rows.append(
                {
                    **case_fields,
                    "readjudication": "not_adjudicable_by_construction",
                    "rationale": (
                        "the cohort panel carries no structure-resolvable vendor identifier, so no "
                        "independent structure exists on that side at all"
                    ),
                    "outside_necs_block": "",
                    "outside_necs_formula": "",
                    "outside_cohort_block": "",
                    "outside_cohort_formula": "",
                }
            )
            continue
        outside_necs = resolver.by_name(str(case["necs_name"]))
        outside_cohort = resolver.by_name(str(case["cohort_name"]))
        outcome, rationale = classify(
            str(case.get("necs_block", "")),
            str(case.get("cohort_block", "")),
            outside_necs,
            outside_cohort,
        )
        rows.append(
            {
                **case_fields,
                "readjudication": outcome,
                "rationale": rationale,
                "outside_necs_block": outside_necs.block or "",
                "outside_necs_formula": outside_necs.formula or "",
                "outside_necs_status": outside_necs.status,
                "outside_cohort_block": outside_cohort.block or "",
                "outside_cohort_formula": outside_cohort.formula or "",
                "outside_cohort_status": outside_cohort.status,
            }
        )
    return pd.DataFrame(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m biomapper.benchmarks.cross_cohort_readjudicate",
        description="Re-adjudicate non-certified cross-cohort cases against an outside source.",
    )
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--cohort", default="arivale")
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Adjudicate at most N adjudicable cases (0 = all). A truncated run says so in the "
        "output rather than looking complete.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    cases_path = args.run_dir / f"certificate_cases_{args.cohort}.csv"
    if not cases_path.exists():
        print(f"[fatal] {cases_path} is missing; run cross_cohort_certify first", file=sys.stderr)
        return 2
    cases = pd.read_csv(cases_path, dtype=str).fillna("")
    if args.limit:
        adjudicable = cases[cases["verdict"] != "certified"]
        keep = set(adjudicable.head(args.limit).index)
        cases = cases.loc[cases.index.isin(keep) | (cases["verdict"] == "certified")]
        print(f"[warn] limited to {args.limit} adjudicable cases; coverage is partial", flush=True)

    resolved = readjudicate(cases, OutsideResolver())
    out = args.run_dir / f"readjudication_{args.cohort}.csv"
    resolved.to_csv(out, index=False)
    tally = dict(Counter(resolved["readjudication"].tolist()))
    summary = {
        "cohort": args.cohort,
        "cases_examined": int(len(resolved)),
        "outcomes": tally,
        "limit_applied": args.limit or None,
        "mass_tolerance_da": MASS_TOLERANCE_DA,
        "outside_source": (
            "PubChem PUG-REST name index (InChIKey, MolecularFormula, MonoisotopicMass)"
        ),
        "caveat": (
            "a first-block comparison over-flags tautomers and charge variants as refuted, and "
            "silently accepts a stereoisomer error when a side has no stereo layer. The second "
            "failure mode is not recoverable by this re-adjudication."
        ),
    }
    (args.run_dir / f"readjudication_{args.cohort}.json").write_text(
        json.dumps(summary, indent=2, default=str)
    )
    for outcome, count in sorted(tally.items(), key=lambda kv: -kv[1]):
        print(f"[readjudicate] {outcome}: {count}", flush=True)
    print(f"[done] {out}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
