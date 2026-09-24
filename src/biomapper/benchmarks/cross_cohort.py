"""Cross-cohort (Monti/NECS) harmonization arm, run against a deployment.

Ported from the biomapper2 engine's ``studies/external_benchmarks/cross_cohort_run.py``
(``origin/dev``, commit ``1ffb571e54fe028ef0ae4e748fc2e7ec093ee603``), rewired from
``biomapper2.mapper.Mapper`` to :class:`biomapper.benchmarks.api_mapper.ApiMapper` exactly as the
11 suite arms were. Running against a deployment is the point: provenance then pins the graph that
actually served the answers rather than a client checkout.

What this measures, and what it does not:

* **Arm M** is BioMapper, names only. Each cohort panel's names go through the mapping API and the
  identifier-only CURIE sets are intersected NECS<->cohort. Linking is set intersection and nothing
  else. Structure-encoding namespaces (``INCHIKEY`` / ``INCHI`` / ``SMILES``) are excluded from the
  linker, because linking on a structure hash would make any downstream structural certificate
  circular and precision 100% by construction.
* **Arm B** is Monti's own published method, reconstructed on the identical row set so the
  comparison is like for like. It is a controlled variable we compute, so the claim that matters is
  against **Monti published**, which we did not compute.
* **Certification is a separate, later judgement** and is deliberately not what forms a link. Three
  of the four cohorts plus NECS ship names only, with no vendor identifiers, so
  ``CohortPanel.certifiable`` is ``False`` for them: their links are countable but never
  structurally certifiable. BLSA is the one to state explicitly in a table, because a blank
  certification column there is a property of the source panel, not a failed certification.

The batch-order hazard the suite already hit applies here unchanged: predictions are joined to
their input rows BY POSITION, so a reordered batch response would attribute one entity's answer to
another. :class:`ApiMapper` escalates that to ``BatchOrderMismatchError`` and aborts rather than
realigning, and :func:`resolve_panel` re-asserts the alignment on the assembled frame. A name is
not unique across a panel in general, so there is nothing safe to realign by.

Usage (resolution is the long pole, so panels are checkpointed and can run as separate processes)::

    python -m biomapper.benchmarks.cross_cohort --out-dir RUN --panel necs
    python -m biomapper.benchmarks.cross_cohort --out-dir RUN --panel arivale
    ...
    python -m biomapper.benchmarks.cross_cohort --out-dir RUN --link-only
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd

from biomapper.benchmarks.adapters.cohort_panel import (
    ARIVALE,
    SPREADSHEET_COHORTS,
    CohortPanel,
    load_cohort_panel,
)
from biomapper.benchmarks.api_mapper import ApiMapper
from biomapper.benchmarks.provenance import (
    DEFAULT_KESTREL_URL,
    build_run_provenance,
    utc_stamp,
)
from biomapper.benchmarks.scorers.arm_b_baseline import (
    MONTI_PUBLISHED,
    MONTI_PUBLISHED_PROVENANCE,
    MONTI_PUBLISHED_SUPERSEDED,
    arm_b_overlap,
)
from biomapper.benchmarks.scorers.cross_cohort_overlap import (
    link_by_intersection,
    row_curie_set,
)
from biomapper.client import DEFAULT_BASE_URL

# Settled: the four cohorts linked to NECS (engine ``cross_cohort_run.py:33``).
COHORTS: tuple[str, ...] = ("arivale", "xuetal", "llfs", "blsa")
PANELS: tuple[str, ...] = ("necs", *COHORTS)

# Default local source artifacts. Each is SHA-pinned into the manifest; a missing one fails loudly
# rather than producing a panel of zero rows, which an adapter would otherwise read as success.
DEFAULT_SPREADSHEET = (
    Path.home()
    / ".claude/uploads/2390692a-d27f-5fa2-bb74-9b032f2d5009/99293a9b-datasets_metabolites.xlsx"
)
DEFAULT_ARIVALE_XLSX = (
    Path.home()
    / "external_benchmark_runs/arivale_public_panel_20260804/watanabe2023_supp_data2_analytes.xlsx"
)
DEFAULT_REFMET_CACHE = (
    Path.home() / "external_benchmark_runs/cohort_panels_20260804/necs_refmet_convert_cache.tsv"
)

# Published panel sizes, for an honest denominator. LLFS is the one that matters: the paper profiled
# 408 metabolites (188 lipid + 220 polar) and the 4-cohort spreadsheet ships only the 364
# RefMet-standardizable subset, so the gap is pre-filtering and is reported. The "345" figure that
# has circulated appears nowhere in the paper and is not used here.
PUBLISHED_PANEL_SIZES: dict[str, int] = {
    "necs": 1213,
    "arivale": 626,
    "xuetal": 821,
    "llfs": 408,
    "blsa": 468,
}

# The request shape the engine driver used, reproduced so Arm M is comparable across runs.
REQUEST = {
    "entity_type": "metabolite",
    "vocab": "CHEBI",
    "annotation_mode": "all",
    "provided_id_columns": [],
}


class PanelAlignmentError(RuntimeError):
    """The mapped frame's query column no longer matches the panel it was built from.

    Same class of failure as ``BatchOrderMismatchError`` and refused for the same reason: every
    downstream count joins a prediction to its input row by position, so a shifted frame silently
    attributes one metabolite's CURIE set to another. Nothing is realigned, because a panel name is
    not guaranteed unique across a cohort.
    """


def sha256_path(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def client_repo_provenance() -> dict[str, str | bool | None]:
    """Commit and cleanliness of the installed client, best effort.

    ``biomapper_version`` comes from installed package metadata, which does not move when a working
    tree does. A run executed from an editable checkout with uncommitted modules would otherwise be
    attributed to whatever version the metadata last recorded, so the commit and the dirty flag are
    pinned alongside it. ``None`` for a wheel install with no repository, which is honest rather
    than invented.
    """
    import subprocess

    repo = Path(__file__).resolve().parents[3]
    try:
        commit = subprocess.run(
            ["git", "-C", str(repo), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=10,
            check=True,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "-C", str(repo), "status", "--porcelain"],
            capture_output=True,
            text=True,
            timeout=10,
            check=True,
        ).stdout.strip()
    except Exception:  # noqa: BLE001 - provenance must never abort a run
        return {"repo": str(repo), "commit": None, "dirty": None}
    return {"repo": str(repo), "commit": commit, "dirty": bool(status)}


def load_panels(spreadsheet: Path, arivale_xlsx: Path) -> dict[str, CohortPanel]:
    """Load all five panels through the ported engine adapter, so exclusions are counted alike."""
    for path in (spreadsheet, arivale_xlsx):
        if not path.exists():
            raise FileNotFoundError(f"cross-cohort source artifact is missing: {path}")
    sheet = pd.read_excel(spreadsheet, sheet_name="Sheet 1", dtype=str).fillna("")
    panels = {
        key: load_cohort_panel(sheet, SPREADSHEET_COHORTS[key])
        for key in ("necs", "xuetal", "llfs", "blsa")
    }
    arivale = pd.read_excel(arivale_xlsx, sheet_name="Arivale_Metabolomics", dtype=str).fillna("")
    panels["arivale"] = load_cohort_panel(arivale, ARIVALE)
    return panels


def resolve_panel(mapper: ApiMapper, panel: CohortPanel, out_dir: Path, label: str) -> pd.DataFrame:
    """Resolve one panel name-only; checkpoint the mapped TSV so a mid-run 5xx loses nothing."""
    dest = out_dir / f"{label}_MAPPED.tsv"
    names = panel.names
    if dest.exists():
        mapped = pd.read_csv(dest, sep="\t", dtype=str).fillna("")
        assert_alignment(mapped, names, label)
        print(f"[resolve] {label}: reusing checkpoint {dest.name} ({len(mapped)} rows)", flush=True)
        return mapped

    print(f"[resolve] {label}: n={len(names)} against {mapper.endpoint}", flush=True)
    out_tsv, stats = mapper.map_dataset_to_kg(
        dataset=pd.DataFrame({"name": names}),
        entity_type=str(REQUEST["entity_type"]),
        name_column="name",
        provided_id_columns=[],
        vocab=str(REQUEST["vocab"]),
        annotation_mode=str(REQUEST["annotation_mode"]),
        output_dir=out_dir,
        output_prefix=label,
    )
    (out_dir / f"{label}_stats.json").write_text(json.dumps(stats, indent=2, default=str))
    mapped = pd.read_csv(out_tsv, sep="\t", dtype=str).fillna("")
    assert_alignment(mapped, names, label)
    return mapped


def assert_alignment(mapped: pd.DataFrame, names: list[str], label: str) -> None:
    """Re-assert position alignment between the mapped frame and the panel that produced it."""
    got = [str(v) for v in mapped["name"].tolist()]
    if got != [str(n) for n in names]:
        first = next(
            (i for i, (a, b) in enumerate(zip(got, names, strict=False)) if str(a) != str(b)),
            min(len(got), len(names)),
        )
        raise PanelAlignmentError(
            f"{label}: mapped frame does not match the panel by position "
            f"(n={len(got)} vs {len(names)}, first divergence at row {first}). Refusing to score: "
            "every count here joins a prediction to its input row positionally."
        )


def errored_names(mapped: pd.DataFrame) -> list[str]:
    """Names whose mapping call itself failed, in panel order.

    An errored row is NOT an unresolved row. Both come back with an empty CURIE set, so folding
    them together would report "this metabolite did not resolve" when the truth is "we never got an
    answer", and it would understate resolution by exactly the number of rows the deployment
    dropped. The deployment returned 5xx under concurrent load during this run, so this is a live
    hazard and not a theoretical one.
    """
    if "mapping_error" not in mapped.columns:
        return []
    return [
        str(row["name"])
        for _, row in mapped.iterrows()
        if str(row.get("mapping_error", "")).strip() not in ("", "nan")
    ]


def repair_errored_rows(
    mapper: ApiMapper, mapped: pd.DataFrame, out_dir: Path, label: str
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Re-map only the errored rows and write them back into their original positions.

    Positional write-back is what keeps this safe: the repaired values land on the same row index
    they came from, so the panel alignment the rest of the module depends on is preserved. A row
    that errors again stays errored and is reported.
    """
    failures = errored_names(mapped)
    if not failures:
        return mapped, {"attempted": 0, "recovered": 0, "still_errored": 0}

    print(f"[repair] {label}: re-mapping {len(failures)} errored rows", flush=True)
    indices = [
        i
        for i, (_, row) in enumerate(mapped.iterrows())
        if str(row.get("mapping_error", "")).strip() not in ("", "nan")
    ]
    retry_frame = pd.DataFrame({"name": [str(mapped.iloc[i]["name"]) for i in indices]})
    out_tsv, _ = mapper.map_dataset_to_kg(
        dataset=retry_frame,
        entity_type=str(REQUEST["entity_type"]),
        name_column="name",
        provided_id_columns=[],
        vocab=str(REQUEST["vocab"]),
        annotation_mode=str(REQUEST["annotation_mode"]),
        output_dir=out_dir,
        output_prefix=f"{label}_repair",
    )
    repaired = pd.read_csv(out_tsv, sep="\t", dtype=str).fillna("")
    if len(repaired) != len(indices):
        raise PanelAlignmentError(
            f"{label}: repair pass returned {len(repaired)} rows for {len(indices)} requests; "
            "refusing to write back a misaligned frame."
        )
    shared_columns = [c for c in mapped.columns if c in repaired.columns]
    # Rebuilt column-wise from plain lists rather than through positional cell assignment: the
    # write-back target is an index into the ORIGINAL panel order, and going via lists keeps that
    # explicit instead of relying on a frame's internal column positions staying put.
    columns = {c: mapped[c].tolist() for c in shared_columns}
    for position, index in enumerate(indices):
        source = repaired.iloc[position]
        if str(source["name"]) != str(mapped.iloc[index]["name"]):
            raise PanelAlignmentError(
                f"{label}: repair row {position} is {source['name']!r} but target row {index} is "
                f"{mapped.iloc[index]['name']!r}; refusing to write back."
            )
        for column in shared_columns:
            columns[column][index] = source[column]
    for column in shared_columns:
        mapped[column] = columns[column]
    mapped.to_csv(out_dir / f"{label}_MAPPED.tsv", sep="\t", index=False)
    still = errored_names(mapped)
    summary = {
        "attempted": len(failures),
        "recovered": len(failures) - len(still),
        "still_errored": len(still),
        "still_errored_names": still[:50],
    }
    print(
        f"[repair] {label}: recovered {summary['recovered']}/{summary['attempted']}, "
        f"{summary['still_errored']} still errored",
        flush=True,
    )
    return mapped, summary


def curies_by_name(mapped: pd.DataFrame, label: str) -> dict[str, frozenset[str]]:
    """``{name: identifier-only CURIE set}`` via the ported cross-cohort scorer."""
    out: dict[str, frozenset[str]] = {}
    for _, row in mapped.iterrows():
        name = str(row.get("name", "")).strip()
        if not name:
            continue
        if name in out:
            raise PanelAlignmentError(
                f"{label}: duplicate panel name {name!r} reached the linker. One of the two rows "
                "would be dropped silently; the adapter is supposed to de-duplicate."
            )
        out[name] = row_curie_set(row)
    return out


def certificate_tally(mapped: pd.DataFrame) -> dict[str, dict[str, int]]:
    """Tallies of the API's resolution-certificate columns for one panel.

    This is the API's ``ResolutionCertificateModel``, which reads the chosen KRAKEN node's own
    InChIKey first and only falls back to an external source. It is therefore NOT fully
    KG-independent and is reported as a resolution property, never as a cross-cohort link
    certificate. The link certificate is
    :func:`biomapper.benchmarks.scorers.link_certificate.certify_link`.
    """
    columns = (
        "certificate_state",
        "certificate_structure_status",
        "certificate_refusal_reason",
        "certificate_tier_b_outcome",
        "certificate_refmet_availability",
        "certificate_lipid_resolution_level",
        "certificate_tier_b_snapshot_version",
        "certificate_independent_source",
        "certificate_independent_of_selection",
        "chosen_kg_id_review",
    )
    tallies: dict[str, dict[str, int]] = {}
    for column in columns:
        if column not in mapped.columns:
            continue
        values = [str(v).strip() or "(blank)" for v in mapped[column].tolist()]
        tallies[column] = dict(sorted(Counter(values).items(), key=lambda kv: -kv[1]))
    return tallies


def chosen_prefix_tally(mapped: pd.DataFrame) -> dict[str, int]:
    """Namespace of the chosen node, per panel. Production answers with ``RM:*`` for most rows."""
    prefixes = [
        str(v).split(":", 1)[0]
        if ":" in str(v)
        else ("(unresolved)" if not str(v).strip() else str(v))
        for v in mapped.get("chosen_kg_id", pd.Series(dtype=str)).tolist()
    ]
    return dict(sorted(Counter(prefixes).items(), key=lambda kv: -kv[1]))


def load_refmet_map(path: Path) -> dict[str, str]:
    if not path.exists():
        raise FileNotFoundError(f"RefMet conversion cache is missing: {path}")
    with path.open() as handle:
        return {row["input"]: row["refmet"] for row in csv.DictReader(handle, delimiter="\t")}


def cross_check_harmonize(
    a_curies: dict[str, frozenset[str]], b_curies: dict[str, frozenset[str]]
) -> dict[str, Any]:
    """Confirm the generic ``biomapper.harmonize`` primitive agrees with the cohort scorer.

    Two implementations of one rule is a real risk: the package ships a pandas-free
    ``harmonize.link_by_intersection`` and this cohort-shaped port. If they ever diverge, a number
    changes for a reason unrelated to the engine, so the agreement is asserted per run rather than
    assumed.
    """
    from biomapper.harmonize import link_by_intersection as generic_link

    engine = link_by_intersection(a_curies, b_curies)
    generic = generic_link(a_curies, b_curies)
    agree = (
        engine.n_links == generic.n_links
        and engine.n_a_linked == generic.n_a_linked
        and engine.n_b_linked == generic.n_b_linked
        and engine.n_a_comparable == generic.n_a_comparable
        and engine.n_b_comparable == generic.n_b_comparable
    )
    return {
        "agree": agree,
        "scorer_n_links": engine.n_links,
        "harmonize_n_links": generic.n_links,
    }


def run_links(
    panels: dict[str, CohortPanel],
    curies: dict[str, dict[str, frozenset[str]]],
    refmet_map: dict[str, str],
    out_dir: Path,
) -> dict[str, Any]:
    necs_names = panels["necs"].names
    results: dict[str, Any] = {}
    for cohort in COHORTS:
        overlap = link_by_intersection(curies["necs"], curies[cohort])
        derived = arm_b_overlap(cohort, necs_names, panels[cohort].names, refmet_map=refmet_map)
        published = MONTI_PUBLISHED[cohort]
        results[cohort] = {
            "certifiable": panels[cohort].certifiable,
            "certifiability_note": (
                "names only, no vendor identifiers: links are countable but never structurally "
                "certifiable"
                if not panels[cohort].certifiable
                else f"vendor id namespaces present: {list(panels[cohort].id_columns)}"
            ),
            "panel_n_loaded": len(panels[cohort].names),
            "panel_n_published": PUBLISHED_PANEL_SIZES[cohort],
            "panel_card": panels[cohort].card,
            "arm_m_links": overlap.n_links,
            "arm_m_necs_linked": overlap.n_a_linked,
            "arm_m_cohort_linked": overlap.n_b_linked,
            "necs_comparable": overlap.n_a_comparable,
            "cohort_comparable": overlap.n_b_comparable,
            "necs_unresolved": len(curies["necs"]) - overlap.n_a_comparable,
            "cohort_unresolved": len(curies[cohort]) - overlap.n_b_comparable,
            "arm_b_rederived": derived.count,
            "arm_b_method": derived.method,
            "monti_published": published,
            "monti_published_provenance": MONTI_PUBLISHED_PROVENANCE[cohort],
            "monti_published_superseded_value": MONTI_PUBLISHED_SUPERSEDED[cohort],
            "arm_b_gap_to_published": derived.count - published,
            "arm_m_vs_arm_b": overlap.n_a_linked - derived.count,
            "arm_m_vs_published": overlap.n_a_linked - published,
            "harmonize_cross_check": cross_check_harmonize(curies["necs"], curies[cohort]),
        }
        print(
            f"[pair] NECS<->{cohort}: Arm-M necs-linked={overlap.n_a_linked} "
            f"cohort-linked={overlap.n_b_linked} links={overlap.n_links} "
            f"Arm-B={derived.count} published={published}",
            flush=True,
        )
        pd.DataFrame(
            [
                {
                    "necs_name": link.a_name,
                    f"{cohort}_name": link.b_name,
                    "shared_curies": "|".join(sorted(link.shared)),
                }
                for link in overlap.links
            ]
        ).to_csv(out_dir / f"links_necs_{cohort}.csv", index=False)
        pd.DataFrame({"name": [n for n, s in curies[cohort].items() if not s]}).to_csv(
            out_dir / f"unresolved_{cohort}.csv", index=False
        )
    pd.DataFrame({"name": [n for n, s in curies["necs"].items() if not s]}).to_csv(
        out_dir / "unresolved_necs.csv", index=False
    )
    return results


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m biomapper.benchmarks.cross_cohort",
        description="Cross-cohort (Monti/NECS) harmonization arm against a deployment.",
    )
    parser.add_argument(
        "--out-dir", type=Path, default=None, help="Run directory (default: timestamped)"
    )
    parser.add_argument("--endpoint", default=DEFAULT_BASE_URL, help="BioMapper API root")
    parser.add_argument(
        "--kestrel-url", default=DEFAULT_KESTREL_URL, help="Kestrel root for /health"
    )
    parser.add_argument("--api-key", default=None, help="Omit for a keyless deployment")
    parser.add_argument("--batch-size", type=int, default=20)
    parser.add_argument("--timeout", type=float, default=1800.0)
    parser.add_argument("--spreadsheet", type=Path, default=DEFAULT_SPREADSHEET)
    parser.add_argument("--arivale-xlsx", type=Path, default=DEFAULT_ARIVALE_XLSX)
    parser.add_argument("--refmet-cache", type=Path, default=DEFAULT_REFMET_CACHE)
    parser.add_argument(
        "--panel",
        choices=PANELS,
        default=None,
        help="Resolve only this panel and exit (checkpointed; run one process per panel)",
    )
    parser.add_argument(
        "--link-only",
        action="store_true",
        help="Skip resolution; link and score from the existing per-panel checkpoints",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    ts = utc_stamp()
    out_dir = args.out_dir or (Path.home() / f"external_benchmark_runs/cross_cohort_{ts}")
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[run] {out_dir}", flush=True)

    panels = load_panels(args.spreadsheet, args.arivale_xlsx)
    for key, panel in panels.items():
        print(
            f"[panel] {key}: n={len(panel.names)} certifiable={panel.certifiable} "
            f"exclusions={panel.card['exclusions']}",
            flush=True,
        )

    if args.panel is not None and not args.link_only:
        mapper = ApiMapper(
            args.endpoint,
            api_key=args.api_key,
            batch_size=args.batch_size,
            timeout=args.timeout,
        )
        resolve_panel(mapper, panels[args.panel], out_dir, args.panel)
        (out_dir / f"{args.panel}_counters.json").write_text(
            json.dumps(mapper.counters.snapshot(), indent=2, default=str)
        )
        print(f"[done] panel {args.panel} -> {out_dir}/{args.panel}_MAPPED.tsv", flush=True)
        return 0

    provenance = build_run_provenance(
        api_endpoint=args.endpoint, kestrel_url=args.kestrel_url, run_id=f"cross_cohort_{ts}"
    )
    if not provenance.pinned:
        print(
            f"[fatal] provenance is not pinned ({provenance.health_error}). Refusing to emit "
            "numbers that cannot be attributed to a graph build.",
            file=sys.stderr,
        )
        return 2
    print(
        f"[kg] kestrel={provenance.kestrel_version} kg={provenance.kg_build.kg_version} "
        f"biolink={provenance.kg_build.biolink_version} "
        f"commit={provenance.kg_build.git_commit[:8]} built={provenance.kg_build.build_timestamp}",
        flush=True,
    )

    mapper = ApiMapper(
        args.endpoint, api_key=args.api_key, batch_size=args.batch_size, timeout=args.timeout
    )
    mapped: dict[str, pd.DataFrame] = {}
    curies: dict[str, dict[str, frozenset[str]]] = {}
    repairs: dict[str, dict[str, Any]] = {}
    for label in PANELS:
        checkpoint = out_dir / f"{label}_MAPPED.tsv"
        if args.link_only and not checkpoint.exists():
            print(f"[fatal] --link-only but {checkpoint} is missing", file=sys.stderr)
            return 2
        frame = resolve_panel(mapper, panels[label], out_dir, label)
        frame, repairs[label] = repair_errored_rows(mapper, frame, out_dir, label)
        mapped[label] = frame
        curies[label] = curies_by_name(frame, label)
        print(
            f"[resolve] {label}: {sum(1 for s in curies[label].values() if s)}/"
            f"{len(curies[label])} with a non-empty identifier-only CURIE set, "
            f"{repairs[label]['still_errored']} errored",
            flush=True,
        )

    unrecovered = {k: v["still_errored"] for k, v in repairs.items() if v["still_errored"]}
    if unrecovered:
        print(
            f"[warn] rows the deployment never answered, after the repair pass: {unrecovered}. "
            "These are recorded as errors, NOT as unresolved metabolites, and the affected panels' "
            "resolution counts are a floor. Re-run those panels before publishing.",
            file=sys.stderr,
        )

    refmet_map = load_refmet_map(args.refmet_cache)
    results = run_links(panels, curies, refmet_map, out_dir)

    manifest = {
        "arm": "M (BioMapper, names only) vs B (Monti method, re-derived) vs Monti published",
        "run_id": provenance.run_id,
        "created_utc": ts,
        "package": {
            "name": "biomapper",
            "version": provenance.biomapper_version,
            "client_repo": client_repo_provenance(),
            "linker": "biomapper.benchmarks.scorers.cross_cohort_overlap.link_by_intersection",
            "linker_generic": "biomapper.harmonize.link_by_intersection (cross-checked per pair)",
            "baseline": "biomapper.benchmarks.scorers.arm_b_baseline.arm_b_overlap",
            "panels": "biomapper.benchmarks.adapters.cohort_panel",
            "mapper": "biomapper.benchmarks.api_mapper.ApiMapper",
            "provenance": "biomapper.benchmarks.provenance.build_run_provenance",
        },
        "ported_from": {
            "repo": "biomapper2",
            "ref": "origin/dev",
            "commit": "1ffb571e54fe028ef0ae4e748fc2e7ec093ee603",
            "modules": [
                "studies/external_benchmarks/adapters/cohort_panel.py",
                "studies/external_benchmarks/scorers/arm_b_baseline.py",
                "studies/external_benchmarks/scorers/cross_cohort_overlap.py",
                "studies/external_benchmarks/scorers/link_certificate.py",
                "studies/external_benchmarks/cross_cohort_run.py (driver, reference)",
            ],
        },
        "provenance": json.loads(provenance.model_dump_json()),
        "request": {
            **REQUEST,
            "batch_size": args.batch_size,
            "endpoint": args.endpoint,
            "api_key_sent": args.api_key is not None,
        },
        "sources": {
            "spreadsheet": {
                "path": str(args.spreadsheet),
                "sha256": sha256_path(args.spreadsheet),
            },
            "arivale_xlsx": {
                "path": str(args.arivale_xlsx),
                "sha256": sha256_path(args.arivale_xlsx),
            },
            "refmet_cache": {
                "path": str(args.refmet_cache),
                "sha256": sha256_path(args.refmet_cache),
            },
            "monti_reference": {
                "doi": "10.1007/s11357-026-02174-2",
                "gold": "MOESM5 (Supplementary Table 5), 1495 metabolites",
                "published_overlaps_read_from": (
                    "Methods, 'Datasets harmonization' and the per-cohort descriptions. "
                    "NOT Table 2, "
                    "which is 'Age-only markers'."
                ),
            },
        },
        "panel_sizes_loaded": {k: len(v.names) for k, v in panels.items()},
        "panel_sizes_published": PUBLISHED_PANEL_SIZES,
        "certifiability": {
            k: {"certifiable": v.certifiable, "id_namespaces": list(v.id_columns)}
            for k, v in panels.items()
        },
        "certificate_tallies": {k: certificate_tally(v) for k, v in mapped.items()},
        "chosen_namespace_tallies": {k: chosen_prefix_tally(v) for k, v in mapped.items()},
        "mapping_errors": {
            "repair_pass": repairs,
            "unrecovered_by_panel": unrecovered,
            "note": (
                "An errored row is a row the deployment never answered. It is counted apart "
                "from an unresolved row, because both carry an empty CURIE set and merging them "
                "would "
                "report a non-resolution that never happened. A panel with unrecovered errors "
                "has a resolution count that is a floor, not a measurement."
            ),
        },
        "circularity_register": {
            "kraken_sources": provenance.kg_build.sources,
            "kraken_source_versions": provenance.kg_build.source_versions,
            "necs_gold_source": (
                "Metabolon CHEMICAL_NAME vendor curation (MOESM5). Not a KRAKEN-ingested "
                "vocabulary, so the NECS name gold is not circular through the graph."
            ),
            "refmet_ingested": "refmet" in {s.lower() for s in provenance.kg_build.sources},
            "refmet_caveat": (
                "RefMet is both ingested into KRAKEN and the resolver's source-weighting target, "
                "and production answers this arm predominantly with RM:* nodes. Any number whose "
                "gold is a RefMet name is coverage, not accuracy."
            ),
        },
        "results": results,
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, default=str))
    print(f"[done] {out_dir}/manifest.json", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
