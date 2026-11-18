#!/usr/bin/env python3
"""Team C runner for Dataset B proteome harmonization via BioMapper 1.9.0.

Runs ``HarmonizeJobs.submit_many(..., entity="proteins", chunk_size=10, capture_max_entities=500)``
with raw unedited ``twinsuk`` strings (default ``TargetFullName``, configurable via
``--name-column``), automatically resumes stalled/paused jobs and retries partial runs until zero
errored rows remain, performs deterministic candidate review via ``protein_review`` +
``SuggestionStore``, and writes ``dataset_b_proteins/choices.tsv`` via
``biomapper.mcp.judgments.write_choices`` alongside ``run_manifest.json`` and executable ``rerun.sh``.
"""

from __future__ import annotations

import argparse
import asyncio
import csv
import hashlib
import json
import os
import shutil
import stat
import subprocess
import time
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import anyio

from biomapper._pins import collect_pins
from biomapper._version import resolve_version
from biomapper.extras.protein_review import (
    find_suspect_twinsuk_rows,
    load_twinsuk_targets,
)
from biomapper.extras.protein_review import (
    decide_review as decide_protein_review,
)
from biomapper.jobs import CohortSpec, HarmonizeJobs
from biomapper.mcp.candidates import RefCandidates
from biomapper.mcp.judgments import (
    CHOICES_MARKER,
    write_choices,
)
from biomapper.mcp.references import RefRecord
from biomapper.models import MappingResult

HARMONIZE_ROOT = Path("/home/trentleslie/agy-team-competition/harmonize")
WORKTREE_DIR = HARMONIZE_ROOT / "worktrees" / "team-c"
DEFAULT_INPUTS_DIR = HARMONIZE_ROOT / "shared" / "inputs"
DEFAULT_SEED_STORE_DIR = HARMONIZE_ROOT / "shared" / "seed_store"
BM_PYTHON = HARMONIZE_ROOT / "bin" / "bm-python"
COHORT_FILES: dict[str, str] = {
    "ukbb": "ukbb_proteins.tsv",
    "arivale": "arivale_proteins.tsv",
    "twinsuk": "twinsuk_somascan_targets.tsv",
}


def sha256_file(path: Path) -> str:
    """Return lowercase hex SHA-256 of ``path``."""
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def snapshot_report_files(subrun_dir: Path) -> dict[str, str]:
    """Hash ``report.json``, JSON files, and main TSVs before sidecar ``choices.tsv`` is written."""
    hashes: dict[str, str] = {}
    for p in sorted(subrun_dir.iterdir()):
        if not p.is_file() or p.name == "choices.tsv":
            continue
        if p.suffix in (".json", ".tsv"):
            hashes[p.name] = sha256_file(p)
    return hashes


def seed_store_if_needed(store_dir: Path, seed_store_dir: Path, *, remap_all: bool) -> bool:
    """Copy Olink seed store into ``store_dir`` unless ``remap_all`` or ``names_only`` exists."""
    store_dir.mkdir(parents=True, exist_ok=True, mode=0o700)
    os.chmod(store_dir, 0o700)
    if remap_all or (store_dir / "names_only").exists():
        return False
    if not seed_store_dir.is_dir():
        return False
    shutil.copytree(seed_store_dir, store_dir, dirs_exist_ok=True, symlinks=False)
    for root, dirs, files in os.walk(store_dir):
        os.chmod(root, 0o700)
        for d in dirs:
            os.chmod(Path(root) / d, 0o700)
        for f in files:
            os.chmod(Path(root) / f, 0o600)
    return True


def _slice_tsv(src: Path, dst: Path, max_rows: int) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    with src.open("r", encoding="utf-8", newline="") as fin:
        with dst.open("w", encoding="utf-8", newline="") as fout:
            for idx, line in enumerate(fin):
                if idx == 0 or idx <= max_rows:
                    fout.write(line)
                else:
                    break


def prepare_input_files(
    *,
    inputs_dir: Path,
    pilot: bool,
    slice_n: int | None,
    staging_dir: Path,
) -> tuple[dict[str, Path], dict[str, str]]:
    """Resolve canonical cohort TSV paths, upstream SHA-256 map, and optional sliced copies."""
    source_dir = inputs_dir / "pilot_20" if pilot else inputs_dir
    canonical_paths = {label: source_dir / fname for label, fname in COHORT_FILES.items()}
    upstream_hashes = {
        fname: sha256_file(canonical_paths[label]) for label, fname in COHORT_FILES.items()
    }
    if slice_n is not None and slice_n > 0:
        sliced_dir = staging_dir / "_sliced_inputs"
        run_paths: dict[str, Path] = {}
        for label, src in canonical_paths.items():
            dst = sliced_dir / src.name
            _slice_tsv(src, dst, slice_n)
            run_paths[label] = dst
        return run_paths, upstream_hashes
    return canonical_paths, upstream_hashes


def render_rerun_script(worktree_dir: Path) -> str:
    """Return the executable ``rerun.sh`` wrapper script contents."""
    python_bin = HARMONIZE_ROOT / "bin" / "bm-python"
    venv_dir = worktree_dir / ".venv"
    script_path = worktree_dir / "scripts" / "run_team_c.py"
    return (
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        f'exec "{python_bin}" --venv "{venv_dir}" "{script_path}" "$@"\n'
    )


def write_executable_rerun(path: Path, worktree_dir: Path) -> None:
    """Write ``rerun.sh`` and ensure it has executable permissions."""
    path.write_text(render_rerun_script(worktree_dir), encoding="utf-8")
    current_mode = path.stat().st_mode
    os.chmod(path, current_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _git_sha(worktree_dir: Path) -> str:
    proc = subprocess.run(
        ["git", "-C", str(worktree_dir), "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        check=False,
    )
    return proc.stdout.strip() if proc.returncode == 0 else "unknown"


def _mapping_result_from_report_row(query_name: str, row: Mapping[str, Any]) -> MappingResult:
    score_raw = row.get("confidence_score")
    score = float(score_raw) if isinstance(score_raw, int | float) else None
    eq_raw = row.get("kg_equivalent_ids") or {}
    eq = {str(k): [str(v) for v in vs] for k, vs in eq_raw.items()}
    chosen = row.get("chosen_kg_id")
    return MappingResult(
        query_name=query_name,
        resolved=bool(row.get("resolved")),
        primary_curie=str(chosen) if chosen else None,
        chosen_kg_id=str(chosen) if chosen else None,
        confidence_score=score,
        kg_equivalent_ids=eq,
        error=str(row["error"]) if row.get("error") else None,
    )


def verify_choices_alignment(
    choices_path: Path, expected_seqids: set[str], pre_hashes: dict[str, str], subrun_dir: Path
) -> dict[str, int]:
    """Verify every ``key`` in ``choices.tsv`` is a ``twinsuk`` ``SeqId`` and main files are untouched."""
    if not choices_path.is_file():
        raise RuntimeError(f"Missing choices file at {choices_path}")
    raw_lines = choices_path.read_text(encoding="utf-8").splitlines()
    if not raw_lines or raw_lines[0].rstrip("\r\n") != CHOICES_MARKER:
        raise RuntimeError(f"{choices_path} missing server-owned CHOICES_MARKER header")
    data_lines = [ln for ln in raw_lines if not ln.startswith("#")]
    counts = {"total": 0, "engine": 0, "candidate": 0, "none": 0}
    seen_keys: set[str] = set()
    for row in csv.DictReader(data_lines, delimiter="\t"):
        key = str(row.get("key") or "").strip()
        if not key or key not in expected_seqids:
            raise RuntimeError(f"choices.tsv key {key!r} is not a twinsuk SeqId in report.json")
        if key in seen_keys:
            raise RuntimeError(f"Duplicate key {key!r} in {choices_path}")
        seen_keys.add(key)
        src = str(row.get("choice_source") or "").strip()
        counts["total"] += 1
        if src in counts:
            counts[src] += 1

    post_hashes = snapshot_report_files(subrun_dir)
    if post_hashes != pre_hashes:
        changed = [k for k in set(pre_hashes) | set(post_hashes) if pre_hashes.get(k) != post_hashes.get(k)]
        raise RuntimeError(f"Report/TSV files modified during review: {changed}")
    return counts


async def run_dataset_b_team_c(
    *,
    output_dir: Path,
    store_dir: Path,
    pilot: bool = False,
    slice_n: int | None = None,
    remap_all: bool = False,
    name_column: str = "TargetFullName",
    inputs_dir: Path = DEFAULT_INPUTS_DIR,
    seed_store_dir: Path = DEFAULT_SEED_STORE_DIR,
    worktree_dir: Path = WORKTREE_DIR,
    base_url: str | None = None,
    kestrel_url: str | None = None,
    api_key: str | None = None,
    score_cutoff: float = 3.2,
    pace_s: float = 1.1,
    poll_s: float = 1.0,
    resume_wait_s: float = 2.0,
    systemic_backoff_s: float = 2.0,
    max_attempts: int = 5,
    chunk_size: int = 10,
    capture_max_entities: int = 500,
    capture_top_n: int = 25,
) -> dict[str, Any]:
    """Execute Dataset B harmonization, weak/suspect candidate review, and manifest placement."""
    output_dir = output_dir.resolve()
    store_dir = store_dir.resolve()
    inputs_dir = inputs_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    subrun_dir = output_dir / "dataset_b_proteins"
    if subrun_dir.exists() and any(subrun_dir.iterdir()):
        raise RuntimeError(
            f"{subrun_dir} already exists and is non-empty; pass a fresh --output-dir."
        )
    jobs_home = output_dir / ".biomapper_home"
    jobs_home.mkdir(parents=True, exist_ok=True, mode=0o700)
    os.chmod(jobs_home, 0o700)

    seeded = seed_store_if_needed(store_dir, seed_store_dir, remap_all=remap_all)
    cohort_paths, upstream_hashes = prepare_input_files(
        inputs_dir=inputs_dir,
        pilot=pilot,
        slice_n=slice_n,
        staging_dir=output_dir,
    )
    cohorts = [
        CohortSpec(
            label="ukbb",
            path=str(cohort_paths["ukbb"]),
            name_column="name",
            key_column="key",
            vocabularies={"UniProtKB": ["UniProtKB"]},
        ),
        CohortSpec(
            label="arivale",
            path=str(cohort_paths["arivale"]),
            name_column="name",
            key_column="key",
            vocabularies={"UniProtKB": ["UniProtKB"]},
        ),
        CohortSpec(
            label="twinsuk",
            path=str(cohort_paths["twinsuk"]),
            name_column=name_column,
            key_column="SeqId",
            vocabularies={},
        ),
    ]

    roots = list(
        dict.fromkeys([output_dir, store_dir, inputs_dir, cohort_paths["twinsuk"].parent])
    )
    started = time.monotonic()
    done_status = None
    attempts_used = 0
    attempts_root = output_dir / ".attempts"
    attempts_root.mkdir(parents=True, exist_ok=True)
    active_subrun_dir: Path | None = None

    async with HarmonizeJobs(
        home=jobs_home,
        roots=roots,
        base_url=base_url,
        kestrel_url=kestrel_url,
        api_key=api_key,
        chunk_size=chunk_size,
        capture_max_entities=capture_max_entities,
        capture_top_n=capture_top_n,
        candidates_default_n=capture_top_n,
        candidates_max_n=max(25, capture_top_n),
        systemic_backoff_s=systemic_backoff_s,
    ) as jobs:
        for attempt in range(max_attempts):
            attempts_used = attempt + 1
            attempt_dir = attempts_root / f"attempt_{attempt}"
            # Distinct attempt_dir guarantees a fresh job_id if a prior attempt finished "partial",
            # so submit_many actually re-runs (reusing store_dir for cached rows).
            submitted = await jobs.submit_many(
                cohorts,
                entity="proteins",
                output_dir=attempt_dir,
                store_dir=store_dir,
                remap_all=(remap_all and attempt == 0),
            )
            while True:
                status = await jobs.wait(submitted.job_id, poll_s=poll_s)
                if status.state in ("complete", "partial", "failed"):
                    done_status = status
                    break
                if status.state in ("stalled", "paused", "expired", "cancelled"):
                    await anyio.sleep(resume_wait_s)
                    await jobs.resume(status.job_id)
                    continue
                done_status = status
                break

            if done_status.state == "complete" and done_status.rows_errored == 0:
                active_subrun_dir = attempt_dir
                break
            if done_status.state == "failed":
                raise RuntimeError(
                    f"Harmonize job {done_status.job_id} failed: {done_status.reason}"
                )
            if attempt + 1 >= max_attempts:
                raise RuntimeError(
                    f"Harmonize job {done_status.job_id} finished {done_status.state} with "
                    f"{done_status.rows_errored} errored rows after {max_attempts} attempts"
                )
            if systemic_backoff_s > 0:
                await anyio.sleep(systemic_backoff_s)

        assert done_status is not None and active_subrun_dir is not None
        pre_hashes = snapshot_report_files(active_subrun_dir)
        report = jobs.load_report(done_status.job_id)
        twinsuk_rows: list[dict[str, Any]] = list(
            ((report.get("arms") or {}).get("names_only") or {})
            .get("mappings", {})
            .get("twinsuk", [])
        )
        row_by_seqid = {str(r["key"]).strip(): r for r in twinsuk_rows}
        targets_by_seqid = load_twinsuk_targets(cohort_paths["twinsuk"])

        reviewed_seqid_refs: list[tuple[str, str]] = []
        weak_twinsuk_names: set[str] = set()
        fetched_extra = 0
        offset: int | None = 0
        weak_review_items: list[Any] = []
        while offset is not None:
            page = jobs.review_items(done_status.job_id, kind="row", offset=offset, limit=50)
            weak_review_items.extend(page.items)
            offset = page.next_offset

        # Stage 1: review job weak_trigger items that touch twinsuk.
        for item in weak_review_items:
            detail = jobs.review_item(done_status.job_id, item.ref)
            assert detail.row is not None
            twinsuk_seqids = [r.key for r in detail.row.rows if r.cohort == "twinsuk"]
            if not twinsuk_seqids:
                continue
            weak_twinsuk_names.add(detail.row.name)
            first_seqid = twinsuk_seqids[0]
            rrow = row_by_seqid.get(first_seqid, {})
            if rrow.get("name"):
                weak_twinsuk_names.add(str(rrow["name"]).strip())

            needs_live_fetch = detail.row.candidates_status not in ("present", "none_searched")
            cand_res = await jobs.fetch_candidates(
                item.ref,
                n=capture_top_n,
                use_captured=True,
                pace_s=pace_s if needs_live_fetch else 0.0,
            )
            if needs_live_fetch:
                fetched_extra += 1
            active_ref = cand_res.ref
            target, full_name = targets_by_seqid.get(first_seqid, ("", detail.row.name))
            score_raw = rrow.get("confidence_score")
            conf_score = float(score_raw) if isinstance(score_raw, int | float) else None
            decision = decide_protein_review(
                target=target,
                target_full_name=full_name,
                weak_trigger=cand_res.weak_trigger,
                engine_choice=cand_res.engine_choice,
                kg_equivalent_ids=rrow.get("kg_equivalent_ids") or {},
                candidates=cand_res.candidates,
                confidence_score=conf_score,
            )
            jobs.record_suggestion(
                active_ref,
                decision.choice,
                decision.rationale,
                update_choices=False,
            )
            for seq_id in twinsuk_seqids:
                reviewed_seqid_refs.append((seq_id, active_ref))

        # Stage 2: review additional resolved suspect rows (collisions, non-HGNC, medium score).
        suspects = find_suspect_twinsuk_rows(
            twinsuk_rows,
            targets_by_seqid,
            already_weak_names=weak_twinsuk_names,
            score_cutoff=score_cutoff,
        )
        suspect_cache: dict[str, tuple[RefRecord, RefCandidates]] = {}
        recorded_refs: set[str] = set()
        for suspect in suspects:
            rrow = row_by_seqid[suspect.seq_id]
            cached = suspect_cache.get(suspect.query_name)
            if cached is None:
                mr = _mapping_result_from_report_row(suspect.query_name, rrow)
                rec = jobs.issue_review_ref(
                    done_status.job_id,
                    name=suspect.query_name,
                    result=mr,
                    weak_trigger=suspect.trigger,
                    row_key=suspect.seq_id,
                )
                cand_res = await jobs.fetch_candidates(
                    rec.ref_id,
                    n=capture_top_n,
                    use_captured=False,
                    pace_s=pace_s if fetched_extra > 0 else 0.0,
                )
                fetched_extra += 1
                suspect_cache[suspect.query_name] = (rec, cand_res)
            else:
                _, cand_res = cached
            active_ref = cand_res.ref
            if active_ref not in recorded_refs:
                decision = decide_protein_review(
                    target=suspect.target,
                    target_full_name=suspect.target_full_name,
                    weak_trigger=suspect.trigger,
                    engine_choice=cand_res.engine_choice,
                    kg_equivalent_ids=rrow.get("kg_equivalent_ids") or {},
                    candidates=cand_res.candidates,
                    confidence_score=suspect.confidence_score,
                )
                jobs.record_suggestion(
                    active_ref,
                    decision.choice,
                    decision.rationale,
                    update_choices=False,
                )
                recorded_refs.add(active_ref)
            reviewed_seqid_refs.append((suspect.seq_id, active_ref))

        choice_rows = jobs.suggestions.cohort_choice_rows(reviewed_seqid_refs)
        choices_path = active_subrun_dir / "choices.tsv"
        write_choices(choices_path, choice_rows, mode=0o600)
        review_counts = verify_choices_alignment(
            choices_path, set(row_by_seqid), pre_hashes, active_subrun_dir
        )
        if subrun_dir.exists():
            subrun_dir.rmdir()
        shutil.move(str(active_subrun_dir), str(subrun_dir))

    wall_clock_s = round(time.monotonic() - started, 3)
    pins_path = jobs_home / "jobs" / done_status.job_id / "pins.json"
    if pins_path.is_file():
        pins_data = json.loads(pins_path.read_text(encoding="utf-8"))
    else:
        pins_data = collect_pins(
            probe=False,
            base_url=base_url or "https://biomapper.expertintheloop.io/api/v1",
            api_key=api_key,
            kestrel_url=kestrel_url or "https://kestrel.expertintheloop.io/api",
        )

    manifest = {
        "schema_version": "dataset_b_run_v1",
        "team": "team-c",
        "mode": "dataset_b",
        "pilot": pilot,
        "slice": slice_n,
        "biomapper_version": resolve_version(),
        "git_sha": _git_sha(worktree_dir),
        "api_health": pins_data.get("api") or {},
        "kestrel_health": pins_data.get("kestrel") or {},
        "upstream_inputs_sha256": upstream_hashes,
        "config": {
            "entity": "proteins",
            "name_column": name_column,
            "cohorts": {
                "ukbb": {"name_column": "name", "key_column": "key", "vocabularies": {"UniProtKB": ["UniProtKB"]}},
                "arivale": {"name_column": "name", "key_column": "key", "vocabularies": {"UniProtKB": ["UniProtKB"]}},
                "twinsuk": {"name_column": name_column, "key_column": "SeqId", "vocabularies": {}},
            },
            "chunk_size": chunk_size,
            "capture_max_entities": capture_max_entities,
            "capture_top_n": capture_top_n,
            "score_cutoff": score_cutoff,
            "store_dir": str(store_dir),
            "seeded_from_seed_store": seeded,
            "remap_all": remap_all,
            "split_complexes": False,
            "annotators": "default",
        },
        "subruns": {
            "dataset_b_proteins": {
                "job_id": done_status.job_id,
                "state": done_status.state,
                "wall_clock_s": wall_clock_s,
                "attempts": attempts_used,
                "rows_errored": done_status.rows_errored,
                "weak_twinsuk_items": len(weak_twinsuk_names),
                "suspect_twinsuk_rows": len(suspects),
                "candidate_extra_fetches": fetched_extra,
                "choices_summary": review_counts,
            }
        },
    }
    (output_dir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_executable_rerun(output_dir / "rerun.sh", worktree_dir)
    return manifest


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--store-dir", type=Path, required=True)
    ap.add_argument("--pilot", action="store_true")
    ap.add_argument("--slice", dest="slice_n", type=int, default=None)
    ap.add_argument("--remap-all", action="store_true")
    ap.add_argument(
        "--name-column",
        choices=("TargetFullName", "Target"),
        default="TargetFullName",
    )
    ap.add_argument("--inputs-dir", type=Path, default=DEFAULT_INPUTS_DIR)
    ap.add_argument("--seed-store-dir", type=Path, default=DEFAULT_SEED_STORE_DIR)
    ap.add_argument("--worktree-dir", type=Path, default=WORKTREE_DIR)
    ap.add_argument("--base-url", type=str, default=None)
    ap.add_argument("--kestrel-url", type=str, default=None)
    ap.add_argument("--api-key", type=str, default=None)
    ap.add_argument("--score-cutoff", type=float, default=3.2)
    ap.add_argument("--pace-s", type=float, default=1.1)
    return ap.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    api_key = args.api_key or os.environ.get("BIOMAPPER_API_KEY")
    manifest = asyncio.run(
        run_dataset_b_team_c(
            output_dir=args.output_dir,
            store_dir=args.store_dir,
            pilot=args.pilot,
            slice_n=args.slice_n,
            remap_all=args.remap_all,
            name_column=args.name_column,
            inputs_dir=args.inputs_dir,
            seed_store_dir=args.seed_store_dir,
            worktree_dir=args.worktree_dir,
            base_url=args.base_url,
            kestrel_url=args.kestrel_url,
            api_key=api_key,
            score_cutoff=args.score_cutoff,
            pace_s=args.pace_s,
        )
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
