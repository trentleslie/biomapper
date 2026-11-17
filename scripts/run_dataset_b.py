"""Dataset B proteome harmonization runner for Team A.

Pre-normalizes SomaScan Target and TargetFullName strings using
somascan_normalizer.py, initializes the local mapping store from the frozen
Olink seed store when applicable, and executes plain
harmonize_many(..., entity="proteins") via HarmonizeJobs with chunk_size=10
and capture_max_entities=0. Automatically resumes stalled/paused jobs and
retries any partial runs until zero errored rows remain, then writes
run_manifest.json, rerun.sh, and dataset_b_proteins/seqid_components.tsv.
"""

from __future__ import annotations

import argparse
import asyncio
import csv
import datetime as dt
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

from biomapper._version import resolve_version
from biomapper.cohorts import HarmonizeProgress
from biomapper.jobs import CohortSpec, HarmonizeJobs

SCRIPT_DIR = Path(__file__).resolve().parent
WORKTREE_DIR = SCRIPT_DIR.parent
HARMONIZE_ROOT = Path("/home/trentleslie/agy-team-competition/harmonize")
SHARED_INPUTS_DIR = HARMONIZE_ROOT / "shared" / "inputs"
SEED_STORE_DIR = HARMONIZE_ROOT / "shared" / "seed_store"
WORKTREE_RERUN_SH = WORKTREE_DIR / "rerun.sh"

if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from somascan_normalizer import write_normalized_inputs  # noqa: E402


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_output(worktree_dir: Path, *args: str) -> str:
    try:
        proc = subprocess.run(
            ["git", "-C", str(worktree_dir), *args],
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError:
        return "unknown"
    if proc.returncode != 0:
        return "unknown"
    return proc.stdout.strip() or "unknown"


def _slice_tsv(src: Path, dst: Path, limit: int) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    with src.open(newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise ValueError(f"Input TSV is empty: {src}")
    header, body = rows[0], rows[1 : 1 + limit]
    with dst.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(body)


def _count_report_errors(report: dict[str, Any]) -> int:
    mapping = report.get("mapping") or {}
    total = 0
    for by_cohort in mapping.values():
        if isinstance(by_cohort, dict):
            for prov in by_cohort.values():
                if isinstance(prov, dict):
                    total += int(prov.get("n_errored") or 0)
    return total


def _init_seed_store(
    store_dir: Path,
    *,
    remap_all: bool,
    seed_store_dir: Path | None = SEED_STORE_DIR,
) -> bool:
    """Copy frozen seed_store into store_dir if names_only is missing and not remap_all."""
    if remap_all or (store_dir / "names_only").exists():
        return False
    if seed_store_dir is None or not seed_store_dir.is_dir():
        return False
    store_dir.mkdir(parents=True, exist_ok=True)
    shutil.copytree(seed_store_dir, store_dir, dirs_exist_ok=True)
    return True


def _archive_dir(path: Path, dest_parent: Path, prefix: str) -> None:
    if not path.exists():
        return
    dest_parent.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now(dt.UTC).strftime("%Y%m%dT%H%M%S%fZ")
    target = dest_parent / f"{prefix}_{stamp}"
    shutil.move(str(path), str(target))


def _write_rerun_script(output_dir: Path) -> None:
    dst = output_dir / "rerun.sh"
    if WORKTREE_RERUN_SH.is_file():
        shutil.copyfile(WORKTREE_RERUN_SH, dst)
    else:
        bm_python = HARMONIZE_ROOT / "bin" / "bm-python"
        venv_dir = WORKTREE_DIR / ".venv"
        runner = SCRIPT_DIR / "run_dataset_b.py"
        dst.write_text(
            "#!/usr/bin/env bash\n"
            "set -euo pipefail\n"
            f'exec "{bm_python}" --venv "{venv_dir}" "{runner}" "$@"\n',
            encoding="utf-8",
        )
    os.chmod(dst, 0o755)


async def run_dataset_b(
    *,
    output_dir: Path,
    store_dir: Path,
    pilot: bool = False,
    slice_n: int | None = None,
    remap_all: bool = False,
    base_url: str | None = None,
    kestrel_url: str | None = None,
    api_key: str | None = None,
    systemic_backoff_s: float | None = None,
    poll_s: float = 1.0,
    max_attempts: int = 5,
    max_resumes: int = 25,
    seed_store_dir: Path | None = SEED_STORE_DIR,
    worktree_dir: Path = WORKTREE_DIR,
) -> dict[str, Any]:
    """Execute the Team A Dataset B harmonization run end-to-end."""
    if slice_n is not None and slice_n < 1:
        raise ValueError(f"--slice must be positive, got {slice_n}")

    output_dir = output_dir.resolve()
    store_dir = store_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    store_dir.mkdir(parents=True, exist_ok=True)

    inputs_dir = (SHARED_INPUTS_DIR / "pilot_20" if pilot else SHARED_INPUTS_DIR).resolve()
    ukbb_upstream = inputs_dir / "ukbb_proteins.tsv"
    arivale_upstream = inputs_dir / "arivale_proteins.tsv"
    twinsuk_upstream = inputs_dir / "twinsuk_somascan_targets.tsv"

    upstream_sha256 = {
        "ukbb_proteins.tsv": _sha256_file(ukbb_upstream),
        "arivale_proteins.tsv": _sha256_file(arivale_upstream),
        "twinsuk_somascan_targets.tsv": _sha256_file(twinsuk_upstream),
    }

    work_dir = output_dir / "_work"
    derived_dir = work_dir / "derived"
    job_home = work_dir / "harmonize_jobs_home"
    attempts_archive = work_dir / "attempts_archive"
    work_dir.mkdir(parents=True, exist_ok=True)
    derived_dir.mkdir(parents=True, exist_ok=True)
    job_home.mkdir(parents=True, exist_ok=True)

    if slice_n is not None:
        sliced_dir = work_dir / "sliced_inputs"
        ukbb_input = sliced_dir / "ukbb_proteins.tsv"
        arivale_input = sliced_dir / "arivale_proteins.tsv"
        _slice_tsv(ukbb_upstream, ukbb_input, slice_n)
        _slice_tsv(arivale_upstream, arivale_input, slice_n)
    else:
        ukbb_input = ukbb_upstream
        arivale_input = arivale_upstream

    components_tsv = derived_dir / "seqid_components.tsv"
    twinsuk_normalized_tsv = derived_dir / "twinsuk_normalized.tsv"
    norm_stats = write_normalized_inputs(
        twinsuk_upstream,
        components_tsv=components_tsv,
        mapped_tsv=twinsuk_normalized_tsv,
        limit=slice_n,
    )

    seed_copied = _init_seed_store(
        store_dir, remap_all=remap_all, seed_store_dir=seed_store_dir
    )

    subrun_dir = output_dir / "dataset_b_proteins"
    if subrun_dir.exists():
        _archive_dir(subrun_dir, attempts_archive, "prior_output")

    cohorts = [
        CohortSpec(
            label="ukbb",
            path=str(ukbb_input),
            name_column="name",
            key_column="key",
            vocabularies={"UniProtKB": ["UniProtKB"]},
        ),
        CohortSpec(
            label="arivale",
            path=str(arivale_input),
            name_column="name",
            key_column="key",
            vocabularies={"UniProtKB": ["UniProtKB"]},
        ),
        CohortSpec(
            label="twinsuk",
            path=str(twinsuk_normalized_tsv),
            name_column="component_query",
            key_column="component_key",
            vocabularies={},
        ),
    ]

    roots = sorted(
        {
            inputs_dir,
            output_dir,
            store_dir,
            work_dir,
        }
    )

    last_printed = {"done": -1, "arm": ""}

    def on_progress(_job_id: str, progress: HarmonizeProgress) -> None:
        done = progress.entities_done
        total = progress.entities_total
        arm = progress.arm
        if arm != last_printed["arm"] or done == total or done - last_printed["done"] >= 50:
            last_printed["done"] = done
            last_printed["arm"] = arm
            print(
                f"[team-a] {arm}: {done}/{total} entities mapped",
                file=sys.stderr,
                flush=True,
            )

    limits: dict[str, Any] = {
        "chunk_size": 10,
        "capture_max_entities": 0,
    }
    if systemic_backoff_s is not None:
        limits["systemic_backoff_s"] = systemic_backoff_s

    started_s = time.monotonic()
    attempts = 0
    resumes = 0
    attempt_remap_all = remap_all
    report: dict[str, Any] = {}

    async with HarmonizeJobs(
        home=job_home,
        roots=roots,
        base_url=base_url,
        api_key=api_key,
        kestrel_url=kestrel_url,
        on_progress=on_progress,
        **limits,
    ) as jobs:
        while True:
            attempts += 1
            if attempts > max_attempts:
                raise RuntimeError(
                    f"Exceeded max_attempts={max_attempts} without achieving 0 errored rows"
                )
            if subrun_dir.exists():
                _archive_dir(subrun_dir, attempts_archive, f"attempt_{attempts - 1}")

            status = await jobs.submit_many(
                cohorts,
                entity="proteins",
                output_dir=subrun_dir,
                store_dir=store_dir,
                remap_all=attempt_remap_all,
            )
            while True:
                status = await jobs.wait(status.job_id, poll_s=poll_s)
                if status.state in ("stalled", "paused", "expired", "cancelled"):
                    if resumes >= max_resumes:
                        raise RuntimeError(
                            f"Job {status.job_id} interrupted {resumes} times "
                            f"({status.state}): {status.reason}"
                        )
                    resumes += 1
                    print(
                        f"[team-a] Resuming job {status.job_id} after {status.state}: "
                        f"{status.reason}",
                        file=sys.stderr,
                        flush=True,
                    )
                    status = await jobs.resume(status.job_id)
                    continue
                break

            if status.state == "failed":
                print(
                    f"[team-a] Attempt {attempts} failed ({status.reason}); retrying.",
                    file=sys.stderr,
                    flush=True,
                )
                # Preserve remap_all only if no superseded backup was created yet.
                has_superseded = any(store_dir.rglob("*.superseded-*"))
                if has_superseded:
                    attempt_remap_all = False
                continue

            if status.state not in ("complete", "partial"):
                raise RuntimeError(
                    f"Unexpected terminal job state {status.state!r}: {status.reason}"
                )

            report = jobs.load_report(status.job_id)
            n_errored = _count_report_errors(report)
            if status.state == "complete" and n_errored == 0 and report.get("status") == "ok":
                break

            print(
                f"[team-a] Attempt {attempts} finished {status.state!r} with "
                f"{n_errored} errored row(s); retrying errored rows using cached mappings.",
                file=sys.stderr,
                flush=True,
            )
            attempt_remap_all = False

    elapsed_s = time.monotonic() - started_s
    shutil.copyfile(components_tsv, subrun_dir / "seqid_components.tsv")
    _write_rerun_script(output_dir)

    pins = report.get("pins") or {}
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "team": "team-a",
        "branch": _git_output(worktree_dir, "rev-parse", "--abbrev-ref", "HEAD"),
        "git_sha": _git_output(worktree_dir, "rev-parse", "HEAD"),
        "created_utc": dt.datetime.now(dt.UTC).isoformat(),
        "biomapper_version": report.get("biomapper_version") or resolve_version(),
        "api_health": pins.get("api") or {},
        "kestrel_health": pins.get("kestrel") or {},
        "upstream_inputs_sha256": upstream_sha256,
        "derived_inputs": {
            "seqid_components.tsv": {
                "path": "dataset_b_proteins/seqid_components.tsv",
                "sha256": _sha256_file(subrun_dir / "seqid_components.tsv"),
                **norm_stats,
            },
            "twinsuk_normalized.tsv": {
                "path": "_work/derived/twinsuk_normalized.tsv",
                "sha256": _sha256_file(twinsuk_normalized_tsv),
            },
        },
        "config": {
            "entity": "proteins",
            "chunk_size": 10,
            "capture_max_entities": 0,
            "pilot": pilot,
            "slice_n": slice_n,
            "remap_all": remap_all,
            "seed_store_copied": seed_copied,
            "store_dir": str(store_dir),
        },
        "subruns": {
            "dataset_b_proteins": {
                "output_dir": "dataset_b_proteins",
                "wall_clock_s": round(elapsed_s, 3),
                "status": report.get("status", "ok"),
                "n_errored": 0,
                "attempts": attempts,
                "resumes": resumes,
            }
        },
    }
    manifest_path = output_dir / "run_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--store-dir", type=Path, required=True)
    ap.add_argument("--pilot", action="store_true")
    ap.add_argument("--slice", type=int, default=None, dest="slice_n")
    ap.add_argument("--remap-all", action="store_true")
    ap.add_argument("--base-url", type=str, default=None)
    ap.add_argument("--kestrel-url", type=str, default=None)
    ap.add_argument("--api-key", type=str, default=None)
    ap.add_argument("--systemic-backoff-s", type=float, default=None)
    ap.add_argument("--poll-s", type=float, default=1.0)
    return ap.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    manifest = asyncio.run(
        run_dataset_b(
            output_dir=args.output_dir,
            store_dir=args.store_dir,
            pilot=args.pilot,
            slice_n=args.slice_n,
            remap_all=args.remap_all,
            base_url=args.base_url,
            kestrel_url=args.kestrel_url,
            api_key=args.api_key,
            systemic_backoff_s=args.systemic_backoff_s,
            poll_s=args.poll_s,
        )
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
