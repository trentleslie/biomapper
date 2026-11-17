#!/usr/bin/env python3
"""Dataset B proteome harmonization runner for Team B (engine-only resolution).

Maps ``ukbb``, ``arivale``, and ``twinsuk`` using unedited strings from the
shared competition inputs. Olink cohorts (``ukbb`` and ``arivale``) reuse the
pre-seeded ``shared/seed_store`` cache, while SomaScan (``twinsuk``) runs the
tiered engine-only cascade over exact ``Target`` and ``TargetFullName`` cells.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
from typing import Any, Callable, Sequence

import httpx

import biomapper
from biomapper import default_protein_engine_policy, harmonize_many
from biomapper._provenance import DEFAULT_KESTREL_URL
from biomapper.client import DEFAULT_BASE_URL
from biomapper.cohorts import BatchedApiMapper

WORKTREE_DIR = Path(__file__).resolve().parent.parent
HARMONIZE_ROOT = WORKTREE_DIR.parent.parent
SHARED_DIR = HARMONIZE_ROOT / "shared"
DEFAULT_INPUTS_DIR = SHARED_DIR / "inputs"
DEFAULT_SEED_STORE_DIR = SHARED_DIR / "seed_store"
DEFAULT_RUNS_DIR = HARMONIZE_ROOT / "runs" / "team-b"
BM_PYTHON_BIN = HARMONIZE_ROOT / "bin" / "bm-python"
VENV_DIR = WORKTREE_DIR / ".venv"

COHORT_SPECS: tuple[tuple[str, str, str, str], ...] = (
    ("ukbb", "ukbb_proteins.tsv", "key", "name"),
    ("arivale", "arivale_proteins.tsv", "key", "name"),
    ("twinsuk", "twinsuk_somascan_targets.tsv", "SeqId", "Target"),
)
TWINSUK_EXTRA_NAME_COLUMNS = ["TargetFullName"]
DEFAULT_SUBRUN_NAME = "dataset_b_proteins"


@dataclass(frozen=True)
class RunnerConfig:
    """Resolved execution settings for a Dataset B Team B run."""

    output_dir: Path
    store_dir: Path
    pilot: bool = False
    slice_n: int | None = None
    remap_all: bool = False
    batch_size: int = 10
    candidate_limit: int = 25
    corroborate_columns: str = "always"
    include_vector_fallback: bool = True
    max_chunk_retries: int = 3
    retry_backoff_s: float = 1.5
    max_resume_passes: int = 10
    resume_backoff_s: float = 2.0
    timeout: float = 300.0
    base_url: str | None = None
    api_key: str | None = None
    kestrel_url: str = DEFAULT_KESTREL_URL
    seed_store_dir: Path = DEFAULT_SEED_STORE_DIR
    inputs_dir: Path | None = None


def sha256_file(path: Path) -> str:
    """Return lowercase hexadecimal SHA-256 digest of ``path``."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(65536)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def resolve_inputs_dir(pilot: bool, inputs_dir: Path | None = None) -> Path:
    """Resolve competition input directory for full or pilot mode."""
    base = inputs_dir if inputs_dir is not None else DEFAULT_INPUTS_DIR
    return base / "pilot_20" if pilot else base


def load_tsv_rows(path: Path, *, limit: int | None = None) -> list[dict[str, str]]:
    """Load exact verbatim TSV rows without whitespace stripping."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV input has no header row: {path}")
        rows: list[dict[str, str]] = []
        for row in reader:
            rows.append({str(k): ("" if v is None else str(v)) for k, v in row.items()})
            if limit is not None and len(rows) >= limit:
                break
    return rows


def prepare_store_dir(
    store_dir: Path,
    seed_store_dir: Path,
    *,
    remap_all: bool,
) -> bool:
    """Copy Olink ``seed_store`` into ``store_dir`` if needed and fix permissions.

    Copies only when ``<store_dir>/names_only`` does not exist yet and
    ``remap_all`` is False. Ensures owner read/write permissions (``0o700`` on
    directories, ``0o600`` on files) so subsequent runs can extend the cache.
    """
    store_dir.mkdir(parents=True, exist_ok=True)
    target_names_only = store_dir / "names_only"
    if remap_all or target_names_only.exists():
        return False
    seed_names_only = seed_store_dir / "names_only"
    if not seed_names_only.exists():
        return False
    shutil.copytree(seed_store_dir, store_dir, dirs_exist_ok=True)
    for root, dirs, files in os.walk(target_names_only):
        root_path = Path(root)
        os.chmod(root_path, 0o700)
        for dirname in dirs:
            os.chmod(root_path / dirname, 0o700)
        for filename in files:
            os.chmod(root_path / filename, 0o600)
    return True


def _git_output(worktree_dir: Path, args: Sequence[str]) -> str:
    proc = subprocess.run(
        ["git", "-C", str(worktree_dir), *args],
        capture_output=True,
        text=True,
        check=True,
    )
    return proc.stdout.strip()


def _probe_health_endpoint(url: str, *, timeout: float = 10.0) -> dict[str, Any]:
    try:
        resp = httpx.get(url, timeout=timeout)
        body: Any
        try:
            body = resp.json()
        except Exception:
            body = resp.text
        return {
            "url": url,
            "status_code": resp.status_code,
            "ok": 200 <= resp.status_code < 300,
            "body": body,
        }
    except Exception as exc:
        return {
            "url": url,
            "status_code": None,
            "ok": False,
            "error": str(exc),
        }


def fetch_service_health(
    *,
    base_url: str | None = None,
    kestrel_url: str = DEFAULT_KESTREL_URL,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Query live BioMapper API and Kestrel health endpoints."""
    api_base = (base_url or DEFAULT_BASE_URL).rstrip("/")
    kestrel_base = kestrel_url.rstrip("/")
    return (
        _probe_health_endpoint(f"{api_base}/health"),
        _probe_health_endpoint(f"{kestrel_base}/health"),
    )


def _extract_coverage_stats(report: Any) -> tuple[dict[str, int], dict[str, dict[str, Any]]]:
    """Extract errored counts and compact coverage stats from a ``ManyCohortReport``."""
    names_arm = report.arms["names_only"]
    provenance_by_cohort = report.mapping.get("names_only", {})
    errored_by_cohort: dict[str, int] = {}
    coverage_by_cohort: dict[str, dict[str, Any]] = {}

    for label in report.labels:
        records = names_arm.records[label]
        results = names_arm.results[label]
        prov = provenance_by_cohort.get(label, {})
        n_rows = len(records)
        n_errored = sum(1 for res in results if res.error)
        n_resolved = sum(1 for res in results if res.resolved and not res.error)
        n_unresolved = n_rows - n_resolved - n_errored
        errored_by_cohort[label] = n_errored
        coverage_by_cohort[label] = {
            "n_rows": n_rows,
            "n_resolved": n_resolved,
            "n_unresolved": n_unresolved,
            "n_errored": n_errored,
            "n_reused": int(prov.get("n_reused", 0) or 0),
            "n_fresh": int(prov.get("n_fresh", 0) or 0),
            "resolution_rate": round((n_resolved / n_rows), 6) if n_rows else 0.0,
        }
    return errored_by_cohort, coverage_by_cohort


def write_rerun_script(path: Path) -> Path:
    """Write an executable single-command ``rerun.sh`` at ``path``."""
    script = (
        "#!/usr/bin/env bash\n"
        f"exec {BM_PYTHON_BIN} --venv {VENV_DIR} {Path(__file__).resolve()} \"$@\"\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(script, encoding="utf-8")
    os.chmod(path, 0o755)
    return path


def _build_cohort_inputs(
    inputs_dir: Path,
    *,
    slice_n: int | None,
) -> tuple[dict[str, Path | list[dict[str, str]]], dict[str, str]]:
    cohorts: dict[str, Path | list[dict[str, str]]] = {}
    upstream_sha256: dict[str, str] = {}
    for cohort_name, filename, _, _ in COHORT_SPECS:
        file_path = inputs_dir / filename
        if not file_path.exists():
            raise FileNotFoundError(f"Required input file not found: {file_path}")
        upstream_sha256[filename] = sha256_file(file_path)
        if slice_n is not None:
            cohorts[cohort_name] = load_tsv_rows(file_path, limit=slice_n)
        else:
            cohorts[cohort_name] = file_path
    return cohorts, upstream_sha256


def _default_mapper_factory(config: RunnerConfig) -> Callable[[], BatchedApiMapper]:
    def _factory() -> BatchedApiMapper:
        return BatchedApiMapper(
            base_url=config.base_url,
            api_key=config.api_key,
            batch_size=config.batch_size,
            timeout=config.timeout,
            max_retries=config.max_chunk_retries,
            retry_backoff_s=config.retry_backoff_s,
        )

    return _factory


def run_dataset_b(
    config: RunnerConfig,
    *,
    mapper_factory: Callable[[], Any] | None = None,
    pins: dict[str, Any] | None = None,
    probe_pins: bool = True,
    health_fetcher: Callable[[], tuple[dict[str, Any], dict[str, Any]]] | None = None,
    worktree_dir: Path = WORKTREE_DIR,
) -> dict[str, Any]:
    """Run Dataset B proteome harmonization with automatic resume on errored chunks.

    Returns the written ``run_manifest.json`` payload.
    """
    inputs_dir = resolve_inputs_dir(config.pilot, config.inputs_dir)
    cohorts, upstream_sha256 = _build_cohort_inputs(inputs_dir, slice_n=config.slice_n)

    config.output_dir.mkdir(parents=True, exist_ok=True)
    subrun_dir = config.output_dir / DEFAULT_SUBRUN_NAME

    seeded = prepare_store_dir(
        config.store_dir,
        config.seed_store_dir,
        remap_all=config.remap_all,
    )

    policy = default_protein_engine_policy(
        extra_name_columns=tuple(TWINSUK_EXTRA_NAME_COLUMNS),
        candidate_limit=config.candidate_limit,
        corroborate_columns=config.corroborate_columns,
        include_vector_fallback=config.include_vector_fallback,
    )

    make_mapper = mapper_factory if mapper_factory is not None else _default_mapper_factory(config)
    key_columns = {name: key_col for name, _, key_col, _ in COHORT_SPECS}
    name_columns = {name: name_col for name, _, _, name_col in COHORT_SPECS}

    start_s = time.monotonic()
    errored_by_cohort: dict[str, int] = {}
    coverage_by_cohort: dict[str, dict[str, Any]] = {}
    attempts_used = 0

    for attempt in range(1, config.max_resume_passes + 1):
        attempts_used = attempt
        mapper = make_mapper()
        effective_remap_all = config.remap_all if attempt == 1 else False
        report = harmonize_many(
            cohorts,
            entity="biolink:Protein",
            name_columns=name_columns,
            key_columns=key_columns,
            extra_name_columns={"twinsuk": list(TWINSUK_EXTRA_NAME_COLUMNS)},
            resolution_policies={"twinsuk": policy},
            strip_names={"ukbb": True, "arivale": True, "twinsuk": False},
            batch_size=config.batch_size,
            timeout=config.timeout,
            progress=True,
            mapper=mapper,
            base_url=config.base_url,
            api_key=config.api_key,
            kestrel_url=config.kestrel_url,
            probe_pins=probe_pins,
            pins=pins,
            store_dir=config.store_dir,
            output_dir=subrun_dir,
            save=True,
            remap_all=effective_remap_all,
            max_retries=config.max_chunk_retries,
        )
        errored_by_cohort, coverage_by_cohort = _extract_coverage_stats(report)
        total_errored = sum(errored_by_cohort.values())
        if total_errored == 0:
            break
        if attempt < config.max_resume_passes and config.resume_backoff_s > 0:
            time.sleep(config.resume_backoff_s * attempt)

    wall_clock_s = round(time.monotonic() - start_s, 3)
    total_errored = sum(errored_by_cohort.values())
    if total_errored > 0:
        raise RuntimeError(
            f"Harmonization finished with {total_errored} errored rows after "
            f"{attempts_used} pass(es): {errored_by_cohort}"
        )

    if health_fetcher is not None:
        api_health, kestrel_health = health_fetcher()
    else:
        api_health, kestrel_health = fetch_service_health(
            base_url=config.base_url,
            kestrel_url=config.kestrel_url,
        )

    git_sha = _git_output(worktree_dir, ["rev-parse", "HEAD"])
    branch = _git_output(worktree_dir, ["rev-parse", "--abbrev-ref", "HEAD"])

    manifest: dict[str, Any] = {
        "schema_version": 1,
        "team": "team-b",
        "branch": branch,
        "git_sha": git_sha,
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "biomapper_version": getattr(biomapper, "__version__", "unknown"),
        "api_health": api_health,
        "kestrel_health": kestrel_health,
        "upstream_inputs_sha256": upstream_sha256,
        "derived_inputs": {},
        "config": {
            "strategy": "engine_only_zero_string_editing",
            "pilot": config.pilot,
            "slice": config.slice_n,
            "remap_all": config.remap_all,
            "seed_store_copied": seeded,
            "store_dir": str(config.store_dir),
            "batch_size": config.batch_size,
            "candidate_limit": config.candidate_limit,
            "corroborate_columns": config.corroborate_columns,
            "include_vector_fallback": config.include_vector_fallback,
            "max_chunk_retries": config.max_chunk_retries,
            "max_resume_passes": config.max_resume_passes,
            "attempts_used": attempts_used,
            "twinsuk_strip_names": False,
            "twinsuk_name_columns": ["Target", *TWINSUK_EXTRA_NAME_COLUMNS],
            "tiers": [tier.name for tier in policy.tiers],
        },
        "subruns": {
            DEFAULT_SUBRUN_NAME: {
                "output_dir": DEFAULT_SUBRUN_NAME,
                "wall_clock_s": wall_clock_s,
                "status": "ok",
                "n_errored": total_errored,
                "cohorts": coverage_by_cohort,
            }
        },
    }

    (config.output_dir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    write_rerun_script(config.output_dir / "rerun.sh")
    return manifest


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse CLI arguments for Team B Dataset B harmonization."""
    default_stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    parser = argparse.ArgumentParser(
        description="Run Team B Dataset B proteome harmonization (engine-only, zero string edits).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_RUNS_DIR / default_stamp,
        help="Run output directory receiving run_manifest.json, rerun.sh, and dataset_b_proteins/.",
    )
    parser.add_argument(
        "--store-dir",
        type=Path,
        default=DEFAULT_RUNS_DIR / "store",
        help="Persistent EntityCache store directory.",
    )
    parser.add_argument(
        "--pilot",
        action="store_true",
        help="Use shared/inputs/pilot_20/ instead of full cohorts.",
    )
    parser.add_argument(
        "--slice",
        dest="slice_n",
        type=int,
        default=None,
        help="Optional row limit per cohort for quick smoke runs.",
    )
    parser.add_argument(
        "--remap-all",
        action="store_true",
        help="Force fresh mapping without copying shared/seed_store/ or reusing cache hits.",
    )
    parser.add_argument("--batch-size", type=int, default=10)
    parser.add_argument("--candidate-limit", type=int, default=25)
    parser.add_argument(
        "--corroborate-columns",
        choices=("always", "selective", "fallback"),
        default="always",
    )
    parser.add_argument("--no-vector-fallback", action="store_true")
    parser.add_argument("--max-chunk-retries", type=int, default=3)
    parser.add_argument("--max-resume-passes", type=int, default=10)
    parser.add_argument("--timeout", type=float, default=300.0)
    parser.add_argument("--base-url", type=str, default=None)
    parser.add_argument("--api-key", type=str, default=None)
    parser.add_argument("--kestrel-url", type=str, default=DEFAULT_KESTREL_URL)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entrypoint."""
    args = parse_args(argv)
    if args.slice_n is not None and args.slice_n <= 0:
        raise ValueError("--slice must be positive")
    config = RunnerConfig(
        output_dir=args.output_dir.resolve(),
        store_dir=args.store_dir.resolve(),
        pilot=args.pilot,
        slice_n=args.slice_n,
        remap_all=args.remap_all,
        batch_size=args.batch_size,
        candidate_limit=args.candidate_limit,
        corroborate_columns=args.corroborate_columns,
        include_vector_fallback=not args.no_vector_fallback,
        max_chunk_retries=args.max_chunk_retries,
        max_resume_passes=args.max_resume_passes,
        timeout=args.timeout,
        base_url=args.base_url,
        api_key=args.api_key,
        kestrel_url=args.kestrel_url,
    )
    manifest = run_dataset_b(config)
    print(json.dumps(manifest["subruns"][DEFAULT_SUBRUN_NAME], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
