"""Offline unit tests for Team B Dataset B runner, seed store logic, and validator contract."""

from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
from typing import Any

from biomapper.cohorts import BatchedApiMapper
from biomapper.models import MappingResult

WORKTREE_DIR = Path(__file__).resolve().parent.parent
HARMONIZE_ROOT = WORKTREE_DIR.parent.parent
SHARED_DIR = HARMONIZE_ROOT / "shared"
RUNNER_PATH = WORKTREE_DIR / "scripts" / "run_team_b_dataset_b.py"
RERUN_SH_PATH = WORKTREE_DIR / "rerun.sh"
VALIDATE_RUN_PATH = SHARED_DIR / "validate_run.py"


def _load_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


runner = _load_module("run_team_b_dataset_b", RUNNER_PATH)
validator = _load_module("shared_validate_run", VALIDATE_RUN_PATH)

SEED_PIN = {
    "api": {
        "base_url": "https://biomapper.expertintheloop.io/api/v1",
        "status": "healthy",
        "error": None,
    },
    "kestrel": {
        "kestrel_version": "0.5.0",
        "kg_build": {
            "kg_version": "2.3.0",
            "git_commit": "3dd08a5b7b17e320ceb0c0091e09b781d4d90242",
        },
        "error": None,
    },
}


def _ok_result(query: str, idx: int) -> MappingResult:
    curie = f"HGNC:{10000 + idx}"
    return MappingResult(
        query_name=query,
        resolved=True,
        primary_curie=curie,
        chosen_kg_id=curie,
        confidence_score=5.1,
        max_annotator_score=5.1,
        kg_equivalent_ids={"HGNC": [curie]},
        curies_by_annotator={"kestrel-hybrid-search": [curie]},
    )


def _offline_health() -> tuple[dict[str, Any], dict[str, Any]]:
    return (
        {"url": "offline://api/health", "status_code": 200, "ok": True, "body": {"status": "ok"}},
        {"url": "offline://kestrel/health", "status_code": 200, "ok": True, "body": {"status": "ok"}},
    )


def test_prepare_store_dir_seed_copy_skip_and_remap_all(tmp_path: Path) -> None:
    seed_dir = tmp_path / "seed_store"
    seed_cache = seed_dir / "names_only" / "biolink_Protein" / "entity_cache.jsonl"
    seed_cache.parent.mkdir(parents=True, exist_ok=True)
    seed_cache.write_text('{"pin": {}}\n', encoding="utf-8")

    # 1) Fresh store without names_only and remap_all=False copies seed_store.
    store_fresh = tmp_path / "store_fresh"
    copied = runner.prepare_store_dir(store_fresh, seed_dir, remap_all=False)
    assert copied is True
    copied_file = store_fresh / "names_only" / "biolink_Protein" / "entity_cache.jsonl"
    assert copied_file.exists()
    assert stat.S_IMODE(copied_file.parent.stat().st_mode) == 0o700
    assert stat.S_IMODE(copied_file.stat().st_mode) == 0o600

    # 2) Existing store with names_only does NOT overwrite existing entries.
    copied_file.write_text('{"pin": {}, "marker": "preserved"}\n', encoding="utf-8")
    copied_again = runner.prepare_store_dir(store_fresh, seed_dir, remap_all=False)
    assert copied_again is False
    assert "preserved" in copied_file.read_text(encoding="utf-8")

    # 3) remap_all=True skips copying seed_store.
    store_remap = tmp_path / "store_remap"
    copied_remap = runner.prepare_store_dir(store_remap, seed_dir, remap_all=True)
    assert copied_remap is False
    assert not (store_remap / "names_only").exists()


def test_parse_args_supports_required_flags(tmp_path: Path) -> None:
    out_dir = tmp_path / "out"
    store_dir = tmp_path / "store"
    args = runner.parse_args(
        [
            "--output-dir",
            str(out_dir),
            "--store-dir",
            str(store_dir),
            "--pilot",
            "--slice",
            "5",
            "--remap-all",
        ]
    )
    assert args.output_dir == out_dir
    assert args.store_dir == store_dir
    assert args.pilot is True
    assert args.slice_n == 5
    assert args.remap_all is True


def test_run_dataset_b_pilot_resumes_transient_errors_and_passes_validator(tmp_path: Path) -> None:
    output_dir = tmp_path / "run_out"
    store_dir = tmp_path / "store"

    # Inject a deterministic transient failure on pass 1 for the first TwinsUK Target,
    # then succeed on pass 2 so the automatic resume loop completes cleanly.
    pilot_twinsuk = runner.load_tsv_rows(SHARED_DIR / "inputs" / "pilot_20" / "twinsuk_somascan_targets.tsv")
    first_target = pilot_twinsuk[0]["Target"]
    state: dict[str, Any] = {"seen_pass1_fail": False, "id_map": {}, "calls": []}

    class OfflineResumingMapper(BatchedApiMapper):
        def __init__(self) -> None:
            super().__init__(batch_size=10)

        def __call__(self, records: list[dict[str, Any]], **kwargs: Any) -> list[MappingResult]:
            names = [str(rec["name"]) for rec in records]
            state["calls"].append(names)
            out: list[MappingResult] = []
            for name in names:
                if name == first_target and not state["seen_pass1_fail"]:
                    state["seen_pass1_fail"] = True
                    out.append(MappingResult(query_name=name, error="simulated 503"))
                    continue
                idx = state["id_map"].setdefault(name, len(state["id_map"]) + 1)
                out.append(_ok_result(name, idx))
            return out

    config = runner.RunnerConfig(
        output_dir=output_dir,
        store_dir=store_dir,
        pilot=True,
        remap_all=False,
        batch_size=10,
        max_chunk_retries=0,
        max_resume_passes=5,
        resume_backoff_s=0.0,
        retry_backoff_s=0.0,
    )

    manifest = runner.run_dataset_b(
        config,
        mapper_factory=OfflineResumingMapper,
        pins=SEED_PIN,
        probe_pins=False,
        health_fetcher=_offline_health,
        worktree_dir=WORKTREE_DIR,
    )

    assert manifest["team"] == "team-b"
    assert manifest["config"]["attempts_used"] == 2
    assert manifest["config"]["seed_store_copied"] is True
    subrun = manifest["subruns"]["dataset_b_proteins"]
    assert subrun["status"] == "ok"
    assert subrun["n_errored"] == 0

    cohorts = subrun["cohorts"]
    assert cohorts["ukbb"]["n_rows"] == 20
    assert cohorts["ukbb"]["n_reused"] == 20
    assert cohorts["ukbb"]["n_fresh"] == 0
    assert cohorts["arivale"]["n_rows"] == 20
    assert cohorts["arivale"]["n_reused"] == 20
    assert cohorts["arivale"]["n_fresh"] == 0
    assert cohorts["twinsuk"]["n_rows"] == 20
    assert cohorts["twinsuk"]["n_resolved"] == 20
    assert cohorts["twinsuk"]["n_errored"] == 0

    rerun_out = output_dir / "rerun.sh"
    assert rerun_out.exists()
    assert os.access(rerun_out, os.X_OK)

    report = validator.validate_run_dir(
        output_dir,
        pilot=True,
        worktree_dir=None,
    )
    assert report["pass"] is True, json.dumps(report, indent=2)


def test_run_dataset_b_slice_and_remap_all_offline(tmp_path: Path) -> None:
    output_dir = tmp_path / "run_slice"
    store_dir = tmp_path / "store_slice"

    seen_queries: dict[str, int] = {}

    class OfflineSliceMapper(BatchedApiMapper):
        def __init__(self) -> None:
            super().__init__(batch_size=5)

        def __call__(self, records: list[dict[str, Any]], **kwargs: Any) -> list[MappingResult]:
            out: list[MappingResult] = []
            for rec in records:
                name = str(rec["name"])
                idx = seen_queries.setdefault(name, len(seen_queries) + 1)
                out.append(_ok_result(name, idx))
            return out

    config = runner.RunnerConfig(
        output_dir=output_dir,
        store_dir=store_dir,
        pilot=True,
        slice_n=5,
        remap_all=True,
        batch_size=5,
        max_chunk_retries=0,
        max_resume_passes=2,
        resume_backoff_s=0.0,
        retry_backoff_s=0.0,
    )
    manifest = runner.run_dataset_b(
        config,
        mapper_factory=OfflineSliceMapper,
        pins=SEED_PIN,
        probe_pins=False,
        health_fetcher=_offline_health,
        worktree_dir=WORKTREE_DIR,
    )

    assert manifest["config"]["slice"] == 5
    assert manifest["config"]["seed_store_copied"] is False
    for cohort_name in ("ukbb", "arivale", "twinsuk"):
        stats = manifest["subruns"]["dataset_b_proteins"]["cohorts"][cohort_name]
        assert stats["n_rows"] == 5
        assert stats["n_resolved"] == 5
        assert stats["n_reused"] == 0


def test_rerun_scripts_executable_and_clean() -> None:
    assert RERUN_SH_PATH.exists()
    assert os.access(RERUN_SH_PATH, os.X_OK)
    assert os.access(RUNNER_PATH, os.X_OK)

    content = RERUN_SH_PATH.read_text(encoding="utf-8")
    assert "bm-python" in content
    assert "--venv" in content
    assert "run_team_b_dataset_b.py" in content
    assert "\u2014" not in content

    help_proc = subprocess.run(
        [str(RERUN_SH_PATH), "--help"],
        capture_output=True,
        text=True,
        check=True,
    )
    for flag in ("--output-dir", "--store-dir", "--pilot", "--slice", "--remap-all"):
        assert flag in help_proc.stdout

    # Verify zero forbidden tokens (R17) and clean editable .venv (R19); uncommitted
    # worktree files are allowed during pre-commit test runs.
    hygiene = [
        issue
        for issue in validator.audit_worktree_hygiene(WORKTREE_DIR)
        if not issue.startswith("Uncommitted changes in worktree")
    ]
    assert hygiene == [], f"Worktree hygiene failures: {hygiene}"
