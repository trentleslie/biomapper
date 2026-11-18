"""Offline tests for Team C SomaScan candidate matching, runner, choices, and validator."""

from __future__ import annotations

import csv
import importlib.util
import json
import os
import stat
from pathlib import Path
from typing import Any

import httpx
import pytest
import respx

from biomapper.extras.protein_review import (
    decide_review,
    find_suspect_twinsuk_rows,
    is_heterocomplex,
)
from biomapper.mcp.judgments import CHOICES_MARKER
from biomapper.mcp.references import CandidateItem
from tests.mcp_server.mcp_helpers import (
    API,
    KG_BUILD,
    KESTREL,
    batch_response,
    kestrel_block,
    kestrel_row,
)

HARMONIZE_ROOT = Path("/home/trentleslie/agy-team-competition/harmonize")
WORKTREE_DIR = HARMONIZE_ROOT / "worktrees" / "team-c"
SHARED_INPUTS = HARMONIZE_ROOT / "shared" / "inputs"


def _load_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _cand(
    curie: str,
    name: str,
    score: float,
    *,
    synonyms: list[str] | None = None,
    categories: list[str] | None = None,
    relation: str | None = None,
) -> CandidateItem:
    return CandidateItem(
        curie=curie,
        name=name,
        synonyms=synonyms or [],
        categories=categories or ["biolink:Gene", "biolink:Protein"],
        score=score,
        endpoints=["hybrid-search"],
        relation=relation,  # type: ignore[arg-type]
    )


def test_protein_review_overrides_ambiguous_tf_tpo_and_urb() -> None:
    # "TF" means Tissue Factor (F3), NOT Transferrin (TF).
    d_tf = decide_review(
        target="TF",
        target_full_name="Tissue Factor",
        weak_trigger="node_collision",
        engine_choice="HGNC:11740",
        kg_equivalent_ids={"HGNC": ["11740"]},
        candidates=[
            _cand("HGNC:11740", "TF", 5.27, synonyms=["Transferrin", "Serotransferrin"], relation="equals"),
            _cand("HGNC:3541", "F3", 3.55, synonyms=["Tissue factor", "TF", "Coagulation factor III"]),
        ],
        confidence_score=5.27,
    )
    assert d_tf.choice == "HGNC:3541"
    assert d_tf.choice_source == "candidate"

    # "Tpo" means Thrombopoietin (THPO), NOT Thyroid peroxidase (TPO).
    d_tpo = decide_review(
        target="Tpo",
        target_full_name="Thrombopoietin",
        weak_trigger="medium_confidence",
        engine_choice="HGNC:12015",
        kg_equivalent_ids={"HGNC": ["12015"]},
        candidates=[
            _cand("HGNC:12015", "TPO", 4.1, synonyms=["Thyroid peroxidase"], relation="equals"),
            _cand("HGNC:11795", "THPO", 3.2, synonyms=["Thrombopoietin", "MPLLG"]),
        ],
        confidence_score=2.6,
    )
    assert d_tpo.choice == "HGNC:11795"
    assert d_tpo.choice_source == "candidate"

    # "URB" + "Steroid-sensitive gene 1" resolves to CCDC80.
    d_urb = decide_review(
        target="URB",
        target_full_name="Steroid-sensitive gene 1",
        weak_trigger="unresolved",
        engine_choice=None,
        kg_equivalent_ids={},
        candidates=[
            _cand("HGNC:30249", "CCDC80", 3.1, synonyms=["URB", "Steroid-sensitive gene 1", "SSG1"]),
        ],
        confidence_score=None,
    )
    assert d_urb.choice == "HGNC:30249"
    assert d_urb.choice_source == "candidate"


def test_protein_review_rejects_receptor_and_paralog_mismatches() -> None:
    # IL-11 RA must pick IL11RA and reject IL11 (missing receptor) and IL11RB (wrong greek).
    d_rec = decide_review(
        target="IL-11 RA",
        target_full_name="Interleukin-11 receptor subunit alpha",
        weak_trigger="outscored",
        engine_choice="HGNC:5966",
        kg_equivalent_ids={"HGNC": ["5966"]},
        candidates=[
            _cand("HGNC:5966", "IL11", 3.9, synonyms=["Interleukin-11"], relation="equals"),
            _cand("HGNC:5968", "IL11RB", 3.4, synonyms=["Interleukin-11 receptor subunit beta"]),
            _cand("HGNC:5967", "IL11RA", 3.1, synonyms=["Interleukin 11 receptor subunit alpha"]),
        ],
        confidence_score=2.5,
    )
    assert d_rec.choice == "HGNC:5967"
    assert d_rec.choice_source == "candidate"

    # Subtype mismatch when only wrong paralog MMP2 exists for MMP-9 -> "none".
    d_paralog = decide_review(
        target="MMP-9",
        target_full_name="Matrix metalloproteinase-9",
        weak_trigger="low_confidence",
        engine_choice="HGNC:7166",
        kg_equivalent_ids={"HGNC": ["7166"]},
        candidates=[
            _cand("HGNC:7166", "MMP2", 1.9, synonyms=["Matrix metalloproteinase-2"], relation="equals"),
        ],
        confidence_score=1.9,
    )
    assert d_paralog.choice == "none"
    assert d_paralog.choice_source == "none"


def test_protein_review_handles_shorthands_complexes_and_sparse_metadata() -> None:
    # Short code "ON" / "Osteonectin" -> SPARC; "41" / "Protein 4.1" -> EPB41.
    d_on = decide_review(
        target="ON",
        target_full_name="Osteonectin",
        weak_trigger="unresolved",
        engine_choice=None,
        kg_equivalent_ids={},
        candidates=[
            _cand("HGNC:7871", "NONO", 2.9, synonyms=["Non-POU domain-containing octamer-binding protein"]),
            _cand("HGNC:11219", "SPARC", 2.6, synonyms=["Osteonectin", "Secreted protein acidic and rich in cysteine"]),
        ],
    )
    assert d_on.choice == "HGNC:11219"

    d_41 = decide_review(
        target="41",
        target_full_name="Protein 4.1",
        weak_trigger="unresolved",
        engine_choice=None,
        kg_equivalent_ids={},
        candidates=[
            _cand("HGNC:3377", "EPB41", 2.5, synonyms=["Protein 4.1", "Erythrocyte membrane protein band 4.1"]),
        ],
    )
    assert d_41.choice == "HGNC:3377"

    # Heterocomplex "C5b, 6 Complex" where candidates only cover single subunits -> "none".
    assert is_heterocomplex("C5b, 6 Complex", "Complement C5b-C6 complex")
    d_complex = decide_review(
        target="C5b, 6 Complex",
        target_full_name="Complement C5b-C6 complex",
        weak_trigger="heterocomplex",
        engine_choice="HGNC:1337",
        kg_equivalent_ids={"HGNC": ["1337"]},
        candidates=[
            _cand("HGNC:1337", "C6", 2.8, synonyms=["Complement component C6"], relation="equals"),
            _cand("HGNC:1331", "C5", 2.6, synonyms=["Complement component C5"]),
        ],
        confidence_score=2.8,
    )
    assert d_complex.choice == "none"

    # Sparse HGNC engine hit at score 2.6 with no contradicting metadata must be KEPT, never erased.
    d_sparse = decide_review(
        target="TIMP-1",
        target_full_name="Metalloproteinase inhibitor 1",
        weak_trigger="medium_confidence",
        engine_choice="HGNC:11820",
        kg_equivalent_ids={"HGNC": ["11820"]},
        candidates=[_cand("HGNC:11820", "TIMP1", 2.6, synonyms=[], relation="equals")],
        confidence_score=2.6,
    )
    assert d_sparse.choice == "HGNC:11820"
    assert d_sparse.choice_source == "engine"

    # Glutamate carboxypeptidase vs Cytosolic non-specific dipeptidase.
    d_cndp2 = decide_review(
        target="Glutamate carboxypeptidase",
        target_full_name="Cytosolic non-specific dipeptidase",
        weak_trigger="node_collision",
        engine_choice="HGNC:3783",
        kg_equivalent_ids={"HGNC": ["3783"]},
        candidates=[
            _cand("HGNC:3783", "FOLH1", 3.1, synonyms=["Glutamate carboxypeptidase 2"], relation="equals"),
            _cand("HGNC:24124", "CNDP2", 3.0, synonyms=["Cytosolic non-specific dipeptidase"]),
        ],
        confidence_score=3.1,
    )
    assert d_cndp2.choice == "HGNC:24124"

    # Viral antigen unresolved with non-human candidate rejected -> "none".
    d_virus = decide_review(
        target="HIV-1 Gag",
        target_full_name="HIV-1 Gag p24 polyprotein",
        weak_trigger="unresolved",
        engine_choice=None,
        kg_equivalent_ids={},
        candidates=[
            _cand("UniProtKB:P12497", "HIV-1 Gag-Pol", 4.1, synonyms=["Human immunodeficiency virus 1"]),
        ],
    )
    assert d_virus.choice == "none"


def test_find_suspect_twinsuk_rows_flags_collisions_non_hgnc_complex_and_score() -> None:
    rows = [
        {"key": "1-1", "name": "Name A", "chosen_kg_id": "HGNC:1", "confidence_score": 5.1},
        {"key": "1-2", "name": "Name B", "chosen_kg_id": "HGNC:1", "confidence_score": 4.9},
        # True duplicate aptamer probe sharing both Target and TargetFullName -> NOT a collision.
        {"key": "2-1", "name": "Same", "chosen_kg_id": "HGNC:2", "confidence_score": 5.3},
        {"key": "2-2", "name": "Same", "chosen_kg_id": "HGNC:2", "confidence_score": 5.3},
        {"key": "3-1", "name": "Foreign", "chosen_kg_id": "UniProtKB:P99999", "confidence_score": 4.0},
        {"key": "4-1", "name": "Fuzzy", "chosen_kg_id": "HGNC:4", "confidence_score": 2.7},
        {"key": "5-1", "name": "AlreadyWeak", "chosen_kg_id": "HGNC:5", "confidence_score": 2.1},
    ]
    targets = {
        "1-1": ("TA", "Name A"),
        "1-2": ("TB", "Name B"),
        "2-1": ("TS", "Same"),
        "2-2": ("TS", "Same"),
        "3-1": ("TF", "Foreign"),
        "4-1": ("TZ", "Fuzzy"),
        "5-1": ("TW", "AlreadyWeak"),
    }
    suspects = find_suspect_twinsuk_rows(
        rows, targets, already_weak_names={"AlreadyWeak"}, score_cutoff=3.2
    )
    by_key = {s.seq_id: s.trigger for s in suspects}
    assert by_key == {
        "1-1": "node_collision",
        "1-2": "node_collision",
        "3-1": "non_hgnc",
        "4-1": "medium_confidence",
    }


def _build_mock_protein_entry(name: str, *, top_n: int | None = None, idx: int = 1) -> dict[str, Any]:
    """Synthesize clean mocked /map/batch responses for pilot_20 and custom test cohorts."""
    curie = f"HGNC:{1000 + idx}"
    score = 5.2
    cert = "not_applicable"
    kestrel = None
    if name == "Interleukin-11 receptor subunit alpha":
        # Simulate a low-score paralog miss rescued via Kestrel candidates.
        curie = "HGNC:5966"
        score = 1.8
        if top_n:
            kestrel = [
                kestrel_block(
                    "hybrid-search",
                    [
                        kestrel_row("HGNC:5966", 1.8, "IL11", synonyms=["Interleukin-11"]),
                        kestrel_row(
                            "HGNC:5967",
                            3.6,
                            "IL11RA",
                            synonyms=["Interleukin 11 receptor subunit alpha"],
                        ),
                    ],
                )
            ]
    elif name == "Ambiguous_Low_Protein":
        curie = "HGNC:7166"
        score = 1.7
        if top_n:
            kestrel = [
                kestrel_block(
                    "hybrid-search",
                    [kestrel_row("HGNC:7166", 1.7, "MMP2", synonyms=["Matrix metalloproteinase-2"])],
                )
            ]
    elif top_n:
        kestrel = [kestrel_block("hybrid-search", [kestrel_row(curie, score, name, synonyms=[name])])]

    return {
        "name": name,
        "curies": [curie],
        "chosen_kg_id": curie,
        "primary_curie": curie,
        "kg_equivalent_ids": {"HGNC": [curie.split(":")[1]]},
        "assigned_ids": {"kestrel-hybrid-search": {"HGNC": {curie.split(":")[1]: {"score": score}}}},
        "resolution_certificate": {"state": cert},
        "kestrel_results": kestrel,
        "error": None,
    }


@pytest.mark.asyncio
async def test_run_team_c_end_to_end_pilot_partial_retry_and_validate_run(tmp_path: Path) -> None:
    runner_mod = _load_module("run_team_c_mod", WORKTREE_DIR / "scripts" / "run_team_c.py")
    val_mod = _load_module("validate_run_mod", HARMONIZE_ROOT / "shared" / "validate_run.py")

    output_dir = tmp_path / "run_out"
    store_dir = tmp_path / "store"
    # Verify seed_store_if_needed skips on remap_all=True and copies cleanly on remap_all=False.
    seed_src = tmp_path / "seed_src"
    (seed_src / "names_only" / "biolink_Protein").mkdir(parents=True)
    (seed_src / "names_only" / "biolink_Protein" / "entity_cache.jsonl").write_text("{}\n", encoding="utf-8")
    assert runner_mod.seed_store_if_needed(tmp_path / "store_remap", seed_src, remap_all=True) is False
    assert not (tmp_path / "store_remap" / "names_only").exists()
    assert runner_mod.seed_store_if_needed(tmp_path / "store_seed", seed_src, remap_all=False) is True
    assert (tmp_path / "store_seed" / "names_only" / "biolink_Protein" / "entity_cache.jsonl").is_file()

    # Load pilot_20 twinsuk SeqIds for exact verification.
    pilot_soma = SHARED_INPUTS / "pilot_20" / "twinsuk_somascan_targets.tsv"
    expected_seqids = {
        r["SeqId"].strip()
        for r in csv.DictReader(pilot_soma.open(encoding="utf-8", newline=""), delimiter="\t")
    }

    # Fail the first 3 non-candidate calls for the first twinsuk protein (14-3-3 protein sigma)
    # so attempt 0 exhausts max_row_attempts=3 and ends partial, then attempt 1 retries only it.
    pilot_soma_rows = list(csv.DictReader(pilot_soma.open(encoding="utf-8", newline=""), delimiter="\t"))
    first_twinsuk_name = pilot_soma_rows[0]["TargetFullName"].strip()
    transient_fail_remaining = 3
    seen_ids: dict[str, int] = {}

    with respx.mock(assert_all_called=False) as router:
        router.get(f"{API}/health").mock(
            return_value=httpx.Response(200, json={"status": "healthy", "version": "0.1.0"})
        )
        router.get(f"{KESTREL}/health").mock(
            return_value=httpx.Response(200, json={"kestrel_version": "0.5.0", "kg_build": KG_BUILD})
        )

        async def mock_batch(request: httpx.Request) -> httpx.Response:
            nonlocal transient_fail_remaining
            body = json.loads(request.content)
            entities = body["entities"]
            # Simulate persistent row failure across attempt 0 so partial retry exercises store reuse.
            out_rows = []
            for ent in entities:
                name = ent["name"]
                top_n = (ent.get("options") or {}).get("kestrel_top_n")
                if name == first_twinsuk_name and not top_n and transient_fail_remaining > 0:
                    transient_fail_remaining -= 1
                    out_rows.append({"name": name, "curies": [], "error": "transient HTTP 503"})
                    continue
                idx = seen_ids.setdefault(name, len(seen_ids) + 1)
                out_rows.append(_build_mock_protein_entry(name, top_n=top_n, idx=idx))
            return httpx.Response(200, json=batch_response(out_rows))

        router.post(f"{API}/map/batch").mock(side_effect=mock_batch)

        manifest = await runner_mod.run_dataset_b_team_c(
            output_dir=output_dir,
            store_dir=store_dir,
            pilot=True,
            slice_n=None,
            remap_all=False,
            inputs_dir=SHARED_INPUTS,
            seed_store_dir=tmp_path / "no_seed_for_mock_kg",
            worktree_dir=WORKTREE_DIR,
            base_url=API,
            kestrel_url=KESTREL,
            pace_s=0.0,
            poll_s=0.02,
            resume_wait_s=0.01,
            systemic_backoff_s=0.0,
            chunk_size=10,
            capture_max_entities=500,
            capture_top_n=25,
        )

    assert manifest["subruns"]["dataset_b_proteins"]["attempts"] == 2
    assert manifest["subruns"]["dataset_b_proteins"]["rows_errored"] == 0
    subrun = output_dir / "dataset_b_proteins"
    choices_path = subrun / "choices.tsv"
    assert choices_path.is_file()
    lines = choices_path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == CHOICES_MARKER
    assert (output_dir / "rerun.sh").is_file()
    assert os.access(output_dir / "rerun.sh", os.X_OK)
    assert stat.S_IMODE((WORKTREE_DIR / "rerun.sh").stat().st_mode) & stat.S_IXUSR

    # Ensure every key in choices.tsv (if any) belongs to twinsuk SeqIds in report.json.
    for row in csv.DictReader([ln for ln in lines if not ln.startswith("#")], delimiter="\t"):
        assert row["key"] in expected_seqids

    val_result = val_mod.validate_run_dir(output_dir, mode="dataset_b", pilot=True, worktree_dir=None)
    assert val_result["pass"] is True, val_result


@pytest.mark.asyncio
async def test_run_team_c_writes_choices_for_duplicate_targets_and_keeps_report_intact(
    tmp_path: Path,
) -> None:
    runner_mod = _load_module("run_team_c_dup_mod", WORKTREE_DIR / "scripts" / "run_team_c.py")

    inputs = tmp_path / "inputs"
    inputs.mkdir()
    for fname in ("ukbb_proteins.tsv", "arivale_proteins.tsv"):
        (inputs / fname).write_text("key\tname\tUniProtKB\nu1\tSPARC\tP09486\n", encoding="utf-8")

    # Duplicate twinsuk TargetFullName entries plus weak IL-11 RA + unmatchable low-score protein.
    (inputs / "twinsuk_somascan_targets.tsv").write_text(
        "SeqId\tTarget\tTargetFullName\n"
        "4931-59\tIL-11 RA\tInterleukin-11 receptor subunit alpha\n"
        "4931-60\tIL-11 RA\tInterleukin-11 receptor subunit alpha\n"
        "7777-1\tMMP-X\tAmbiguous_Low_Protein\n"
        "8888-1\tON\tSPARC\n",
        encoding="utf-8",
    )

    output_dir = tmp_path / "dup_run"
    store_dir = tmp_path / "dup_store"
    seen_ids: dict[str, int] = {}

    with respx.mock(assert_all_called=False) as router:
        router.get(f"{API}/health").mock(
            return_value=httpx.Response(200, json={"status": "healthy", "version": "0.1.0"})
        )
        router.get(f"{KESTREL}/health").mock(
            return_value=httpx.Response(200, json={"kestrel_version": "0.5.0", "kg_build": KG_BUILD})
        )

        async def mock_batch(request: httpx.Request) -> httpx.Response:
            body = json.loads(request.content)
            out_rows = []
            for ent in body["entities"]:
                name = ent["name"]
                top_n = (ent.get("options") or {}).get("kestrel_top_n")
                idx = seen_ids.setdefault(name, len(seen_ids) + 1)
                out_rows.append(_build_mock_protein_entry(name, top_n=top_n, idx=idx))
            return httpx.Response(200, json=batch_response(out_rows))

        router.post(f"{API}/map/batch").mock(side_effect=mock_batch)

        await runner_mod.run_dataset_b_team_c(
            output_dir=output_dir,
            store_dir=store_dir,
            pilot=False,
            inputs_dir=inputs,
            seed_store_dir=tmp_path / "none",
            worktree_dir=WORKTREE_DIR,
            base_url=API,
            kestrel_url=KESTREL,
            pace_s=0.0,
            poll_s=0.02,
            resume_wait_s=0.01,
            systemic_backoff_s=0.0,
            chunk_size=10,
            capture_max_entities=500,
            capture_top_n=25,
        )

    subrun = output_dir / "dataset_b_proteins"
    lines = (subrun / "choices.tsv").read_text(encoding="utf-8").splitlines()
    parsed = {
        r["key"]: r
        for r in csv.DictReader([ln for ln in lines if not ln.startswith("#")], delimiter="\t")
    }
    assert parsed["4931-59"]["suggested_choice"] == "HGNC:5967"
    assert parsed["4931-60"]["suggested_choice"] == "HGNC:5967"
    assert parsed["7777-1"]["suggested_choice"] == "none"
    # Engine report remains untouched with original weak resolution (HGNC:5966).
    report = json.loads((subrun / "report.json").read_text(encoding="utf-8"))
    tw_rows = {r["key"]: r for r in report["arms"]["names_only"]["mappings"]["twinsuk"]}
    assert tw_rows["4931-59"]["chosen_kg_id"] == "HGNC:5966"
    assert tw_rows["4931-60"]["chosen_kg_id"] == "HGNC:5966"
    assert set(parsed).issubset(set(tw_rows))


def test_prepare_input_files_slice_and_hygiene_patterns(tmp_path: Path) -> None:
    runner_mod = _load_module("run_team_c_slice_mod", WORKTREE_DIR / "scripts" / "run_team_c.py")
    val_mod = _load_module("validate_run_hyg_mod", HARMONIZE_ROOT / "shared" / "validate_run.py")

    paths, hashes = runner_mod.prepare_input_files(
        inputs_dir=SHARED_INPUTS,
        pilot=True,
        slice_n=4,
        staging_dir=tmp_path,
    )
    for label, p in paths.items():
        with p.open(encoding="utf-8", newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        assert len(rows) == 4, label
    assert set(hashes) == {"ukbb_proteins.tsv", "arivale_proteins.tsv", "twinsuk_somascan_targets.tsv"}

    # Verify zero forbidden tokens across Team C touched files and zero em dashes in new files.
    new_files = [
        WORKTREE_DIR / "src" / "biomapper" / "extras" / "protein_review.py",
        WORKTREE_DIR / "scripts" / "run_team_c.py",
        WORKTREE_DIR / "rerun.sh",
        WORKTREE_DIR / "tests" / "mcp_server" / "test_protein_weak_review.py",
        WORKTREE_DIR / "tests" / "test_team_c_protein_review.py",
    ]
    touched_files = [
        *new_files,
        WORKTREE_DIR / "src" / "biomapper" / "jobs.py",
        WORKTREE_DIR / "src" / "biomapper" / "mcp" / "tools.py",
        WORKTREE_DIR / "src" / "biomapper" / "mcp" / "candidates.py",
        WORKTREE_DIR / "src" / "biomapper" / "mcp" / "judgments.py",
        WORKTREE_DIR / "src" / "biomapper" / "mcp" / "jobs.py",
        WORKTREE_DIR / "src" / "biomapper" / "mcp" / "review_entries.py",
    ]
    for path in new_files:
        assert "\u2014" not in path.read_text(encoding="utf-8"), f"em dash found in {path}"
    for path in touched_files:
        if path.parts[-2] == "tests" or "tests" in path.parts:
            continue
        text = path.read_text(encoding="utf-8")
        for tok in val_mod.FORBIDDEN_CODE_PATTERNS:
            assert tok not in text, f"forbidden token {tok!r} in {path}"
