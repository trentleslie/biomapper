"""Offline unit and integration tests for Team A SomaScan normalization and runner."""

from __future__ import annotations

import csv
import json
import os
import re
import stat
import subprocess
import sys
from pathlib import Path
from typing import Any

import httpx
import pytest
import respx

HARMONIZE_ROOT = Path("/home/trentleslie/agy-team-competition/harmonize")
WORKTREE_DIR = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = WORKTREE_DIR / "scripts"
SHARED_DIR = HARMONIZE_ROOT / "shared"
TWINSUK_TSV = SHARED_DIR / "inputs" / "twinsuk_somascan_targets.tsv"

if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))
if str(SHARED_DIR) not in sys.path:
    sys.path.insert(0, str(SHARED_DIR))

from biomapper._provenance import DEFAULT_KESTREL_URL  # noqa: E402
from biomapper.client import DEFAULT_BASE_URL  # noqa: E402
from run_dataset_b import run_dataset_b  # noqa: E402
from somascan_normalizer import (  # noqa: E402
    normalize_somascan_row,
    normalize_somascan_table,
    write_normalized_inputs,
)
from validate_run import load_seqid_components, validate_run_dir  # noqa: E402

SEED_KG_BUILD = {
    "kg_version": "2.3.0",
    "git_commit": "3dd08a5b7b17e320ceb0c0091e09b781d4d90242",
    "build_timestamp": "2026-10-01",
}


def _fake_hgnc_entry(name: str, idx: int) -> dict[str, Any]:
    curie = f"HGNC:{100000 + idx}"
    return {
        "name": name,
        "curies": [curie],
        "chosen_kg_id": curie,
        "primary_curie": curie,
        "kg_equivalent_ids": {"HGNC": [str(100000 + idx)], "UniProtKB": [f"P{idx:05d}"]},
        "assigned_ids": {
            "kestrel-hybrid-search": {"HGNC": {str(100000 + idx): {"score": 5.0}}}
        },
        "error": None,
    }


def _mount_mock_api(
    router: respx.MockRouter,
    *,
    base_url: str,
    kestrel_url: str,
    kg_build: dict[str, str],
    fail_chunk_calls: int = 0,
    error_once_names: set[str] | None = None,
) -> dict[str, Any]:
    state: dict[str, Any] = {
        "batch_calls": 0,
        "sent_names": [],
        "batch_sizes": [],
        "fail_chunk_calls": fail_chunk_calls,
        "error_once_names": set(error_once_names or ()),
        "errored_once": set(),
    }

    router.get(f"{base_url}/health").mock(
        return_value=httpx.Response(200, json={"status": "healthy", "version": "0.1.0"})
    )
    router.get(f"{kestrel_url}/health").mock(
        return_value=httpx.Response(
            200, json={"kestrel_version": "0.5.0", "kg_build": kg_build}
        )
    )

    name_ids: dict[str, int] = {}

    def handle_batch(request: httpx.Request) -> httpx.Response:
        state["batch_calls"] += 1
        body = json.loads(request.content)
        entities = body.get("entities", [])
        state["batch_sizes"].append(len(entities))

        if state["fail_chunk_calls"] > 0:
            state["fail_chunk_calls"] -= 1
            return httpx.Response(500, json={"detail": "transient server error"})

        results = []
        for idx, ent in enumerate(entities):
            qname = ent["name"]
            state["sent_names"].append(qname)
            if qname in state["error_once_names"] and qname not in state["errored_once"]:
                state["errored_once"].add(qname)
                results.append(
                    {
                        "name": qname,
                        "curies": [],
                        "chosen_kg_id": None,
                        "kg_equivalent_ids": {},
                        "assigned_ids": {},
                        "error": f"transient row timeout for {qname}",
                    }
                )
                continue
            eid = name_ids.setdefault(qname, len(name_ids) + 1 + idx * 0)
            results.append(_fake_hgnc_entry(qname, eid))

        return httpx.Response(
            200,
            json={
                "results": results,
                "metadata": {"request_id": "test-req", "processing_time_ms": 1.0},
            },
        )

    router.post(f"{base_url}/map/batch").mock(side_effect=handle_batch)
    return state


# ---------------------------------------------------------------------------
# Unit tests: complex splitting vs single-protein guards
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("seq_id", "target", "full_name", "expected_queries"),
    [
        (
            "3422-4",
            "CDK1/cyclin B",
            "Cyclin-dependent kinase 1:G2/mitotic-specific cyclin-B1 complex",
            ["CDK1", "CCNB1"],
        ),
        (
            "4482-66",
            "C5b, 6 Complex",
            "Complement C5b-C6 complex",
            ["C5", "C6"],
        ),
        (
            "2765-4",
            "GDF-11/8",
            "Growth/differentiation factor 11/8",
            ["GDF11", "MSTN"],
        ),
        (
            "5183-53",
            "AMPK a1b1g1",
            "AMP Kinase (alpha1beta1gamma1)",
            ["PRKAA1", "PRKAB1", "PRKAG1"],
        ),
        (
            "5245-40",
            "AMPK a2b2g1",
            "AMP Kinase (alpha2beta2gamma1)",
            ["PRKAA2", "PRKAB2", "PRKAG1"],
        ),
        ("2455-17", "IL-12", "Interleukin-12", ["IL12A", "IL12B"]),
        ("4494-63", "IL-23", "Interleukin-23", ["IL23A", "IL12B"]),
        ("2829-19", "IL-27", "Interleukin-27", ["IL27", "EBI3"]),
        (
            "3503-4",
            "Integrin a1b1",
            "Integrin alpha-I: beta-1 complex",
            ["ITGA1", "ITGB1"],
        ),
        (
            "4917-62",
            "Integrin aVb5",
            "Integrin alpha-V: beta-5 complex",
            ["ITGAV", "ITGB5"],
        ),
        (
            "0000-1",
            "Integrin aVb3",
            "Integrin alpha-V: beta-3 complex",
            ["ITGAV", "ITGB3"],
        ),
        ("0000-2", "Integrin aVb3", "Integrin aVb3", ["ITGAV", "ITGB3"]),
        (
            "3714-49",
            "CK-MB",
            "Creatine kinase M-type:Creatine kinase B-type heterodimer",
            ["CKM", "CKB"],
        ),
        (
            "5225-50",
            "CK2-A1:B",
            "Casein kinase II 2-alpha:2-beta heterotetramer",
            ["CSNK2A1", "CSNK2B"],
        ),
        (
            "5226-36",
            "CK2-A2:B",
            "Casein kinase II 2-alpha':2-beta heterotetramer",
            ["CSNK2A2", "CSNK2B"],
        ),
        (
            "3505-6",
            "Lymphotoxin a1/b2",
            "Lymphotoxin alpha1:beta2",
            ["LTA", "LTB"],
        ),
        ("3390-72", "PIK3CA/PIK3R1", "PIK3CA/PIK3R1", ["PIK3CA", "PIK3R1"]),
        (
            "3392-68",
            "PKB a/b/g",
            "RAC-alpha/beta/gamma serine/threonine-protein kinase",
            ["AKT1", "AKT2", "AKT3"],
        ),
        ("3441-64", "GSK-3 alpha/beta", "Glycogen synthase kinase-3 alpha/beta", ["GSK3A", "GSK3B"]),
        ("2986-49", "Gro-b/g", "Gro-beta/gamma", ["CXCL2", "CXCL3"]),
        ("2625-53", "HSP 90a/b", "Heat shock protein HSP 90-alpha/beta", ["HSP90AA1", "HSP90AB1"]),
        ("3032-11", "FSH", "Follicle stimulating hormone", ["CGA", "FSHB"]),
        ("2753-2", "C1q", "Complement C1q subcomponent", ["C1QA", "C1QB", "C1QC"]),
        ("2429-27", "C8", "Complement component C8", ["C8A", "C8B", "C8G"]),
        ("2796-62", "Fibrinogen", "Fibrinogen", ["FGA", "FGB", "FGG"]),
        ("5934-1", "Ferritin", "Ferritin", ["FTH1", "FTL"]),
        ("4903-72", "Calcineurin", "Calcineurin", ["PPP3CA", "PPP3R1"]),
        ("2668-70", "Calpain I", "Calpain I", ["CAPN1", "CAPNS1"]),
        ("4915-64", "Hemoglobin", "Hemoglobin", ["HBA1", "HBB"]),
        ("2728-62", "Laminin", "Laminin", ["LAMA1", "LAMB1", "LAMC1"]),
    ],
)
def test_split_complex_components(
    seq_id: str, target: str, full_name: str, expected_queries: list[str]
) -> None:
    comps = normalize_somascan_row(seq_id, target, full_name)
    assert [c.component_key for c in comps] == [
        f"{seq_id}__{i}" for i in range(1, len(expected_queries) + 1)
    ]
    assert [c.component_query for c in comps] == expected_queries


@pytest.mark.parametrize(
    ("seq_id", "target", "full_name", "expected_query"),
    [
        ("14156-33", "14-3-3 protein beta/alpha", "14-3-3 protein beta/alpha", "YWHAB"),
        ("5858-6", "14-3-3 protein zeta/delta", "14-3-3 protein zeta/delta", "YWHAZ"),
        (
            "10370-21",
            "STAT1",
            "Signal transducer and activator of transcription 1-alpha/beta",
            "STAT1",
        ),
        ("9183-7", "IFN-a/b R1", "Interferon alpha/beta receptor 1", "IFNAR1"),
        (
            "4450-26",
            "hnRNP A/B",
            "Heterogeneous nuclear ribonucleoprotein A/B",
            "HNRNPAB",
        ),
        (
            "5351-52",
            "hnRNP A2/B1",
            "Heterogeneous nuclear ribonucleoproteins A2/B1",
            "HNRNPA2B1",
        ),
        ("2827-23", "Fractalkine/CX3CL-1", "Fractalkine", "CX3CL1"),
        ("4866-59", "TrkB", "BDNF/NT-3 growth factors receptor", "NTRK2"),
        ("13117-232", "CHKB", "Choline/ethanolamine kinase", "CHKB"),
        (
            "3684-78",
            "Macrophage scavenger receptor",
            "Macrophage scavenger receptor types I and II",
            "MSR1",
        ),
        (
            "4153-11",
            "alpha-1-antichymotrypsin complex",
            "Alpha-1-antichymotrypsin complex",
            "SERPINA3",
        ),
        ("4179-57", "14-3-3", "14-3-3 protein family", "14-3-3 protein family"),
        ("11089-7", "IgA", "Immunoglobulin A", "Immunoglobulin A"),
        ("4907-56", "D-dimer", "D-dimer", "D-dimer"),
    ],
)
def test_single_protein_slash_and_family_guards_not_split(
    seq_id: str, target: str, full_name: str, expected_query: str
) -> None:
    comps = normalize_somascan_row(seq_id, target, full_name)
    assert len(comps) == 1
    assert comps[0].component_key == seq_id
    assert comps[0].component_query == expected_query


@pytest.mark.parametrize(
    ("seq_id", "target", "full_name", "expected_query"),
    [
        ("4931-59", "TF", "Tissue Factor", "F3"),
        ("3873-51", "TMA", "Thyroid peroxidase", "TPO"),
        ("4706-17", "41", "Protein 4.1", "EPB41"),
        ("4995-16", "HPG-", "15-hydroxyprostaglandin dehydrogenase [NAD(+)]", "HPGD"),
        ("4330-4", "PSA", "Prostate-specific antigen", "KLK3"),
        ("5478-50", "PSMA", "Glutamate carboxypeptidase 2", "FOLH1"),
        ("3192-3", "Glutamate carboxypeptidase", "Cytosolic non-specific dipeptidase", "CNDP2"),
        ("5843-60", "NG36", "Histone-lysine N-methyltransferase EHMT2", "EHMT2"),
        (
            "4304-18",
            "transcription factor MLR1, isoform CRA_b",
            "Ligand-dependent nuclear receptor corepressor-like protein",
            "LCORL",
        ),
        ("5861-78", "3HAO", "3-hydroxyanthranilate 3,4-dioxygenase", "HAAO"),
        ("9202-309", "ATPO", "ATP synthase subunit O, mitochondrial", "ATP5PO"),
        ("4910-21", "GIB", "Phospholipase A2", "PLA2G1B"),
        ("2449-1", "GV", "Calcium-dependent phospholipase A2", "PLA2G5"),
        ("10356-21", "c-Jun", "Transcription factor AP-1", "JUN"),
        ("2846-24", "Ubiquitin+1", "Ubiquitin+1, truncated mutation for UbB", "UBB"),
        ("5328-33", "EGFRvIII", "Epidermal growth factor receptor variant III", "EGFR"),
        ("5312-49", "Apo E2", "Apolipoprotein E (isoform E2)", "APOE"),
        ("4394-71", "FGF-8A", "Fibroblast growth factor 8 isoform A", "FGF8"),
        ("4498-62", "NCAM-120", "Neural cell adhesion molecule 1, 120 kDa isoform", "NCAM1"),
        ("4867-15", "VEGF121", "Vascular endothelial growth factor A, isoform 121", "VEGFA"),
        (
            "9182-3",
            "sLRP1",
            "Low-density lipoprotein receptor-related protein 1, soluble",
            "LRP1",
        ),
        (
            "3052-8",
            "Fas ligand, soluble",
            "Tumor necrosis factor ligand superfamily member 6, soluble form",
            "FASLG",
        ),
        ("3710-49", "Angiostatin", "Angiostatin", "PLG"),
        ("2201-17", "Endostatin", "Endostatin", "COL18A1"),
        ("4900-8", "C3a", "C3a anaphylatoxin", "C3"),
        ("2755-8", "C3adesArg", "C3a anaphylatoxin des Arginine", "C3"),
        ("5803-24", "C3d", "Complement C3d fragment", "C3"),
        ("2683-1", "iC3b", "Complement C3b, inactivated", "C3"),
        ("4890-10", "ACTH", "Corticotropin", "POMC"),
        ("2558-51", "b-Endorphin", "Beta-endorphin", "POMC"),
        ("4544-4", "CTAP-III", "Connective tissue-activating peptide III", "PPBP"),
        ("2790-54", "NAP-2", "Neutrophil-activating peptide 2", "PPBP"),
        ("3434-34", "FN1.3", "Fibronectin Fragment 3", "FN1"),
        ("7655-11", "N-terminal pro-BNP", "N-terminal pro-BNP", "NPPB"),
        ("3723-1", "BNP-32", "Brain natriuretic peptide 32", "NPPB"),
        (
            "8446-4",
            "PACAP-27",
            "Pituitary adenylate cyclase-activating polypeptide 27",
            "ADCYAP1",
        ),
        ("5957-30", "Somatostatin-28", "Somatostatin-28", "SST"),
        ("4157-2", "Thrombin", "Thrombin", "F2"),
        ("3758-68", "Activated Protein C", "Activated Protein C", "PROC"),
        ("5307-12", "Coagulation Factor IXab", "Coagulation factor IXab", "F9"),
        ("3077-66", "Coagulation Factor Xa", "Coagulation factor Xa", "F10"),
        ("4923-79", "MIS", "Muellerian-inhibiting factor", "AMH"),
        ("3642-4", "SLAF5", "SLAM family member 5", "CD84"),
        ("11094-104", "LPPL", "Galectin-10", "CLC"),
        ("3083-71", "XEDAR", "Tumor necrosis factor receptor superfamily member 27", "EDA2R"),
        ("4314-12", "XTP3A", "dCTP pyrophosphatase 1", "DCTPP1"),
        ("3388-58", "PAK7", "Serine/threonine-protein kinase PAK 7", "PAK5"),
        ("3347-9", "BARK1", "beta-adrenergic receptor kinase 1", "GRK2"),
        ("11096-57", "HEMK2", "HemK methyltransferase family member 2", "N6AMT1"),
        ("9999-1", "NovelTarget", "Novel protein, soluble form", "Novel protein"),
        ("9999-2", "NOVEL1", "Unseen hypothetical factor (isoform B)", "NOVEL1"),
    ],
)
def test_nomenclature_normalization_and_disambiguation(
    seq_id: str, target: str, full_name: str, expected_query: str
) -> None:
    comps = normalize_somascan_row(seq_id, target, full_name)
    assert len(comps) == 1
    assert comps[0].component_query == expected_query


# ---------------------------------------------------------------------------
# Full 1,364-row seqid_components.tsv validation
# ---------------------------------------------------------------------------


def test_full_1364_twinsuk_normalization_and_seqid_components(tmp_path: Path) -> None:
    components_tsv = tmp_path / "seqid_components.tsv"
    mapped_tsv = tmp_path / "twinsuk_normalized.tsv"
    stats = write_normalized_inputs(
        TWINSUK_TSV,
        components_tsv=components_tsv,
        mapped_tsv=mapped_tsv,
    )

    assert stats["n_input_seqids"] == 1364
    assert stats["n_split_seqids"] == 45
    assert stats["n_single_seqids"] == 1319
    assert stats["n_total_components"] > 1364

    with TWINSUK_TSV.open(newline="", encoding="utf-8") as fh:
        expected_seqids = {str(r["SeqId"]).strip() for r in csv.DictReader(fh, delimiter="\t")}
    assert len(expected_seqids) == 1364

    comp_map, comp_errs = load_seqid_components(tmp_path, expected_seqids)
    assert comp_errs == []
    assert set(comp_map) == expected_seqids

    # Ruling R17: never emit CURIEs or database IDs in component_query.
    with components_tsv.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            q = row["component_query"]
            assert q.strip() != ""
            assert not re.search(
                r"\b(?:HGNC|UniProtKB|NCBIGene|PR|ENSEMBL|CHEBI|HMDB):", q, re.I
            )


def test_normalize_somascan_table_limit() -> None:
    comps, stats = normalize_somascan_table(TWINSUK_TSV, limit=20)
    assert stats["n_input_seqids"] == 20
    assert len(comps) >= 20


# ---------------------------------------------------------------------------
# End-to-end runner, retry/resume, seed_store reuse, and validator tests
# ---------------------------------------------------------------------------


def test_runner_pilot_offline_end_to_end(tmp_path: Path) -> None:
    api = "https://api.test/api/v1"
    kestrel = "https://kestrel.test/api"
    out_dir = tmp_path / "run_pilot"
    store_dir = tmp_path / "store_pilot"

    with respx.mock(assert_all_called=False) as router:
        mock_state = _mount_mock_api(
            router, base_url=api, kestrel_url=kestrel, kg_build=SEED_KG_BUILD
        )
        manifest = __import__("asyncio").run(
            run_dataset_b(
                output_dir=out_dir,
                store_dir=store_dir,
                pilot=True,
                remap_all=False,
                base_url=api,
                kestrel_url=kestrel,
                systemic_backoff_s=0.0,
                poll_s=0.01,
                seed_store_dir=None,
            )
        )

    assert manifest["schema_version"] == 1
    assert manifest["config"]["chunk_size"] == 10
    assert manifest["config"]["capture_max_entities"] == 0
    assert manifest["subruns"]["dataset_b_proteins"]["n_errored"] == 0
    assert all(sz <= 10 for sz in mock_state["batch_sizes"])

    rerun_path = out_dir / "rerun.sh"
    assert rerun_path.is_file()
    assert bool(rerun_path.stat().st_mode & stat.S_IXUSR)

    # Confirm strict names-only on twinsuk and UniProtKB declared on ukbb/arivale.
    settings = json.loads((out_dir / "dataset_b_proteins" / "settings.json").read_text())
    assert settings["vocabularies"]["declared"]["twinsuk"] == []
    assert settings["vocabularies"]["declared"]["ukbb"] == ["UNIPROTKB"]
    assert settings["vocabularies"]["declared"]["arivale"] == ["UNIPROTKB"]

    val = validate_run_dir(out_dir, mode="dataset_b", pilot=True, worktree_dir=None)
    assert val["pass"] is True, val["failures"]


def test_runner_full_dataset_reuses_seed_store_and_passes_validator(tmp_path: Path) -> None:
    out_dir = tmp_path / "run_full"
    store_dir = tmp_path / "store_full"

    with respx.mock(assert_all_called=False) as router:
        mock_state = _mount_mock_api(
            router,
            base_url=DEFAULT_BASE_URL,
            kestrel_url=DEFAULT_KESTREL_URL,
            kg_build=SEED_KG_BUILD,
        )
        manifest = __import__("asyncio").run(
            run_dataset_b(
                output_dir=out_dir,
                store_dir=store_dir,
                pilot=False,
                remap_all=False,
                systemic_backoff_s=0.0,
                poll_s=0.01,
            )
        )

    assert manifest["config"]["seed_store_copied"] is True
    assert (store_dir / "names_only").is_dir()

    report = json.loads((out_dir / "dataset_b_proteins" / "report.json").read_text())
    prov = report["mapping"]["names_only"]
    assert prov["ukbb"]["n_fresh"] == 0
    assert prov["ukbb"]["n_reused"] == 2923
    assert prov["arivale"]["n_fresh"] == 0
    assert prov["arivale"]["n_reused"] == 1197
    # Only TwinsUK queries absent from the Olink seed store hit the API, deduplicated once each.
    assert prov["twinsuk"]["n_fresh"] > 0
    assert prov["twinsuk"]["n_reused"] > 0
    assert len(mock_state["sent_names"]) == len(set(mock_state["sent_names"]))
    assert len(mock_state["sent_names"]) <= prov["twinsuk"]["n_fresh"]
    assert all(sz <= 10 for sz in mock_state["batch_sizes"])

    val = validate_run_dir(out_dir, mode="dataset_b", pilot=False, worktree_dir=None)
    assert val["pass"] is True, val["failures"]
    cohorts = val["subruns"]["dataset_b_proteins"]["cohorts"]
    for label in ("ukbb", "arivale", "twinsuk"):
        assert cohorts[label]["errored"] == 0
        assert cohorts[label]["grouped_coverage"] >= 0.90


def test_runner_resumes_stalled_and_retries_partial_until_zero_errors(tmp_path: Path) -> None:
    api = "https://api.test/api/v1"
    kestrel = "https://kestrel.test/api"
    out_dir = tmp_path / "run_resilient"
    store_dir = tmp_path / "store_resilient"

    # Fail 3 chunks in a row (triggering stalled -> resume) and fail one
    # TwinsUK entity on its first attempt (even across internal 3-row retries,
    # if we fail it 4 times total it turns partial on attempt 1 and succeeds
    # on attempt 2 via the saved-mapping store!).
    with respx.mock(assert_all_called=False) as router:
        state = _mount_mock_api(
            router,
            base_url=api,
            kestrel_url=kestrel,
            kg_build=SEED_KG_BUILD,
            fail_chunk_calls=3,
        )

        fail_budget = {"YWHAB": 4}
        orig_handler = router.routes[-1].side_effect

        def wrapped_batch(request: httpx.Request) -> httpx.Response:
            body = json.loads(request.content)
            entities = body.get("entities", [])
            if state["fail_chunk_calls"] > 0:
                return orig_handler(request)
            names = [e["name"] for e in entities]
            if any(fail_budget.get(n, 0) > 0 for n in names):
                results = []
                for n in names:
                    if fail_budget.get(n, 0) > 0:
                        fail_budget[n] -= 1
                        results.append(
                            {
                                "name": n,
                                "curies": [],
                                "chosen_kg_id": None,
                                "kg_equivalent_ids": {},
                                "assigned_ids": {},
                                "error": "forced row failure",
                            }
                        )
                    else:
                        results.append(_fake_hgnc_entry(n, abs(hash(n)) % 90000 + 1000))
                return httpx.Response(
                    200,
                    json={
                        "results": results,
                        "metadata": {"request_id": "r", "processing_time_ms": 1.0},
                    },
                )
            return orig_handler(request)

        router.post(f"{api}/map/batch").mock(side_effect=wrapped_batch)

        manifest = __import__("asyncio").run(
            run_dataset_b(
                output_dir=out_dir,
                store_dir=store_dir,
                pilot=True,
                remap_all=True,
                base_url=api,
                kestrel_url=kestrel,
                systemic_backoff_s=0.0,
                poll_s=0.01,
                seed_store_dir=None,
            )
        )

    subrun_meta = manifest["subruns"]["dataset_b_proteins"]
    assert subrun_meta["resumes"] >= 1
    assert subrun_meta["attempts"] >= 2
    assert subrun_meta["n_errored"] == 0
    assert manifest["config"]["remap_all"] is True
    assert manifest["config"]["seed_store_copied"] is False

    val = validate_run_dir(out_dir, mode="dataset_b", pilot=True, worktree_dir=None)
    assert val["pass"] is True, val["failures"]


def test_rerun_sh_help_invocation() -> None:
    rerun_sh = WORKTREE_DIR / "rerun.sh"
    assert os.access(rerun_sh, os.X_OK)
    proc = subprocess.run(
        [str(rerun_sh), "--help"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0
    assert "--output-dir" in proc.stdout
    assert "--store-dir" in proc.stdout
    assert "--pilot" in proc.stdout
    assert "--slice" in proc.stdout
    assert "--remap-all" in proc.stdout


def test_runner_slice_option(tmp_path: Path) -> None:
    api = "https://api.test/api/v1"
    kestrel = "https://kestrel.test/api"
    out_dir = tmp_path / "run_slice"
    store_dir = tmp_path / "store_slice"

    with respx.mock(assert_all_called=False) as router:
        _mount_mock_api(
            router, base_url=api, kestrel_url=kestrel, kg_build=SEED_KG_BUILD
        )
        manifest = __import__("asyncio").run(
            run_dataset_b(
                output_dir=out_dir,
                store_dir=store_dir,
                pilot=False,
                slice_n=10,
                remap_all=True,
                base_url=api,
                kestrel_url=kestrel,
                systemic_backoff_s=0.0,
                poll_s=0.01,
                seed_store_dir=None,
            )
        )

    assert manifest["config"]["slice_n"] == 10
    assert manifest["subruns"]["dataset_b_proteins"]["n_errored"] == 0
    report = json.loads((out_dir / "dataset_b_proteins" / "report.json").read_text())
    assert report["mapping"]["names_only"]["ukbb"]["n_rows"] == 10
    assert report["mapping"]["names_only"]["arivale"]["n_rows"] == 10
    assert report["mapping"]["names_only"]["twinsuk"]["n_rows"] >= 10
