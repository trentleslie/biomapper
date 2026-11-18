"""Offline tests for 1.9.0 weak protein review, candidate merging, 429 retry, and choices."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import httpx
import pytest
import respx
from mcp import Client

from biomapper.jobs import HarmonizeJobs
from biomapper.mcp import candidates as cands
from biomapper.mcp import judgments, tools
from biomapper.mcp.references import RefEngine, RefStore
from biomapper.mcp.server import build_server
from biomapper.models import MappingResult
from tests.mcp_server.harmonize_helpers import (
    UNRESOLVED,
    Scripted,
    cfg,
    spec,
    wait_h,
    write_cohort,
)
from tests.mcp_server.mcp_helpers import (
    API,
    batch_response,
    entity_response,
    healthy_routes,
    kestrel_block,
    kestrel_row,
)

PROT = "biolink:Protein"
SM = "biolink:SmallMolecule"


@pytest.fixture()
def router():
    with respx.mock(assert_all_called=False) as r:
        healthy_routes(r)
        yield r


def _mr(name: str = "x", **kw: Any) -> MappingResult:
    return MappingResult.from_api_response(entity_response(name, **kw), query_name=name)


def test_weak_trigger_flags_no_node_outscored_and_low_confidence_for_proteins() -> None:
    no_node = _mr(chosen="INCHIKEY:ABCDEFGHIJKLMNOP", equivalents={"INCHIKEY": ["ABC"]}, cert="not_applicable")
    assert tools.weak_trigger(no_node, PROT) == "no_node"

    empty_chosen = MappingResult(
        query_name="p",
        resolved=True,
        primary_curie="HGNC:1",
        chosen_kg_id=None,
        kg_equivalent_ids={},
    )
    assert tools.weak_trigger(empty_chosen, PROT) == "no_node"

    outscored = MappingResult.from_api_response(
        {
            "result": {
                "name": "p",
                "chosen_kg_id": "HGNC:1",
                "primary_curie": "HGNC:1",
                "kg_equivalent_ids": {"HGNC": ["1"]},
                "assigned_ids": {
                    "kestrel-hybrid-search": {
                        "HGNC": {"1": {"score": 2.4}, "2": {"score": 3.9}}
                    }
                },
                "resolution_certificate": {"state": "not_applicable"},
                "error": None,
            },
            "metadata": {},
        },
        query_name="p",
    )
    assert tools.weak_trigger(outscored, PROT) == "outscored"

    low_prot = MappingResult.from_api_response(
        {
            "result": {
                "name": "p",
                "chosen_kg_id": "HGNC:1",
                "primary_curie": "HGNC:1",
                "kg_equivalent_ids": {"HGNC": ["1"]},
                "assigned_ids": {"kestrel-hybrid-search": {"HGNC": {"1": {"score": 1.65}}}},
                "resolution_certificate": {"state": "not_applicable"},
                "error": None,
            },
            "metadata": {},
        },
        query_name="p",
    )
    assert tools.weak_trigger(low_prot, PROT) == "low_confidence"
    # Small molecules use certificate semantics, not the non-small-molecule < 2.0 threshold.
    assert tools.weak_trigger(low_prot, SM) is None

    ok_prot = MappingResult.from_api_response(
        {
            "result": {
                "name": "p",
                "chosen_kg_id": "HGNC:1",
                "primary_curie": "HGNC:1",
                "kg_equivalent_ids": {"HGNC": ["1"]},
                "assigned_ids": {"kestrel-hybrid-search": {"HGNC": {"1": {"score": 2.4}}}},
                "resolution_certificate": {"state": "not_applicable"},
                "error": None,
            },
            "metadata": {},
        },
        query_name="p",
    )
    assert tools.weak_trigger(ok_prot, PROT) is None


def test_compact_merges_up_to_12_synonyms_and_metadata_across_blocks() -> None:
    assert cands.MAX_SYNONYMS == 12
    syns_a = [f"syn-{i}" for i in range(1, 9)]
    syns_b = [f"syn-{i}" for i in range(7, 16)]
    block_a = kestrel_block(
        "hybrid-search",
        [{"id": "HGNC:3541", "score": 3.1, "name": "", "synonyms": syns_a, "categories": ["biolink:Gene"]}],
    )
    block_b = kestrel_block(
        "text-search",
        [
            {
                "id": "HGNC:3541",
                "score": 4.4,
                "name": "F3",
                "synonyms": syns_b,
                "categories": ["biolink:Protein"],
            }
        ],
    )
    result = _mr(chosen="HGNC:11740", cert="not_applicable", kestrel=[block_a, block_b])
    status, reason, items, _ = cands.compact(result, RefEngine.of(result), n=10)
    assert status == "present" and reason is None
    assert len(items) == 1
    item = items[0]
    assert item.curie == "HGNC:3541"
    assert item.name == "F3"
    assert item.score == 4.4
    assert item.endpoints == ["hybrid-search", "text-search"]
    assert item.categories == ["biolink:Gene", "biolink:Protein"]
    assert len(item.synonyms) == 12
    assert "syn-12" in item.synonyms


def test_judgments_write_choices_and_cohort_choice_rows(tmp_path: Path) -> None:
    refs = RefStore(tmp_path / "home", wall_clock=lambda: 1_700_000_000.0)
    store = judgments.SuggestionStore(refs, wall_clock=lambda: 1_700_000_100.0)
    rec = refs.issue(
        source="job",
        job_id="a" * 16,
        row_key="twinsuk:4931-59 (+1 more rows)",
        query={"name": "Tissue Factor", "entity_type": PROT, "identifiers": {}, "annotation_mode": "missing"},
        result=_mr("Tissue Factor", chosen="HGNC:11740", cert="not_applicable"),
        weak_trigger="medium_confidence",
        kg_build="2.3.0 (3dd08a5b)",
    )
    store.record(rec, "HGNC:11740", "engine", "verified")
    expanded = store.cohort_choice_rows([("4931-59", rec.ref_id), ("4931-60", rec.ref_id)])
    choices_path = tmp_path / "choices.tsv"
    judgments.write_choices(choices_path, expanded, mode=0o600)

    lines = choices_path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == judgments.CHOICES_MARKER
    parsed = list(csv.DictReader([ln for ln in lines if not ln.startswith("#")], delimiter="\t"))
    assert [(r["key"], r["suggested_choice"], r["ref"]) for r in parsed] == [
        ("4931-59", "HGNC:11740", rec.ref_id),
        ("4931-60", "HGNC:11740", rec.ref_id),
    ]
    assert (choices_path.stat().st_mode & 0o777) == 0o600

    # Refuses existing user-owned file and symlinks.
    foreign = tmp_path / "foreign.tsv"
    foreign.write_text("user notes\n", encoding="utf-8")
    with pytest.raises(judgments.JudgmentError, match="left alone"):
        judgments.write_choices(foreign, expanded)
    link_path = tmp_path / "link.tsv"
    link_path.symlink_to(choices_path)
    with pytest.raises(judgments.JudgmentError, match="link"):
        judgments.write_choices(link_path, expanded)


async def test_capture_retries_http_429_outside_gate_without_failing(
    router: respx.MockRouter, tmp_path: Path
) -> None:
    data = tmp_path / "data"
    data.mkdir()
    c1 = write_cohort(data / "c1.tsv", ["glucose", "weakprot"])
    c2 = write_cohort(data / "c2.tsv", ["glucose", "weakprot"])

    table = {
        "glucose": {"chosen": "HGNC:1", "cert": "not_applicable", "equivalents": {"HGNC": ["1"]}},
        "weakprot": UNRESOLVED,
    }
    base = Scripted(table)
    cand_429_left = 2

    async def gated_batch(request: httpx.Request) -> httpx.Response:
        nonlocal cand_429_left
        body = json.loads(request.content)
        is_cand = any((e.get("options") or {}).get("kestrel_top_n") for e in body["entities"])
        if is_cand and cand_429_left > 0:
            cand_429_left -= 1
            return httpx.Response(429, headers={"Retry-After": "0"})
        return await base.batch(request)

    router.post(f"{API}/map/batch").mock(side_effect=gated_batch)
    config = cfg(tmp_path, data, capture_max_entities=10, candidate_budget=3)
    async with Client(build_server(config)) as client:
        res = await client.call_tool(
            "harmonize_many",
            {"cohorts": [spec(c1), spec(c2)], "entity": "proteins"},
        )
        done = await wait_h(client, res.structured_content["job_id"], "complete")
        assert done["state"] == "complete"
        assert done["review_items"]["candidates_captured"] == 1
        assert cand_429_left == 0


async def test_harmonize_jobs_public_review_candidate_and_suggestion_apis(
    router: respx.MockRouter, tmp_path: Path
) -> None:
    data = tmp_path / "data"
    data.mkdir()
    c1 = write_cohort(data / "c1.tsv", ["StrongP", "AmbigTF"])
    c2 = write_cohort(data / "c2.tsv", ["StrongP", "AmbigTF"])

    async def batch_handler(request: httpx.Request) -> httpx.Response:
        body = json.loads(request.content)
        entries = []
        for ent in body["entities"]:
            name = ent["name"]
            top_n = (ent.get("options") or {}).get("kestrel_top_n")
            if name == "StrongP":
                entries.append(
                    {
                        "name": name,
                        "curies": ["HGNC:100"],
                        "chosen_kg_id": "HGNC:100",
                        "primary_curie": "HGNC:100",
                        "kg_equivalent_ids": {"HGNC": ["100"]},
                        "assigned_ids": {"kestrel-hybrid-search": {"HGNC": {"100": {"score": 5.1}}}},
                        "resolution_certificate": {"state": "not_applicable"},
                        "error": None,
                    }
                )
            else:
                kestrel = None
                if top_n:
                    kestrel = [
                        kestrel_block(
                            "hybrid-search",
                            [
                                kestrel_row("HGNC:11740", 5.2, "TF", synonyms=["Transferrin"]),
                                kestrel_row("HGNC:3541", 3.4, "F3", synonyms=["Tissue Factor", "TF"]),
                            ],
                        )
                    ]
                entries.append(
                    {
                        "name": name,
                        "curies": ["HGNC:11740"],
                        "chosen_kg_id": "HGNC:11740",
                        "primary_curie": "HGNC:11740",
                        "kg_equivalent_ids": {"HGNC": ["11740"]},
                        "assigned_ids": {"kestrel-hybrid-search": {"HGNC": {"11740": {"score": 1.8}}}},
                        "resolution_certificate": {"state": "not_applicable"},
                        "kestrel_results": kestrel,
                        "error": None,
                    }
                )
        return httpx.Response(200, json=batch_response(entries))

    router.post(f"{API}/map/batch").mock(side_effect=batch_handler)

    out_dir = tmp_path / "out"
    out_dir.mkdir()
    store_dir = tmp_path / "store"
    store_dir.mkdir()

    async with HarmonizeJobs(
        home=tmp_path / "home",
        roots=[data, out_dir, store_dir],
        base_url=API,
        kestrel_url="https://kestrel.test/api",
        chunk_size=10,
        capture_max_entities=500,
        capture_top_n=25,
        candidates_default_n=25,
        candidates_max_n=25,
        systemic_backoff_s=0.0,
    ) as jobs:
        job = await jobs.submit_many(
            [spec(c1), spec(c2)],
            entity="proteins",
            output_dir=out_dir,
            store_dir=store_dir,
        )
        done = await jobs.wait(job.job_id, poll_s=0.02)
        assert done.state == "complete"
        page = jobs.review_items(done.job_id, kind="row")
        assert page.total == 1
        ref = page.items[0].ref
        cands_out = await jobs.fetch_candidates(ref, n=25, use_captured=True)
        assert cands_out.status == "present"
        assert [c.curie for c in cands_out.candidates] == ["HGNC:11740", "HGNC:3541"]

        sugg = jobs.record_suggestion(ref, "HGNC:3541", "Tissue Factor matches F3")
        assert sugg.choice == "HGNC:3541" and sugg.choice_source == "candidate"
        extra = jobs.issue_review_ref(
            done.job_id,
            name="StrongP",
            result=MappingResult(query_name="StrongP", resolved=True, chosen_kg_id="HGNC:100"),
            weak_trigger="medium_confidence",
            row_key="100-1",
        )
        assert jobs.load_ref(extra.ref_id).row_key == "100-1"
