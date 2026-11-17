"""Offline unit tests for Team B engine-only resolution and client extensions."""

from __future__ import annotations

import asyncio
import builtins
import json
from pathlib import Path
import sys
from typing import Any

import httpx
import pytest
import respx

import biomapper
from biomapper import (
    CohortResolutionPolicy,
    EngineSearchOptions,
    ResolutionTier,
    default_protein_engine_policy,
    harmonize_many,
)
from biomapper._checkpoint import EntityCache, cache_key
from biomapper.client import BioMapperClient
from biomapper.cohorts import BatchedApiMapper
from biomapper.mcp.config import ServerConfig
from biomapper.mcp.harmonize_jobs import CohortSpec, HarmonizeRunner
from biomapper.mcp.jobs import JobManager, RequestGate
from biomapper.models import MappingResult
from biomapper.resolution import PolicyOutcome, call_mapper_with_options

TEST_PIN = {
    "api_endpoint": "https://biomapper.expertintheloop.io/api/v1",
    "kestrel_version": "0.5.0",
    "kg_version": "2.3.0",
    "kg_git_commit": "3dd08a5b7b17e320ceb0c0091e09b781d4d90242",
    "entity_type": "biolink:Protein",
    "annotation_mode": "missing",
}

TEST_RAW_PINS = {
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


def _mapping_result(
    query: str,
    *,
    curie: str | None = None,
    score: float = 5.0,
    annotator: str = "kestrel-hybrid-search",
    equivalent_ids: dict[str, list[str]] | None = None,
    error: str | None = None,
) -> MappingResult:
    if error is not None:
        return MappingResult(query_name=query, error=error)
    if curie is None:
        return MappingResult(query_name=query, resolved=False)
    prefix = curie.split(":", 1)[0]
    eq = dict(equivalent_ids) if equivalent_ids is not None else {prefix: [curie]}
    return MappingResult(
        query_name=query,
        resolved=True,
        primary_curie=curie,
        chosen_kg_id=curie,
        confidence_score=score,
        max_annotator_score=score,
        kg_equivalent_ids=eq,
        curies_by_annotator={annotator: [curie]},
    )


class ScriptedEngineMapper(BatchedApiMapper):
    """Offline deterministic BatchedApiMapper recording exact query strings and options."""

    def __init__(self, responder: Any) -> None:
        super().__init__(batch_size=10)
        self.responder = responder
        self.calls: list[tuple[list[str], dict[str, Any]]] = []

    def __call__(self, records: list[dict[str, Any]], **kwargs: Any) -> list[MappingResult]:
        # Record verbatim strings without mutation to verify zero-string-editing invariant.
        names = [str(rec["name"]) for rec in records]
        self.calls.append((names, dict(kwargs)))
        return [self.responder(name, kwargs) for name in names]


def _write_seed_store(store_dir: Path, entries: dict[str, tuple[str, float]]) -> None:
    cache_dir = store_dir / "names_only" / "biolink_Protein"
    cache = EntityCache(cache_dir, TEST_PIN)
    base_options = {"annotation_mode": "missing", "entity_type": "biolink:Protein"}
    keyed: list[tuple[str, MappingResult]] = []
    for raw_name, (kg_id, score) in entries.items():
        key = cache_key(raw_name, None, {}, base_options)
        keyed.append((key, _mapping_result(raw_name, curie=kg_id, score=score)))
    cache.commit_batch(vocab=None, batch_index=0, keyed_results=keyed)


def _batch_api_entry(name: str, curie: str = "HGNC:6018", score: float = 5.3) -> dict[str, Any]:
    prefix = curie.split(":", 1)[0]
    return {
        "name": name,
        "curies": [curie],
        "primary_curie": curie,
        "chosen_kg_id": curie,
        "kg_equivalent_ids": {prefix: [curie]},
        "curies_assigned": {"kestrel-hybrid-search": [curie]},
        "kestrel_results": [
            {
                "endpoint": "kestrel-hybrid-search",
                "rows": [
                    {
                        "id": curie,
                        "curie": curie,
                        "name": name,
                        "score": score,
                        "category": "biolink:Protein",
                    }
                ],
            }
        ],
    }


@respx.mock
def test_batched_api_mapper_forwards_engine_options_and_retries() -> None:
    captured_payloads: list[dict[str, Any]] = []
    attempt_counter = {"n": 0}

    def _handle_map(request: httpx.Request) -> httpx.Response:
        attempt_counter["n"] += 1
        payload = json.loads(request.content.decode("utf-8"))
        captured_payloads.append(payload)
        if attempt_counter["n"] == 1:
            return httpx.Response(503, json={"detail": "transient unavailable"})
        entries = [_batch_api_entry(item["name"]) for item in payload["entities"]]
        return httpx.Response(200, json={"results": entries, "metadata": {}, "summary": {}})

    respx.post("https://api.example.test/api/v1/map/batch").mock(side_effect=_handle_map)

    mapper = BatchedApiMapper(
        base_url="https://api.example.test/api/v1",
        api_key="test-key",
        batch_size=10,
        max_retries=2,
        retry_backoff_s=0.0,
    )
    results = mapper(
        [{"name": "IL6", "identifiers": {}}],
        entity_type="biolink:Gene",
        annotation_mode="missing",
        prefer_human=True,
        prefer_canonical=True,
        candidate_limit=25,
        annotators=["kestrel-text-search", "kestrel-hybrid-search"],
        kestrel_top_n=30,
        vocab="HGNC",
        array_delimiters=[";"],
    )

    assert len(results) == 1
    assert results[0].chosen_kg_id == "HGNC:6018"
    assert attempt_counter["n"] == 2
    sent = captured_payloads[-1]["entities"][0]
    assert sent["name"] == "IL6"
    assert sent["entity_type"] == "biolink:Gene"
    opts = sent["options"]
    assert opts["prefer_human"] is True
    assert opts["prefer_canonical"] is True
    assert opts["candidate_limit"] == 25
    assert opts["annotators"] == ["kestrel-text-search", "kestrel-hybrid-search"]
    assert opts["kestrel_top_n"] == 30
    assert opts["vocab"] == "HGNC"
    assert opts["array_delimiters"] == [";"]

    # Default max_retries=0 preserves strict error reporting without extra retries.
    attempt_counter["n"] = 0
    strict_mapper = BatchedApiMapper(
        base_url="https://api.example.test/api/v1",
        api_key="test-key",
        batch_size=10,
        max_retries=0,
    )
    err_results = strict_mapper(
        [{"name": "IL6", "identifiers": {}}],
        entity_type="biolink:Protein",
        annotation_mode="missing",
    )
    assert attempt_counter["n"] == 1
    assert err_results[0].error is not None


def test_engine_search_options_scope_cache_keys_and_preserve_seed_store(tmp_path: Path) -> None:
    store_dir = tmp_path / "store"
    _write_seed_store(
        store_dir,
        {
            "IL6": ("HGNC:6018", 5.32),
            "APOE": ("HGNC:613", 5.29),
        },
    )

    def _responder(name: str, kwargs: dict[str, Any]) -> MappingResult:
        if name == "IL6":
            return _mapping_result(name, curie="HGNC:6018", score=5.3)
        if name == "Interleukin-6":
            return _mapping_result(name, curie="HGNC:6018", score=5.1)
        return _mapping_result(name)

    mapper = ScriptedEngineMapper(_responder)
    cohorts = {
        "ukbb": [{"key": "u1", "name": "IL6"}, {"key": "u2", "name": "APOE"}],
        "arivale": [{"key": "a1", "name": "IL6"}],
        "twinsuk": [
            {"SeqId": "1000-1", "Target": "IL6", "TargetFullName": "Interleukin-6"},
        ],
    }

    report1 = harmonize_many(
        cohorts,
        entity="biolink:Protein",
        name_columns={"ukbb": "name", "arivale": "name", "twinsuk": "Target"},
        key_columns={"ukbb": "key", "arivale": "key", "twinsuk": "SeqId"},
        extra_name_columns={"twinsuk": ["TargetFullName"]},
        resolution_policies={
            "twinsuk": default_protein_engine_policy(extra_name_columns=("TargetFullName",))
        },
        strip_names={"ukbb": True, "arivale": True, "twinsuk": False},
        mapper=mapper,
        probe_pins=False,
        pins=TEST_RAW_PINS,
        store_dir=store_dir,
        output_dir=tmp_path / "out1",
        save=True,
    )

    # ukbb and arivale reused 100% from seed_store, twinsuk ran 100% live despite sharing "IL6".
    prov1 = report1.mapping["names_only"]
    assert prov1["ukbb"]["n_reused"] == 2
    assert prov1["ukbb"]["n_fresh"] == 0
    assert prov1["arivale"]["n_reused"] == 1
    assert prov1["arivale"]["n_fresh"] == 0
    assert prov1["twinsuk"]["n_reused"] == 0
    assert prov1["twinsuk"]["n_fresh"] == 1
    called_queries = [q for names, _ in mapper.calls for q in names]
    assert "APOE" not in called_queries
    assert "IL6" in called_queries
    assert "Interleukin-6" in called_queries

    # Second pass reuses both Olink seed hits and newly persisted TwinsUK engine hits.
    mapper2 = ScriptedEngineMapper(_responder)
    report2 = harmonize_many(
        cohorts,
        entity="biolink:Protein",
        name_columns={"ukbb": "name", "arivale": "name", "twinsuk": "Target"},
        key_columns={"ukbb": "key", "arivale": "key", "twinsuk": "SeqId"},
        extra_name_columns={"twinsuk": ["TargetFullName"]},
        resolution_policies={
            "twinsuk": default_protein_engine_policy(extra_name_columns=("TargetFullName",))
        },
        strip_names={"ukbb": True, "arivale": True, "twinsuk": False},
        mapper=mapper2,
        probe_pins=False,
        pins=TEST_RAW_PINS,
        store_dir=store_dir,
        output_dir=tmp_path / "out2",
        save=True,
    )
    assert report2.mapping["names_only"]["twinsuk"]["n_reused"] == 1
    assert report2.mapping["names_only"]["twinsuk"]["n_fresh"] == 0
    assert mapper2.calls == []


def test_strict_zero_string_editing_invariant(tmp_path: Path) -> None:
    raw_target = "  C3b / iC3b  "
    raw_full = "  Complement C3b/iC3b fragment (alpha' chain)  "

    def _responder(name: str, kwargs: dict[str, Any]) -> MappingResult:
        if name == raw_full:
            return _mapping_result(name, curie="HGNC:1318", score=4.9)
        return _mapping_result(name)

    mapper = ScriptedEngineMapper(_responder)
    cohorts = {
        "ukbb": [{"key": "u1", "name": "C3"}],
        "twinsuk": [{"SeqId": "2000-1", "Target": raw_target, "TargetFullName": raw_full}],
    }
    report = harmonize_many(
        cohorts,
        entity="biolink:Protein",
        name_columns={"ukbb": "name", "twinsuk": "Target"},
        key_columns={"ukbb": "key", "twinsuk": "SeqId"},
        extra_name_columns={"twinsuk": ["TargetFullName"]},
        resolution_policies={
            "twinsuk": default_protein_engine_policy(extra_name_columns=("TargetFullName",))
        },
        strip_names={"ukbb": True, "twinsuk": False},
        mapper=mapper,
        probe_pins=False,
        pins=TEST_RAW_PINS,
        store_dir=tmp_path / "store",
        save=False,
    )

    twinsuk_queries = [
        q for names, _ in mapper.calls for q in names if q != "C3"
    ]
    assert raw_target in twinsuk_queries
    assert raw_full in twinsuk_queries
    assert raw_target.strip() not in twinsuk_queries
    assert raw_full.strip() not in twinsuk_queries
    assert report.arms["names_only"].results["twinsuk"][0].chosen_kg_id == "HGNC:1318"


def test_tiered_engine_cascade_fallback_and_abstention(tmp_path: Path) -> None:
    def _responder(name: str, kwargs: dict[str, Any]) -> MappingResult:
        annotators = kwargs.get("annotators")
        etype = kwargs.get("entity_type")
        # 1) Rescue via TargetFullName when Target ("3HAO") misses.
        if name == "3-hydroxyanthranilate 3,4-dioxygenase" and etype == "biolink:Protein":
            return _mapping_result(name, curie="HGNC:4796", score=4.8)
        # 2) Rescue via kestrel-text-search on tier 2.
        if name == "LEXICAL_ONLY" and list(annotators or ()) == ["kestrel-text-search"]:
            return _mapping_result(
                name,
                curie="HGNC:1001",
                score=3.4,
                annotator="kestrel-text-search",
            )
        # 3) Rescue via biolink:Gene on tier 3.
        if name == "GENE_ONLY" and etype == "biolink:Gene":
            return _mapping_result(name, curie="HGNC:1002", score=4.2)
        # 4) Rescue via biolink:GeneOrGeneProduct on tier 4.
        if name == "GOGP_ONLY" and etype == "biolink:GeneOrGeneProduct":
            return _mapping_result(name, curie="HGNC:1003", score=3.8)
        # 5) Rescue via kestrel-vector-search on tier 5 when score >= 0.75.
        if name == "VECTOR_OK" and list(annotators or ()) == ["kestrel-vector-search"]:
            return _mapping_result(
                name,
                curie="HGNC:1004",
                score=1.1,
                annotator="kestrel-vector-search",
            )
        # 6) Weak vector hit (< 0.75) must abstain rather than accept noise.
        if name == "VECTOR_WEAK" and list(annotators or ()) == ["kestrel-vector-search"]:
            return _mapping_result(
                name,
                curie="HGNC:9999",
                score=0.55,
                annotator="kestrel-vector-search",
            )
        return _mapping_result(name)

    mapper = ScriptedEngineMapper(_responder)
    cohorts = {
        "ukbb": [{"key": "u0", "name": "IL6"}],
        "twinsuk": [
            {
                "SeqId": "s1",
                "Target": "3HAO",
                "TargetFullName": "3-hydroxyanthranilate 3,4-dioxygenase",
            },
            {"SeqId": "s2", "Target": "LEXICAL_ONLY", "TargetFullName": "LEXICAL_ONLY"},
            {"SeqId": "s3", "Target": "GENE_ONLY", "TargetFullName": "GENE_ONLY"},
            {"SeqId": "s4", "Target": "GOGP_ONLY", "TargetFullName": "GOGP_ONLY"},
            {"SeqId": "s5", "Target": "VECTOR_OK", "TargetFullName": "VECTOR_OK"},
            {"SeqId": "s6", "Target": "VECTOR_WEAK", "TargetFullName": "VECTOR_WEAK"},
        ],
    }
    report = harmonize_many(
        cohorts,
        entity="biolink:Protein",
        name_columns={"ukbb": "name", "twinsuk": "Target"},
        key_columns={"ukbb": "key", "twinsuk": "SeqId"},
        extra_name_columns={"twinsuk": ["TargetFullName"]},
        resolution_policies={
            "twinsuk": default_protein_engine_policy(extra_name_columns=("TargetFullName",))
        },
        strip_names={"ukbb": True, "twinsuk": False},
        mapper=mapper,
        probe_pins=False,
        pins=TEST_RAW_PINS,
        store_dir=tmp_path / "store",
        save=False,
    )

    arm = report.arms["names_only"]
    by_id = {
        rec.key: res
        for rec, res in zip(arm.records["twinsuk"], arm.results["twinsuk"], strict=True)
    }
    assert by_id["s1"].chosen_kg_id == "HGNC:4796"
    assert by_id["s2"].chosen_kg_id == "HGNC:1001"
    assert by_id["s3"].chosen_kg_id == "HGNC:1002"
    assert by_id["s4"].chosen_kg_id == "HGNC:1003"
    assert by_id["s5"].chosen_kg_id == "HGNC:1004"
    assert by_id["s6"].chosen_kg_id is None


def test_human_preference_and_dual_column_conflict_resolution(tmp_path: Path) -> None:
    def _responder(name: str, kwargs: dict[str, Any]) -> MappingResult:
        # Primary hits a non-human MGI node even with high score, while secondary
        # (or lexical tier) resolves an HGNC human node -> human node wins.
        if name == "HUMAN_PREF_TARGET":
            return _mapping_result(
                name,
                curie="MGI:1098",
                score=5.4,
                equivalent_ids={"MGI": ["MGI:1098"]},
            )
        if name == "Human Interleukin 6":
            return _mapping_result(name, curie="HGNC:6018", score=4.9)
        # Dual-column tie / near-tie: TargetFullName ("Thrombopoietin") wins specificity tie.
        if name == "Tpo":
            return _mapping_result(name, curie="HGNC:12015", score=4.1)
        if name == "Thrombopoietin":
            return _mapping_result(name, curie="HGNC:11785", score=4.3)
        # Large Target lead >= conflict_score_delta keeps Target.
        if name == "STRONG_TARGET":
            return _mapping_result(name, curie="HGNC:5001", score=5.2)
        if name == "WEAK_FULL":
            return _mapping_result(name, curie="HGNC:5002", score=3.6)
        return _mapping_result(name)

    mapper = ScriptedEngineMapper(_responder)
    cohorts = {
        "ukbb": [{"key": "u0", "name": "IL6"}],
        "twinsuk": [
            {
                "SeqId": "h1",
                "Target": "HUMAN_PREF_TARGET",
                "TargetFullName": "Human Interleukin 6",
            },
            {"SeqId": "t1", "Target": "Tpo", "TargetFullName": "Thrombopoietin"},
            {"SeqId": "c1", "Target": "STRONG_TARGET", "TargetFullName": "WEAK_FULL"},
        ],
    }
    report = harmonize_many(
        cohorts,
        entity="biolink:Protein",
        name_columns={"ukbb": "name", "twinsuk": "Target"},
        key_columns={"ukbb": "key", "twinsuk": "SeqId"},
        extra_name_columns={"twinsuk": ["TargetFullName"]},
        resolution_policies={
            "twinsuk": default_protein_engine_policy(extra_name_columns=("TargetFullName",))
        },
        strip_names={"ukbb": True, "twinsuk": False},
        mapper=mapper,
        probe_pins=False,
        pins=TEST_RAW_PINS,
        store_dir=tmp_path / "store",
        save=False,
    )

    arm = report.arms["names_only"]
    by_id = {
        rec.key: res
        for rec, res in zip(arm.records["twinsuk"], arm.results["twinsuk"], strict=True)
    }
    assert by_id["h1"].chosen_kg_id == "HGNC:6018"
    assert by_id["t1"].chosen_kg_id == "HGNC:11785"
    assert by_id["c1"].chosen_kg_id == "HGNC:5001"


def test_subrequest_error_propagates_and_never_becomes_unresolved(tmp_path: Path) -> None:
    store_dir = tmp_path / "store"

    def _responder(name: str, kwargs: dict[str, Any]) -> MappingResult:
        if name == "FAIL_ME":
            return _mapping_result(name, error="HTTP 502 upstream error")
        return _mapping_result(name, curie="HGNC:6018", score=5.2)

    mapper = ScriptedEngineMapper(_responder)
    report = harmonize_many(
        {
            "ukbb": [{"key": "u0", "name": "IL6"}],
            "twinsuk": [
                {"SeqId": "s_ok", "Target": "IL6", "TargetFullName": "Interleukin-6"},
                {"SeqId": "s_err", "Target": "FAIL_ME", "TargetFullName": "Fail Protein"},
            ],
        },
        entity="biolink:Protein",
        name_columns={"ukbb": "name", "twinsuk": "Target"},
        key_columns={"ukbb": "key", "twinsuk": "SeqId"},
        extra_name_columns={"twinsuk": ["TargetFullName"]},
        resolution_policies={
            "twinsuk": default_protein_engine_policy(extra_name_columns=("TargetFullName",))
        },
        strip_names={"ukbb": True, "twinsuk": False},
        mapper=mapper,
        probe_pins=False,
        pins=TEST_RAW_PINS,
        store_dir=store_dir,
        output_dir=tmp_path / "out",
        save=True,
    )

    prov = report.mapping["names_only"]["twinsuk"]
    assert prov["n_errored"] == 1
    arm = report.arms["names_only"]
    by_id = {
        rec.key: res
        for rec, res in zip(arm.records["twinsuk"], arm.results["twinsuk"], strict=True)
    }
    assert by_id["s_ok"].chosen_kg_id == "HGNC:6018"
    assert by_id["s_err"].error == "HTTP 502 upstream error"

    # Errored item must not be written into EntityCache as a finished miss.
    cache_files = list(store_dir.rglob("entity_cache.jsonl"))
    all_lines = "\n".join(p.read_text(encoding="utf-8") for p in cache_files)
    assert "FAIL_ME" not in all_lines


@pytest.mark.asyncio
async def test_harmonize_runner_forwards_cohort_engine_options_and_resumes(tmp_path: Path) -> None:
    input_dir = tmp_path / "inputs"
    input_dir.mkdir(parents=True, exist_ok=True)
    ukbb_tsv = input_dir / "ukbb.tsv"
    ukbb_tsv.write_text("key\tname\nu1\tIL6\n", encoding="utf-8")
    twinsuk_tsv = input_dir / "twinsuk.tsv"
    twinsuk_tsv.write_text(
        "SeqId\tTarget\tTargetFullName\n1000-1\tIL6\tInterleukin-6\n",
        encoding="utf-8",
    )

    seen_bodies: list[dict[str, Any]] = []

    with respx.mock(assert_all_called=False) as router:
        router.get("https://api.example.test/api/v1/health").mock(
            return_value=httpx.Response(
                200,
                json={"status": "healthy", "version": "1.9.0"},
            )
        )
        router.get("https://kestrel.example.test/api/health").mock(
            return_value=httpx.Response(
                200,
                json={
                    "status": "ok",
                    "kestrel_version": TEST_PIN["kestrel_version"],
                    "kg_build": {
                        "kg_version": TEST_PIN["kg_version"],
                        "kraken_package_version": TEST_PIN["kg_version"],
                        "git_commit": TEST_PIN["kg_git_commit"],
                    },
                },
            )
        )

        def _batch(request: httpx.Request) -> httpx.Response:
            body = json.loads(request.content.decode("utf-8"))
            seen_bodies.append(body)
            entries = [_batch_api_entry(item["name"]) for item in body["entities"]]
            return httpx.Response(200, json={"results": entries, "metadata": {}, "summary": {}})

        router.post("https://api.example.test/api/v1/map/batch").mock(side_effect=_batch)

        config = ServerConfig(
            base_url="https://api.example.test/api/v1",
            kestrel_url="https://kestrel.example.test/api",
            home=tmp_path / "mcp_home",
            env_roots=[tmp_path],
        )
        async with BioMapperClient(
            base_url=config.base_url,
            api_key=None,
            timeout=10.0,
        ) as client:
            manager = JobManager(config, gate=RequestGate(), client=client)
            runner = HarmonizeRunner(manager)
            policy = CohortResolutionPolicy(
                tiers=(
                    ResolutionTier(
                        name="protein_hybrid",
                        entity_type="biolink:Protein",
                        annotators=("kestrel-hybrid-search",),
                        prefer_human=True,
                        prefer_canonical=True,
                        candidate_limit=20,
                    ),
                ),
                name_columns=("TargetFullName",),
                corroborate_columns="selective",
                strip_names=False,
            )
            search_opts = EngineSearchOptions(
                entity_type="biolink:Protein",
                annotators=("kestrel-hybrid-search",),
                prefer_human=True,
                prefer_canonical=True,
                candidate_limit=20,
            )
            specs = [
                CohortSpec(
                    label="ukbb",
                    path=str(ukbb_tsv),
                    name_column="name",
                    key_column="key",
                ),
                CohortSpec(
                    label="twinsuk",
                    path=str(twinsuk_tsv),
                    name_column="Target",
                    key_column="SeqId",
                    extra_name_columns=["TargetFullName"],
                    search_options=search_opts.to_dict(),
                    resolution_policy=policy.to_dict(),
                    strip_names=False,
                ),
            ]

            store_dir = str(tmp_path / "mcp_store")
            worker = asyncio.create_task(manager.run_worker())
            try:
                status = await runner.submit(
                    kind="harmonize_many",
                    specs=specs,
                    entity="biolink:Protein",
                    category=None,
                    output_dir=str(tmp_path / "mcp_out1"),
                    store_dir=store_dir,
                    remap_all=False,
                    client_roots=[tmp_path],
                )

                for _ in range(250):
                    listing = await manager.status(status.job_id, client_roots=[tmp_path])
                    cur = listing.harmonize_jobs[0]
                    if cur.state in ("complete", "partial", "failed", "stalled"):
                        status = cur
                        break
                    await asyncio.sleep(0.02)

                assert status.state == "complete", f"unexpected state {status.state}: {status.reason}"
                twinsuk_requests = [
                    ent
                    for body in seen_bodies
                    for ent in body["entities"]
                    if (ent.get("options") or {}).get("candidate_limit") == 20
                ]
                assert twinsuk_requests
                first_opts = twinsuk_requests[0]["options"]
                assert first_opts["prefer_human"] is True
                assert first_opts["prefer_canonical"] is True
                assert first_opts["candidate_limit"] == 20

                # Second submission with same store_dir reuses cached sub-requests (0 new HTTP batch calls).
                calls_before = len(seen_bodies)
                status2 = await runner.submit(
                    kind="harmonize_many",
                    specs=specs,
                    entity="biolink:Protein",
                    category=None,
                    output_dir=str(tmp_path / "mcp_out2"),
                    store_dir=store_dir,
                    remap_all=False,
                    client_roots=[tmp_path],
                )
                for _ in range(250):
                    listing2 = await manager.status(status2.job_id, client_roots=[tmp_path])
                    cur2 = listing2.harmonize_jobs[0]
                    if cur2.state in ("complete", "partial", "failed", "stalled"):
                        status2 = cur2
                        break
                    await asyncio.sleep(0.02)
                assert status2.state == "complete", f"unexpected state {status2.state}: {status2.reason}"
                assert len(seen_bodies) == calls_before
            finally:
                manager.begin_shutdown()
                worker.cancel()
                try:
                    await worker
                except asyncio.CancelledError:
                    pass
                manager.close()


def test_call_mapper_with_options_ignores_unsupported_kwargs() -> None:
    def _legacy_mapper(
        records: list[dict[str, Any]],
        *,
        entity_type: str = "biolink:Protein",
        annotation_mode: str = "missing",
    ) -> list[MappingResult]:
        return [_mapping_result(str(rec["name"]), curie="HGNC:6018") for rec in records]

    opts = EngineSearchOptions(
        entity_type="biolink:Gene",
        annotators=("kestrel-hybrid-search",),
        prefer_human=True,
        candidate_limit=15,
    )
    results = call_mapper_with_options(
        _legacy_mapper,
        [{"name": "IL6", "identifiers": {}}],
        entity_type="biolink:Protein",
        annotation_mode="missing",
        search_options=opts,
    )
    assert len(results) == 1
    assert results[0].chosen_kg_id == "HGNC:6018"
    outcome = PolicyOutcome(results=results, reused=[False], n_errored=0, diagnostics={})
    assert outcome.n_errored == 0


def test_import_biomapper_resolution_without_pandas(monkeypatch: pytest.MonkeyPatch) -> None:
    real_import = builtins.__import__
    saved_modules = {
        name: mod
        for name, mod in sys.modules.items()
        if name == "biomapper" or name.startswith("biomapper.") or name == "pandas"
    }

    def _no_pandas_import(name: str, *args: Any, **kwargs: Any) -> Any:
        if name == "pandas" or name.startswith("pandas."):
            raise ImportError("pandas blocked for no-pandas test")
        return real_import(name, *args, **kwargs)

    try:
        monkeypatch.setattr(builtins, "__import__", _no_pandas_import)
        for mod in list(saved_modules):
            sys.modules.pop(mod, None)

        resolution_mod = builtins.__import__(
            "biomapper.resolution",
            globals(),
            locals(),
            ["CohortResolutionPolicy", "default_protein_engine_policy"],
        )
        assert hasattr(resolution_mod, "CohortResolutionPolicy")
        assert hasattr(resolution_mod, "default_protein_engine_policy")
    finally:
        for mod in list(sys.modules):
            if mod == "biomapper" or mod.startswith("biomapper."):
                sys.modules.pop(mod, None)
        sys.modules.update(saved_modules)
