"""Cross-cohort certification and refusal adjudication (offline; no PubChem, no KG).

The load-bearing assertions here are about refusal accounting, because refusal is the outcome the
write-up has to be able to explain: a names-only panel refuses by construction, a failed lookup is a
transient run artifact that must not be reported as an absent structure, and a first-block
comparison never claims to have checked stereochemistry.
"""

from __future__ import annotations

import pandas as pd

from biomapper.benchmarks.cross_cohort_certify import (
    adjudicate_cases,
    first_block,
    read_links,
)
from biomapper.benchmarks.scorers.cross_cohort_overlap import Link
from biomapper.benchmarks.scorers.independent_inchikey import ProvidedBlock

GLUCOSE_BLOCK = "WQZGKKKJIJFFOK"
OTHER_BLOCK = "ZZZZZZZZZZZZZZ"


def test_first_block_extracts_connectivity_layer_only():
    assert first_block("WQZGKKKJIJFFOK-GASJEMHNSA-N") == GLUCOSE_BLOCK
    assert first_block("WQZGKKKJIJFFOK") == GLUCOSE_BLOCK
    assert first_block("") is None and first_block(None) is None


def _necs(block: str | None, name: str = "glucose") -> ProvidedBlock:
    return ProvidedBlock(block, "gold-necs-moesm5", "success", record_id=f"moesm5:{name}")


def _cohort(
    block: str | None,
    status: str = "success",
    source: str = "provided-pubchem",
    name: str = "glucose",
) -> ProvidedBlock:
    return ProvidedBlock(block, source, status, record_id=f"arivale:{name}")


def _link(a: str = "glucose", b: str = "glucose") -> Link:
    return Link(a_name=a, b_name=b, shared=frozenset({"CHEBI:17234"}))


def test_agreeing_blocks_certify_but_never_claim_stereo():
    cases = adjudicate_cases(
        [_link()], {"glucose": _necs(GLUCOSE_BLOCK)}, {"glucose": _cohort(GLUCOSE_BLOCK)}
    )
    row = cases.iloc[0]
    assert row["verdict"] == "certified"
    # Both sides are first blocks, so the stereo layer was not compared. Saying otherwise would
    # pass a stereoisomer error silently.
    assert bool(row["stereo_checked"]) is False


def test_connectivity_disagreement_is_refuted_with_both_blocks_recorded():
    cases = adjudicate_cases(
        [_link()], {"glucose": _necs(GLUCOSE_BLOCK)}, {"glucose": _cohort(OTHER_BLOCK)}
    )
    row = cases.iloc[0]
    assert row["verdict"] == "refuted"
    assert row["necs_block"] == GLUCOSE_BLOCK and row["cohort_block"] == OTHER_BLOCK
    assert row["refusal_class"] == ""


def test_missing_necs_gold_is_refused_and_classified():
    cases = adjudicate_cases([_link()], {}, {"glucose": _cohort(GLUCOSE_BLOCK)})
    row = cases.iloc[0]
    assert row["verdict"] == "refused"
    assert row["refusal_class"] == "necs_gold_has_no_curated_inchikey"


def test_transient_lookup_failure_is_not_reported_as_absent_structure():
    cases = adjudicate_cases(
        [_link()],
        {"glucose": _necs(GLUCOSE_BLOCK)},
        {"glucose": _cohort(None, status="lookup_failed")},
    )
    row = cases.iloc[0]
    assert row["verdict"] == "refused"
    assert row["refusal_class"] == "cohort_lookup_failed_transient"


def test_clean_miss_is_distinguished_from_a_failure():
    cases = adjudicate_cases(
        [_link()],
        {"glucose": _necs(GLUCOSE_BLOCK)},
        {"glucose": _cohort(None, status="clean_miss")},
    )
    assert cases.iloc[0]["refusal_class"] == "cohort_lookup_clean_miss"


def test_panel_without_structure_resolvable_id_is_classified_separately():
    cases = adjudicate_cases(
        [_link()],
        {"glucose": _necs(GLUCOSE_BLOCK)},
        {"glucose": _cohort(None, status="no_structure_resolvable_id", source="none")},
    )
    assert cases.iloc[0]["refusal_class"] == "cohort_panel_has_no_structure_resolvable_id"


def test_both_sides_missing_is_its_own_class():
    cases = adjudicate_cases([_link()], {}, {})
    assert cases.iloc[0]["refusal_class"] == "no_independent_structure_either_side"


def test_same_curator_record_refuses_rather_than_self_certifying():
    shared = ProvidedBlock(GLUCOSE_BLOCK, "provided-hmdb", "success", record_id="same:glucose")
    cases = adjudicate_cases([_link()], {"glucose": shared}, {"glucose": shared})
    row = cases.iloc[0]
    assert row["verdict"] == "refused" and "same curator record" in row["reason"]


def test_read_links_round_trips_the_driver_output(tmp_path):
    path = tmp_path / "links_necs_arivale.csv"
    pd.DataFrame(
        [
            {
                "necs_name": "glucose",
                "arivale_name": "glucose",
                "shared_curies": "CHEBI:17234|KEGG:C00031",
            }
        ]
    ).to_csv(path, index=False)
    links = read_links(path, "arivale")
    assert len(links) == 1
    assert links[0].a_name == "glucose" and links[0].b_name == "glucose"
    assert links[0].shared == frozenset({"CHEBI:17234", "KEGG:C00031"})


def test_read_links_handles_an_empty_link_file(tmp_path):
    path = tmp_path / "links_necs_llfs.csv"
    pd.DataFrame(columns=["necs_name", "llfs_name", "shared_curies"]).to_csv(path, index=False)
    assert read_links(path, "llfs") == []
