"""Re-adjudication of non-certified cross-cohort cases (offline; a stub outside source).

The assertions that matter are the ones that stop a first-block mismatch from being reported as a
wrong molecule, and the ones that keep "we could not check" apart from "we checked and it is fine".
"""

from __future__ import annotations

import pandas as pd

from biomapper.benchmarks.cross_cohort_readjudicate import (
    OutsideRecord,
    OutsideResolver,
    classify,
    formulas_differ_only_in_hydrogen,
    readjudicate,
)

GLUCOSE = "WQZGKKKJIJFFOK"
OTHER = "ZZZZZZZZZZZZZZ"
THIRD = "QQQQQQQQQQQQQQ"


def _hit(block: str, formula: str, mass: float) -> OutsideRecord:
    return OutsideRecord(block=block, formula=formula, mass=mass, status="success")


_MISS = OutsideRecord(None, None, None, "clean_miss")
_FAILED = OutsideRecord(None, None, None, "lookup_failed")


# ==================================================================================================
# Formula comparison
# ==================================================================================================


def test_hydrogen_only_difference_is_recognized():
    # A protonation or tautomer shift. The compounds are the same for harmonization.
    assert formulas_differ_only_in_hydrogen("C7H9N2O", "C7H8N2O") is True
    assert formulas_differ_only_in_hydrogen("C6H12O6", "C6H12O6") is True


def test_real_composition_difference_is_not_excused():
    assert formulas_differ_only_in_hydrogen("C6H12O6", "C6H12O5") is False
    assert formulas_differ_only_in_hydrogen("C6H12O6", "C5H12O6") is False


def test_charge_suffix_does_not_defeat_the_comparison():
    assert formulas_differ_only_in_hydrogen("C7H9N2O+", "C7H8N2O") is True


def test_unparseable_formula_never_buys_a_free_pass():
    assert formulas_differ_only_in_hydrogen(None, "C6H12O6") is False
    assert formulas_differ_only_in_hydrogen("", "C6H12O6") is False
    assert formulas_differ_only_in_hydrogen("not a formula", "C6H12O6") is False


# ==================================================================================================
# Classification
# ==================================================================================================


def test_same_composition_with_different_blocks_is_an_artifact_not_an_error():
    outcome, rationale = classify(
        GLUCOSE,
        OTHER,
        _hit(GLUCOSE, "C6H12O6", 180.0634),
        _hit(OTHER, "C6H13O6", 181.0712),
    )
    assert outcome == "tautomer_or_charge_or_salt_artifact"
    assert "not invariant" in rationale


def test_outside_source_backing_the_cohort_implicates_the_gold():
    outcome, _ = classify(
        OTHER,  # what the curated gold says
        GLUCOSE,  # what the cohort's vendor id says
        _hit(GLUCOSE, "C6H12O6", 180.0634),  # outside source, NECS name -> glucose
        _hit(GLUCOSE, "C6H12O6", 180.0634),
    )
    assert outcome == "necs_gold_suspect"


def test_outside_source_backing_the_gold_implicates_the_cohort_annotation():
    outcome, _ = classify(
        GLUCOSE,
        OTHER,
        _hit(GLUCOSE, "C6H12O6", 180.0634),
        _hit(GLUCOSE, "C6H12O6", 180.0634),
    )
    assert outcome == "cohort_id_suspect"


def test_neither_side_matching_is_escalated_not_guessed():
    outcome, _ = classify(
        GLUCOSE,
        OTHER,
        _hit(THIRD, "C9H8O4", 180.0423),
        _hit(THIRD, "C9H8O4", 180.0423),
    )
    assert outcome == "both_sides_disagree_with_outside"


def test_both_sides_corroborated_at_different_structures_is_a_real_wrong_link():
    outcome, rationale = classify(
        GLUCOSE,
        OTHER,
        _hit(GLUCOSE, "C6H12O6", 180.0634),
        _hit(OTHER, "C9H8O4", 180.0423),
    )
    assert outcome == "genuine_structural_disagreement"
    assert "genuinely measured different molecules" in rationale


def test_unresolvable_outside_source_concludes_nothing():
    assert classify(GLUCOSE, OTHER, _MISS, _MISS)[0] == "outside_source_unresolved"
    assert classify(GLUCOSE, OTHER, _FAILED, _FAILED)[0] == "outside_source_unresolved"
    # One side only is also not enough to place the defect.
    assert (
        classify(GLUCOSE, OTHER, _hit(GLUCOSE, "C6H12O6", 180.0634), _MISS)[0]
        == "outside_source_unresolved"
    )


def test_mass_disagreement_blocks_the_artifact_call():
    # Same formula string but masses far apart: not one compound, so no free pass.
    outcome, _ = classify(
        GLUCOSE,
        OTHER,
        _hit(GLUCOSE, "C6H12O6", 180.0634),
        _hit(OTHER, "C6H12O6", 200.0),
    )
    assert outcome != "tautomer_or_charge_or_salt_artifact"


# ==================================================================================================
# Driving it over a case table
# ==================================================================================================


class _StubResolver:
    def __init__(self, table: dict[str, OutsideRecord]) -> None:
        self.table = table
        self.asked: list[str] = []

    def by_name(self, name: str) -> OutsideRecord:
        self.asked.append(name)
        return self.table.get(name, _MISS)


def _cases() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "necs_name": "certified_one",
                "cohort_name": "certified_one",
                "verdict": "certified",
                "refusal_class": "",
                "necs_block": GLUCOSE,
                "cohort_block": GLUCOSE,
            },
            {
                "necs_name": "refuted_one",
                "cohort_name": "refuted_other",
                "verdict": "refuted",
                "refusal_class": "",
                "necs_block": GLUCOSE,
                "cohort_block": OTHER,
            },
            {
                "necs_name": "names_only",
                "cohort_name": "names_only",
                "verdict": "refused",
                "refusal_class": "cohort_panel_has_no_structure_resolvable_id",
                "necs_block": GLUCOSE,
                "cohort_block": "",
            },
        ]
    )


def test_readjudicate_skips_certified_and_never_calls_out_for_construction_refusals():
    resolver = _StubResolver(
        {
            "refuted_one": _hit(GLUCOSE, "C6H12O6", 180.0634),
            "refuted_other": _hit(GLUCOSE, "C6H12O6", 180.0634),
        }
    )
    result = readjudicate(_cases(), resolver)  # type: ignore[arg-type]
    assert len(result) == 2  # the certified case is not re-examined
    outcomes = dict(zip(result["necs_name"], result["readjudication"], strict=True))
    assert outcomes["refuted_one"] == "cohort_id_suspect"
    assert outcomes["names_only"] == "not_adjudicable_by_construction"
    # A panel with no structure-resolvable identifier is never sent to the outside source: there is
    # nothing to adjudicate, and asking would invent a verdict the certificate never made.
    assert "names_only" not in resolver.asked


def test_readjudicate_records_the_outside_evidence_for_review():
    resolver = _StubResolver(
        {
            "refuted_one": _hit(THIRD, "C9H8O4", 180.0423),
            "refuted_other": _hit(OTHER, "C6H12O6", 180.0634),
        }
    )
    result = readjudicate(_cases(), resolver)  # type: ignore[arg-type]
    row = result[result["necs_name"] == "refuted_one"].iloc[0]
    assert row["outside_necs_block"] == THIRD
    assert row["outside_cohort_block"] == OTHER
    assert row["outside_necs_formula"] == "C9H8O4"
    assert row["rationale"]


def test_outside_resolver_treats_an_ambiguous_name_as_a_miss():
    body = (
        '{"PropertyTable": {"Properties": ['
        '{"InChIKey": "WQZGKKKJIJFFOK-GASJEMHNSA-N", "MolecularFormula": "C6H12O6",'
        ' "MonoisotopicMass": "180.06"},'
        '{"InChIKey": "ZZZZZZZZZZZZZZ-YYYYYYYYYY-N", "MolecularFormula": "C6H12O6",'
        ' "MonoisotopicMass": "180.06"}]}}'
    )
    # Two different structures behind one name: adjudicating off PubChem's rank 1 would be a guess.
    assert OutsideResolver._parse(body).status == "clean_miss"


def test_outside_resolver_parses_a_single_hit():
    body = (
        '{"PropertyTable": {"Properties": [{"InChIKey": "WQZGKKKJIJFFOK-GASJEMHNSA-N",'
        ' "MolecularFormula": "C6H12O6", "MonoisotopicMass": "180.063388"}]}}'
    )
    record = OutsideResolver._parse(body)
    assert record.resolved is True
    assert record.block == GLUCOSE and record.formula == "C6H12O6"
    assert record.mass is not None and abs(record.mass - 180.063388) < 1e-6


def test_outside_resolver_reports_a_malformed_body_as_a_failure_not_a_miss():
    assert OutsideResolver._parse("not json").status == "lookup_failed"
    assert OutsideResolver._parse('{"PropertyTable": {"Properties": []}}').status == "clean_miss"
