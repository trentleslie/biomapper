"""Edge cases of the LMSD capability/regression gate not covered by test_scorers.py."""

from __future__ import annotations

import pytest

from biomapper.benchmarks.scorers.regression import (
    assert_capability_floor,
    capability_resolvability,
)


def _regime(name: str, total: int, fraction: float) -> dict:
    return {"by_name_source_regime": {name: {"coverage": {"total": total, "fraction": fraction}}}}


@pytest.mark.parametrize(
    "result",
    [
        {},
        {"by_name_source_regime": None},
        {"by_name_source_regime": {}},
        {
            "by_name_source_regime": {
                "common_systematic": {"coverage": {"total": 5, "fraction": 1.0}}
            }
        },
        {"by_name_source_regime": {"shorthand": {}}},
        {"by_name_source_regime": {"shorthand": None}},
    ],
    ids=["key-missing", "none", "empty-dict", "regime-absent", "regime-empty", "regime-none"],
)
def test_resolvability_is_none_when_regime_missing(result):
    assert capability_resolvability(result) is None


def test_resolvability_is_none_when_regime_total_is_zero():
    assert capability_resolvability(_regime("shorthand", 0, 1.0)) is None


def test_resolvability_is_none_when_coverage_lacks_total():
    result = {"by_name_source_regime": {"shorthand": {"coverage": {"fraction": 0.9}}}}
    assert capability_resolvability(result) is None


def test_resolvability_reads_a_non_default_regime():
    result = {
        "by_name_source_regime": {
            "shorthand": {"coverage": {"total": 10, "fraction": 0.1}},
            "common_systematic": {"coverage": {"total": 20, "fraction": 0.75}},
        }
    }
    assert capability_resolvability(result, regime="common_systematic") == pytest.approx(0.75)
    assert capability_resolvability(result) == pytest.approx(0.1)


def test_resolvability_never_falls_back_to_blended_coverage():
    result = {
        "coverage": {"total": 100, "fraction": 0.99},
        "by_name_source_regime": {
            "common_systematic": {"coverage": {"total": 100, "fraction": 0.99}}
        },
    }
    assert capability_resolvability(result, regime="shorthand") is None
    with pytest.raises(ValueError, match="no 'shorthand' observations"):
        assert_capability_floor(result, 0.5, regime="shorthand")


@pytest.mark.parametrize(
    "result",
    [{}, {"by_name_source_regime": None}, _regime("shorthand", 0, 1.0)],
    ids=["absent", "none", "zero-total"],
)
def test_floor_raises_no_observations_when_regime_absent_or_empty(result):
    with pytest.raises(ValueError, match="no 'shorthand' observations"):
        assert_capability_floor(result, 0.0)


def test_floor_no_observations_message_names_the_requested_regime():
    with pytest.raises(ValueError, match="no 'common_systematic' observations"):
        assert_capability_floor(_regime("shorthand", 10, 1.0), 0.5, regime="common_systematic")


def test_floor_raises_floor_not_met_when_below():
    with pytest.raises(ValueError, match="floor not met") as excinfo:
        assert_capability_floor(_regime("shorthand", 100, 0.899), 0.90)
    assert "0.899" in str(excinfo.value)
    assert "0.900" in str(excinfo.value)


def test_floor_passes_when_value_equals_floor():
    assert assert_capability_floor(_regime("shorthand", 100, 0.90), 0.90) is None


def test_floor_passes_at_equality_for_a_non_default_regime():
    assert_capability_floor(_regime("common_systematic", 4, 0.5), 0.5, regime="common_systematic")
