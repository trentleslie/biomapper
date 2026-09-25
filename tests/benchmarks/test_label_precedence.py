"""Label precedence: the weaker claim wins, because overstating is the error that misleads."""

from __future__ import annotations

import pytest

from biomapper.benchmarks.suite import weakest_claim


@pytest.mark.parametrize(
    "declared,circularity,expected,why",
    [
        # RefMet: role defaults to "accuracy", but refmet is an ingested KRAKEN source.
        ("accuracy", "coverage", "coverage", "a config default must not outrank a live verdict"),
        # MetaboliteAnnotator: coverage BY CONSTRUCTION (name-hit rate). Circularity cannot see
        # that, because the gold source is not in the graph, so it says accuracy_candidate.
        (
            "coverage",
            "accuracy_candidate",
            "coverage",
            "an arm that is coverage by construction stays coverage",
        ),
        ("partly_circular", "accuracy_candidate", "partly_circular", "partly-circular is weaker"),
        ("capability_regression", "coverage", "capability_regression", "both weakest, stable"),
        ("accuracy", "accuracy_candidate", "accuracy_candidate", "candidate is the weaker claim"),
        ("", "coverage", "coverage", "one source only"),
        ("accuracy", "", "accuracy", "one source only"),
        ("", "", "", "no source"),
        ("something_new", "accuracy", "something_new", "an unknown label is not strong evidence"),
    ],
)
def test_weakest_claim(declared: str, circularity: str, expected: str, why: str) -> None:
    assert weakest_claim(declared, circularity) == expected, why


def test_never_upgrades_a_coverage_arm() -> None:
    """The property that matters: coverage in, coverage out, whatever the other source says."""
    for other in ("accuracy", "accuracy_candidate", "partly_circular", "coverage", ""):
        assert weakest_claim("coverage", other) == "coverage"
