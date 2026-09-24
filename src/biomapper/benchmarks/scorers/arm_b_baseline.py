"""Unit 3 — Arm-B baseline reconstruction (Monti et al. 2026's published per-pair method).

Ported from the biomapper2 engine at ``origin/dev``
(``studies/external_benchmarks/scorers/arm_b_baseline.py``, commit
``1ffb571e54fe028ef0ae4e748fc2e7ec093ee603``). The reconstruction logic is unchanged.

**The published-overlap table IS changed, and that is the point of this docstring.** The engine
constant read ``{"arivale": 615, "xuetal": 432, "llfs": 163, "blsa": 99}`` and attributed all four
to "Monti Table 2". Re-read against the paper PDF (DOI 10.1007/s11357-026-02174-2), two of the four
values are wrong and the citation points at the wrong table. Both the corrected and the superseded
values are kept below so a number can never move without the reason being visible. See
:data:`MONTI_PUBLISHED_PROVENANCE`.

Arm B is the number BioMapper is measured against. It is reconstructed here on the IDENTICAL row
set the M/M+ID arms use, so the comparison is like-for-like (R21a). Two methods, per the paper's
own ``Datasets harmonization`` text:

  - NECS↔Arivale, NECS↔Xu : Metabolon ``CHEMICAL_NAME`` string match (Arivale case-insensitive,
    Xu case-sensitive — the settings that reproduced the paper in the 2026-08-20 replication).
  - NECS↔LLFS, NECS↔BLSA  : RefMet standardized-name join, with ``drop_na(refmet_name)`` applied
    BEFORE the join (a name that does not RefMet-standardize cannot match).

Because WE build this baseline, it is a controlled variable: it must be **frozen** (locked as a
characterization test at a recorded commit) before any BioMapper live run, and the recovery claim
must beat **Monti-published** — the number we did NOT compute — by more than the per-pair
``|re-derived − published|`` reconstruction gap. This module exposes that gap on every result.

Fully offline: pure set logic over name lists and a precomputed RefMet map (a dict, loaded by the
caller from the persisted cache — never a network call here).
"""

from __future__ import annotations

from dataclasses import dataclass

# Monti et al. 2026 published NECS<->cohort overlaps (the un-gamed comparator; we did not compute
# these). Corrected 2026-09-23 by re-reading the paper. The values live in the Methods section
# "Datasets harmonization" and in the per-cohort descriptions, NOT in Table 2 (Table 2 is "Age-only
# markers"). Two entries in the engine's table were wrong:
#
#   xuetal: the engine had 432, taken from the Xu cohort description ("432 metabolites were matched
#           to metabolites in the NECS dataset"). The Methods harmonization sentence says 385:
#           "when comparing NECS with Arivale and with Xu et al., we used the Metabolon-provided
#           CHEMICAL_NAME's, which yielded an overlap of 615 and 385 metabolites, respectively."
#           The paper contradicts itself on this one pair. 385 is used here because it comes from
#           the sentence describing the harmonization PROCEDURE this module reconstructs, and that
#           same sentence's Arivale value (615) agrees with the Arivale cohort description. 432 is
#           retained in MONTI_PUBLISHED_PROVENANCE so the conflict is reportable, never hidden.
#
#   blsa:   the engine had 99, which is NOT a NECS overlap. The BLSA cohort description says
#           "Ninety-nine metabolites were in common with the LLFS using RefMet standardized names"
#           - that is BLSA<->LLFS. Scoring a NECS<->BLSA reconstruction against it compares two
#           different pairs. The NECS<->BLSA value is 188, from the same Methods sentence as LLFS:
#           "When comparing NECS with LLFS and BLSA, we mapped all metabolites to RefMet
#           identifiers, which yielded an overlap of 163 and 188 metabolites, respectively."
#
# Caveat worth carrying into any write-up: 188 is also the paper's count of LIPID metabolites in the
# LLFS panel ("408 metabolites (188 lipid and 220 polar)"). The coincidence is not evidence of a
# transcription error in the paper, and the harmonization sentence is explicit, so 188 is taken at
# face value. Flag it rather than silently preferring one reading.
MONTI_PUBLISHED: dict[str, int] = {"arivale": 615, "xuetal": 385, "llfs": 163, "blsa": 188}

# Superseded engine values, kept so a changed number is always traceable to a reason.
MONTI_PUBLISHED_SUPERSEDED: dict[str, int] = {
    "arivale": 615,
    "xuetal": 432,
    "llfs": 163,
    "blsa": 99,
}

# Per-pair citation for every published value, so a manifest can carry the quote rather than a bare
# integer. ``conflict`` is populated only where the paper itself disagrees with itself.
MONTI_PUBLISHED_PROVENANCE: dict[str, dict[str, object]] = {
    "arivale": {
        "published": 615,
        "quote": "which yielded an overlap of 615 and 385 metabolites, respectively",
        "section": "Methods, 'Datasets harmonization'",
        "corroborated_by": "Arivale cohort description: 'After name curation, 615 metabolites were matched to metabolites in the NECS dataset'",
        "conflict": None,
    },
    "xuetal": {
        "published": 385,
        "quote": "which yielded an overlap of 615 and 385 metabolites, respectively",
        "section": "Methods, 'Datasets harmonization'",
        "corroborated_by": None,
        "conflict": {
            "alternate": 432,
            "quote": "A total of 821 lipid and polar metabolites were generated using non-targeted metabolomic analysis, and 432 metabolites were matched to metabolites in the NECS dataset",
            "section": "Xu et al. cohort description",
            "note": "the paper states both 385 and 432 for NECS<->Xu; unresolved in the source",
        },
    },
    "llfs": {
        "published": 163,
        "quote": "which yielded an overlap of 163 and 188 metabolites, respectively",
        "section": "Methods, 'Datasets harmonization'",
        "corroborated_by": "LLFS cohort description: 'Of these, 163 metabolites were matched to the NECS metabolites based on RefMet standardized names'",
        "conflict": None,
    },
    "blsa": {
        "published": 188,
        "quote": "which yielded an overlap of 163 and 188 metabolites, respectively",
        "section": "Methods, 'Datasets harmonization'",
        "corroborated_by": None,
        "conflict": {
            "alternate": 99,
            "quote": "Ninety-nine metabolites were in common with the LLFS using RefMet standardized names",
            "section": "BLSA cohort description",
            "note": "99 is the BLSA<->LLFS overlap, a different pair; it was the engine's value for NECS<->BLSA and is wrong for this comparison",
        },
    },
}

# Per-pair method: ("name", case_sensitive) or ("refmet",). Unknown cohort → fail loud.
# The union type is spelled out rather than left as ``tuple[str, ...]`` (the engine's annotation) so
# the ``case_sensitive`` element keeps its bool type through to ``name_match_overlap``.
PAIR_METHOD: dict[str, tuple[str, bool] | tuple[str]] = {
    "arivale": ("name", False),
    "xuetal": ("name", True),
    "llfs": ("refmet",),
    "blsa": ("refmet",),
}


@dataclass(frozen=True)
class ArmBResult:
    cohort: str
    method: str
    count: int  # re-derived overlap on the identical row set
    published: int  # Monti Table 2
    gap: int  # count - published; the error bar recovery must exceed


def _name_key(name: str, case_sensitive: bool) -> str:
    s = name.strip()
    return s if case_sensitive else s.lower()


def name_match_overlap(names_a: list[str], names_b: list[str], *, case_sensitive: bool) -> set[str]:
    """Overlap by exact CHEMICAL_NAME string match (Arm B for same-vendor pairs)."""
    a = {_name_key(n, case_sensitive) for n in names_a if n.strip()}
    b = {_name_key(n, case_sensitive) for n in names_b if n.strip()}
    return a & b


def refmet_join_overlap(
    names_a: list[str],
    names_b: list[str],
    refmet_map: dict[str, str],
) -> set[str]:
    """Overlap by RefMet standardized-name join, drop_na(refmet_name) BEFORE the join.

    ``refmet_map`` maps a raw name to its RefMet name; a name absent from the map or mapping to
    an empty string does not standardize and is dropped before the intersection (exactly Monti's
    ``drop_na`` step — a non-standardizing name can never match).
    """

    def refmet_set(names: list[str]) -> set[str]:
        out: set[str] = set()
        for n in names:
            r = refmet_map.get(n.strip(), "").strip()
            if r:
                out.add(r.lower())
        return out

    return refmet_set(names_a) & refmet_set(names_b)


def arm_b_overlap(
    cohort: str,
    necs_names: list[str],
    cohort_names: list[str],
    *,
    refmet_map: dict[str, str] | None = None,
) -> ArmBResult:
    """Compute Arm B for one NECS↔cohort pair using the cohort's published method. Fails loud on
    an unknown cohort (never defaults to name-match)."""
    if cohort not in PAIR_METHOD:
        raise ValueError(
            f"no Arm-B method registered for cohort {cohort!r}; known: {sorted(PAIR_METHOD)}"
        )
    method_name, *rest = PAIR_METHOD[cohort]
    if method_name == "name":
        case_sensitive = bool(rest[0]) if rest else False
        matched = name_match_overlap(necs_names, cohort_names, case_sensitive=case_sensitive)
        method = f"CHEMICAL_NAME ({'case-sensitive' if case_sensitive else 'case-insensitive'})"
    else:
        if refmet_map is None:
            raise ValueError(
                f"cohort {cohort!r} uses the RefMet method but no refmet_map was provided"
            )
        matched = refmet_join_overlap(necs_names, cohort_names, refmet_map)
        method = "RefMet join (drop_na before join)"
    published = MONTI_PUBLISHED[cohort]
    return ArmBResult(
        cohort=cohort,
        method=method,
        count=len(matched),
        published=published,
        gap=len(matched) - published,
    )
