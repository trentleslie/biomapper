"""Cross-dataset equivalence linking over two already-resolved datasets.

``link_by_intersection``, :class:`Link` and the overlap counters are ported from the engine's
``studies/external_benchmarks/scorers/cross_cohort_overlap.py`` (read at ``origin/dev``). The
refusal accounting — naming the entities that did not resolve, and keeping an API error apart
from a non-resolution — is written for this package, because a client-side harmonization report
has to be able to say which inputs it could not judge.

Fully offline: it consumes CURIE sets that are already resolved (real ones from a live mapping
run, or literals in tests) and never calls the API or the knowledge graph.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from typing import Any

from biomapper.harmonize.curies import curie_set
from biomapper.models import MappingResult


@dataclass(frozen=True)
class Link:
    """One harmonized pair. ``shared`` carries the CURIE(s) that formed it, for audit."""

    a_key: str
    b_key: str
    shared: frozenset[str]


@dataclass(frozen=True)
class OverlapResult:
    """Outcome of intersecting two sides' CURIE sets."""

    links: tuple[Link, ...]
    n_links: int  # distinct (a, b) linked pairs
    n_a_linked: int  # distinct A-side keys appearing in a link
    n_b_linked: int  # distinct B-side keys appearing in a link
    n_a_comparable: int  # A-side keys with a non-empty CURIE set (the shared denominator)
    n_b_comparable: int  # B-side keys with a non-empty CURIE set
    a_unresolved: tuple[str, ...]  # A-side keys whose CURIE set was EMPTY — refusal candidates
    b_unresolved: tuple[str, ...]  # B-side keys whose CURIE set was EMPTY — refusal candidates


def link_by_intersection(
    a_curies: dict[str, frozenset[str]],
    b_curies: dict[str, frozenset[str]],
) -> OverlapResult:
    """Link A<->B entities whose normalized CURIE sets intersect.

    Uses an inverted CURIE->keys index so cost is O(total CURIEs) rather than O(|A|.|B|). A pair
    sharing several CURIEs yields ONE link carrying all the shared CURIEs. The comparable
    denominator is the count of keys with a non-empty CURIE set on each side, so an entity that
    never resolved is never scored as a failure to link — it is reported as a refusal candidate in
    ``a_unresolved`` / ``b_unresolved`` instead of being silently dropped.
    """
    idx_b: dict[str, set[str]] = defaultdict(set)
    for b_key, curies in b_curies.items():
        for c in curies:
            idx_b[c].add(b_key)

    pair_shared: dict[tuple[str, str], set[str]] = defaultdict(set)
    for a_key, curies in a_curies.items():
        for c in curies:
            for b_key in idx_b.get(c, ()):
                pair_shared[(a_key, b_key)].add(c)

    links = tuple(
        Link(a_key=a, b_key=b, shared=frozenset(shared))
        for (a, b), shared in sorted(pair_shared.items())
    )
    return OverlapResult(
        links=links,
        n_links=len(links),
        n_a_linked=len({lk.a_key for lk in links}),
        n_b_linked=len({lk.b_key for lk in links}),
        n_a_comparable=sum(1 for s in a_curies.values() if s),
        n_b_comparable=sum(1 for s in b_curies.values() if s),
        a_unresolved=tuple(k for k, s in a_curies.items() if not s),
        b_unresolved=tuple(k for k, s in b_curies.items() if not s),
    )


def curie_sets_from_results(
    results: Iterable[MappingResult],
    key: Callable[[MappingResult, int], str] | None = None,
) -> dict[str, frozenset[str]]:
    """Build the linker's input from mapping results: ``{key: identifier-only CURIE set}``.

    ``key`` defaults to the result's ``query_name``. A duplicate key raises ``ValueError`` rather
    than collapsing two distinct entities onto one entry, which would silently drop one of them;
    pass a ``key`` callable (it receives the result and its index) when a cohort genuinely has
    repeated names.
    """
    out: dict[str, frozenset[str]] = {}
    for i, r in enumerate(results):
        k = key(r, i) if key is not None else r.query_name
        if k in out:
            raise ValueError(
                f"duplicate key {k!r} in the input results. Two entities would collapse onto one "
                "entry and one of them would be dropped. Pass key=... to disambiguate."
            )
        out[k] = curie_set(r.chosen_kg_id, r.kg_equivalent_ids)
    return out


@dataclass(frozen=True)
class HarmonizationResult:
    """Cross-dataset equivalence between two cohorts, plus what could not be judged.

    Every input entity lands in exactly one bucket per side:

    - **comparable** — resolved to at least one identifier, so it could link (whether it did or not)
    - **unresolved** — resolved to nothing. A refusal candidate, never a link, never a miss.
    - **errored** — the mapping call itself failed. "We do not know", which is a different claim
      from "it did not resolve", so the two are counted apart.
    """

    overlap: OverlapResult
    a_label: str
    b_label: str
    a_errors: tuple[str, ...]
    b_errors: tuple[str, ...]
    n_a_total: int
    n_b_total: int

    # -- linking ---------------------------------------------------------

    @property
    def links(self) -> tuple[Link, ...]:
        return self.overlap.links

    @property
    def n_links(self) -> int:
        return self.overlap.n_links

    @property
    def n_a_linked(self) -> int:
        return self.overlap.n_a_linked

    @property
    def n_b_linked(self) -> int:
        return self.overlap.n_b_linked

    # -- what could link -------------------------------------------------

    @property
    def n_a_comparable(self) -> int:
        return self.overlap.n_a_comparable

    @property
    def n_b_comparable(self) -> int:
        return self.overlap.n_b_comparable

    # -- what could not ---------------------------------------------------

    @property
    def a_unresolved(self) -> tuple[str, ...]:
        """A-side entities that resolved to no identifier. Refusal candidates."""
        return self.overlap.a_unresolved

    @property
    def b_unresolved(self) -> tuple[str, ...]:
        """B-side entities that resolved to no identifier. Refusal candidates."""
        return self.overlap.b_unresolved

    @property
    def n_a_unresolved(self) -> int:
        return len(self.a_unresolved)

    @property
    def n_b_unresolved(self) -> int:
        return len(self.b_unresolved)

    @property
    def n_a_errors(self) -> int:
        return len(self.a_errors)

    @property
    def n_b_errors(self) -> int:
        return len(self.b_errors)

    # -- rates -------------------------------------------------------------

    @property
    def a_link_rate(self) -> float | None:
        """Linked A-side entities over COMPARABLE A-side entities; ``None`` if none comparable.

        The denominator excludes unresolved entities deliberately: an entity that resolved to
        nothing had no opportunity to link, so counting it as a miss would understate the rate and
        conflate non-resolution with non-equivalence.
        """
        return (self.n_a_linked / self.n_a_comparable) if self.n_a_comparable else None

    @property
    def b_link_rate(self) -> float | None:
        """Linked B-side entities over COMPARABLE B-side entities; ``None`` if none comparable."""
        return (self.n_b_linked / self.n_b_comparable) if self.n_b_comparable else None

    def summary(self) -> dict[str, Any]:
        """Counts-only summary, keyed by the two cohort labels. Safe to log or serialize."""
        return {
            "n_links": self.n_links,
            self.a_label: {
                "total": self.n_a_total,
                "comparable": self.n_a_comparable,
                "unresolved": self.n_a_unresolved,
                "errors": self.n_a_errors,
                "linked": self.n_a_linked,
                "link_rate": self.a_link_rate,
            },
            self.b_label: {
                "total": self.n_b_total,
                "comparable": self.n_b_comparable,
                "unresolved": self.n_b_unresolved,
                "errors": self.n_b_errors,
                "linked": self.n_b_linked,
                "link_rate": self.b_link_rate,
            },
        }


def _split_errors(
    results: Sequence[MappingResult],
    key: Callable[[MappingResult, int], str] | None,
) -> tuple[list[MappingResult], tuple[str, ...]]:
    """Partition results into (mappable, errored-keys). An error is not a non-resolution."""
    mappable: list[MappingResult] = []
    errored: list[str] = []
    for i, r in enumerate(results):
        k = key(r, i) if key is not None else r.query_name
        if r.error:
            errored.append(k)
        else:
            mappable.append(r)
    return mappable, tuple(errored)


def harmonize(
    a_results: Sequence[MappingResult],
    b_results: Sequence[MappingResult],
    *,
    a_label: str = "a",
    b_label: str = "b",
    key: Callable[[MappingResult, int], str] | None = None,
) -> HarmonizationResult:
    """Harmonize two already-resolved datasets by cross-dataset equivalence.

    Two entities are equivalent when they resolve to the same canonical KRAKEN node, detected as a
    non-empty intersection of their identifier-only CURIE sets. This is a pure set operation over
    results you already have; it issues no requests.

    Args:
        a_results: Mapping results for cohort A (e.g. from :func:`biomapper.map_entities`).
        b_results: Mapping results for cohort B.
        a_label:   Name for cohort A in :meth:`HarmonizationResult.summary`.
        b_label:   Name for cohort B in :meth:`HarmonizationResult.summary`.
        key:       Optional ``(result, index) -> str`` key. Defaults to ``query_name``; supply one
                   when a cohort has repeated names (a duplicate key otherwise raises).

    Returns:
        A :class:`HarmonizationResult`. Entities that resolved to nothing are reported in
        ``a_unresolved`` / ``b_unresolved`` as refusal candidates, and entities whose mapping call
        errored in ``a_errors`` / ``b_errors``. Neither is silently dropped.

    Raises:
        ValueError: If either side has duplicate keys (see ``key``).
    """
    a_mappable, a_errors = _split_errors(a_results, key)
    b_mappable, b_errors = _split_errors(b_results, key)

    # The index passed to `key` is the position within the ERROR-FREE subset, so a custom key must
    # not assume it lines up with the caller's original list. It stays unique, which is all the
    # duplicate-key guard needs.
    overlap = link_by_intersection(
        curie_sets_from_results(a_mappable, key=key),
        curie_sets_from_results(b_mappable, key=key),
    )
    return HarmonizationResult(
        overlap=overlap,
        a_label=a_label,
        b_label=b_label,
        a_errors=a_errors,
        b_errors=b_errors,
        n_a_total=len(a_results),
        n_b_total=len(b_results),
    )
