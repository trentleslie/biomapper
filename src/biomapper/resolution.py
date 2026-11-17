"""Engine-only multi-tier resolution over unedited cohort columns.

Resolves cohorts whose primary name column may contain aliases, abbreviations, or
multi-word descriptions alongside secondary unedited name columns (such as
``Target`` and ``TargetFullName`` in SomaScan panels) without editing, stripping,
splitting, or rewriting any input string.

All resolution decisions come strictly from BioMapper engine calls across:

* engine search flags (``prefer_human``, ``prefer_canonical``, ``candidate_limit``),
* annotator subsets (``kestrel-hybrid-search``, ``kestrel-text-search``,
  ``kestrel-vector-search``),
* Biolink entity types (``biolink:Protein``, ``biolink:Gene``,
  ``biolink:GeneOrGeneProduct``), and
* dual-column corroboration and confidence/specificity adjudication when two unedited
  columns disagree.

Every individual sub-request is deduplicated and persisted through
:class:`biomapper._checkpoint.EntityCache` keyed by exact query string, identifiers,
and request options, preserving build-pinning, resumability, and zero interference
with standard single-tier cache entries.
"""

from __future__ import annotations

import inspect
from collections import Counter, defaultdict
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Literal

from biomapper._checkpoint import EntityCache, cache_key
from biomapper.client import BioMapperClient
from biomapper.harmonize.curies import curie_set, is_structure_curie
from biomapper.models import MappingResult

if TYPE_CHECKING:
    from biomapper._mapping_store import ArmStore
    from biomapper.cohorts import Cohort, CohortRecord, Mapper

__all__ = [
    "CohortResolutionPolicy",
    "EngineSearchOptions",
    "PolicyOutcome",
    "ResolutionTier",
    "default_protein_engine_policy",
    "is_human_node",
    "is_usable_node",
    "nodes_agree",
    "resolve_cohort_with_policy",
]

CorroborationMode = Literal["always", "selective", "fallback"]

_DEFAULT_SCORE_THRESHOLD = 3.0
_DEFAULT_FALLBACK_SCORE = 1.5
_DEFAULT_STRONG_SCORE = 4.5
_DEFAULT_CONFLICT_DELTA = 0.8
_DEFAULT_VECTOR_SCORE = 0.75


def _as_str_tuple(value: str | Sequence[str] | None) -> tuple[str, ...] | None:
    if value is None:
        return None
    if isinstance(value, str):
        return (value,)
    return tuple(str(v) for v in value)


@dataclass(frozen=True)
class EngineSearchOptions:
    """Engine search options forwarded to ``BioMapperClient.map_entities``.

    Any option left as ``None`` is omitted from both the API request and the
    cache key so default callers stay byte-identical to pre-existing behavior.
    """

    entity_type: str | None = None
    annotators: tuple[str, ...] | None = None
    prefer_human: bool | None = None
    prefer_canonical: bool | None = None
    candidate_limit: int | None = None
    kestrel_top_n: int | None = None
    vocab: str | tuple[str, ...] | None = None
    array_delimiters: tuple[str, ...] | None = None

    def __post_init__(self) -> None:
        BioMapperClient._check_bounds("candidate_limit", self.candidate_limit)
        BioMapperClient._check_bounds("kestrel_top_n", self.kestrel_top_n)

    def to_cache_options(
        self, default_entity_type: str, annotation_mode: str
    ) -> dict[str, Any]:
        """Return the ``options`` dict embedded into ``EntityCache`` keys."""
        out: dict[str, Any] = {
            "entity_type": self.entity_type or default_entity_type,
            "annotation_mode": annotation_mode,
        }
        if self.annotators is not None:
            out["annotators"] = list(self.annotators)
        if self.prefer_human is not None:
            out["prefer_human"] = self.prefer_human
        if self.prefer_canonical is not None:
            out["prefer_canonical"] = self.prefer_canonical
        if self.candidate_limit is not None:
            out["candidate_limit"] = self.candidate_limit
        if self.kestrel_top_n is not None:
            out["kestrel_top_n"] = self.kestrel_top_n
        if self.vocab is not None:
            out["vocab"] = list(self.vocab) if isinstance(self.vocab, tuple) else self.vocab
        if self.array_delimiters is not None:
            out["array_delimiters"] = list(self.array_delimiters)
        return out

    def to_mapper_kwargs(self) -> dict[str, Any]:
        """Return non-``None`` request options for ``BatchedApiMapper`` / ``map_entities``."""
        kwargs: dict[str, Any] = {}
        if self.annotators is not None:
            kwargs["annotators"] = list(self.annotators)
        if self.prefer_human is not None:
            kwargs["prefer_human"] = self.prefer_human
        if self.prefer_canonical is not None:
            kwargs["prefer_canonical"] = self.prefer_canonical
        if self.candidate_limit is not None:
            kwargs["candidate_limit"] = self.candidate_limit
        if self.kestrel_top_n is not None:
            kwargs["kestrel_top_n"] = self.kestrel_top_n
        if self.vocab is not None:
            kwargs["vocab"] = list(self.vocab) if isinstance(self.vocab, tuple) else self.vocab
        if self.array_delimiters is not None:
            kwargs["array_delimiters"] = list(self.array_delimiters)
        return kwargs

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-safe representation."""
        out: dict[str, Any] = {}
        if self.entity_type is not None:
            out["entity_type"] = self.entity_type
        if self.annotators is not None:
            out["annotators"] = list(self.annotators)
        if self.prefer_human is not None:
            out["prefer_human"] = self.prefer_human
        if self.prefer_canonical is not None:
            out["prefer_canonical"] = self.prefer_canonical
        if self.candidate_limit is not None:
            out["candidate_limit"] = self.candidate_limit
        if self.kestrel_top_n is not None:
            out["kestrel_top_n"] = self.kestrel_top_n
        if self.vocab is not None:
            out["vocab"] = list(self.vocab) if isinstance(self.vocab, tuple) else self.vocab
        if self.array_delimiters is not None:
            out["array_delimiters"] = list(self.array_delimiters)
        return out

    @classmethod
    def from_value(
        cls, value: EngineSearchOptions | Mapping[str, Any]
    ) -> EngineSearchOptions:
        """Coerce a mapping or ``EngineSearchOptions`` into ``EngineSearchOptions``."""
        if isinstance(value, cls):
            return value
        if not isinstance(value, Mapping):
            raise TypeError(
                f"expected EngineSearchOptions or mapping, got {type(value).__name__}"
            )
        vocab_raw = value.get("vocab")
        vocab: str | tuple[str, ...] | None
        if vocab_raw is None or isinstance(vocab_raw, str):
            vocab = vocab_raw
        else:
            vocab = tuple(str(v) for v in vocab_raw)
        return cls(
            entity_type=value.get("entity_type"),
            annotators=_as_str_tuple(value.get("annotators")),
            prefer_human=value.get("prefer_human"),
            prefer_canonical=value.get("prefer_canonical"),
            candidate_limit=value.get("candidate_limit"),
            kestrel_top_n=value.get("kestrel_top_n"),
            vocab=vocab,
            array_delimiters=_as_str_tuple(value.get("array_delimiters")),
        )


@dataclass(frozen=True)
class ResolutionTier:
    """One search pass in a tiered engine-only resolution cascade."""

    name: str
    entity_type: str = "biolink:Protein"
    annotators: tuple[str, ...] | None = None
    prefer_human: bool | None = True
    prefer_canonical: bool | None = True
    candidate_limit: int | None = 25
    kestrel_top_n: int | None = None
    vocab: str | tuple[str, ...] | None = None
    min_accept_score: float | None = _DEFAULT_SCORE_THRESHOLD
    min_fallback_score: float | None = _DEFAULT_FALLBACK_SCORE
    require_human: bool = True
    allow_non_human_fallback: bool = True
    is_vector_only: bool = False

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("ResolutionTier.name must be non-empty")
        if not self.entity_type.startswith("biolink:"):
            raise ValueError(
                f"ResolutionTier.entity_type must be a Biolink category, got {self.entity_type!r}"
            )
        BioMapperClient._check_bounds("candidate_limit", self.candidate_limit)
        BioMapperClient._check_bounds("kestrel_top_n", self.kestrel_top_n)

    def search_options(self) -> EngineSearchOptions:
        """Return :class:`EngineSearchOptions` for this tier."""
        return EngineSearchOptions(
            entity_type=self.entity_type,
            annotators=self.annotators,
            prefer_human=self.prefer_human,
            prefer_canonical=self.prefer_canonical,
            candidate_limit=self.candidate_limit,
            kestrel_top_n=self.kestrel_top_n,
            vocab=self.vocab,
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-safe representation."""
        return {
            "name": self.name,
            "entity_type": self.entity_type,
            "annotators": list(self.annotators) if self.annotators is not None else None,
            "prefer_human": self.prefer_human,
            "prefer_canonical": self.prefer_canonical,
            "candidate_limit": self.candidate_limit,
            "kestrel_top_n": self.kestrel_top_n,
            "vocab": list(self.vocab) if isinstance(self.vocab, tuple) else self.vocab,
            "min_accept_score": self.min_accept_score,
            "min_fallback_score": self.min_fallback_score,
            "require_human": self.require_human,
            "allow_non_human_fallback": self.allow_non_human_fallback,
            "is_vector_only": self.is_vector_only,
        }

    @classmethod
    def from_value(cls, value: ResolutionTier | Mapping[str, Any]) -> ResolutionTier:
        """Coerce a mapping or ``ResolutionTier`` into ``ResolutionTier``."""
        if isinstance(value, cls):
            return value
        if not isinstance(value, Mapping):
            raise TypeError(f"expected ResolutionTier or mapping, got {type(value).__name__}")
        vocab_raw = value.get("vocab")
        vocab: str | tuple[str, ...] | None
        if vocab_raw is None or isinstance(vocab_raw, str):
            vocab = vocab_raw
        else:
            vocab = tuple(str(v) for v in vocab_raw)
        return cls(
            name=str(value["name"]),
            entity_type=str(value.get("entity_type", "biolink:Protein")),
            annotators=_as_str_tuple(value.get("annotators")),
            prefer_human=value.get("prefer_human", True),
            prefer_canonical=value.get("prefer_canonical", True),
            candidate_limit=value.get("candidate_limit", 25),
            kestrel_top_n=value.get("kestrel_top_n"),
            vocab=vocab,
            min_accept_score=value.get("min_accept_score", _DEFAULT_SCORE_THRESHOLD),
            min_fallback_score=value.get("min_fallback_score", _DEFAULT_FALLBACK_SCORE),
            require_human=bool(value.get("require_human", True)),
            allow_non_human_fallback=bool(value.get("allow_non_human_fallback", True)),
            is_vector_only=bool(value.get("is_vector_only", False)),
        )


@dataclass(frozen=True)
class CohortResolutionPolicy:
    """Principled engine-only resolution policy for one cohort.

    ``name_columns`` names unedited secondary columns (such as
    ``("TargetFullName",)``) queried alongside the cohort's primary
    ``name_column``. Every query string is sent verbatim.
    """

    tiers: tuple[ResolutionTier, ...]
    name_columns: tuple[str, ...] = ()
    corroborate_columns: CorroborationMode = "always"
    strong_score_threshold: float = _DEFAULT_STRONG_SCORE
    conflict_score_delta: float = _DEFAULT_CONFLICT_DELTA
    prefer_specific_column_on_tie: bool = True
    strip_names: bool = False

    def __post_init__(self) -> None:
        if not self.tiers:
            raise ValueError("CohortResolutionPolicy.tiers must contain at least one tier")
        if self.corroborate_columns not in ("always", "selective", "fallback"):
            raise ValueError(
                f"unknown corroborate_columns mode {self.corroborate_columns!r}; expected "
                "'always', 'selective', or 'fallback'"
            )

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-safe representation."""
        return {
            "name_columns": list(self.name_columns),
            "corroborate_columns": self.corroborate_columns,
            "strong_score_threshold": self.strong_score_threshold,
            "conflict_score_delta": self.conflict_score_delta,
            "prefer_specific_column_on_tie": self.prefer_specific_column_on_tie,
            "strip_names": self.strip_names,
            "tiers": [t.to_dict() for t in self.tiers],
        }

    @classmethod
    def from_value(
        cls, value: CohortResolutionPolicy | Mapping[str, Any]
    ) -> CohortResolutionPolicy:
        """Coerce a mapping or ``CohortResolutionPolicy`` into ``CohortResolutionPolicy``."""
        if isinstance(value, cls):
            return value
        if not isinstance(value, Mapping):
            raise TypeError(
                f"expected CohortResolutionPolicy or mapping, got {type(value).__name__}"
            )
        tiers = tuple(ResolutionTier.from_value(t) for t in value["tiers"])
        return cls(
            tiers=tiers,
            name_columns=tuple(str(c) for c in (value.get("name_columns") or ())),
            corroborate_columns=value.get("corroborate_columns", "always"),
            strong_score_threshold=float(
                value.get("strong_score_threshold", _DEFAULT_STRONG_SCORE)
            ),
            conflict_score_delta=float(
                value.get("conflict_score_delta", _DEFAULT_CONFLICT_DELTA)
            ),
            prefer_specific_column_on_tie=bool(
                value.get("prefer_specific_column_on_tie", True)
            ),
            strip_names=bool(value.get("strip_names", False)),
        )


def default_protein_engine_policy(
    extra_name_columns: Sequence[str] = ("TargetFullName",),
    *,
    candidate_limit: int = 25,
    corroborate_columns: CorroborationMode = "always",
    include_vector_fallback: bool = True,
) -> CohortResolutionPolicy:
    """Build the 5-tier engine-only protein harmonization policy.

    Tiers execute in precision-first order:
    1. ``protein_hybrid``: ``biolink:Protein`` + ``kestrel-hybrid-search``.
    2. ``protein_text``: ``biolink:Protein`` + ``kestrel-text-search`` (strict lexical).
    3. ``gene_hybrid``: ``biolink:Gene`` + ``kestrel-hybrid-search``.
    4. ``gene_or_product_hybrid``: ``biolink:GeneOrGeneProduct`` + ``kestrel-hybrid-search``.
    5. ``protein_vector`` (optional): ``biolink:Protein`` + ``kestrel-vector-search``.
    """
    tiers: list[ResolutionTier] = [
        ResolutionTier(
            name="protein_hybrid",
            entity_type="biolink:Protein",
            annotators=("kestrel-hybrid-search",),
            prefer_human=True,
            prefer_canonical=True,
            candidate_limit=candidate_limit,
            min_accept_score=_DEFAULT_SCORE_THRESHOLD,
            min_fallback_score=_DEFAULT_FALLBACK_SCORE,
            require_human=True,
            allow_non_human_fallback=True,
        ),
        ResolutionTier(
            name="protein_text",
            entity_type="biolink:Protein",
            annotators=("kestrel-text-search",),
            prefer_human=True,
            prefer_canonical=True,
            candidate_limit=candidate_limit,
            min_accept_score=None,
            min_fallback_score=None,
            require_human=True,
            allow_non_human_fallback=True,
        ),
        ResolutionTier(
            name="gene_hybrid",
            entity_type="biolink:Gene",
            annotators=("kestrel-hybrid-search",),
            prefer_human=True,
            prefer_canonical=True,
            candidate_limit=candidate_limit,
            min_accept_score=_DEFAULT_SCORE_THRESHOLD,
            min_fallback_score=_DEFAULT_FALLBACK_SCORE,
            require_human=True,
            allow_non_human_fallback=True,
        ),
        ResolutionTier(
            name="gene_or_product_hybrid",
            entity_type="biolink:GeneOrGeneProduct",
            annotators=("kestrel-hybrid-search",),
            prefer_human=True,
            prefer_canonical=True,
            candidate_limit=candidate_limit,
            min_accept_score=_DEFAULT_SCORE_THRESHOLD,
            min_fallback_score=_DEFAULT_FALLBACK_SCORE,
            require_human=True,
            allow_non_human_fallback=True,
        ),
    ]
    if include_vector_fallback:
        tiers.append(
            ResolutionTier(
                name="protein_vector",
                entity_type="biolink:Protein",
                annotators=("kestrel-vector-search",),
                prefer_human=True,
                prefer_canonical=True,
                candidate_limit=candidate_limit,
                min_accept_score=_DEFAULT_VECTOR_SCORE,
                min_fallback_score=_DEFAULT_VECTOR_SCORE,
                require_human=True,
                allow_non_human_fallback=False,
                is_vector_only=True,
            )
        )
    return CohortResolutionPolicy(
        tiers=tuple(tiers),
        name_columns=tuple(extra_name_columns),
        corroborate_columns=corroborate_columns,
        strong_score_threshold=_DEFAULT_STRONG_SCORE,
        conflict_score_delta=_DEFAULT_CONFLICT_DELTA,
        prefer_specific_column_on_tie=True,
        strip_names=False,
    )


# ---------------------------------------------------------------------------
# Node / agreement predicates
# ---------------------------------------------------------------------------


def is_usable_node(result: MappingResult | None) -> bool:
    """True when ``result`` resolved cleanly to a non-structure KG node."""
    if result is None or result.error:
        return False
    if not result.chosen_kg_id or is_structure_curie(result.chosen_kg_id):
        return False
    return bool(curie_set(result.chosen_kg_id, result.kg_equivalent_ids))


def is_human_node(result: MappingResult | None) -> bool:
    """True when ``result`` is usable and bears an HGNC human identifier."""
    if not is_usable_node(result) or result is None:
        return False
    chosen = (result.chosen_kg_id or "").upper()
    if chosen.startswith("HGNC:"):
        return True
    eq_hgnc = (result.kg_equivalent_ids or {}).get("HGNC")
    if eq_hgnc:
        return True
    id_hgnc = (result.identifiers or {}).get("HGNC")
    return bool(id_hgnc)


def nodes_agree(a: MappingResult | None, b: MappingResult | None) -> bool:
    """True when two usable mapping results describe the same KG node or CURIE set."""
    if not (is_usable_node(a) and is_usable_node(b)) or a is None or b is None:
        return False
    if a.chosen_kg_id == b.chosen_kg_id:
        return True
    set_a = curie_set(a.chosen_kg_id, a.kg_equivalent_ids)
    set_b = curie_set(b.chosen_kg_id, b.kg_equivalent_ids)
    return bool(set_a & set_b)


def _node_rank_key(result: MappingResult) -> tuple[int, int, float]:
    chosen_hgnc = 1 if (result.chosen_kg_id or "").upper().startswith("HGNC:") else 0
    human = 1 if is_human_node(result) else 0
    score = result.confidence_score if result.confidence_score is not None else 0.0
    return (chosen_hgnc, human, score)


def _pick_preferred_agreed(a: MappingResult, b: MappingResult) -> MappingResult:
    return a if _node_rank_key(a) >= _node_rank_key(b) else b


def _meets_score(score: float | None, threshold: float | None) -> bool:
    if threshold is None:
        return True
    if score is None:
        return False
    return score >= threshold


def _is_lexical_result(tier: ResolutionTier, result: MappingResult) -> bool:
    if tier.annotators and "kestrel-text-search" in tier.annotators:
        return True
    return "kestrel-text-search" in (result.curies_by_annotator or {})


# ---------------------------------------------------------------------------
# Mapper invocation & cache helpers
# ---------------------------------------------------------------------------


def _accepts_option_kwargs(mapper: Mapper) -> set[str] | None:
    """Return supported optional kwargs, or ``None`` if mapper accepts ``**kwargs``."""
    try:
        sig = inspect.signature(mapper)
    except (TypeError, ValueError):
        return set()
    params = sig.parameters
    if any(p.kind == inspect.Parameter.VAR_KEYWORD for p in params.values()):
        return None
    return set(params)


def call_mapper_with_options(
    mapper: Mapper,
    requests: list[dict[str, Any]],
    *,
    entity_type: str,
    annotation_mode: str,
    search_options: EngineSearchOptions | None = None,
) -> list[MappingResult]:
    """Call ``mapper``, forwarding non-``None`` engine options when accepted."""
    kwargs: dict[str, Any] = {
        "entity_type": search_options.entity_type or entity_type if search_options else entity_type,
        "annotation_mode": annotation_mode,
    }
    if search_options is not None:
        extra = search_options.to_mapper_kwargs()
        supported = _accepts_option_kwargs(mapper)
        if supported is None:
            kwargs.update(extra)
        else:
            for k, v in extra.items():
                if k in supported:
                    kwargs[k] = v
    results = list(mapper(requests, **kwargs))
    if len(results) != len(requests):
        raise RuntimeError(
            f"mapper returned {len(results)} results for {len(requests)} requests; a mapper must "
            "return one result per request, in order"
        )
    return results


def _request_key(
    name: str,
    identifiers: Mapping[str, Sequence[str]],
    options: Mapping[str, Any],
) -> str:
    flat = {v: "|".join(sorted(str(c) for c in codes)) for v, codes in identifiers.items()}
    return cache_key(name, None, flat, dict(options))


def _execute_subrequests(
    mapper: Mapper,
    items: list[tuple[str, dict[str, list[str]]]],
    *,
    cache: EntityCache | None,
    arm_name: str,
    next_batch_index: Callable[[], int],
    entity_type: str,
    annotation_mode: str,
    search_options: EngineSearchOptions,
    batch_size: int,
) -> tuple[dict[str, MappingResult], set[str]]:
    """Map unique ``(name, identifiers)`` queries for one tier, returning ``(by_key, disk_hits)``."""
    cache_opts = search_options.to_cache_options(entity_type, annotation_mode)
    found: dict[str, MappingResult] = {}
    disk_hits: set[str] = set()
    order: list[str] = []
    first: dict[str, tuple[str, dict[str, list[str]]]] = {}

    for name, ids in items:
        k = _request_key(name, ids, cache_opts)
        if k in found or k in first:
            continue
        hit = cache.get(k) if cache is not None else None
        if hit is not None:
            found[k] = hit
            disk_hits.add(k)
        else:
            first[k] = (name, dict(ids))
            order.append(k)

    step = max(batch_size, 1) if cache is not None else max(len(order), 1)
    for start in range(0, len(order), step):
        batch_keys = order[start : start + step]
        batch = [{"name": first[k][0], "identifiers": dict(first[k][1])} for k in batch_keys]
        results = call_mapper_with_options(
            mapper,
            batch,
            entity_type=entity_type,
            annotation_mode=annotation_mode,
            search_options=search_options,
        )
        keyed = list(zip(batch_keys, results, strict=True))
        if cache is not None:
            cache.commit_batch(vocab=arm_name, batch_index=next_batch_index(), keyed_results=keyed)
        found.update(keyed)

    return found, disk_hits


# ---------------------------------------------------------------------------
# Multi-tier / multi-column adjudication
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class _Observation:
    tier_index: int
    tier: ResolutionTier
    column_role: str  # "primary" or auxiliary column name
    query_text: str
    result: MappingResult


@dataclass
class _RowState:
    record: CohortRecord
    identifiers: dict[str, list[str]]
    secondary_queries: tuple[tuple[str, str], ...]  # ((column_name, verbatim_value), ...)
    observations: list[_Observation] = field(default_factory=list)
    final_result: MappingResult | None = None
    winner_tier: str | None = None
    winner_column: str | None = None
    errored: str | None = None
    reused_from_disk: bool = True
    queried_live: bool = False


def _secondary_queries_for(
    record: CohortRecord, policy: CohortResolutionPolicy
) -> tuple[tuple[str, str], ...]:
    cols = policy.name_columns or tuple(record.extra_names.keys())
    seen: set[str] = {record.name}
    out: list[tuple[str, str]] = []
    for col in cols:
        val = record.extra_names.get(col)
        if val is None or val in seen:
            continue
        seen.add(val)
        out.append((col, val))
    return tuple(out)


def _needs_corroboration_in_tier1(
    primary_res: MappingResult,
    tier: ResolutionTier,
    policy: CohortResolutionPolicy,
) -> bool:
    if policy.corroborate_columns == "always":
        return True
    if primary_res.error is not None:
        return True
    if not is_human_node(primary_res):
        return True
    if not _meets_score(primary_res.confidence_score, tier.min_accept_score):
        return True
    if policy.corroborate_columns == "fallback":
        return False
    # selective mode: corroborate unless primary hit strongly as a long uppercase symbol
    name = primary_res.query_name
    looks_long_upper = len(name) >= 5 and all(
        ch.isupper() or ch.isdigit() or ch == "-" for ch in name
    )
    is_strong = _meets_score(primary_res.confidence_score, policy.strong_score_threshold)
    return not (looks_long_upper and is_strong)


def _try_early_accept(
    state: _RowState,
    tier: ResolutionTier,
    primary_res: MappingResult,
    secondary_obs: list[tuple[str, MappingResult]],
    policy: CohortResolutionPolicy,
) -> tuple[MappingResult, str, str] | None:
    """Return ``(winner_result, tier_name, column_role)`` if row is settled at ``tier``."""
    if state.errored is not None:
        return None

    usable_sec = [(col, res) for col, res in secondary_obs if is_usable_node(res)]
    primary_usable = is_usable_node(primary_res)

    # 1. Dual-column agreement: primary and any secondary column resolve to the same node.
    if primary_usable:
        for col, sec_res in usable_sec:
            if nodes_agree(primary_res, sec_res) and (
                is_human_node(primary_res) or is_human_node(sec_res)
            ):
                winner = _pick_preferred_agreed(primary_res, sec_res)
                if not tier.is_vector_only or _meets_score(
                    max(
                        primary_res.confidence_score or 0.0,
                        sec_res.confidence_score or 0.0,
                    ),
                    tier.min_accept_score,
                ):
                    win_col = "primary" if winner is primary_res else col
                    return winner, tier.name, win_col

    # 2. No secondary column values (e.g., Target == TargetFullName): accept strong primary.
    if not state.secondary_queries:
        if is_human_node(primary_res) and _meets_score(
            primary_res.confidence_score, tier.min_accept_score
        ):
            return primary_res, tier.name, "primary"
        return None

    # 3. Only one side produced a human node and the other side did not conflict with a human node.
    human_sec = [(col, res) for col, res in usable_sec if is_human_node(res)]
    primary_human = is_human_node(primary_res)

    if primary_human and not human_sec:
        if _meets_score(primary_res.confidence_score, tier.min_accept_score):
            return primary_res, tier.name, "primary"
        return None

    if human_sec and not primary_human:
        # Take the strongest human secondary result if it meets acceptance score and no sec conflict
        if len(human_sec) == 1 or all(
            nodes_agree(human_sec[0][1], other) for _, other in human_sec[1:]
        ):
            best_col, best_res = max(
                human_sec, key=lambda item: _node_rank_key(item[1])
            )
            if _meets_score(best_res.confidence_score, tier.min_accept_score):
                return best_res, tier.name, best_col
        return None

    # 4. Both sides returned human nodes that DISAGREE.
    if primary_human and human_sec:
        best_sec_col, best_sec = max(human_sec, key=lambda item: _node_rank_key(item[1]))
        st = primary_res.confidence_score
        sf = best_sec.confidence_score
        if (
            not tier.is_vector_only
            and st is not None
            and sf is not None
            and tier.min_accept_score is not None
        ):
            if (
                st >= policy.strong_score_threshold
                and sf < tier.min_accept_score
                and (st - sf) >= policy.conflict_score_delta
            ):
                return primary_res, tier.name, "primary"
            if (
                sf >= policy.strong_score_threshold
                and st < tier.min_accept_score
                and (sf - st) >= policy.conflict_score_delta
            ):
                return best_sec, tier.name, best_sec_col
        # Genuine conflict: defer to subsequent tiers and final cluster adjudication.
        return None

    return None


def _adjudicate_observations(
    state: _RowState, policy: CohortResolutionPolicy
) -> tuple[MappingResult, str | None, str | None]:
    """Select the highest-evidence node across all gathered observations, or abstain."""
    if state.errored is not None:
        return MappingResult(query_name=state.record.name, error=state.errored), None, None

    eligible: list[_Observation] = []
    for obs in state.observations:
        res = obs.result
        if not is_usable_node(res):
            continue
        if obs.tier.require_human and not is_human_node(res) and not (
            obs.tier.allow_non_human_fallback
            and _meets_score(res.confidence_score, obs.tier.min_accept_score)
        ):
            continue
        eligible.append(obs)

    if not eligible:
        return MappingResult(query_name=state.record.name), None, None

    # Partition by human-ness first if any human observation exists.
    human_obs = [o for o in eligible if is_human_node(o.result)]
    pool = human_obs if human_obs else eligible

    # Build equivalence clusters across all (tier, column) observations.
    clusters: list[list[_Observation]] = []
    for obs in pool:
        placed = False
        for cluster in clusters:
            if any(nodes_agree(obs.result, existing.result) for existing in cluster):
                cluster.append(obs)
                placed = True
                break
        if not placed:
            clusters.append([obs])

    # Evaluate each cluster's evidence and filter uncorroborated sub-floor noise.
    candidates: list[dict[str, Any]] = []
    for cluster in clusters:
        cols_in_cluster = {o.column_role for o in cluster}
        dual_column_support = "primary" in cols_in_cluster and len(cols_in_cluster) > 1
        non_vector_obs = [o for o in cluster if not o.tier.is_vector_only]
        has_non_vector = bool(non_vector_obs)

        # Drop single-source noise below its tier's min_fallback_score unless corroborated.
        passes_floor = dual_column_support or any(
            _meets_score(o.result.confidence_score, o.tier.min_fallback_score)
            for o in cluster
        )
        if not passes_floor:
            continue

        # In vector-only clusters, require either dual-column agreement or vector accept score.
        if not has_non_vector and not dual_column_support:
            vector_ok = any(
                _meets_score(o.result.confidence_score, o.tier.min_accept_score)
                for o in cluster
            )
            if not vector_ok:
                continue

        has_lexical = any(_is_lexical_result(o.tier, o.result) for o in cluster)
        best_obs = max(
            cluster,
            key=lambda o: (
                0 if o.tier.is_vector_only else 1,
                1 if _is_lexical_result(o.tier, o.result) else 0,
                *_node_rank_key(o.result),
                -o.tier_index,
            ),
        )
        best_non_vec_score = max(
            (o.result.confidence_score or 0.0) for o in (non_vector_obs or cluster)
        )
        distinct_tiers = len({(o.tier.name, o.column_role) for o in cluster})
        has_secondary = any(o.column_role != "primary" for o in cluster)
        specificity_bonus = 0
        if policy.prefer_specific_column_on_tie and has_secondary:
            sec_texts = [o.query_text for o in cluster if o.column_role != "primary"]
            if sec_texts and max(len(t) for t in sec_texts) > len(state.record.name):
                specificity_bonus = 1

        candidates.append(
            {
                "has_non_vector": 1 if has_non_vector else 0,
                "dual_column_support": 1 if dual_column_support else 0,
                "has_lexical": 1 if has_lexical else 0,
                "distinct_tiers": distinct_tiers,
                "best_score": best_non_vec_score,
                "specificity_bonus": specificity_bonus,
                "best_obs": best_obs,
            }
        )

    if not candidates:
        return MappingResult(query_name=state.record.name), None, None

    # Guard against two conflicting pure vector-only hits when neither has dual support.
    if len(candidates) > 1 and all(not c["has_non_vector"] for c in candidates):
        by_vec_score = sorted(candidates, key=lambda c: c["best_score"], reverse=True)
        if not by_vec_score[0]["dual_column_support"]:
            gap = by_vec_score[0]["best_score"] - by_vec_score[1]["best_score"]
            if gap < policy.conflict_score_delta:
                return MappingResult(query_name=state.record.name), None, None

    # Narrow to top evidence class (non-vector, dual-column support), then within
    # conflict_score_delta apply lexical/tier/specificity tie-breaking cleanly without bin seams.
    best_class = max((c["has_non_vector"], c["dual_column_support"]) for c in candidates)
    top_class = [
        c for c in candidates if (c["has_non_vector"], c["dual_column_support"]) == best_class
    ]
    max_score = max(c["best_score"] for c in top_class)
    near_top = [
        c for c in top_class if (max_score - c["best_score"]) < policy.conflict_score_delta
    ]
    winner_entry = max(
        near_top,
        key=lambda c: (
            c["has_lexical"],
            c["specificity_bonus"],
            c["distinct_tiers"],
            c["best_score"],
            *_node_rank_key(c["best_obs"].result),
            -c["best_obs"].tier_index,
        ),
    )
    winner = winner_entry["best_obs"]
    return winner.result, winner.tier.name, winner.column_role


@dataclass(frozen=True)
class PolicyOutcome:
    """Per-row results, reuse accounting, and tier/column diagnostics."""

    results: list[MappingResult]
    reused: list[bool]
    n_errored: int
    diagnostics: dict[str, Any]


def resolve_cohort_with_policy(
    cohort: Cohort,
    sent_identifiers: Sequence[Mapping[str, Sequence[str]]],
    *,
    mapper: Mapper,
    store: ArmStore | None,
    policy: CohortResolutionPolicy,
    default_entity_type: str,
    annotation_mode: str = "missing",
    batch_size: int = 10,
    on_batch: Callable[[int, int], None] | None = None,
) -> PolicyOutcome:
    """Resolve every record of ``cohort`` through ``policy`` with cache resumability."""
    records = cohort.records
    if len(records) != len(sent_identifiers):
        raise ValueError(
            f"{cohort.label}: {len(records)} records vs {len(sent_identifiers)} identifier rows"
        )

    cache = store.cache if store is not None else None
    arm_name = store.arm if store is not None else "names_only"
    next_batch = store.next_batch_index if store is not None else (lambda: 0)
    effective_batch = store.batch_size if store is not None else batch_size

    states: list[_RowState] = [
        _RowState(
            record=rec,
            identifiers={v: list(codes) for v, codes in ids.items()},
            secondary_queries=_secondary_queries_for(rec, policy),
        )
        for rec, ids in zip(records, sent_identifiers, strict=True)
    ]

    pending: list[int] = list(range(len(states)))
    total_rows = len(states)
    if on_batch is not None:
        on_batch(0, total_rows)

    for tier_idx, tier in enumerate(policy.tiers):
        if not pending:
            break
        search_opts = tier.search_options()
        cache_opts = search_opts.to_cache_options(default_entity_type, annotation_mode)

        # 1. Primary column queries for all pending rows.
        primary_items = [
            (states[idx].record.name, states[idx].identifiers) for idx in pending
        ]
        primary_results, primary_disk_hits = _execute_subrequests(
            mapper,
            primary_items,
            cache=cache,
            arm_name=arm_name,
            next_batch_index=next_batch,
            entity_type=tier.entity_type,
            annotation_mode=annotation_mode,
            search_options=search_opts,
            batch_size=effective_batch,
        )

        # 2. Decide which pending rows need secondary column queries in this tier.
        sec_needed: list[int] = []
        primary_by_idx: dict[int, MappingResult] = {}
        for idx in pending:
            st = states[idx]
            pkey = _request_key(st.record.name, st.identifiers, cache_opts)
            pres = primary_results[pkey]
            primary_by_idx[idx] = pres
            if pkey not in primary_disk_hits:
                st.reused_from_disk = False
                st.queried_live = True
            if pres.error is not None:
                st.errored = st.errored or pres.error
            else:
                st.observations.append(
                    _Observation(
                        tier_index=tier_idx,
                        tier=tier,
                        column_role="primary",
                        query_text=st.record.name,
                        result=pres,
                    )
                )
            if st.secondary_queries and st.errored is None and (
                tier_idx > 0 or _needs_corroboration_in_tier1(pres, tier, policy)
            ):
                sec_needed.append(idx)

        # 3. Secondary column queries (e.g. TargetFullName), verbatim and deduplicated.
        sec_results: dict[str, MappingResult] = {}
        sec_disk_hits: set[str] = set()
        if sec_needed:
            sec_items = [
                (qtext, states[idx].identifiers)
                for idx in sec_needed
                for _col, qtext in states[idx].secondary_queries
            ]
            sec_results, sec_disk_hits = _execute_subrequests(
                mapper,
                sec_items,
                cache=cache,
                arm_name=arm_name,
                next_batch_index=next_batch,
                entity_type=tier.entity_type,
                annotation_mode=annotation_mode,
                search_options=search_opts,
                batch_size=effective_batch,
            )

        # 4. Evaluate early acceptance vs deferring to subsequent tiers.
        still_pending: list[int] = []
        for idx in pending:
            st = states[idx]
            sec_obs: list[tuple[str, MappingResult]] = []
            if idx in sec_needed:
                for col, qtext in st.secondary_queries:
                    skey = _request_key(qtext, st.identifiers, cache_opts)
                    sres = sec_results[skey]
                    if skey not in sec_disk_hits:
                        st.reused_from_disk = False
                        st.queried_live = True
                    if sres.error is not None:
                        st.errored = st.errored or sres.error
                    else:
                        sec_obs.append((col, sres))
                        st.observations.append(
                            _Observation(
                                tier_index=tier_idx,
                                tier=tier,
                                column_role=col,
                                query_text=qtext,
                                result=sres,
                            )
                        )
            accepted = _try_early_accept(st, tier, primary_by_idx[idx], sec_obs, policy)
            if accepted is not None:
                winner_res, win_tier, win_col = accepted
                st.final_result = winner_res
                st.winner_tier = win_tier
                st.winner_column = win_col
            elif st.errored is None:
                still_pending.append(idx)

        pending = still_pending
        if on_batch is not None:
            on_batch(total_rows - len(pending), total_rows)

    # Adjudicate any rows that did not early-accept.
    tier_counts: Counter[str] = Counter()
    column_counts: Counter[str] = Counter()
    final_results: list[MappingResult] = []
    reused_flags: list[bool] = []
    n_errored = 0

    for st in states:
        if st.final_result is None:
            winner_res, win_tier, win_col = _adjudicate_observations(st, policy)
            st.final_result = winner_res
            st.winner_tier = win_tier
            st.winner_column = win_col

        final = st.final_result.model_copy(update={"query_name": st.record.name})
        if final.error:
            n_errored += 1
        if st.winner_tier is not None:
            tier_counts[st.winner_tier] += 1
        elif final.error:
            tier_counts["errored"] += 1
        else:
            tier_counts["unresolved"] += 1
        if st.winner_column is not None:
            column_counts[st.winner_column] += 1

        final_results.append(final)
        reused_flags.append(st.reused_from_disk and not st.queried_live)

    if on_batch is not None and total_rows:
        on_batch(total_rows, total_rows)

    diagnostics = {
        "tiers_configured": [t.name for t in policy.tiers],
        "corroborate_columns": policy.corroborate_columns,
        "extra_name_columns": list(policy.name_columns),
        "by_winning_tier": dict(sorted(tier_counts.items())),
        "by_winning_column": dict(sorted(column_counts.items())),
        "n_grouped_candidates": sum(1 for r in final_results if is_usable_node(r)),
        "n_human_nodes": sum(1 for r in final_results if is_human_node(r)),
        "n_errored": n_errored,
    }
    return PolicyOutcome(
        results=final_results,
        reused=reused_flags,
        n_errored=n_errored,
        diagnostics=diagnostics,
    )
