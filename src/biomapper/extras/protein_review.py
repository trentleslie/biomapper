"""Deterministic weak-result and suspect-row candidate review for SomaScan human proteins.

Uses untrusted cohort context (``Target`` and ``TargetFullName``) together with BioMapper's
engine result and Kestrel candidate metadata (``curie``, ``name``, ``synonyms``, ``categories``,
``score``, ``relation``) to decide whether to keep the engine's choice, override with a better
Kestrel candidate CURIE, or record ``"none"`` when no candidate matches the human target.
"""

from __future__ import annotations

import csv
import re
import unicodedata
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal

from biomapper.harmonize.curies import curie_set, is_structure_curie, normalize_curie
from biomapper.mcp.references import CandidateItem

ChoiceSource = Literal["engine", "candidate", "none"]

_HUMAN_PRIMARY_PREFIXES: frozenset[str] = frozenset({"HGNC"})
_PROTEIN_GENE_PREFIXES: frozenset[str] = frozenset(
    {"HGNC", "UNIPROTKB", "UNIPROT", "NCBIGENE", "PR", "ENSEMBL"}
)
_PROTEIN_GENE_CATEGORIES: frozenset[str] = frozenset(
    {
        "biolink:protein",
        "biolink:gene",
        "biolink:polypeptide",
        "biolink:geneorprotein",
        "biolink:genomicentity",
    }
)
_REJECTED_CATEGORIES: frozenset[str] = frozenset(
    {
        "biolink:smallmolecule",
        "biolink:chemicalentity",
        "biolink:disease",
        "biolink:phenotypicfeature",
        "biolink:organismtaxon",
        "biolink:molecularactivity",
        "biolink:biologicalprocess",
        "biolink:anatomicalentity",
        "biolink:cellularcomponent",
    }
)
_NON_HUMAN_MARKERS: frozenset[str] = frozenset(
    {
        "mouse",
        "mus musculus",
        "rat",
        "rattus",
        "bovine",
        "porcine",
        "canine",
        "feline",
        "equine",
        "gallus",
        "chicken",
        "xenopus",
        "zebrafish",
        "drosophila",
        "yeast",
        "saccharomyces",
        "escherichia",
        "e coli",
        "hiv",
        "hiv-1",
        "hiv-2",
        "hepatitis c",
        "hepatitis b",
        "sars",
        "influenza",
        "dengue",
        "vaccinia",
        "botulinum",
        "tetanus",
        "cholera",
        "diphtheria",
        "pertussis",
        "plasmodium",
        "malaria",
        "staphylococcus",
        "streptococcus",
        "bacillus",
        "anthrax",
    }
)

_GREEK_TO_WORD: dict[str, str] = {
    "\u03b1": "alpha",
    "\u0391": "alpha",
    "\u03b2": "beta",
    "\u0392": "beta",
    "\u03b3": "gamma",
    "\u0393": "gamma",
    "\u03b4": "delta",
    "\u0394": "delta",
    "\u03b5": "epsilon",
    "\u0395": "epsilon",
    "\u03b6": "zeta",
    "\u0396": "zeta",
    "\u03ba": "kappa",
    "\u039a": "kappa",
    "\u03bb": "lambda",
    "\u039b": "lambda",
    "\u03bc": "mu",
    "\u039c": "mu",
}

_GREEK_WORDS: frozenset[str] = frozenset(
    {"alpha", "beta", "gamma", "delta", "epsilon", "zeta", "kappa", "lambda", "mu"}
)
_GREEK_SHORT: dict[str, str] = {
    "a": "alpha",
    "b": "beta",
    "g": "gamma",
    "d": "delta",
    "e": "epsilon",
    "z": "zeta",
    "k": "kappa",
    "l": "lambda",
    "m": "mu",
}

_ROMAN_TO_INT: dict[str, str] = {
    "i": "1",
    "ii": "2",
    "iii": "3",
    "iv": "4",
    "v": "5",
    "vi": "6",
    "vii": "7",
    "viii": "8",
    "ix": "9",
    "x": "10",
    "xi": "11",
    "xii": "12",
}

_STOPWORDS: frozenset[str] = frozenset(
    {
        "a",
        "an",
        "and",
        "at",
        "by",
        "chain",
        "domain",
        "family",
        "for",
        "from",
        "gene",
        "homolog",
        "in",
        "isoform",
        "member",
        "of",
        "or",
        "precursor",
        "protein",
        "soluble",
        "subunit",
        "the",
        "to",
        "type",
        "with",
    }
)

_QUALIFIER_PATTERNS: tuple[re.Pattern[str], ...] = (
    re.compile(r"\[cleaved into[^\]]*\]", re.IGNORECASE),
    re.compile(
        r",\s*(?:mitochondrial|cytoplasmic|nuclear|cytosolic|secreted|extracellular|"
        r"soluble|plasma|liver|bone|cardiac muscle|skeletal muscle|muscle|brain|kidney|"
        r"testis|placenta|erythrocyte|isoform\s+\w+|tissue-nonspecific isozyme|"
        r"membrane[- ]associated)\b.*$",
        re.IGNORECASE,
    ),
    re.compile(r"\s+\((?:fragment|precursor|soluble)\)\s*$", re.IGNORECASE),
)

_ROLE_GROUPS: dict[str, tuple[str, ...]] = {
    "receptor": ("receptor", "receptors"),
    "ligand": ("ligand", "ligands"),
    "inhibitor": ("inhibitor", "inhibitors", "antagonist", "antagonists"),
    "activator": ("activator", "activators", "agonist", "agonists"),
    "binding": ("binding",),
    "converting": ("converting",),
    "related": ("related",),
    "like": ("like",),
    "kinase": ("kinase", "kinases"),
    "phosphatase": ("phosphatase", "phosphatases"),
}

_DIRECT_ROLE_WORDS: frozenset[str] = frozenset(
    {"receptor", "inhibitor", "activator", "binding", "converting", "related", "like"}
)

_CHEM_POS_RE = re.compile(r"\b\d+(?:,\d+)+-\w+")
_TOKEN_RE = re.compile(r"[a-z0-9']+")
_SUFFIX_NUMBER_RE = re.compile(r"([a-z]{2,})(\d+[a-z]?)$")
_NUMERIC_TOKEN_RE = re.compile(r"^(\d+)([a-z]?)$")
_COLON_COMPLEX_RE = re.compile(r"[^:]{3,}:[^:]{3,}")
_DEFAULT_SCORE_REVIEW_CUTOFF = 3.2


@dataclass(frozen=True)
class SuspectProteinRow:
    """One twinsuk protein row selected for supplemental candidate review."""

    seq_id: str
    target: str
    target_full_name: str
    query_name: str
    trigger: str
    chosen_kg_id: str | None
    confidence_score: float | None


@dataclass(frozen=True)
class CandidateAssessment:
    """Evaluation of one Kestrel candidate against ``(Target, TargetFullName)``."""

    curie: str
    name: str
    tier: int
    semantic_score: float
    kestrel_score: float
    is_hgnc: bool
    is_engine: bool
    eligible: bool
    contradicts_target: bool
    reasons: tuple[str, ...] = field(default_factory=tuple)
    conflicts: tuple[str, ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class ReviewDecision:
    """Deterministic review outcome for one reviewed twinsuk target."""

    choice: str
    choice_source: ChoiceSource
    rationale: str
    matched_tier: int
    engine_contradicted: bool
    assessed: tuple[CandidateAssessment, ...] = field(default_factory=tuple)


def load_twinsuk_targets(path: Path) -> dict[str, tuple[str, str]]:
    """Read ``SeqId -> (Target, TargetFullName)`` from a SomaScan targets TSV."""
    out: dict[str, tuple[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            seq_id = str(row.get("SeqId") or "").strip()
            target = str(row.get("Target") or "").strip()
            full_name = str(row.get("TargetFullName") or "").strip()
            if seq_id:
                out[seq_id] = (target, full_name)
    return out


def _replace_greek(text: str) -> str:
    return "".join(_GREEK_TO_WORD.get(ch, ch) for ch in text)


def _strip_qualifiers(text: str) -> str:
    cleaned = text.strip()
    for pat in _QUALIFIER_PATTERNS:
        cleaned = pat.sub("", cleaned).strip()
    return cleaned.rstrip(",; ")


def _normalize_phrase(text: str) -> str:
    nfkc = unicodedata.normalize("NFKC", text or "")
    expanded = _replace_greek(nfkc).lower()
    expanded = re.sub(r"(\w)-\s*(\w)", r"\1 \2", expanded)
    tokens = _TOKEN_RE.findall(expanded)
    return " ".join(tokens)


def _squash_alnum(text: str) -> str:
    nfkc = unicodedata.normalize("NFKC", text or "")
    expanded = _replace_greek(nfkc).lower()
    return re.sub(r"[^a-z0-9]+", "", expanded)


def _phrase_variants(text: str, *, allow_soluble_prefix: bool = False) -> list[str]:
    base = (text or "").strip()
    if not base:
        return []
    candidates = [base]
    stripped = _strip_qualifiers(base)
    if stripped and stripped != base:
        candidates.append(stripped)
    if allow_soluble_prefix and len(base) >= 4 and base.startswith("s") and base[1].isupper():
        candidates.append(base[1:])
    seen: set[str] = set()
    out: list[str] = []
    for raw in candidates:
        norm = _normalize_phrase(raw)
        if norm and norm not in seen:
            seen.add(norm)
            out.append(norm)
    return out


def _informative_tokens(text: str) -> set[str]:
    norm = _normalize_phrase(_strip_qualifiers(text))
    out: set[str] = set()
    for tok in norm.split():
        clean = tok.strip("'")
        if not clean or clean in _STOPWORDS:
            continue
        if clean in _ROMAN_TO_INT and len(norm.split()) > 1:
            clean = _ROMAN_TO_INT[clean]
        out.add(clean)
        m = _SUFFIX_NUMBER_RE.match(clean)
        if m:
            out.add(m.group(1))
            out.add(m.group(2))
    return out


def _role_flags(text: str) -> set[str]:
    norm = _normalize_phrase(text)
    toks = set(norm.split())
    roles: set[str] = set()
    for role, aliases in _ROLE_GROUPS.items():
        if any(a in toks for a in aliases):
            roles.add(role)
    # Recognize common receptor abbreviations such as "IL-11 RA" or "VEGF sR2".
    if re.search(r"\b(?:s?r\d+|ra|rb|rc)\b", norm):
        roles.add("receptor")
    return roles


def _subtype_markers(text: str) -> tuple[set[str], set[str]]:
    """Extract numeric and Greek subtype markers, excluding chemical positions."""
    cleaned = _strip_qualifiers(text)
    cleaned = _CHEM_POS_RE.sub(" ", cleaned)
    cleaned = re.sub(r"\b14-3-3\b", " fourteen_three_three ", cleaned)
    cleaned = re.sub(r"\b\d+['\u2032]-?\w*", " ", cleaned)
    norm = _normalize_phrase(cleaned)
    raw_toks = norm.split()
    nums: set[str] = set()
    greeks: set[str] = set()
    for i, tok in enumerate(raw_toks):
        if tok in _GREEK_WORDS:
            greeks.add(tok)
            continue
        if tok in _ROMAN_TO_INT and i > 0:
            nums.add(_ROMAN_TO_INT[tok])
            continue
        num_m = _NUMERIC_TOKEN_RE.match(tok)
        if num_m and i > 0 and raw_toks[i - 1] not in {"band", "protein", "cd"}:
            nums.add(num_m.group(1))
            continue
        suf_m = _SUFFIX_NUMBER_RE.match(tok)
        if suf_m:
            digits = re.match(r"^(\d+)", suf_m.group(2))
            if digits:
                nums.add(digits.group(1))
    # Expand shorthand "ra"/"rb" at token boundary when "receptor" context exists.
    if re.search(r"\bra\b", norm) and "receptor" in _role_flags(text):
        greeks.add("alpha")
    if re.search(r"\brb\b", norm) and "receptor" in _role_flags(text):
        greeks.add("beta")
    for short_g, full_g in _GREEK_SHORT.items():
        if re.search(rf"\b[a-z]+\s+{short_g}\b", norm) and short_g in ("a", "b", "g", "d"):
            if len(raw_toks) >= 2 and raw_toks[-1] == short_g:
                greeks.add(full_g)
    return nums, greeks


def is_heterocomplex(target: str, target_full_name: str) -> bool:
    """True when SomaScan describes a multi-subunit heteromeric complex (e.g. ``A:B complex``)."""
    return len(_complex_parts(target, target_full_name)) >= 2


def _complex_parts(target: str, target_full_name: str) -> list[str]:
    full = (target_full_name or "").strip()
    short = (target or "").strip()
    if _COLON_COMPLEX_RE.search(full):
        body = re.sub(r"\bcomplex\b.*$", "", full, flags=re.IGNORECASE).strip()
        parts = [p.strip() for p in body.split(":") if p.strip()]
        if len(parts) >= 2:
            return parts
    if "complex" in full.lower():
        # Handle hyphen/slash subunit pairs such as "Complement C5b-C6 complex".
        m = re.search(r"\b(c\d+[a-z]?)[-/,]\s*(c?\d+[a-z]?)\b", full, flags=re.IGNORECASE)
        if m:
            p1, p2 = m.group(1), m.group(2)
            if p2[0].isdigit():
                p2 = p1[0] + p2
            return [p1, p2]
    if "complex" in short.lower() and ("," in short or "/" in short):
        body = re.sub(r"\bcomplex\b.*$", "", short, flags=re.IGNORECASE).strip()
        raw_parts = [p.strip() for p in re.split(r"[,/]+", body) if p.strip()]
        if len(raw_parts) >= 2:
            out: list[str] = [raw_parts[0]]
            lead = re.match(r"^([A-Za-z]+)", raw_parts[0])
            for part in raw_parts[1:]:
                if part[0].isdigit() and lead is not None:
                    out.append(lead.group(1) + part)
                else:
                    out.append(part)
            return out
    return []


def _non_human_reasons(item: CandidateItem) -> list[str]:
    reasons: list[str] = []
    curie = normalize_curie(item.curie) or item.curie
    prefix = curie.split(":", 1)[0]
    if prefix in _HUMAN_PRIMARY_PREFIXES:
        return []
    blob = " ".join([item.name, *item.synonyms]).lower()
    for marker in _NON_HUMAN_MARKERS:
        if re.search(rf"\b{re.escape(marker)}\b", blob):
            reasons.append(f"non_human:{marker}")
            break
    return reasons


def _category_ok(item: CandidateItem) -> bool:
    curie = normalize_curie(item.curie) or item.curie
    prefix = curie.split(":", 1)[0]
    cats = {c.strip().lower() for c in item.categories if c and c.strip()}
    if not cats:
        return prefix in _PROTEIN_GENE_PREFIXES
    if cats & _PROTEIN_GENE_CATEGORIES:
        return True
    if cats & _REJECTED_CATEGORIES:
        return False
    return prefix in _PROTEIN_GENE_PREFIXES


def _has_descriptive_text(phrases: Sequence[str]) -> bool:
    for p in phrases:
        stripped = (p or "").strip()
        if not stripped:
            continue
        words = stripped.split()
        if len(words) >= 2:
            return True
        if stripped.isalpha() and len(stripped) >= 9 and not stripped.isupper():
            return True
    return False


def _symbol_matches_initials(symbol_squash: str, full_name: str) -> bool:
    if len(symbol_squash) < 2:
        return False
    toks = [t for t in _normalize_phrase(full_name).split() if t and t not in _STOPWORDS]
    if not toks:
        return False
    initials = "".join(t[0] for t in toks)
    if initials == symbol_squash or (len(initials) >= 2 and initials in symbol_squash):
        return True
    # Also accept contiguous sub-consonant / stem matches like "thpo" vs "thrombopoietin".
    first = toks[0]
    if len(symbol_squash) >= 3 and len(first) >= 4 and symbol_squash[0] == first[0]:
        pos = 0
        for ch in symbol_squash:
            idx = first.find(ch, pos)
            if idx == -1:
                pos = -1
                break
            pos = idx + 1
        if pos != -1:
            return True
    return False


def _best_phrase_pair_score(
    query_variants: Sequence[str],
    query_tokens: set[str],
    cand_phrases: Sequence[str],
) -> tuple[bool, float]:
    exact = False
    best_overlap = 0.0
    q_squash = {_squash_alnum(q) for q in query_variants if _squash_alnum(q)}
    for cp in cand_phrases:
        for cv in _phrase_variants(cp):
            if cv in query_variants:
                exact = True
                best_overlap = 1.0
            c_sq = _squash_alnum(cv)
            if c_sq and len(c_sq) >= 2 and c_sq in q_squash:
                exact = True
                best_overlap = 1.0
            c_toks = _informative_tokens(cv)
            if query_tokens and c_toks:
                shared = len(query_tokens & c_toks)
                if shared:
                    contain_q = shared / len(query_tokens)
                    contain_c = shared / len(c_toks)
                    jaccard = shared / len(query_tokens | c_toks)
                    sim = 0.55 * min(contain_q, contain_c) + 0.45 * jaccard
                    if contain_q == 1.0 and contain_c >= 0.7:
                        sim = max(sim, 0.90)
                    best_overlap = max(best_overlap, sim)
    return exact, best_overlap


def assess_candidate(
    item: CandidateItem,
    *,
    target: str,
    target_full_name: str,
    engine_curies: set[str],
) -> CandidateAssessment:
    """Assess how cleanly one Kestrel candidate matches ``(Target, TargetFullName)``."""
    norm_curie = normalize_curie(item.curie) or item.curie
    prefix = norm_curie.split(":", 1)[0]
    is_hgnc = prefix in _HUMAN_PRIMARY_PREFIXES
    is_engine = (item.relation in ("equals", "equivalent")) or (norm_curie in engine_curies)
    kestrel_score = float(item.score or 0.0)

    cand_phrases = [p for p in [item.name, *item.synonyms] if p and p.strip()]
    cand_desc_phrases = [p for p in cand_phrases if _has_descriptive_text([p])]
    full_variants = _phrase_variants(target_full_name)
    full_tokens = _informative_tokens(target_full_name)
    short_variants = _phrase_variants(target, allow_soluble_prefix=True)
    short_tokens = _informative_tokens(target)
    short_squash = _squash_alnum(target)

    reasons: list[str] = []
    conflicts: list[str] = []

    if is_structure_curie(norm_curie):
        conflicts.append("structure_curie")
    if not _category_ok(item):
        conflicts.append("non_protein_category")
    conflicts.extend(_non_human_reasons(item))

    # Role symmetry check across full target description vs candidate labels.
    t_roles = _role_flags(f"{target} {target_full_name}")
    c_blob = " ".join(cand_phrases)
    c_roles = _role_flags(c_blob)
    for role in sorted((t_roles ^ c_roles) & _DIRECT_ROLE_WORDS):
        if role in t_roles:
            # Allow HGNC symbol suffix like IL11RA or PLAUR, and do not penalize sparse
            # symbol-only candidates lacking descriptive synonyms when the target symbol matches.
            if role == "receptor" and any(
                re.search(r"(?:ra|rb|rc|r\d*|rec)$", _squash_alnum(p)) for p in cand_phrases
            ):
                continue
            exact_short_only = any(_squash_alnum(p) == short_squash for p in cand_phrases if short_squash)
            if cand_desc_phrases or (role == "receptor" and not exact_short_only):
                conflicts.append(f"missing_role:{role}")
        else:
            conflicts.append(f"extra_role:{role}")

    # Subtype number / Greek letter paralog check when descriptive text is present.
    t_nums, t_greeks = _subtype_markers(f"{target} {target_full_name}")
    c_nums, c_greeks = _subtype_markers(c_blob)
    if t_nums and c_nums and not (t_nums & c_nums):
        conflicts.append(f"subtype_number_mismatch:{sorted(t_nums)}!={sorted(c_nums)}")
    if t_greeks and c_greeks and not (t_greeks & c_greeks):
        conflicts.append(f"subtype_greek_mismatch:{sorted(t_greeks)}!={sorted(c_greeks)}")

    # Heteromeric complex check: reject candidates that fail to cover all complex subunits.
    complex_parts = _complex_parts(target, target_full_name)
    if complex_parts and len(complex_parts) >= 2:
        covered_parts = 0
        cand_squashed = {_squash_alnum(cp) for cp in cand_phrases if _squash_alnum(cp)}
        for part in complex_parts:
            p_vars = _phrase_variants(part)
            p_toks = _informative_tokens(part)
            p_sq = _squash_alnum(part)
            p_exact, p_ov = _best_phrase_pair_score(p_vars, p_toks, cand_phrases)
            if p_exact or p_ov >= 0.72 or (p_sq and p_sq in cand_squashed):
                covered_parts += 1
        if covered_parts < len(complex_parts):
            conflicts.append("partial_complex_only")

    exact_full, full_sim = _best_phrase_pair_score(full_variants, full_tokens, cand_phrases)
    exact_short, short_sim = _best_phrase_pair_score(short_variants, short_tokens, cand_phrases)

    # Homonym collision check: e.g. "TF" -> "Transferrin" when TargetFullName is "Tissue Factor".
    desc_overlap = 0.0
    if cand_desc_phrases and full_tokens:
        _, desc_overlap = _best_phrase_pair_score(full_variants, full_tokens, cand_desc_phrases)
        if desc_overlap == 0.0 and not exact_full:
            conflicts.append("homonym_contradicts_full_name")

    # Symbol-to-initials fallback check when Kestrel candidate only carries symbol name.
    initials_hit = False
    if not exact_full and item.name and len(cand_phrases) == 1 and not cand_desc_phrases:
        c_sq = _squash_alnum(item.name)
        if _symbol_matches_initials(c_sq, target_full_name) or (
            short_squash and len(short_squash) >= 2 and c_sq == short_squash
        ):
            initials_hit = True

    tier = 0
    sem_score = 0.0
    if not conflicts:
        if exact_full:
            tier = 4
            sem_score = 1.0
            reasons.append("exact_full_name")
        elif full_sim >= 0.72:
            tier = 3 if (exact_short or short_sim >= 0.7) else 2
            sem_score = full_sim
            reasons.append("full_name_token_match")
        elif exact_short and (len(short_squash) >= 3 or not cand_desc_phrases or desc_overlap > 0):
            if cand_desc_phrases and desc_overlap > 0:
                tier = 3
                sem_score = max(0.75, desc_overlap)
                reasons.append("symbol_with_full_support")
            elif len(short_squash) >= 3 and not full_tokens:
                tier = 2
                sem_score = 0.70
                reasons.append("symbol_only_no_full")
            elif initials_hit:
                tier = 2
                sem_score = 0.65
                reasons.append("symbol_and_initials")
            elif len(short_squash) >= 4 and not cand_desc_phrases:
                tier = 1
                sem_score = 0.55
                reasons.append("symbol_only_sparse_metadata")
        elif initials_hit:
            tier = 1
            sem_score = 0.55
            reasons.append("initials_match_sparse_metadata")
    else:
        sem_score = max(0.0, full_sim - 0.5)

    eligible = (not conflicts) and (prefix in _PROTEIN_GENE_PREFIXES) and _category_ok(item)
    contradicts_target = bool(
        {"homonym_contradicts_full_name", "partial_complex_only", "non_protein_category"}
        & set(conflicts)
        or any(
            c.startswith(("non_human:", "missing_role:", "extra_role:", "subtype_"))
            for c in conflicts
        )
    )

    return CandidateAssessment(
        curie=norm_curie,
        name=item.name,
        tier=tier,
        semantic_score=round(sem_score, 4),
        kestrel_score=kestrel_score,
        is_hgnc=is_hgnc,
        is_engine=is_engine,
        eligible=eligible,
        contradicts_target=contradicts_target,
        reasons=tuple(reasons),
        conflicts=tuple(conflicts),
    )


def decide_review(
    *,
    target: str,
    target_full_name: str,
    weak_trigger: str | None,
    engine_choice: str | None,
    kg_equivalent_ids: Mapping[str, Sequence[str]] | None,
    candidates: Sequence[CandidateItem],
    confidence_score: float | None = None,
) -> ReviewDecision:
    """Decide ``engine``, overriding ``candidate``, or ``"none"`` for one twinsuk row."""
    raw_engine = (engine_choice or "").strip() or None
    norm_engine = normalize_curie(raw_engine) if raw_engine else None
    engine_curies = set(curie_set(raw_engine, dict(kg_equivalent_ids or {})))
    if norm_engine:
        engine_curies.add(norm_engine)

    assessments = tuple(
        assess_candidate(
            c,
            target=target,
            target_full_name=target_full_name,
            engine_curies=engine_curies,
        )
        for c in candidates
    )

    engine_assessments = [a for a in assessments if a.is_engine]
    best_engine = max(
        engine_assessments,
        key=lambda a: (a.tier, int(a.is_hgnc), a.semantic_score, a.kestrel_score),
        default=None,
    )

    valid_candidates = sorted(
        (a for a in assessments if a.eligible and a.tier >= 1),
        key=lambda a: (a.tier, int(a.is_hgnc), a.semantic_score, a.kestrel_score),
        reverse=True,
    )
    best_overall = valid_candidates[0] if valid_candidates else None
    strong_candidates = [a for a in valid_candidates if a.tier >= 2]
    best_strong = strong_candidates[0] if strong_candidates else None

    engine_is_missing_or_invalid = (
        raw_engine is None
        or is_structure_curie(raw_engine)
        or (weak_trigger in ("unresolved", "refused", "no_node"))
    )
    engine_contradicted = engine_is_missing_or_invalid or any(
        a.contradicts_target for a in engine_assessments
    )

    # Reject completely ungrounded low-confidence guesses (< 2.0) when sparse metadata shows zero overlap.
    if (
        not engine_is_missing_or_invalid
        and not engine_contradicted
        and confidence_score is not None
        and confidence_score < 2.0
        and best_engine is not None
        and best_engine.tier == 0
    ):
        engine_contradicted = True

    # Case 1: Engine chose something missing, non-node, or explicitly contradicted.
    if engine_is_missing_or_invalid or engine_contradicted:
        if best_strong is not None:
            if best_strong.is_engine and norm_engine and best_strong.is_hgnc:
                return ReviewDecision(
                    choice=raw_engine or best_strong.curie,
                    choice_source="engine",
                    rationale=(
                        f"Kept engine {raw_engine} (equivalent candidate {best_strong.curie} "
                        f"verified: {','.join(best_strong.reasons)})"
                    ),
                    matched_tier=best_strong.tier,
                    engine_contradicted=False,
                    assessed=assessments,
                )
            return ReviewDecision(
                choice=best_strong.curie,
                choice_source="candidate",
                rationale=(
                    f"Overrode engine ({raw_engine or 'unresolved'}, trigger={weak_trigger}) with "
                    f"{best_strong.curie} ({best_strong.name}): {','.join(best_strong.reasons)}"
                ),
                matched_tier=best_strong.tier,
                engine_contradicted=engine_contradicted,
                assessed=assessments,
            )
        why = (
            "; ".join(best_engine.conflicts)
            if best_engine and best_engine.conflicts
            else (weak_trigger or "no_matching_human_candidate")
        )
        return ReviewDecision(
            choice="none",
            choice_source="none",
            rationale=f"No Kestrel candidate verified TargetFullName={target_full_name!r} ({why})",
            matched_tier=0,
            engine_contradicted=engine_contradicted,
            assessed=assessments,
        )

    # Case 2: Engine has a valid node, check if another HGNC candidate strictly beats it cleanly.
    if (
        best_strong is not None
        and not best_strong.is_engine
        and (
            best_engine is None
            or best_strong.tier > best_engine.tier + (0 if best_engine.tier <= 2 else 1)
            or (
                best_strong.tier >= 3
                and best_engine.tier <= 2
                and best_strong.semantic_score > best_engine.semantic_score + 0.15
            )
        )
    ):
        return ReviewDecision(
            choice=best_strong.curie,
            choice_source="candidate",
            rationale=(
                f"Overrode engine {raw_engine} with stronger Kestrel match "
                f"{best_strong.curie} ({best_strong.name}): {','.join(best_strong.reasons)}"
            ),
            matched_tier=best_strong.tier,
            engine_contradicted=False,
            assessed=assessments,
        )

    # Case 3: Keep engine if HGNC (or upgrade equivalent non-HGNC CURIE to an equivalent HGNC candidate).
    if raw_engine is not None and norm_engine is not None:
        if not norm_engine.startswith("HGNC:"):
            equiv_hgnc = next(
                (a for a in valid_candidates if a.is_engine and a.is_hgnc),
                None,
            )
            if equiv_hgnc is not None:
                return ReviewDecision(
                    choice=equiv_hgnc.curie,
                    choice_source="candidate",
                    rationale=(
                        f"Canonicalized equivalent non-HGNC engine choice {raw_engine} "
                        f"to {equiv_hgnc.curie} ({equiv_hgnc.name})"
                    ),
                    matched_tier=equiv_hgnc.tier,
                    engine_contradicted=False,
                    assessed=assessments,
                )
        tier = best_engine.tier if best_engine else (best_overall.tier if best_overall else 0)
        return ReviewDecision(
            choice=raw_engine,
            choice_source="engine",
            rationale=(
                f"Verified engine choice {raw_engine} for {target!r} / {target_full_name!r} "
                f"(trigger={weak_trigger or 'suspect'}, score={confidence_score})"
            ),
            matched_tier=tier,
            engine_contradicted=False,
            assessed=assessments,
        )

    return ReviewDecision(
        choice="none",
        choice_source="none",
        rationale=f"Unresolved row with no matching human protein candidate for {target_full_name!r}",
        matched_tier=0,
        engine_contradicted=True,
        assessed=assessments,
    )


def find_suspect_twinsuk_rows(
    twinsuk_rows: Sequence[Mapping[str, Any]],
    targets_by_seqid: Mapping[str, tuple[str, str]],
    *,
    already_weak_names: set[str] | frozenset[str] = frozenset(),
    score_cutoff: float = _DEFAULT_SCORE_REVIEW_CUTOFF,
) -> list[SuspectProteinRow]:
    """Identify resolved twinsuk rows not flagged by ``weak_trigger`` that warrant review.

    Detects:
    - ``node_collision``: distinct ``(Target, TargetFullName)`` pairs sharing one ``chosen_kg_id``;
    - ``non_hgnc``: resolved CURIE without an ``HGNC:`` primary ID;
    - ``heterocomplex``: multi-subunit complex descriptions that resolved to a single node;
    - ``medium_confidence``: resolved rows scoring below ``score_cutoff``.
    """
    node_pairs: dict[str, set[tuple[str, str]]] = {}
    for row in twinsuk_rows:
        if row.get("error"):
            continue
        chosen = normalize_curie(row.get("chosen_kg_id"))
        if not chosen or is_structure_curie(chosen):
            continue
        seq_id = str(row.get("key") or "").strip()
        t_pair = targets_by_seqid.get(seq_id)
        if t_pair is not None:
            node_pairs.setdefault(chosen, set()).add((t_pair[0].lower(), t_pair[1].lower()))

    colliding_nodes = {node for node, pairs in node_pairs.items() if len(pairs) > 1}
    suspects: list[SuspectProteinRow] = []

    for row in twinsuk_rows:
        if row.get("error"):
            continue
        seq_id = str(row.get("key") or "").strip()
        query_name = str(row.get("name") or "").strip()
        if not seq_id or not query_name or query_name in already_weak_names:
            continue
        raw_chosen = row.get("chosen_kg_id")
        chosen = normalize_curie(raw_chosen)
        if not chosen or is_structure_curie(chosen):
            continue
        target, full_name = targets_by_seqid.get(seq_id, ("", query_name))
        score_raw = row.get("confidence_score")
        score = float(score_raw) if isinstance(score_raw, int | float) else None

        reasons: list[str] = []
        if chosen in colliding_nodes:
            reasons.append("node_collision")
        if not chosen.startswith("HGNC:"):
            reasons.append("non_hgnc")
        if is_heterocomplex(target, full_name):
            reasons.append("heterocomplex")
        if score is not None and score < score_cutoff:
            reasons.append("medium_confidence")
        if reasons:
            suspects.append(
                SuspectProteinRow(
                    seq_id=seq_id,
                    target=target,
                    target_full_name=full_name,
                    query_name=query_name,
                    trigger=reasons[0],
                    chosen_kg_id=str(raw_chosen),
                    confidence_score=score,
                )
            )
    return suspects
