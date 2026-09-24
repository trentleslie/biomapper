"""biomapper — Python client for the BioMapper2 API.

Quick start::

    from biomapper import map_entity, map_entities, BioMapperClient

    # Single lookup (synchronous)
    result = map_entity("L-Histidine")
    print(result.primary_curie)      # RM:0129894
    print(result.confidence_tier)    # high

    # Batch (synchronous, with progress bar)
    results = map_entities(
        [{"name": "L-Histidine"}, {"name": "Glucose"}],
        progress=True,
    )

    # File-based (synchronous, with progress bar)
    from pathlib import Path
    from biomapper import map_dataset_file_sync

    result = map_dataset_file_sync(
        Path("compounds.tsv"),
        name_column="name",
        provided_id_columns=["hmdb_id"],
        progress=True,
    )
    result.raise_for_error()  # opt-in: raise if the stream truncated
    print(f"resolved {sum(1 for r in result.results if r.resolved)}")

    # Async (in an async context)
    async with BioMapperClient() as client:
        result = await client.map_entity("L-Histidine")

    # Harmonize two already-resolved cohorts (local, offline, no extra requests).
    # `harmonize` is deliberately NOT re-exported at the package root: the name would
    # shadow the `biomapper.harmonize` submodule and break `biomapper.harmonize.curie_set`.
    from biomapper.harmonize import harmonize

    report = harmonize(ukbb_results, arivale_results, a_label="ukbb", b_label="arivale")
    print(report.n_links, report.a_unresolved)
"""

from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _dist_version

from biomapper.client import BioMapperClient
from biomapper.dataset import map_dataset_file_sync
from biomapper.exceptions import (
    BioMapperAuthError,
    BioMapperConfigError,
    BioMapperError,
    BioMapperRateLimitError,
    BioMapperServerError,
    BioMapperTimeoutError,
)
from biomapper.mapper import (
    list_annotators,
    list_entity_types,
    list_vocabularies,
    map_entities,
    map_entity,
    summarize,
)
from biomapper.models import (
    AnnotatorInfo,
    DatasetMappingResult,
    EntityTypeInfo,
    KestrelRequestParams,
    KestrelSearchResult,
    LipidResolution,
    MappingResult,
    MappingSummary,
    ResolutionCertificate,
    VocabularyInfo,
)

# Single-sourced from the installed distribution metadata, which Poetry builds from
# ``pyproject.toml``. It is deliberately NOT a literal here: a second literal is a second source of
# truth, and the two drifted (pyproject 1.5.1 against a hardcoded 1.4.0) for long enough that no run
# manifest could name its own package version unambiguously. ``pyproject.toml`` is the one source;
# ``tests/test_version.py`` fails the build if this module ever reintroduces a literal.
#
# The fallback fires only for a source tree with no installed metadata at all (a bare checkout that
# was never installed). It is a loud sentinel rather than a plausible-looking number, because a
# provenance field that reads "1.4.0" when nothing is installed is worse than one that reads
# unknown: only the second tells a reader the value cannot be trusted.
try:  # pragma: no cover - exercised by tests/test_version.py in both branches
    __version__ = _dist_version("biomapper")
except PackageNotFoundError:
    __version__ = "0.0.0+unknown"

__all__ = [
    # Client
    "BioMapperClient",
    # Sync helpers
    "map_entity",
    "map_entities",
    "map_dataset_file_sync",
    "list_entity_types",
    "list_annotators",
    "list_vocabularies",
    "summarize",
    # Models
    "MappingResult",
    "MappingSummary",
    "DatasetMappingResult",
    "EntityTypeInfo",
    "AnnotatorInfo",
    "VocabularyInfo",
    "ResolutionCertificate",
    "LipidResolution",
    "KestrelSearchResult",
    "KestrelRequestParams",
    # Exceptions
    "BioMapperError",
    "BioMapperAuthError",
    "BioMapperConfigError",
    "BioMapperRateLimitError",
    "BioMapperServerError",
    "BioMapperTimeoutError",
]
