# Changelog

All notable changes to the `biomapper` Python client are recorded here.
This project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.5.0] - 2026-09-24

### Added

- **`biomapper.benchmarks.cross_cohort` and `cross_cohort_certify` - the Monti/NECS cross-cohort
  arm, finished rather than stubbed.** Phase 2b migrated the 11 external benchmark suite arms but
  not the cross-cohort machinery, so the arm still had to be run from an engine checkout, which
  would have pinned its numbers to a different backend than the suite's. These modules close that
  gap: the whole arm now runs through the package against a deployment, so one provenance story
  covers the suite, the cohort analyses, and the cross-cohort comparison.

  Ported from the biomapper2 engine at `origin/dev`, commit
  `1ffb571e54fe028ef0ae4e748fc2e7ec093ee603`:

  - `benchmarks.adapters.cohort_panel` - the arivale / xuetal / llfs / blsa panels, with the
    exclusion accounting (blank name, unnamed vendor feature, de-duplication) and the
    `certifiable` determination unchanged.
  - `benchmarks.scorers.cross_cohort_overlap` - the cohort-shaped linker. It deliberately
    EXCLUDES `INCHIKEY`, `INCHI` and `SMILES` from the linker: linking on a structure hash would
    make the downstream structural certificate circular and precision 100% by construction. The
    comments explaining that are preserved verbatim.
  - `benchmarks.scorers.link_certificate` and
    `benchmarks.scorers.independent_link_certificate_overlap` - the KG-independent link
    certificate, which refuses rather than certifying off a KG-derived structure.
  - `benchmarks.scorers.arm_b_baseline` - Monti's own published method, reconstructed on the
    identical row set.

  The driver is rewired from `biomapper2.mapper.Mapper` to
  `benchmarks.api_mapper.ApiMapper`, exactly as the suite arms were, so provenance names the graph
  that served the answers instead of a client checkout.

### Changed

- **The Monti published-overlap table is corrected.** `arm_b_baseline.MONTI_PUBLISHED` read
  `{"arivale": 615, "xuetal": 432, "llfs": 163, "blsa": 99}` and attributed all four values to
  "Monti Table 2". Re-read against the paper (DOI `10.1007/s11357-026-02174-2`), the citation and
  two of the four values were wrong:

  - Table 2 is "Age-only markers", not an overlap table. The overlaps are in the Methods section
    "Datasets harmonization" and in the per-cohort descriptions.
  - NECS to Xu is **385** in the harmonization methods. The Xu cohort description separately says
    432, so the paper contradicts itself on this one pair. 385 is used because it comes from the
    sentence describing the procedure Arm B reconstructs, and that sentence's Arivale value agrees
    with the Arivale cohort description.
  - NECS to BLSA is **188**. The old 99 was the BLSA to LLFS overlap ("Ninety-nine metabolites
    were in common with the LLFS"), a different pair entirely.

  Both the corrected and the superseded tables ship, as `MONTI_PUBLISHED` and
  `MONTI_PUBLISHED_SUPERSEDED`, with a per-pair quote and section in
  `MONTI_PUBLISHED_PROVENANCE`, so a changed number is always traceable to a reason.

### Fixed

- **An errored row is no longer indistinguishable from an unresolved one in the cross-cohort
  path.** Both come back with an empty CURIE set, so folding them together reports a
  non-resolution that never happened and understates resolution by however many rows the
  deployment dropped. `cross_cohort.errored_names` counts them apart, `repair_errored_rows`
  re-maps only the failed rows and writes them back into their original positions, and any row
  that still errors is reported in the manifest with the affected panel's counts labelled a floor.
  The deployment did return 5xx under concurrent load during this work, so this is a live hazard.

### Fixed (review round 2)

- **Checkpoints now carry the backend that answered them.** Panels resolve as separate processes and
  are combined by `--link-only`, which previously fetched `/health` once at finalization and stamped
  it on every checkpoint. A deployment or graph-build change between panel runs would then present a
  mixed-backend result as one pinned run. Each panel writes a provenance sidecar at resolve time, and
  finalization compares `endpoint`, `kestrel_version`, `kg_version`, `biolink_version`,
  `build_timestamp` and `git_commit` per panel. A mismatch, or a checkpoint with no sidecar at all,
  raises `BackendDriftError`. `--allow-unpinned-checkpoints` is an explicit escape that records the
  gap in the manifest instead of letting the run read as pinned.
- **A missing per-pair link artifact no longer reads as zero links.** `manifest.json` existing does
  not imply the artifacts do. `read_links` raises `MissingLinkArtifactError`, and each link file's
  row count is cross-checked against the manifest's `arm_m_links`.
- **A transient PubChem failure is no longer overwritten by a clean fallback miss.** A `lookup_failed`
  on the CID route followed by a `clean_miss` on the HMDB route previously reported a retryable run
  artifact as a genuine coverage gap. The failure is sticky unless a fallback actually resolves.

### Added (review round 2)

- **`biomapper.benchmarks.cross_cohort_readjudicate`.** Every non-certified case is checked against
  PubChem's name index, which is a different lookup than either side of the certificate used. A
  first-block mismatch is classified as a tautomer, charge or salt artifact only when the outside
  source corroborates BOTH sides and their heavy-atom compositions match with a mass gap consistent
  with the hydrogen-count difference. Outcomes separate `necs_gold_suspect` from `cohort_id_suspect`
  from `genuine_structural_disagreement` from `outside_source_unresolved`, so a refusal that could
  not be checked is never folded into one that was.

### Notes

- Cohorts that ship names only (NECS, Xu, LLFS, BLSA have no vendor identifier column) are
  `certifiable=False`. Their links are countable but never structurally certifiable, and
  `cross_cohort_certify` reports that as refused **by construction**, not as an uncertified
  failure. Only Arivale carries structure-resolvable vendor identifiers in this arm.
- Every certificate verdict here is an InChIKey **first block** comparison, which is neither
  tautomer- nor charge-invariant: it over-flags tautomers as refuted and silently accepts a
  stereoisomer error whenever a side lacks the stereo layer. `stereo_checked` records which
  happened. The NECS curated gold also carries roughly 5% InChIKey errors, so a refuted verdict is
  a candidate for re-adjudication against an outside source, not a finding on its own.

## [1.4.0] - 2026-09-23

### Added

- **`biomapper.harmonize` — local, client-side cross-dataset harmonization.**
  Two entities from two cohorts are equivalent when they resolve to the same canonical KRAKEN
  node. The linker is an identifier-set intersection, never string matching, and it runs entirely
  on the client over results you already have. It issues no requests and needs no knowledge-graph
  access, so it is fully testable offline.

  - `harmonize(a_results, b_results)` -> `HarmonizationResult` over two lists of `MappingResult`.
  - `link_by_intersection(a_curies, b_curies)` -> `OverlapResult` for raw CURIE-set inputs.
  - `curie_set()`, `normalize_curie()`, `canonical_prefix()`, `predicted_curies()`.

  Two invariants are carried over from the engine's own linking code, along with the comments
  that explain them:

  1. **The linker is identifier-only.** `INCHIKEY`, `INCHI` and `SMILES` are excluded. Linking on
     a structure hash would make any downstream structural certificate circular and would make
     precision 100% by construction.
  2. **CURIE prefix synonyms normalize.** `KEGG.COMPOUND:C00031` equals `KEGG:C00031`, while
     genuinely different identifier spaces such as `KEGG.GLYCAN` stay distinct.

  An entity that resolved to no identifier is a **refusal candidate, not a link**. It is reported
  by name in `a_unresolved` / `b_unresolved` and counted in `summary()`, never silently dropped,
  and it is excluded from the link-rate denominator so non-resolution is never scored as
  non-equivalence. An entity whose mapping call errored is counted separately again, in
  `a_errors` / `b_errors`: "we do not know" is a different claim from "it did not resolve".

  Keys are derived once per input row, against its position in the **original** list, so an
  index-based custom `key` cannot hand the same string to an errored row and a later resolved
  one. Uniqueness is enforced across errored and resolved rows together. Cohort labels key
  `summary()`, so labels that are equal to each other or to a reserved field are rejected rather
  than allowed to silently overwrite a sibling entry.

  `harmonize` is deliberately **not** re-exported at the package root. Binding that name on
  `biomapper` would shadow the `biomapper.harmonize` submodule and break
  `biomapper.harmonize.curie_set` for anyone who reaches for it the obvious way. Import it as
  `from biomapper.harmonize import harmonize`.

- **Mapping options the client could not previously send.** `map_entity`, `map_entities` (async
  and sync) and `map_dataset_file_iter` now accept `vocab`, `prefer_human`, `prefer_canonical`,
  `candidate_limit` and `kestrel_top_n`; the entity and batch paths also accept
  `array_delimiters` (the dataset routes do not take it). Any option left as `None` is omitted
  from the request, so the server's own default applies and the payload is byte-unchanged for
  callers who do not use them. `candidate_limit` and `kestrel_top_n` are validated against the
  API's 1..100 bound locally, turning a wasted round trip into an immediate error.

- **Response fields the client was discarding.** `MappingResult` and `RawApiResult` now carry the
  full `resolution_certificate` (as `MappingResult.certificate`), `lipid_resolution`,
  `chosen_kg_id_lipid_hint`, `refmet_availability`, `refmet_source`, `refmet_snapshot_version`,
  `tier_b_snapshot_version` and the opt-in `kestrel_results` passthrough rows. New models
  `ResolutionCertificate`, `LipidResolution`, `KestrelSearchResult` and `KestrelRequestParams`
  are exported from the package root. `MappingResult.refusal_reason` reads the certificate's
  refusal reason, which is distinct from `error`: an error means the call failed, a refusal means
  the pipeline ran and declined to assert an answer.

- `identifiers` now accepts a list per vocabulary (`{"KEGG": ["C00031", "C00267"]}`), matching
  the API, which types it as `dict[str, str | list[str]]`.

### Fixed

- **`map_entities` failed against the current API.** `BatchMappingResponse.summary` was typed
  `dict[str, int]`, but the API returns a nested `refmet_source_counts` map alongside the scalar
  tallies (it types the field as `{str: int | {str: int}}`). Validation rejected the whole batch
  response, and the batch loop's chunk-level error handling then converted that into a per-entity
  error for **every record in the chunk**. The repo's own test fixture omitted
  `refmet_source_counts`, which is why the suite stayed green while live calls did not.

### Notes

- **Not added: `POST /api/v1/map/dataset`** (the non-streaming dataset route). It returns a
  path on the server's filesystem rather than per-entity JSON, so a remote client cannot read its
  output. `map_dataset_file_iter` / `map_dataset_file_sync` (the streaming route) remain the way
  to map a file. Recorded here so the omission is a decision rather than an oversight.
- No new runtime dependencies. `biomapper.harmonize` is pure Python and does **not** require
  pandas, which is an optional extra of this package.

## [1.3.0]

- Surface `chosen_kg_id_review` from the mapping API.
- Support the new `EntityType[]` response shape with `defaultPrefixes`.

## [1.1.0]

- Add `kg_equivalent_ids` support.

## [1.0.0]

- Initial release of `biomapper` following the rename from `ddharmon`.
