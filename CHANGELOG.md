# Changelog

All notable changes to the `biomapper` Python client are recorded here.
This project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
