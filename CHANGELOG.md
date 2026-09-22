# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
uses [semantic versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-09-22

First release. Extracts a temporal-clustering method from a published MCF7
ATAC-seq analysis into a reusable, validated Snakemake workflow.

### Added

- Snakemake workflow from counts matrix to HTML report, fanning out per
  treatment arm, with region-specific rules that activate on coordinates and
  are skipped silently otherwise.
- Multi-arm timecourses sharing one unstimulated baseline, expressed by a
  sentinel value in the samplesheet's `group` column.
- Two-gate dynamic selection: a significance test plus an effect-size floor on
  the range of per-timepoint means.
- Three differential methods: per-arm LRT, a union of per-timepoint Wald
  contrasts, and a model-free variance ranking for single-replicate designs.
- Automated trajectory-class assignment, replacing the hand-written cluster
  map in the source analysis, with a dendrogram and a silhouette diagnostic
  emitted whatever the method.
- An optional second stage that splits one broad class and names the pieces
  from their centroids.
- Six figures per arm, one HTML report, and a `run_manifest.json` recording
  the resolved config, package versions, every seed, `sessionInfo()` and input
  checksums.
- Every threshold in `config/config.yaml`, validated by a JSON Schema at parse
  time. Presets for the MCF7 and T cell parameter sets.
- Three test tiers: static checks, an end-to-end run with golden tables, and a
  correctness tier that recovers known trajectory shapes from simulated
  counts. Plus an opt-in fidelity check against a published assignment.
- A bundled demo derived from a published CD8+ T cell ATAC timecourse
  (GEO GSE228171).

### Changed, relative to the source analysis

- `lrt_lfc` is renamed `min_range`. It is a range on the transformed scale and
  never was a log2 fold change; the old name causes it to be set an order of
  magnitude wrong.
- The supercluster map is derived rather than typed. `superclusters.overrides`
  reproduces a historical hand assignment when needed.
- The default assignment labels each cluster by the shape of its own centroid
  rather than cutting the dendrogram. Cutting at k=3 on the bundled demo
  produces a class holding one Decreasing, one Increasing and one Transient
  cluster together.
- The stage-2 labeller's hard-coded timepoints are derived from the grid. The
  mid timepoint is the observed time nearest a fraction of the last, which
  reproduces the source's value exactly on its own grid and, unlike the
  middle-index reading, transfers to unevenly spaced grids.
- `cluster.max_features` defaults to 6000 rather than 30000. `degPatterns`
  cost is cubic, not quadratic: measured 72 s at 3,251 features, extrapolating
  to roughly 18 hours at 30,000.
- Features dropped by `minc` are reported as `Unassigned` rather than
  disappearing from the output.
- Every stochastic step is seeded, including `degPatterns`, which the source
  did not seed despite its methods section saying otherwise.

### Known limitations

- `cluster.degpatterns.cutoff` has no effect in any released DEGreport
  version. The key is kept so configs written against the original notebooks
  parse, and chromadyn warns when it is set. See `docs/method.md`.
- `degPatterns` chooses its own number of clusters and merges shapes whose
  z-scored profiles correlate strongly. `n_clusters` overrides the cut, but
  `minc` is applied afterwards and can reduce the result below what you asked
  for.
- `vignettes/public-geo-timecourse.Rmd` has a verified download but its
  analysis has not been run end to end here, and is marked as such.
- Annotation packages are not bundled; `annotate.txdb` requires you to install
  the one your genome needs.
