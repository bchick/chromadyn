# timecourse-patterns

**A Snakemake workflow for finding the features that change over a
timecourse and grouping them by the shape of their change.** 

![Trajectory classes recovered from an ATAC-seq timecourse of antiviral CD8+ T cells](docs/img/supercluster_ribbon.png)

*Trajectory classes the workflow recovers from ATAC-seq of antiviral CD8+ T
cells in mice infected with LCMV Armstrong, sampled from naive (day 0) to day
8 post infection. Each panel is one class: the line is the class mean
accessibility (z-score) and the bands show its spread across peaks. Data from
McDonald, Chick et al., Immunity 2023 ([details](#demo-data)).*

## Pipeline

![The workflow as a metro map: validate, prefilter, transform and QC over all samples; then, per treatment arm, a DESeq2 test over time, two-gate selection, clustering, trajectory classes and figures; then a cross-arm comparison, the report and the manifest](docs/img/pipeline.svg)

1. **Validate** the samplesheet against the counts: exact ID match, numeric
   timepoints, enough of them per arm, non-negative integer counts.
2. **Prefilter** features too sparse to have been measured.
3. **Transform** with DESeq2's variance-stabilizing transform (or rlog), blind
   to the design. **QC**: PCA and sample and replicate correlations.
4. **Test over time**, per treatment arm, with a DESeq2 likelihood-ratio test
   (or per-timepoint Wald contrasts). Single-replicate designs rank by
   variance instead.
5. **Select dynamic features** that pass both a significance gate and a floor
   on the range of their per-timepoint means.
6. **Cluster** their z-scored trajectories with DEGreport `degPatterns`, as in
   the source analysis, or with k-means.
7. **Name trajectory classes** (Decreasing, Transient, Late Increasing, ...)
   from each cluster's mean profile, with an optional second stage that
   splits one broad class.
8. **Report**: six figures per arm, a cross-arm comparison, BED files per
   class and optional genomic annotation in region mode, an HTML report and
   `run_manifest.json`.

## Quickstart

```bash
git clone https://github.com/bchick/timecourse-patterns
cd timecourse-patterns
pixi install
pixi run demo          # about 2 minutes, writes results/report.html
```

That runs the bundled demo, and the figure above is its actual output.

### Demo data

The demo is drawn from a published ATAC-seq timecourse of CD8+ T cells
responding to acute viral infection. Mice were infected with lymphocytic
choriomeningitis virus (LCMV) Armstrong, a strain the immune system clears
within about a week. Accessibility was profiled in naive CD8+ T cells (day 0),
in antiviral CD8+ T cells at days 3 and 5 post infection, and in sorted
terminal effector cells at day 8, the peak of the response.

It is a subset of the full samples, so that it runs in about two minutes:

- **Samples**: 9 of the study's 62 ATAC libraries: two wild-type replicates
  each at days 0, 3 and 5, and three at day 8. A 48-hour timepoint was left
  out because it failed QC.
- **Peaks**: 5,000 of the 129,076 consensus peaks on standard chromosomes,
  sampled so that every trajectory shape is represented, along with static
  and low-count peaks, so every step of the pipeline has real work to do.
- **Counts**: unchanged. Every retained peak keeps its full read count from
  each library; nothing is downsampled.

The [T cell example](#example-trajectory-classes-carry-biology) below runs the
same nine libraries on every peak.

> McDonald BD\*, Chick BY\*, Ahmed NU, et al. Canonical BAF complex activity
> shapes the enhancer landscape that licenses CD8+ T cell effector and memory
> fates. *Immunity* 56(6):1303-1319.e5 (2023).
> [doi:10.1016/j.immuni.2023.05.005](https://doi.org/10.1016/j.immuni.2023.05.005).
> GEO [GSE228381](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228381)
> (ATAC sub-series [GSE228171](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228171)).

Which libraries were used, and why, is documented in
[`demo/PROVENANCE.md`](demo/PROVENANCE.md).

## Scope

The workflow addresses one question: **given counts over a timecourse, which
features change over time, and what distinct shapes do those changes take?**

It starts at a counts matrix and never touches FASTQ, BAM or peak calling, so
it runs downstream of nf-core/atacseq, nf-core/rnaseq, DiffBind,
featureCounts, salmon or any other quantifier. Features may be genomic regions
(ATAC, CUT&RUN, ChIP) or genes (RNA-seq); region-specific steps switch on
when you supply coordinates and are skipped silently when you do not.

The step it adds beyond the tools it wraps is the naming: turning a few dozen
clusters into a handful of named trajectory classes by a stated rule, where
the source analysis typed the mapping in by hand.

## Input

### `samplesheet.tsv`

| Column | Required | Meaning |
|---|---|---|
| `sample` | yes | Unique ID. **Must match a counts column name exactly.** |
| `group` | yes | Treatment arm, for example `DrugA` or `DrugB`. The sentinel value (default `shared`) joins that row to every arm. |
| `time` | yes | Numeric timepoint. Units are declared in the config, never parsed from the value. |
| `replicate` | yes | Replicate label within group and time. |
| `batch` | no | Used only when `batch_correct` is enabled. |
| anything else | no | Carried through to `colData` and available to custom designs. |

The `shared` sentinel is how you express a **multi-arm timecourse with one
unstimulated baseline**: the same t=0 libraries serve as the zero timepoint in
every arm, and each arm is fitted separately against them.

```
sample              group   time  replicate
vehicle_r1          shared  0     r1
vehicle_r2          shared  0     r2
drugA_2h_r1         DrugA   2     r1
drugB_2h_r1         DrugB   2     r1
```

### Counts

A `.tsv`, `.csv`, `.rds` matrix or a serialized `SummarizedExperiment`. Rows are
features, columns are samples, raw integer counts by default. A mismatch
between samplesheet IDs and counts columns is a hard error listing both sides,
never a silent subset.

### `features.bed` (optional)

BED3 or BED6 mapping to the counts row names. Supplying it switches on region
mode: BED export per trajectory class, and genomic annotation if you configure
a TxDb. Without it, timecourse-patterns runs in gene mode.

## Configuration

Everything lives in `config/config.yaml`, which ships with every key present,
its default, and a comment explaining why the key exists. A JSON Schema
validates it at parse time, so a typo fails immediately rather than three rules
later. The keys you are most likely to touch:

| Key | Default | What it controls |
|---|---|---|
| `filter.min_mean_counts` | 10 | Drops features too sparse for the assay to have measured |
| `differential.method` | `lrt` | `lrt`, `wald_union`, or `none` for single-replicate designs |
| `differential.fdr` | 0.01 | Significance gate |
| `differential.min_range` | 0.5 | Effect-size gate, on the **transformed** scale. Not a fold change. |
| `cluster.method` | `degpatterns` | `degpatterns` as in the source analysis, or `kmeans` to cluster every feature without subsampling |
| `cluster.degpatterns.minc` | 50 | Minimum cluster size; `auto` scales it to your feature count |
| `cluster.max_features` | 6000 | Above this, subsample and assign the rest by correlation |
| `superclusters.method` | `shape` | How clusters become named classes |
| `superclusters.split` | `Increasing`, k=3 | Optional second stage that splits one broad class |

Ready-made presets are in `config/presets/`.

## Output

```
results/
├── qc/           validation report, library sizes, PCA, correlations
├── differential/ <arm>_results.tsv.gz  per-feature statistics and the dynamic flag
├── clusters/     <arm>_clusters.tsv    feature to cluster to trajectory class
│                 <arm>_k_diagnostics.tsv, bed/, cross_arm_crosstab.tsv
├── figures/      six panels per arm, in pdf, svg and png
├── objects/      cached intermediates, so re-running a late rule is cheap
├── report.html
└── run_manifest.json
```

`run_manifest.json` holds the fully resolved config, package versions, every
seed, `sessionInfo()` and SHA-256 checksums of the inputs. A result that cannot
be traced to its parameters is not a result.

## Example: trajectory classes carry biology

Classes are defined from counts and timepoints alone, so a fair test is
whether they also differ in things the clustering never saw. This example runs
the workflow on every peak of the timecourse behind the demo: ATAC-seq of CD8+ T
cells in LCMV Armstrong infected mice, naive through day 8 post infection
(129,076 peaks, 54,493 of them dynamic). Each class has its own genomic
context and its own motif signature compared with static peaks.

![Trajectory classes, genomic context and motif enrichment](docs/img/tcell_classes_overview.png)

*A: the five trajectory classes, with the number of peaks in each. B: where
each class sits in the genome, with static peaks as the reference; the dashed
line marks the promoter fraction of static peaks. C: known motifs (JASPAR2024)
enriched in each class over static peaks, by MEME-suite SEA; a dot marks
q < 1e-5.*

Peaks that close after activation (Decreasing) carry motifs of TCF7 and LEF1,
factors of the naive and memory T cell state. Peaks that open transiently,
around day 3, carry motifs of the AP-1 and BATF factors induced by T cell
activation. Peaks that open late carry motifs of ETS factors. The clustering
was never told any of this.

Scripts, tables and full reproduction steps are in
[`examples/tcell_motifs/`](examples/tcell_motifs/README.md).

## Relationship to other tools

**DESeq2** and **DEGreport** do the statistical work: DESeq2 decides whether
a feature changes over time, and `degPatterns` groups the ones that do by
correlation of their trajectories. Run by hand, they leave you with numbered
clusters, a set of thresholds spread across a notebook, and the job of deciding
what each cluster means. This workflow runs them with every threshold in one
validated config, and makes that last step a stated rule.

**Mfuzz** (soft clustering) and **TCseq** are alternative ways to cluster
timecourse profiles. They could stand in for the clustering step here; the
naming step works on cluster mean profiles and does not care which clusterer
produced them.

**ImpulseDE2** and **maSigPro** fit parametric models of expression over time,
impulse and polynomial respectively. They are more powerful when your
trajectories really do follow that form, and answer a different question:
model fit and significance rather than a taxonomy of shapes.

**Tempora** and other trajectory-inference methods order *cells* along a
pseudotime. This workflow is for bulk timecourses with real sampled
timepoints; you tell it the times.

What the workflow contributes on top of the tools it wraps: multi-arm designs
that share one baseline, the effect-size gate, automated naming of trajectory
classes, guards against silent failures in the underlying tools, and a
reproducibility record.

## Documentation

- [`docs/method.md`](docs/method.md) why the method is what it is
- [`docs/parameters.md`](docs/parameters.md) how to choose the thresholds
- [`docs/troubleshooting.md`](docs/troubleshooting.md) symptoms and fixes
- [`docs/faq.md`](docs/faq.md)
- [`vignettes/getting-started.Rmd`](vignettes/getting-started.Rmd)

## Installation

pixi is the primary route and installs Snakemake itself, so nothing needs to be
on your PATH first:

```bash
pixi install && pixi run demo
```

Two other routes to the same pinned stack: `snakemake --use-conda` picks up
`workflow/envs/r.yaml` if you already run conda, and the container published to
GHCR carries the whole environment. `pixi.lock` is committed and pins every
package to an exact build for linux-64 and osx-arm64.

## Testing

```bash
bash tests/test_static.sh        # tier 1, seconds, no execution
bash tests/test_endtoend.sh      # tier 2, minutes, full run on the demo
Rscript tests/test_correctness.R # tier 3, recovery of known trajectory shapes
```

Tier 3 simulates counts from five known shapes and checks what comes back.

## Citation

To cite the workflow, see [`CITATION.cff`](CITATION.cff). Please also cite
the tools that do the statistical work:

- **DESeq2**: Love MI, Huber W, Anders S. Moderated estimation of fold change
  and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15:550
  (2014). [doi:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)
- **DEGreport**: Pantano L. DEGreport: Report of DEG analysis. R package,
  Bioconductor. [doi:10.18129/B9.bioc.DEGreport](https://doi.org/10.18129/B9.bioc.DEGreport)
- **Snakemake**: Mölder F, Jablonski KP, Letcher B, et al. Sustainable data
  analysis with Snakemake. *F1000Research* 10:33 (2021).
  [doi:10.12688/f1000research.29032.2](https://doi.org/10.12688/f1000research.29032.2)

If you use the bundled demo data or the T cell example, cite the paper the
data come from, not this repository:

> McDonald BD\*, Chick BY\*, Ahmed NU, et al. Canonical BAF complex activity
> shapes the enhancer landscape that licenses CD8+ T cell effector and memory
> fates. *Immunity* 56(6):1303-1319.e5 (2023).
> [doi:10.1016/j.immuni.2023.05.005](https://doi.org/10.1016/j.immuni.2023.05.005)

\* Equal contribution.

## License

MIT. See [`LICENSE`](LICENSE).
