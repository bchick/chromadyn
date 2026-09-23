# Demo data provenance

The bundled demo is a 5,000-feature subset of a published ATAC-seq timecourse.
It is real data, not simulated, and it is already public.

## Source

| | |
|---|---|
| Publication | McDonald BD\*, Chick BY\*, Ahmed NU, et al. Canonical BAF complex activity shapes the enhancer landscape that licenses CD8+ T cell effector and memory fates. *Immunity* 56(6):1303-1319.e5 (2023) |
| DOI | [10.1016/j.immuni.2023.05.005](https://doi.org/10.1016/j.immuni.2023.05.005) |
| GEO SuperSeries | [GSE228381](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228381) |
| ATAC sub-series | [GSE228171](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228171) |
| Genome | GRCm39 (mm39), GENCODE vM35 |
| Quantification | nf-core/atacseq 2.1.2, consensus peaks, featureCounts v2.0.1 in SAF mode with `-O --fracOverlap 0.2 -p` |
| Source matrix | `consensus_peaks.mRp.clN.featureCounts.txt`, 129,314 peaks x 62 libraries |

The source matrix is not committed. `build_demo.R` regenerates everything here
from it. It takes the matrix path as its first argument and otherwise reads
`demo/_source_featureCounts.txt`, a git-ignored location, so a copy can sit
there without being committed.

## Design

Nine libraries, one arm, four timepoints, in days post infection.

| sample | group | time | replicate | reads in peaks |
|---|---|---|---|---|
| `Naive_WT_REP1` | WT | 0 | r1 | 12.59 M |
| `Naive_WT_REP2` | WT | 0 | r2 | 16.91 M |
| `D3_WT_REP1` | WT | 3 | r1 | 6.12 M |
| `D3_WT_REP2` | WT | 3 | r2 | 3.83 M |
| `D5_WT_REP1` | WT | 5 | r1 | 12.45 M |
| `D5_WT_REP2` | WT | 5 | r2 | 12.79 M |
| `D8_WT_TE_Exp2_REP1` | WT | 8 | r1 | 6.55 M |
| `D8_WT_TE_Exp2_REP2` | WT | 8 | r2 | 7.13 M |
| `D8_WT_TE_Exp2_REP3` | WT | 8 | r3 | 9.00 M |

Three choices worth stating, because each is a judgement call someone might
reasonably make differently.

**D8 is the three terminal-effector (TE) libraries from experiment 2.** D8 in
the source has 15 WT libraries, across three sorted subsets (TE, EEC, MP) and
two experiments, while Naive, D3 and D5 are unsorted total CD8. The source
project's own temporal analysis pseudobulks across subsets by averaging
variance-stabilized values. chromadyn starts from raw counts and does its own
transform, so it cannot consume averaged values; rather than pool at count
level, the demo takes the TE libraries as they are. They are the deepest TE
set and they sit within a single experiment. D8 therefore has three replicates
against two at the other timepoints, which is a realistic and perfectly
analysable imbalance.

**The 48h timepoint is excluded.** FRiP is 0.115 for `48h_WT_REP2` and 0.314
for `48h_WT_REP1`, against 0.38 to 0.65 for every other library, and the two
carry 0.87 M and 1.73 M reads in peaks against 3.8 M to 16.9 M. Every temporal
analysis in the source project drops it as a QC failure.

**There is no `batch` column.** The only available batch variable, experiment 1
versus experiment 2, is perfectly confounded with time here: every D8 library
is from experiment 2 and no earlier timepoint carries an experiment label at
all. Shipping the column would invite a model that cannot be fit.

## Feature selection

Seed 42, in `build_demo.R`, which prints every count below when it runs.

1. Restrict to standard chromosomes, `chr1` to `chr19`, `chrX`, `chrY`. This
   drops 238 peaks on unplaced scaffolds. It is done by chromosome and never
   by row order: contig naming in this matrix mixes UCSC and Ensembl styles,
   scaffolds sort first, and `Interval_1` is on `GL456210.1`, so taking the
   first N rows would yield a demo made entirely of scaffold junk.
2. Split into four pools on the range of log2 CPM per-timepoint means, among
   features with mean count >= 10: **dynamic** (range in the top 20%),
   **flat** (bottom 15%), **middle** (everything between), plus a
   **low-count** pool of features with mean count between 1 and 10.
3. Sample 2,500 dynamic, 1,500 flat, 500 middle, 500 low-count.

The dynamic pool is **shape-stratified**: its z-scored profiles are k-means
clustered into 6 groups and sampled in proportion, so no trajectory shape is
lost to the subsample. The flat features exist so that the dynamic-versus-static
figure has something to contrast and so the false-positive rate of
`select_dynamic` is measurable. The middle band keeps the demo from being
artificially bimodal. The low-count features exist so the pipeline's
`filter.min_mean_counts` prefilter has real work to do.

An earlier version drew the flat pool from the bottom **half** of the range
distribution. On ATAC data this deep that is not flat at all: 89% of the demo
came back dynamic and the dynamic-versus-static figure had nothing to show.
Hence the bottom 15%.

## What the demo produces

With shipped defaults, as a reference for anyone changing this file:

| stage | result |
|---|---|
| prefilter, `min_mean_counts: 10` | 4,500 of 5,000 kept, 500 dropped |
| LRT, `padj < 0.01` | 3,466 |
| range gate, `min_range >= 0.5` | 3,253 |
| both (dynamic) | 3,251 of 4,500, 72% |
| degPatterns, `minc: 50` | 8 clusters, 3,251 features, about 70 s |

The dynamic fraction is not a sign that the subsample is enriched for change.
On the full matrix, 54,493 of the 73,528 peaks that pass the same prefilter
are dynamic (74%), so the demo is representative in this respect. Most peaks
in the full matrix that do not change are removed by the prefilter.

## Redistribution

These counts derive from data already deposited in GEO under GSE228381 and
published in *Immunity*. The subset is committed here for the demo under the
repository's MIT license; the underlying data remain subject to the terms of
the GEO deposition. Cite the paper, not this repository, when using the data
itself.
