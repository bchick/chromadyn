# Method

Why chromadyn does what it does. For *how* to set the thresholds, see
[parameters.md](parameters.md).

The method is not new here. It was extracted from an MCF7 ATAC-seq analysis
whose notebooks will be released with its manuscript; see
[`provenance/`](provenance/). File and line references below, such as
`01_degpatterns_clustering.Rmd:380`, point into those notebooks. This document
explains the reasoning behind each step and records where chromadyn departs
from the original, and why.

## The question

Given counts over a timecourse: which features change over time, and what
distinct shapes do those changes take?

That second half is the part most tools leave to you. A differential test
returns a list. A clustering returns numbered groups. Neither tells you that
your data contains five temporal programs and which features belong to each.

## Two gates, not one

A feature is dynamic when it passes **both**:

1. a significance test over time, and
2. a floor on the **range of its per-timepoint means**, on the transformed scale.

The second gate is the one people leave out, and it is the one that matters at
scale. With enough replication the LRT returns thousands of features whose
total excursion across the entire timecourse is a few percent. They are real,
in the sense that the change is larger than the noise. They are also useless
as trajectories, and clustering them produces confident-looking classes built
on nothing. The range gate is an effect-size floor that says: changed, and
changed by enough to have a shape worth naming.

chromadyn reports the two gates separately, and says so explicitly when the
range gate removes nothing, rather than letting you believe an effect-size
filter was applied when it was not.

### `min_range` is not a fold change

In the source analysis this parameter was called `lrt_lfc`. It is not a log2
fold change. It is the range, maximum minus minimum, of the per-timepoint
replicate means **on the variance-stabilized scale**.

The names are close enough that anyone reading the original notebook will set
it an order of magnitude wrong. chromadyn renames it `min_range` and states
the scale in the config comment. If you change `transform.method`, the scale
changes and this threshold must be re-chosen.

## One model per arm, sharing a baseline

A timecourse with two treatment arms and a single unstimulated baseline is
fitted as **two separate models**, each using the same t=0 libraries, with no
arm term anywhere:

```r
dds_hrg <- DESeqDataSetFromMatrix(counts[, hrg_and_shared], design = ~ time)
dds_hrg <- DESeq(dds_hrg, test = "LRT", reduced = ~ 1)
```

The obvious alternative is one model with `~ arm * time` and an interaction
term. That answers a different question: *does the arm modify the time
response*, a comparison between arms. The question here is *what shape does
each arm's response have*, which is a description of each arm on its own
terms. Fitting them jointly would also make each arm's dispersion estimate
depend on the other's, so adding a third arm would change the first arm's
results.

Expressing this needs nothing more than a sentinel value in the `group`
column, `shared` by default. Rows carrying it join every arm. In the source
project this was a line of recoding buried in a preprocessing notebook
(`00_preprocessing.Rmd:244`); making it a samplesheet value is the single
change that most improves reusability, because it is the first thing a new
user with this design needs and will not find in other tools.

## Cluster first, then name

chromadyn clusters into however many groups the data supports, then collapses
those into a few named classes. It does not cluster directly to k.

Clustering directly to five groups forces the answer. Clustering to the
data's own granularity and then merging keeps the evidence: you can look at
the dendrogram and the silhouette diagnostic and see whether the merge you
took was defensible. Both are emitted on every run whatever the assignment
method is, so you can see the cut you did not make.

### The half that had to be automated

In the source analysis, the merge was done by eye. Fourteen cluster numbers
were typed into a `tribble()` by hand (`01_degpatterns_clustering.Rmd:380-403`).
That is the single reason the original assignment could not be reproduced
from the code, and automating it is most of the point of this repository.

### Why shape labelling rather than cutting the tree

The obvious automation is `cutree(hc, k)`. chromadyn's default instead labels
each cluster from the shape of its own mean profile, and lets clusters sharing
a label form a class. `hclust` remains selectable.

The reason is empirical. On the bundled T cell demo, cutting the dendrogram at
k=3 produces a class containing one Decreasing, one Increasing and one
Transient cluster together, because those three happen to sit in the same
branch. Labelling each cluster by its own centroid separates all three. The
reference classes are a per-feature statement about trajectory shape, not a
property of dendrogram structure, so reading them off the structure is the
wrong operation. You can see this directly in
`results/figures/<arm>_dendrogram.png`, where the tip colours are the shape
labels and the tree is the structure they cut across.

## Naming the shapes

Two stages, because one broad class usually contains several real shapes: an
"Increasing" group is often a mixture of transient, sustained and late rises
that only separate when split again.

Stage 1 asks whether a cluster rises, falls, or spikes and returns. Stage 2
takes one class, splits it by k-means on per-feature profiles, and names the
pieces. The predicates are the source's, with its literal timepoints
generalized.

### Generalizing the hard-coded timepoints

The source's rules refer to the literal strings `"60"` and `"240"`, so they
work on one experiment and nothing else. `240` is straightforwardly the last
timepoint. `60` took more care.

Reading it as the **middle post-baseline index** reproduces 60 on the
`0/30/60/120/240` grid, and is still wrong. That grid is roughly log spaced,
so its middle index falls a quarter of the way through in *time*. On a
near-linear grid such as the demo's `0/3/5/8` days, the middle index falls at
62% of the way through, by which point a genuinely late-rising trajectory has
already risen, and the late-versus-sustained test silently stops working.

chromadyn takes the observed timepoint nearest `mid_fraction` of the last
timepoint, default 0.25. That is exactly 60 on the source grid and 3 rather
than 5 on the demo's.

Stage 1 and stage 2 have separate threshold blocks. Sharing them is tempting
and wrong: a `cross` margin tuned to separate a late rise from a sustained one
also stops stage 1 recognising a transient, because the same comparison means
something different when every candidate is already rising.

### When the rules cannot name a shape

If two of the pieces of a split would get the same name, chromadyn numbers
them all instead, prints each centroid's statistics, and names the parameter
to adjust. This is the source's behaviour preserved. Numbered labels are a
signal, not a failure: they mean your grid or your data does not support the
distinction the vocabulary is trying to draw. `config/demo.yaml` carries a
worked example of tuning it.

## What the transform is for

Everything downstream works on the variance-stabilized matrix: the range gate,
the clustering, the heatmaps. That is why `min_range` is expressed against it
rather than as a fold change, and why `transform.blind` defaults to `TRUE`:
the same matrix feeds both selection and clustering, so a transform informed
by the comparison it will be used to test is a subtle form of double dipping.

Turning on `batch_correct` deliberately breaks that, which is why it is off by
default and recorded in the manifest when used. chromadyn refuses to correct a
batch that is confounded with time, checking whether the model matrix for
`~ batch + time` is full rank rather than counting timepoints per batch: a
batch spanning t=0,3,5 against another holding only t=8 passes the naive check
while being perfectly confounded.

## Notes on DEGreport

chromadyn uses `DEGreport::degPatterns` for the clustering. Three things about
it are worth knowing, all verified in versions 1.36.0 and 1.42.0 rather than
taken from the documentation.

**`cutoff` does nothing.** It is documented everywhere as a correlation cutoff
for merging clusters. `degPatterns` passes it to an internal function that
declares the argument and never reads it, and the reduce step does not accept
it at all. Setting it changes nothing, and the original analyses that set it
were unaffected by it. chromadyn keeps the key so old configs still parse, and
warns if you set it. Use `n_clusters` to control cluster count.

**`reduce` is not a cluster merger.** It walks each cluster and drops any
feature that is a boxplot outlier in any single column. It shrinks clusters
and never combines them.

**`minc` drops features rather than merging them**, and the comparison is
strictly greater than. Features in undersized clusters are removed from the
result entirely. chromadyn reports them as `Unassigned` with an
`assign_method`, rather than letting them disappear from the output.

`degPatterns` also overwrites its input rownames with `make.names()`, so
`chr1:7401731-7402231` comes back as `chr1.7401731.7402231`. That is the
origin of the corrupted-coordinate failure described below.

## Guardrails

Each of these is a real silent failure from the source analyses, implemented
as an assertion rather than a warning.

**Numeric time ordering.** `factor(as.character(c(0,15,30,60,120,240)))` sorts
lexicographically to `0, 120, 15, 240, 30, 60`. Every trajectory becomes noise
and nothing errors. The timepoint grid is derived once from the samplesheet
and asserted numeric and increasing, at parse time and again at every point a
factor is built.

**Metadata alignment.** `identical(rownames(metadata), colnames(matrix))`, not
`setequal`. Order matters, and a reordered `merge()` result attached
positionally is a well-known way to produce a confident wrong answer.

**Coordinates are never round-tripped through a feature ID.** A variant of the
source analysis parsed coordinates back out of `peak_id` with a regex that
stopped matching once `make.names()` had replaced the separators, and wrote
thousands of BED rows reading `chr1.7401731.7402231  NA  NA`, silently.
chromadyn carries `chr`, `start` and `end` as real columns from the input BED
to the output, keeps an explicit key table to undo the mangling, and asserts
the coordinates are valid immediately before writing each file.

**Everything is seeded.** The source seeded its k-means but not `degPatterns`,
while its methods section said all modalities were seeded; that clustering was
reproducible only through a cached RDS. chromadyn seeds before every
stochastic step and records every seed in `run_manifest.json`.

**Supercluster map coverage.** A join against an incomplete map silently drops
the unmapped clusters. chromadyn fails, naming them.

**Scale.** `degPatterns` cost is cubic in feature count. See
[parameters.md](parameters.md) for the measured curve and what chromadyn does
about it.

## What the samplesheet replaced

The source project had no samplesheet. It parsed condition, time and replicate
out of BAM filenames with a regex, including a hard-coded correction for a
mislabelled timepoint:

```r
time_min <- ifelse(time_min == 5, 15, time_min)
```

That correction is invisible to anyone reading the results. The samplesheet
exists to make that class of hidden fix impossible: if a timepoint is wrong,
you fix it in a file that is an input to the analysis and travels with it.
