# Troubleshooting

Organized by the symptom you see, not by the code that produced it.

## The run stops at validation

That is the intended behaviour. `01_validate.R` runs 16 checks before anything
expensive and fails rather than warning, because the failures it guards
against produce a confident wrong answer rather than a crash. Every check is
written to `results/qc/validation_report.tsv` with a status and a message.

### "samplesheet and counts columns disagree"

Sample IDs must match counts column names **exactly**. timecourse-patterns prints which
IDs are on each side and will not guess a mapping or quietly analyse a subset.

Most often this is a suffix the quantifier left on. nf-core/atacseq, for
example, leaves `.mLb.clN.sorted.bam` on the column names of its consensus
matrix. Strip it in your samplesheet or in the matrix, not by hand in one
place.

### "arm 'X' has one library per timepoint"

`DESeq(test = "LRT")` over `~ time` cannot estimate dispersion without
replication. Use:

```yaml
differential:
  method: none
  top_n: 5000
```

See [parameters.md](parameters.md#single-replicate-timecourses).

### "cluster.degpatterns.minc is N but only M features survive the prefilter"

`minc` scales with the number of features reaching the clusterer, not with
your input. Use `minc: auto`, or the value the message suggests.

### "batch column is confounded with time"

Every level of your batch variable occurs at one timepoint, or the model
matrix for `~ batch + time` is otherwise rank deficient. Removing that batch
would remove the trajectory along with it. Either drop the batch correction or
get a design where batches span timepoints.

The check is on matrix rank, not on counting timepoints per batch, because a
batch spanning t=0,3,5 against another holding only t=8 fails the naive check
while being perfectly confounded.

### "time must be numeric"

Put `4h` in `time_unit`, not in `time`. The column is a number; units are
declared once in the config and never parsed from values.

## Every trajectory looks like noise

Almost always the time-ordering trap. If your time column is character,
`factor()` sorts it lexicographically:

```
0, 120, 15, 240, 30, 60
```

Every profile is then reordered into nonsense and nothing errors. timecourse-patterns
asserts numeric, increasing factor levels at validation and again wherever a
time factor is built, so if you are seeing this with timecourse-patterns, check that the
figure you are looking at came from this pipeline.

Second possibility: replicates that do not agree. Look at
`results/qc/pca.png`. Replicates should sit on top of each other and
timepoints should separate. If they do not, nothing downstream is trustworthy.

## Almost everything is called dynamic

Expected on a small or deliberately enriched feature set; the bundled demo is
72% dynamic because it was built that way.

On a whole-genome run it means a gate is not biting. The `select_dynamic` log
reports each gate separately:

```
[select_dynamic] arm 'WT': padj < 0.01 -> 3466 | range >= 0.5 -> 3253 | both -> 3251
```

If the range count is nearly your whole matrix, raise `min_range`. timecourse-patterns
prints an explicit note when the range gate removes nothing at all.

## Nothing is called dynamic

The run fails with the count each gate passed. Loosen `differential.fdr` or
`differential.min_range`. If the significance gate passed nothing, check the
PCA first: a timecourse with no structure has nothing to find.

## The clustering takes forever

`degPatterns` cost is **cubic** in feature count: about 72 s at 3,251
features, an hour at 12,000, and roughly 18 hours at 30,000. The default
`cluster.max_features: 6000` caps it; above that timecourse-patterns subsamples and
assigns the remainder by correlation. Lower the cap while iterating. See
[parameters.md](parameters.md#clustermax_features-and-runtime).

## I get "Increasing-1", "Increasing-2" instead of names

The shape rules could not tell the pieces of the split apart, so they were
numbered. That is a signal, not a failure. The log prints each centroid's
statistics and names the parameter to adjust; usually
`superclusters.split.cross`. There is a worked example in
[parameters.md](parameters.md#naming-thresholds).

## Many features come back "Unassigned"

Two possible causes, distinguished by `assign_method` in the cluster table and
by `clusters/<arm>_cluster_assignment_qc.tsv`.

`unassigned_minc`: they were in clusters smaller than `minc`, which DEGreport
drops rather than merges. Lower `minc`.

`unassigned_low_cor`: they were held out by the `max_features` cap and did not
correlate with any cluster centroid above `cluster.assign_min_cor`. Raise
`max_features` if you can afford the runtime, or lower `assign_min_cor` if you
are willing to place weakly matching features.

## The BED files have NA coordinates

They should not, and timecourse-patterns asserts against exactly this immediately before
writing each file, because the source analysis shipped thousands of rows of
`chr1.7401731.7402231  NA  NA` without noticing.

If you hit it, your features BED does not cover every feature in the counts
matrix. Validation reports the count and the first few missing IDs.

Never reconstruct coordinates from a feature ID. `degPatterns` runs
`make.names()` over rownames, so `chr1:100-200` becomes `chr1.100.200` and any
regex expecting `:` and `-` silently produces NA.

## `library(DEGreport)` fails with "no package called GenomeInfoDbData"

pixi skips conda post-link scripts by default, and that package installs its
payload in one. Run `pixi run postinstall`, which the `demo` and `test` tasks
depend on anyway.

## A rule reruns when nothing changed

Snakemake reruns on changes to code, parameters and inputs, not only on
timestamps. After editing a script, that is correct. To see why:

```bash
pixi run --frozen snakemake --configfile config/demo.yaml -n -r
```

To rerun on timestamps only, add `--rerun-triggers mtime`.

## The log files are empty

They should not be. Snakemake only redirects output into `log:` for `shell:`
directives, not `script:`, so timecourse-patterns tees messages itself. If a log is
empty the script probably failed before `cd_init()`; the error will be on the
console and in the Snakemake log under `.snakemake/log/`.

## I want to rerun one step without redoing the clustering

The expensive objects are cached in `results/objects/`. Deleting a downstream
output reruns only that rule:

```bash
rm results/clusters/WT_clusters.tsv
pixi run --frozen snakemake --configfile config/demo.yaml -j 4
```

Any script can also be run on its own for debugging, outside Snakemake
entirely. Each begins with a block that builds an equivalent `snakemake`
object:

```bash
pixi run --frozen Rscript workflow/scripts/08_superclusters.R
```

Edit the `wildcards` and `configfile` in that block to point at your run.
