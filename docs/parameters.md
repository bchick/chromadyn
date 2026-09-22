# Choosing parameters

Every threshold is in `config/config.yaml`. This is how to pick the ones that
matter. For *why* they exist, see [method.md](method.md).

## Three real parameter sets

From analyses this method has actually been run on, as worked examples of the
range that works in practice.

| | Large ATAC timecourse | CUT&RUN timecourse | T cell ATAC demo |
|---|---|---|---|
| features in | ~50,000 peaks | ~30,000 peaks | 5,000 peaks |
| features clustered | ~14,000 | ~2,000 | 3,251 |
| timepoints | 5 | 5 | 4 (0-8 days) |
| `min_mean_counts` | 10 | n/a | 10 |
| `fdr` | 0.01 | 0.05 | 0.01 |
| `min_range` | 0.5 | 0.5 | 0.5 |
| `minc` | 50 | 30 | 50 |

The pattern: `fdr` and `min_range` are stable across datasets, `minc` scales
with how many features reach the clusterer.

## `minc`, minimum cluster size

**Scale it with the number of features being clustered, not with the input.**
50 is a sensible floor when ~14,000 features reach the clusterer. On 2,000 it
means a cluster must hold 2.5% of your data to survive, which is far too
aggressive, and chromadyn fails at validation rather than let you find out
from a suspiciously short cluster list.

`minc: auto` sets `max(15, round(nrow / 100))`, roughly 1% of the clustered
features. Use it unless you have a reason not to.

Remember that DEGreport **drops** the members of undersized clusters; it does
not merge them elsewhere. chromadyn reports them as `Unassigned` in the
cluster table, with `assign_method` saying why. If that count is large, `minc`
is too high.

## `min_range`, the effect-size gate

On the **transformed** scale, not a fold change. 0.5 on a variance-stabilized
matrix is a reasonable starting point and has held across every dataset above.

How to check it is doing something: the `select_dynamic` log prints how many
features each gate passes separately.

```
[select_dynamic] arm 'WT': padj < 0.01 -> 3466 | range >= 0.5 -> 3253 | both -> 3251
```

If the range count is much larger than the significance count, the gate is not
constraining anything and chromadyn says so explicitly. Raise it until it
starts removing features that the test called significant but whose
trajectories are flat when you plot them.

If you change `transform.method`, re-choose this. The scale changes.

## `cutoff`

Leave it alone. It has no effect in any released DEGreport version; see
[method.md](method.md#notes-on-degreport). To control cluster count use
`cluster.degpatterns.n_clusters`.

## `n_clusters` and how many clusters you get

By default DEGreport cuts its tree at the diana divisive coefficient, which is
data-determined: you do not choose the number of clusters and you cannot
predict it. Setting `n_clusters` cuts at a fixed k instead.

One caveat: `minc` is applied *after* the cut, and drops groups smaller than
it. Asking for 5 can therefore give you 3. If that happens, lower `minc`.

## `k`, the number of trajectory classes

Only used when `superclusters.method` is `hclust` or `kmeans`. The default
`shape` method does not need it.

Every run writes `clusters/<arm>_k_diagnostics.tsv` whatever the method:

```
k  silhouette  wss
2  0.354       10.16
3  0.413        4.82
4  0.410        2.12
5  0.348        1.02
```

Higher silhouette is better. Within-cluster sum of squares always falls with
k, so look for where it stops falling steeply. Here both point at k=3.

Read it as evidence, not as an answer. It describes cutting the dendrogram of
cluster profiles, and the default assignment method does not cut that
dendrogram at all.

## `cluster.max_features` and runtime

**`degPatterns` cost is cubic in feature count.** Measured on 9 libraries,
R 4.4.3, DEGreport 1.42.0:

| features | time |
|---|---|
| 1,500 | 7 s |
| 3,251 | 72 s |
| 4,004 | 144 s |

That fits `n^3.06`, which extrapolates to:

| features | extrapolated |
|---|---|
| 6,000 | ~8 min |
| 8,000 | ~20 min |
| 12,000 | ~1 hour |
| 30,000 | ~18 hours |

The default cap is 6,000. Above it, chromadyn clusters a stratified subsample
with the recorded seed and assigns the remainder by correlation to cluster
centroids, reporting both counts in
`clusters/<arm>_cluster_assignment_qc.tsv`.

The subsample path is not a large loss. Capping the demo at 1,500 recovered
the same 8 clusters and assigned 1,750 of 1,751 held-out features. Raise the
cap if you would rather wait; lower it if you are iterating.

`cluster.assign_min_cor` (default 0.6) is the floor for a held-out feature to
join a centroid. Below it the feature is left `Unassigned` rather than forced
somewhere.

## Naming thresholds

`superclusters.shape` controls stage 1, `superclusters.split` stage 2. They
are separate on purpose.

**`cross`** is the value a centroid must be below, at the mid timepoint, to
count as a late riser. Profiles are z-scores, so 0 is the feature's own mean
and a negative value is a margin below it. The default 0.0 is the source
analysis's value, kept for fidelity, and it is often too permissive: a
sustained riser that is merely *at* its mean early looks late.

Worked example, from the bundled demo. On a 0/3/5/8 day grid the split
reported:

```
[classify] shape labelling was ambiguous (2 of 3 centroids share or lack a label)
[classify]   reference timepoints: first=0 mid=3 last=8
[classify]   1   peak_time=8  mid_value=-0.02  last_value=+0.85  late_ratio=+1.00
[classify]   2   peak_time=3  mid_value=+1.24  last_value=+0.02  late_ratio=+0.02
[classify]   3   peak_time=8  mid_value=-0.91  last_value=+1.20  late_ratio=+1.00
```

Centroids 1 and 3 both sit below their mean at day 3, so both satisfy
"below the mean at mid" and collide. Centroid 1 is at -0.02, essentially at
its mean; centroid 3 is at -0.91, clearly still down. Setting `cross: -0.5`
separates them into Sustained Increasing and Late Increasing.

**`late_ratio_cutoff`** (0.5) is the centroid's last value over its maximum.
Below it, the class peaked and came back down: Transient Increasing.

**`mid_fraction`** (0.25) picks the mid timepoint as the observed time nearest
that fraction of the last. On `0/30/60/120/240` it gives exactly 60, the value
the source hard-coded. Override with `mid_time` if your sampling is very
uneven.

**Numbered labels are a signal.** If you get `Increasing-1`, `Increasing-2`,
the rules could not tell the pieces apart. The log prints each centroid's
statistics; usually `cross` is the knob.

## Single-replicate timecourses

`DESeq(test = "LRT")` over `~ time` needs replication. chromadyn detects a
design with one library per timepoint at validation and fails the `lrt` path
with a message pointing at the supported alternative:

```yaml
differential:
  method: none
  top_n: 5000
```

That ranks features by the variance of their timepoint means and takes the top
N. The `min_range` gate still applies, so `top_n` is a ceiling and not a quota.
No significance is claimed, and none should be reported.

## Number of timepoints

Three is the minimum for a shape to exist, and chromadyn warns at three: the
transient and late classes cannot separate cleanly. Four works, as the demo
shows. Five or more is comfortable.

Uneven spacing is fine and is handled explicitly; see `mid_fraction` above.
