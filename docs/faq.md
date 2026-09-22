# FAQ

**Does it work on RNA-seq?**

Yes. chromadyn is assay-agnostic; it starts at a counts matrix and does not
care what the rows are. Leave `input.features` null for gene mode, and the
region-only rules are skipped. Genes are usually fewer than peaks, so use
`minc: auto`.

**Can I use already-normalized data?**

Set `transform.method: none`. The differential step is then limited to
`method: none` as well, since DESeq2 needs raw counts. Remember `min_range` is
expressed on whatever scale you supply.

**Do I need the same timepoints in every arm?**

Yes. Validation fails if an arm is missing a timepoint the others have,
because the shapes would not be comparable.

**Can arms have different numbers of replicates?**

Yes, and different numbers per timepoint within an arm. The bundled demo has
two replicates at three timepoints and three at the fourth. Validation warns
so the imbalance is on the record.

**What if each arm has its own baseline rather than a shared one?**

Set `input.baseline_group` to a value no row uses, and give each arm its own
t=0 rows with that arm's group label.

**How many timepoints do I need?**

Three minimum, four works, five or more is comfortable. Uneven spacing is
handled explicitly.

**Why are my class names not the five canonical ones?**

The five names are a vocabulary, not a promise. You get the ones the data
supports. If you get numbered labels instead, the rules could not separate the
shapes; see [troubleshooting.md](troubleshooting.md).

**Can I reproduce an old hand-made cluster assignment?**

Yes, with `superclusters.method: overrides` and a `{cluster_id: name}` map.
chromadyn fails if the map misses any cluster, naming the ones it missed.

One caveat that matters: cluster IDs are meaningful only within a single run
over a single feature set. A map from one run cannot be applied to another
with different features, and doing so silently means something else.

**Is it deterministic?**

Yes. Every stochastic step is seeded and every seed is recorded in
`run_manifest.json`. The tier 2 golden tables are diffed across independent
runs, which is what demonstrates it.

**Can I run it on a cluster?**

It is an ordinary Snakemake workflow, so the usual executor plugins apply.
`--use-conda` picks up `workflow/envs/r.yaml` if you would rather not use
pixi.

**Why is `cutoff` in the config if it does nothing?**

So configs written against the original notebooks still parse. chromadyn warns
if you set it. See [method.md](method.md#notes-on-degreport).

**Can I add my own trajectory vocabulary?**

Not through the config today. The vocabularies live in `CD_VOCAB` in
`workflow/scripts/lib/classify.R` and the predicates beside them; adding one
is a small, self-contained edit.

**How do I cite it?**

See `CITATION.cff`. If you use the bundled demo data, cite the paper it comes
from rather than this repository: `demo/PROVENANCE.md` has the details.
