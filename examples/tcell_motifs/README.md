# Example: are trajectory classes biology?

chromadyn defines trajectory classes from counts and timepoints alone. This
example asks whether those classes also differ in things the clustering never
saw: where in the genome they sit, and which transcription factor motifs they
carry. If they do, the classes reflect regulatory programs rather than an
artefact of the clustering.

![Trajectory classes, genomic context and motif enrichment](tcell_classes_overview.png)

## Data

The full wild-type CD8+ T cell ATAC timecourse behind the bundled demo
([Immunity 2023](https://doi.org/10.1016/j.immuni.2023.05.005), GEO
GSE228171): the same nine libraries and design as `demo/`, at 0, 3, 5 and 8
days, but every consensus peak on a standard chromosome (129,076) instead of
a 5,000-peak subsample. chromadyn runs with the demo settings unchanged, so
any difference from the demo comes from the data.

It calls 54,282 peaks dynamic and 19,035 static. The `Increasing` class is
split in three, and the split names its pieces numerically. They are renamed
here by the shape they actually have:

| chromadyn label | Shape | Peaks |
|---|---|---|
| Decreasing | Decreasing | 16,580 |
| Transient | Transient | 10,685 |
| Increasing-1 | Transient Increasing: rises by day 3, then slowly declines | 7,963 |
| Increasing-2 | Late Increasing: flat to day 3, then rises | 10,078 |
| Increasing-3 | Dip and recovery: falls at day 3, recovers above baseline | 8,976 |

## Method

Static peaks (accessible regions the same run judged not to change) are the
reference throughout. Against shuffled or genomic sequence, every class would
look enriched for the motifs of whatever keeps chromatin open. Against static
peaks, enrichment means the motif is specific to that temporal program.

- **Genomic context** (`annotate_classes.R`): ChIPseeker against GENCODE
  vM35, collapsed to five categories: promoter (within 1 kb of a TSS),
  promoter flank (1 to 3 kb), exon/UTR, intronic and distal intergenic. No
  category is called "enhancer"; that would need histone marks or perturbation
  data this analysis does not use.
- **Motifs** (`run_motifs.sh`): MEME-suite SEA 5.5.4 per class against 20,000
  static peaks sampled with a fixed seed. Every region is trimmed to 200 bp
  around the peak centre, so differences in peak width between classes cannot
  pass for enrichment. Motifs are JASPAR2024 CORE vertebrates non-redundant
  (879 motifs):

  ```bash
  curl -L -o examples/tcell_motifs/work/motifs/jaspar2024_vert_nr.meme \
    https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt
  # sha256 dd494278d356a4e170908c74d6d8eb746a2c7504d6210cd51017141f21a65b18
  ```

## Results

- **Decreasing**: TCF7, LEF1, TCF7L2 and HNF1A motifs, the TCF/LEF family
  associated with the naive and memory T cell state, which is closing as the
  cells activate.
- **Transient** and **Transient Increasing**: dominated by AP-1 and BATF
  (FOS, JUN, JUNB, FOSL1, BATF, BATF3, BACH1), the activation-induced factors,
  with enrichment ratios above 4 in the Transient class. Both classes are
  depleted of promoters relative to static peaks, which fits distal elements
  opened by activation.
- **Late Increasing**: ETS factors (ETS1, ETS2, ERG, GABPA, ETV1) and IKZF3.
- **Dip and recovery**: 56% promoter peaks, against 30% of static peaks. Its
  top motifs (KLF15, ZBED4, ZNF610, EGR4, TFDP1) are GC-rich, which may reflect
  the base composition of CpG-island promoters as much as specific factor
  binding. Treat these hits with that caveat in mind.

Full tables: [`motif_enrichment.tsv`](motif_enrichment.tsv) (every SEA result
per class) and [`genomic_annotation.tsv`](genomic_annotation.tsv).

## Reproducing

Run from the repository root. Beyond the chromadyn environment, the example
needs bedtools, the MEME suite, and the Bioconductor packages ChIPseeker and
txdbmaker, plus a GRCm39 FASTA with a `.fai` index and the GENCODE vM35 GTF.
Each script takes its input paths as arguments; the defaults point at the
author's local copies.

```bash
Rscript examples/tcell_motifs/build_inputs.R <featureCounts matrix> examples/tcell_motifs/work/inputs
pixi run --frozen snakemake -j 8 \
    --configfile config/demo.yaml examples/tcell_motifs/config.yaml
bash examples/tcell_motifs/run_motifs.sh examples/tcell_motifs/work/results <genome.fa> <motifs.meme>
Rscript examples/tcell_motifs/annotate_classes.R examples/tcell_motifs/work/results <gencode.vM35.gtf>
pixi run --frozen Rscript examples/tcell_motifs/plot_motifs.R
pixi run --frozen Rscript examples/tcell_motifs/make_figure.R
```

Pass both config files after a single `--configfile` flag. Repeating the flag
replaces the first file instead of stacking on it.

Intermediates go to `work/` (about 2.5 GB, not committed).
