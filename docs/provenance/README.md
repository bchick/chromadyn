# Provenance

The timecourse-patterns method was extracted from an MCF7 ATAC-seq timecourse analysis.
Its notebooks are not included in this repository yet: they will be released
together with the manuscript they belong to.

References elsewhere in this repository to files such as
`00_preprocessing.Rmd` or `01_degpatterns_clustering.Rmd:380-403` point into
those notebooks, and record which part of the original each piece of timecourse-patterns
was ported from.

| File | What it is |
|---|---|
| `00_preprocessing.Rmd` | Consensus peaks, prefilter, VST, QC |
| `01_degpatterns_clustering.Rmd` | The method proper: per-arm LRT, the range gate, degPatterns, the hand-written supercluster map, the k-means sub-split |
| `run_chunks_1to6.R` | A plain-Rscript re-run of the first six chunks |
| a CUT&RUN clustering notebook | The same method on CUT&RUN data, which shows which parts are method and which are dataset |

`docs/method.md` explains what timecourse-patterns changed and why. The clearest example
is the hand-typed cluster map at `01_degpatterns_clustering.Rmd:380-403`:
automating it is most of the point of this repository.
