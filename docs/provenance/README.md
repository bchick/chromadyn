# Provenance

The notebooks chromadyn's method was extracted from, kept verbatim for
reference. They are not run by anything and are not maintained here.

| File | What it is |
|---|---|
| `00_preprocessing.Rmd` | MCF7 ATAC preprocessing: consensus peaks, prefilter, VST, QC |
| `01_degpatterns_clustering.Rmd` | The method proper: per-arm LRT, the range gate, degPatterns, the hand-written supercluster map, the k-means sub-split |
| `run_chunks_1to6.R` | A plain-Rscript re-run of the first six chunks, the closest thing in the original to a pipeline step |
| `02_arid1a_degpatterns_clustering.Rmd` | The same method on CUT&RUN data, which is what shows which parts are method and which are dataset |
| `mcf7_methods.md` | The manuscript methods section |

Read `docs/method.md` for what chromadyn changed and why. Where that document
and these notebooks disagree about what the method *is*, the notebooks win.
Where they disagree about what the pipeline should *do*, chromadyn wins: the
hand-typed cluster map in `01_degpatterns_clustering.Rmd:380-403` is the
clearest example, and automating it is most of the point of this repository.
