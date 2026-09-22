# chromadyn

Temporal clustering of omics timecourses. Give it a counts matrix and a
samplesheet; get back the features that change over time, grouped into named
trajectory classes, with publication figures and an HTML report.

> Status: under construction. The quickstart below is not yet runnable.

```bash
git clone https://github.com/bchick/chromadyn
cd chromadyn
pixi install
pixi run demo
```

chromadyn starts at a counts matrix. It never touches FASTQ, BAM, or peak
calling, so it composes downstream of nf-core/atacseq, nf-core/rnaseq, DiffBind,
featureCounts, salmon, or any other quantifier, and it stays assay-agnostic:
features may be genomic regions (ATAC, CUT&RUN, ChIP) or genes (RNA-seq).

Full documentation lands with v0.1.0.
