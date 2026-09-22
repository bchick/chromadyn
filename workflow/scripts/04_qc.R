# ---------------------------------------------------------------------------
# 04_qc.R: look at the data before modelling it.
#
# PCA, sample correlation and within-condition replicate correlation. None of
# these gate the pipeline; they exist so that a design problem is visible
# before you spend an afternoon interpreting clusters built on it.
# ---------------------------------------------------------------------------

if (!exists("snakemake")) {
  source("workflow/scripts/lib/standalone.R")
  snakemake <- cd_standalone(
    rule       = "qc",
    configfile = "config/demo.yaml",
    input      = list(dds = "obj_dds", transformed = "obj_transformed"),
    output     = list(rep_cor_table = "replicate_cor_tsv"))
}
source(file.path(snakemake@scriptdir, "lib", "common.R"))
source(file.path(snakemake@scriptdir, "lib", "namespace_fix.R"))
source(file.path(snakemake@scriptdir, "lib", "theme.R"))
source(file.path(snakemake@scriptdir, "lib", "plots.R"))
cd_init(snakemake)

suppressPackageStartupMessages(library(DESeq2))

dds <- load_obj(snakemake@input$dds)
mat <- load_obj(snakemake@input$transformed)
coldata <- as.data.frame(SummarizedExperiment::colData(dds))
coldata$time <- as.numeric(as.character(coldata$time))
assert_aligned(coldata, mat, "colData against transformed matrix")

formats <- cfg("figures.formats")
dpi <- cfg("figures.dpi")

pca <- run_pca(mat, coldata)
save_figure(pca$plot, cd_path("pca", fmt = formats[1]),
            width = FIG_W_1COL, height = FIG_W_1COL, formats = formats, dpi = dpi)
message(sprintf("[qc] PCA: PC1 %.1f%%, PC2 %.1f%%", pca$var_exp[1], pca$var_exp[2]))

save_figure(plot_sample_correlation(mat, coldata),
            cd_path("sample_cor", fmt = formats[1]),
            width = FIG_W_1_5COL, height = FIG_W_1_5COL, formats = formats, dpi = dpi)

rep_cor <- compute_rep_cor(mat, coldata)
save_figure(plot_replicate_correlation(rep_cor),
            cd_path("replicate_cor", fmt = formats[1]),
            width = FIG_W_1COL, height = FIG_PANEL, formats = formats, dpi = dpi)
write_tsv_atomic(rep_cor, snakemake@output$rep_cor_table)

if (nrow(rep_cor)) {
  message(sprintf("[qc] replicate correlation: %d pairs, r %.3f to %.3f (median %.3f)",
                  nrow(rep_cor), min(rep_cor$r), max(rep_cor$r), stats::median(rep_cor$r)))
  cd_note("replicate_cor_min", min(rep_cor$r))
  cd_note("replicate_cor_median", stats::median(rep_cor$r))
} else {
  message("[qc] no replicated conditions, replicate correlation skipped")
}
cd_note("pca_var_explained", as.numeric(pca$var_exp))

cd_finish()
