# ---------------------------------------------------------------------------
# 03_transform.R: variance-stabilize the counts, optionally remove a batch.
#
# Everything downstream works on this matrix: the range gate in
# 06_select_dynamic, the clustering in 07, the heatmaps in 10. It is the
# common scale, which is why min_range is expressed on it and not as a fold
# change.
# ---------------------------------------------------------------------------

if (!exists("snakemake")) {
  source("workflow/scripts/lib/standalone.R")
  snakemake <- cd_standalone(
    rule       = "transform",
    configfile = "config/demo.yaml",
    input      = list(dds = "obj_dds"),
    output     = list(transformed = "obj_transformed"))
}
source(file.path(snakemake@scriptdir, "lib", "common.R"))
source(file.path(snakemake@scriptdir, "lib", "namespace_fix.R"))
cd_init(snakemake)

suppressPackageStartupMessages(library(DESeq2))

dds <- load_obj(snakemake@input$dds)
method <- cfg_req("transform.method")
blind <- cfg("transform.blind", TRUE)

mat <- switch(
  method,
  vst  = SummarizedExperiment::assay(vst(dds, blind = blind)),
  rlog = SummarizedExperiment::assay(rlog(dds, blind = blind)),
  none = {
    message("[transform] method: none, using counts as supplied")
    as.matrix(counts(dds))
  },
  stop(sprintf("unknown transform.method '%s'", method), call. = FALSE)
)
message(sprintf("[transform] %s%s: %d features x %d samples",
                method, if (method == "none") "" else sprintf(" (blind = %s)", blind),
                nrow(mat), ncol(mat)))

# --- optional batch correction ---------------------------------------------
# Off by default. It makes the matrix design-aware, which is exactly what
# transform.blind = TRUE was avoiding, so it is a deliberate choice and is
# recorded in the manifest.
bc <- cfg("batch_correct.method", "none")
if (bc != "none") {
  bcol <- cfg_req("batch_correct.column")
  cd <- as.data.frame(SummarizedExperiment::colData(dds))
  cd_assert(bcol %in% colnames(cd),
            "batch_correct.column '%s' is not a colData column. Available: %s",
            bcol, paste(colnames(cd), collapse = ", "))
  batch <- factor(cd[[bcol]])
  cd_assert(bc == "limma", "unknown batch_correct.method '%s'", bc)
  cd_assert(requireNamespace("limma", quietly = TRUE),
            "batch_correct.method is limma but the limma package is not installed.")
  # Preserve the time effect explicitly, or removeBatchEffect will regress out
  # part of the signal along with the batch.
  design <- stats::model.matrix(~ cd$time)
  mat <- limma::removeBatchEffect(mat, batch = batch, design = design)
  message(sprintf("[transform] limma::removeBatchEffect on '%s' (%d levels), time preserved",
                  bcol, nlevels(batch)))
  cd_note("batch_correct", list(method = bc, column = bcol, levels = levels(batch)))
} else {
  cd_note("batch_correct", "none")
}

cd_note("transform_method", method)
cd_note("transform_blind", blind)
cd_note("n_features", nrow(mat))

save_obj(mat, snakemake@output$transformed)
cd_finish()
