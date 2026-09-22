# ---------------------------------------------------------------------------
# 05_differential.R: does each feature change over time, within one arm?
#
# One model per arm, fitted against the baseline libraries the arm shares with
# every other arm. There is no treatment or arm term anywhere: see docs/method.md
# for why a shared-baseline multi-arm timecourse is analysed as separate
# per-arm fits rather than as one model with an interaction.
#
# This script produces the statistics only. The two gates that decide which
# features are "dynamic" are applied in 06_select_dynamic.R, so that you can
# change a threshold without refitting the model.
# ---------------------------------------------------------------------------

if (!exists("snakemake")) {
  source("workflow/scripts/lib/standalone.R")
  snakemake <- cd_standalone(
    rule       = "differential",
    wildcards  = list(arm = "WT"),
    configfile = "config/demo.yaml",
    input      = list(dds = "obj_dds", transformed = "obj_transformed"),
    output     = list(obj = "obj_differential"),
    params     = list(samples = c("Naive_WT_REP1", "Naive_WT_REP2",
                                  "D3_WT_REP1", "D3_WT_REP2",
                                  "D5_WT_REP1", "D5_WT_REP2",
                                  "D8_WT_TE_Exp2_REP1", "D8_WT_TE_Exp2_REP2",
                                  "D8_WT_TE_Exp2_REP3")))
}
source(file.path(snakemake@scriptdir, "lib", "common.R"))
source(file.path(snakemake@scriptdir, "lib", "namespace_fix.R"))
cd_init(snakemake)

suppressPackageStartupMessages(library(DESeq2))

arm <- snakemake@wildcards$arm
samples <- unlist(snakemake@params$samples)
dds_all <- load_obj(snakemake@input$dds)
vst_all <- load_obj(snakemake@input$transformed)

cd_assert(all(samples %in% colnames(dds_all)),
          "arm '%s': sample(s) not in the counts matrix: %s",
          arm, paste(setdiff(samples, colnames(dds_all)), collapse = ", "))

# Column order is whatever arm_samples() decided, once, in the Snakefile. Both
# matrices are indexed by the same vector so they cannot drift apart.
cts <- counts(dds_all)[, samples, drop = FALSE]
vst_arm <- vst_all[, samples, drop = FALSE]
coldata <- as.data.frame(SummarizedExperiment::colData(dds_all))[samples, , drop = FALSE]

# Rebuild the time factor on the arm's own timepoints. Taking the whole
# experiment's levels would leave empty levels in the model matrix and DESeq2
# would fail on a rank-deficient design.
coldata$time <- factor(as.numeric(as.character(coldata$time)),
                       levels = sort(unique(as.numeric(as.character(coldata$time)))))
assert_time_levels(coldata$time, sprintf("arm '%s' time", arm))
assert_aligned(coldata, cts, sprintf("arm '%s' colData against counts", arm))

times <- as.numeric(levels(coldata$time))
message(sprintf("[differential] arm '%s': %d libraries, %d timepoints (%s)",
                arm, length(samples), length(times), paste(times, collapse = ", ")))

# --- per-timepoint means and the range statistic ---------------------------
# Computed on the TRANSFORMED scale, which is what differential.min_range is
# expressed against. Reused downstream by the sub-split and the heatmaps, so
# it is stored rather than recomputed.
tp_means <- timepoint_means(vst_arm, as.numeric(as.character(coldata$time)))
rng <- apply(tp_means, 1, function(x) max(x) - min(x))

method <- cfg_req("differential.method")
res <- NULL

if (method == "lrt") {
  # ---- likelihood ratio test: any difference across time at all -----------
  design_f <- stats::as.formula(cfg_req("differential.design"))
  reduced_f <- stats::as.formula(cfg_req("differential.reduced"))
  message(sprintf("[differential] LRT, full %s against reduced %s",
                  deparse(design_f), deparse(reduced_f)))

  dds <- DESeqDataSetFromMatrix(cts, colData = coldata, design = design_f)
  dds <- DESeq(dds, test = "LRT", reduced = reduced_f, quiet = TRUE)
  r <- as.data.frame(results(dds))
  res <- data.frame(feature_id = rownames(r), baseMean = r$baseMean,
                    stat = r$stat, pvalue = r$pvalue, padj = r$padj,
                    stringsAsFactors = FALSE)

} else if (method == "wald_union") {
  # ---- union of per-timepoint contrasts against t=0 -----------------------
  # For users whose question is "when did it change", not "did it change".
  # min_contrasts raises the bar from any single timepoint to a sustained one.
  dds <- DESeqDataSetFromMatrix(cts, colData = coldata, design = ~ time)
  dds <- DESeq(dds, quiet = TRUE)
  ref <- levels(coldata$time)[1]
  others <- levels(coldata$time)[-1]
  fdr <- cfg_req("differential.fdr")

  per <- lapply(others, function(lv) {
    r <- as.data.frame(results(dds, contrast = c("time", lv, ref)))
    data.frame(feature_id = rownames(r), timepoint = lv,
               padj = r$padj, stat = r$stat, baseMean = r$baseMean,
               stringsAsFactors = FALSE)
  })
  per <- do.call(rbind, per)
  per$sig <- !is.na(per$padj) & per$padj < fdr
  n_sig <- tapply(per$sig, per$feature_id, sum)
  best <- tapply(per$padj, per$feature_id,
                 function(p) if (all(is.na(p))) NA_real_ else min(p, na.rm = TRUE))
  maxstat <- tapply(abs(per$stat), per$feature_id,
                    function(s) if (all(is.na(s))) NA_real_ else max(s, na.rm = TRUE))
  ids <- rownames(cts)
  res <- data.frame(
    feature_id = ids,
    baseMean = rowMeans(counts(dds, normalized = TRUE))[ids],
    stat = as.numeric(maxstat[ids]),
    pvalue = NA_real_,
    padj = as.numeric(best[ids]),
    n_contrasts = as.integer(n_sig[ids]),
    stringsAsFactors = FALSE
  )
  message(sprintf("[differential] wald_union: %d contrasts against t=%s", length(others), ref))

} else if (method == "none") {
  # ---- no model: rank by variance of the timepoint means ------------------
  # The supported path for single-replicate timecourses, where no dispersion
  # can be estimated. 06_select_dynamic takes differential.top_n of these.
  res <- data.frame(
    feature_id = rownames(cts),
    baseMean = rowMeans(cts),
    stat = apply(tp_means, 1, stats::var),
    pvalue = NA_real_,
    padj = NA_real_,
    stringsAsFactors = FALSE
  )
  message("[differential] method: none, ranking by variance of timepoint means")

} else {
  stop(sprintf("unknown differential.method '%s'", method), call. = FALSE)
}

cd_assert(identical(res$feature_id, rownames(tp_means)),
          "arm '%s': result rows and timepoint means are not aligned.", arm)
res$range <- as.numeric(rng[res$feature_id])

if (method %in% c("lrt", "wald_union")) {
  n_sig <- sum(!is.na(res$padj) & res$padj < cfg_req("differential.fdr"))
  message(sprintf("[differential] padj < %s: %d of %d features",
                  format(cfg_req("differential.fdr")), n_sig, nrow(res)))
  cd_note("n_significant", n_sig)
}
cd_note("method", method)
cd_note("n_features", nrow(res))
cd_note("timepoints", times)

save_obj(list(res = res, timepoint_means = tp_means, times = times,
              samples = samples, method = method, arm = arm),
         snakemake@output$obj)
cd_finish()
