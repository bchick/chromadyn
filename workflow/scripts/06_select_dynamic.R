# ---------------------------------------------------------------------------
# 06_select_dynamic.R: two gates, not one.
#
# A feature is dynamic when it passes BOTH a significance test and an
# effect-size floor on the range of its per-timepoint means. With enough
# replication the test alone returns thousands of features whose total
# excursion over the whole timecourse is a few percent, and clustering those
# produces confident-looking classes built on noise.
#
# The floor is on the TRANSFORMED scale, not a log2 fold change. See
# docs/method.md for why the source parameter name was misleading.
# ---------------------------------------------------------------------------

if (!exists("snakemake")) {
  source("workflow/scripts/lib/standalone.R")
  snakemake <- cd_standalone(
    rule       = "select_dynamic",
    wildcards  = list(arm = "WT"),
    configfile = "config/demo.yaml",
    input      = list(obj = "obj_differential"),
    output     = list(results = "differential"))
}
source(file.path(snakemake@scriptdir, "lib", "common.R"))
source(file.path(snakemake@scriptdir, "lib", "namespace_fix.R"))
cd_init(snakemake)

arm <- snakemake@wildcards$arm
d <- load_obj(snakemake@input$obj)
res <- d$res
method <- d$method

min_range <- cfg_req("differential.min_range")
fdr <- cfg_req("differential.fdr")

if (method == "lrt") {
  gate_sig <- !is.na(res$padj) & res$padj < fdr
  sig_label <- sprintf("padj < %s", format(fdr))

} else if (method == "wald_union") {
  min_contrasts <- cfg_req("differential.min_contrasts")
  gate_sig <- !is.na(res$n_contrasts) & res$n_contrasts >= min_contrasts
  sig_label <- sprintf("significant in >= %d contrast(s)", min_contrasts)

} else {
  # method: none. Rank by the variance statistic and take the top N; the range
  # gate still applies on top, so top_n is a ceiling and not a quota.
  top_n <- cfg_req("differential.top_n")
  ord <- order(res$stat, decreasing = TRUE, na.last = TRUE)
  gate_sig <- logical(nrow(res))
  gate_sig[utils::head(ord, top_n)] <- TRUE
  sig_label <- sprintf("top %d by variance", top_n)
}

gate_range <- !is.na(res$range) & res$range >= min_range
res$dynamic <- gate_sig & gate_range

n <- nrow(res)
message(sprintf(
  "[select_dynamic] arm '%s': %s -> %d | range >= %s -> %d | both -> %d of %d (%.1f%%)",
  arm, sig_label, sum(gate_sig), format(min_range), sum(gate_range),
  sum(res$dynamic), n, 100 * sum(res$dynamic) / n))

# A gate that removes nothing is not a gate. Say so rather than letting the
# user believe an effect-size filter was applied when it was not.
only_sig <- sum(gate_sig & !gate_range)
if (sum(gate_sig) > 0 && only_sig == 0) {
  message(sprintf(paste0("[select_dynamic] note: the range gate removed nothing. Every ",
                         "feature passing %s already had range >= %s, so min_range is ",
                         "not currently constraining the result."),
                  sig_label, format(min_range)))
}
cd_assert(sum(res$dynamic) > 0,
          paste0("arm '%s': no features passed both gates (%s gave %d, range >= %s gave %d). ",
                 "Loosen differential.fdr or differential.min_range."),
          arm, sig_label, sum(gate_sig), format(min_range), sum(gate_range))

cd_note("n_gate_significant", sum(gate_sig))
cd_note("n_gate_range", sum(gate_range))
cd_note("n_dynamic", sum(res$dynamic))
cd_note("min_range", min_range)

cols <- c("feature_id", "baseMean", "stat", "pvalue", "padj",
          if ("n_contrasts" %in% colnames(res)) "n_contrasts", "range", "dynamic")
write_tsv_atomic(res[, cols], snakemake@output$results)
cd_finish()
