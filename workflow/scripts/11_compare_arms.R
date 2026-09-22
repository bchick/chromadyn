# ---------------------------------------------------------------------------
# 11_compare_arms.R: which class does a feature land in under each arm?
#
# In the source project this question ("which features that are Transient
# under one treatment are Increasing under the other") was answered by reading
# two tables side by side and typing the counts into a third. It is a join, so
# it is a rule.
#
# Degrades to a single-arm summary rather than failing, because a one-arm
# timecourse is a perfectly ordinary input and `rule all` should not depend
# on having something to compare.
# ---------------------------------------------------------------------------

if (!exists("snakemake")) {
  source("workflow/scripts/lib/standalone.R")
  snakemake <- cd_standalone(
    rule       = "compare_arms",
    configfile = "config/demo.yaml",
    input      = list(clusters = "results/clusters/WT_clusters.tsv"),
    output     = list(crosstab = "crosstab"))
}
source(file.path(snakemake@scriptdir, "lib", "common.R"))
source(file.path(snakemake@scriptdir, "lib", "namespace_fix.R"))
cd_init(snakemake)

paths <- unlist(snakemake@input$clusters)
arms <- unlist(snakemake@params$arms)
if (is.null(arms)) arms <- sub("_clusters\\.tsv$", "", basename(paths))

tabs <- stats::setNames(lapply(paths, function(p) {
  read_tsv_strict(p, c("feature_id", "supercluster_label"))
}), arms)

if (length(arms) == 1L) {
  # One arm: emit the class sizes in the same shape a crosstab would take, so
  # that a downstream consumer does not need to special-case it.
  t1 <- tabs[[1]]
  out <- as.data.frame(table(t1$supercluster_label), stringsAsFactors = FALSE)
  colnames(out) <- c(arms[1], "n")
  out <- out[order(-out$n), ]
  message(sprintf("[compare_arms] single arm '%s': %d classes, no cross-arm comparison possible",
                  arms[1], nrow(out)))
  cd_note("n_arms", 1L)
} else {
  # Pairwise over every arm combination, long format so any number of arms
  # works without reshaping.
  combos <- utils::combn(arms, 2L, simplify = FALSE)
  rows <- lapply(combos, function(pair) {
    a <- tabs[[pair[1]]]
    b <- tabs[[pair[2]]]
    m <- merge(a[, c("feature_id", "supercluster_label")],
               b[, c("feature_id", "supercluster_label")],
               by = "feature_id", suffixes = c("_a", "_b"))
    if (!nrow(m)) return(NULL)
    tt <- as.data.frame(table(m$supercluster_label_a, m$supercluster_label_b),
                        stringsAsFactors = FALSE)
    colnames(tt) <- c("class_a", "class_b", "n")
    tt$arm_a <- pair[1]
    tt$arm_b <- pair[2]
    tt[tt$n > 0, c("arm_a", "class_a", "arm_b", "class_b", "n")]
  })
  out <- do.call(rbind, rows)
  cd_assert(!is.null(out) && nrow(out) > 0,
            "no features are shared between arms, so there is nothing to cross-tabulate.")
  out <- out[order(out$arm_a, out$arm_b, -out$n), ]
  message(sprintf("[compare_arms] %d arms, %d pairwise comparison(s), %d non-empty cells",
                  length(arms), length(combos), nrow(out)))
  cd_note("n_arms", length(arms))
}

write_tsv_atomic(out, snakemake@output$crosstab)
cd_finish()
