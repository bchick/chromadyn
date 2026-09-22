# ---------------------------------------------------------------------------
# 07_cluster.R: group the dynamic features by the shape of their trajectory.
#
# The core of the method, and the step with the most ways to be silently
# wrong. Four of them are handled explicitly here:
#
#   Metadata alignment. degPatterns attaches metadata to columns positionally.
#   rownames(metadata) must equal colnames(matrix) in the same ORDER, not just
#   as a set, or every feature is clustered against the wrong timepoints and
#   nothing errors. assert_aligned() checks identity, not setequal.
#
#   Feature ID mangling. degPatterns runs make.names() over the rownames, so
#   "chr1:7401731-7402231" comes back as "chr1.7401731.7402231". The source
#   analysis parsed IDs back out with a regex that did not match afterwards
#   and wrote thousands of BED rows reading "chr1.7401731.7402231 NA NA". We
#   keep an explicit key table and rejoin through it; nothing downstream ever
#   parses an ID.
#
#   minc drops, it does not merge. DEGreport keeps clusters with size strictly
#   greater than minc and REMOVES the members of the rest. Those features are
#   reported as Unassigned rather than quietly vanishing from the output.
#
#   Scale. degPatterns builds an all-pairs correlation matrix, so cost grows
#   quadratically. Above cluster.max_features we cluster a stratified
#   subsample and assign the remainder by correlation to cluster centroids,
#   reporting both counts.
# ---------------------------------------------------------------------------

if (!exists("snakemake")) {
  source("workflow/scripts/lib/standalone.R")
  snakemake <- cd_standalone(
    rule       = "cluster",
    wildcards  = list(arm = "WT"),
    configfile = "config/demo.yaml",
    input      = list(dds = "obj_dds", transformed = "obj_transformed",
                      diff_obj = "obj_differential", results = "differential"),
    output     = list(obj = "obj_degpatterns", qc = "assignment_qc"))
}
source(file.path(snakemake@scriptdir, "lib", "common.R"))
source(file.path(snakemake@scriptdir, "lib", "namespace_fix.R"))
cd_init(snakemake)

suppressPackageStartupMessages({
  library(DESeq2)
  library(DEGreport)
})

arm <- snakemake@wildcards$arm
d <- load_obj(snakemake@input$diff_obj)
vst_all <- load_obj(snakemake@input$transformed)
dds_all <- load_obj(snakemake@input$dds)
res <- read_tsv_strict(snakemake@input$results, c("feature_id", "dynamic"))

dynamic_ids <- res$feature_id[as.logical(res$dynamic)]
cd_assert(length(dynamic_ids) > 0, "arm '%s': nothing to cluster.", arm)

samples <- d$samples
mat <- vst_all[dynamic_ids, samples, drop = FALSE]
coldata <- as.data.frame(SummarizedExperiment::colData(dds_all))[samples, , drop = FALSE]
times_num <- as.numeric(as.character(coldata$time))

# --- guardrail 9.3: zero-variance rows -------------------------------------
# degPatterns correlates every pair of features; a constant row gives NA and
# poisons the distance matrix.
fvar <- apply(mat, 1, stats::var, na.rm = TRUE)
n_zero_var <- sum(!is.finite(fvar) | fvar <= 0)
if (n_zero_var > 0) {
  message(sprintf("[cluster] dropping %d zero-variance feature(s)", n_zero_var))
  mat <- mat[is.finite(fvar) & fvar > 0, , drop = FALSE]
}
cd_assert(nrow(mat) > 0, "arm '%s': every dynamic feature had zero variance.", arm)

# --- metadata, built the way degPatterns needs it --------------------------
# Time must be a factor whose levels are in numeric order, and the sample-ID
# column must not be present: degPatterns would treat it as a grouping
# variable. Rownames carry the sample identity instead.
meta <- data.frame(
  time = factor(times_num, levels = sort(unique(times_num))),
  replicate = factor(coldata$replicate),
  row.names = rownames(coldata)
)
assert_time_levels(meta$time, sprintf("arm '%s' clustering metadata", arm))
assert_aligned(meta, mat, sprintf("arm '%s' degPatterns input", arm))

# --- minc ------------------------------------------------------------------
minc_cfg <- cfg_req("cluster.degpatterns.minc")
minc <- if (identical(minc_cfg, "auto")) max(15L, round(nrow(mat) / 100)) else as.integer(minc_cfg)
if (identical(minc_cfg, "auto")) {
  message(sprintf("[cluster] minc: auto resolves to %d for %d features", minc, nrow(mat)))
}
# Guardrail 9.8 again, now against the number actually reaching the clusterer.
cd_assert(minc <= nrow(mat) / 20,
          paste0("arm '%s': cluster.degpatterns.minc is %d but only %d features reach ",
                 "clustering, so a cluster would need %.0f%% of them. Try minc: %d, ",
                 "or minc: auto."),
          arm, minc, nrow(mat), 100 * minc / nrow(mat), max(15L, round(nrow(mat) / 100)))

# `cutoff` is inert in DEGreport: it is passed to .select_genes(), which
# declares the argument and never reads it, and .reduce() does not accept it.
# Verified in 1.36.0 and 1.42.0. Say so rather than let the user believe a
# threshold they set is doing something.
cutoff <- cfg("cluster.degpatterns.cutoff")
if (!isTRUE(all.equal(cutoff, 0.5))) {
  warning(sprintf(paste0("cluster.degpatterns.cutoff is set to %s, but this parameter has no ",
                         "effect in DEGreport %s: degPatterns passes it to an internal ",
                         "function that ignores it. Use cluster.degpatterns.n_clusters to ",
                         "control cluster count. See docs/parameters.md."),
                  format(cutoff), as.character(utils::packageVersion("DEGreport"))),
          call. = FALSE)
}

# --- profiles used for centroids and remainder assignment ------------------
# The same representation degPatterns clusters on: summarize by timepoint,
# then z-score each feature across timepoints.
profiles <- zscore_rows(timepoint_means(mat, times_num))

# --- guardrail 9.9: scale --------------------------------------------------
max_features <- cfg_req("cluster.max_features")
if (nrow(mat) > max_features) {
  cd_set_seed("cluster_subsample")
  rng <- apply(profiles, 1, function(x) max(x) - min(x))
  # Stratify by effect size so a rare shape is not lost to a random draw.
  strata <- cut(rng, stats::quantile(rng, probs = seq(0, 1, 0.1)), include.lowest = TRUE)
  per <- split(rownames(mat), strata)
  want <- round(max_features * lengths(per) / nrow(mat))
  core <- unlist(Map(function(g, k) sample(g, min(k, length(g))), per, want), use.names = FALSE)
  rest <- setdiff(rownames(mat), core)
  minc_core <- max(15L, round(minc * length(core) / nrow(mat)))
  message(sprintf(paste0("[cluster] %d features exceeds cluster.max_features (%d). Clustering a ",
                         "stratified subsample of %d (minc scaled to %d); the remaining %d will ",
                         "be assigned by correlation to cluster centroids."),
                  nrow(mat), max_features, length(core), minc_core, length(rest)))
} else {
  core <- rownames(mat)
  rest <- character(0)
  minc_core <- minc
}

# --- cluster ----------------------------------------------------------------
# Guardrail 9.6: seed immediately before the stochastic call, and record it.
# The source project seeded its k-means but not degPatterns, while claiming in
# its methods that all modalities were seeded.
seed <- cd_set_seed("degpatterns")
message(sprintf("[cluster] degPatterns on %d features, minc = %d, seed = %d",
                length(core), minc_core, seed))

deg_args <- list(
  ma = mat[core, , drop = FALSE],
  metadata = meta,
  time = "time",
  col = cfg("cluster.degpatterns.col"),
  minc = minc_core,
  summarize = cfg("cluster.degpatterns.summarize"),
  consensusCluster = cfg("cluster.degpatterns.consensus_cluster"),
  reduce = cfg("cluster.degpatterns.reduce"),
  scale = cfg("cluster.degpatterns.scale"),
  cutoff = cutoff,
  plot = TRUE
)
# nClusters and skipDendrogram exist from DEGreport 1.40. Pass them only when
# the installed version accepts them, so the pipeline still runs on 1.36.
known <- names(formals(DEGreport::degPatterns))
n_clusters <- cfg("cluster.degpatterns.n_clusters")
if (!is.null(n_clusters)) {
  cd_assert("nClusters" %in% known,
            paste0("cluster.degpatterns.n_clusters is set but DEGreport %s does not support ",
                   "it (added in 1.40). Leave it null or upgrade."),
            as.character(utils::packageVersion("DEGreport")))
  deg_args$nClusters <- as.integer(n_clusters)
}
if ("skipDendrogram" %in% known) {
  deg_args$skipDendrogram <- cfg("cluster.degpatterns.skip_dendrogram", TRUE)
}

t0 <- Sys.time()
deg <- do.call(DEGreport::degPatterns, deg_args)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
message(sprintf(
  "[cluster] degPatterns returned %d cluster(s) over %d features in %.1fs",
  length(unique(deg$df$cluster)), nrow(deg$df), elapsed))

# --- undo the make.names() mangling ----------------------------------------
mangled <- make.names(rownames(mat))
cd_assert(!anyDuplicated(mangled),
          paste0("arm '%s': make.names() maps two different feature IDs onto the same ",
                 "string, so degPatterns results cannot be mapped back unambiguously. ",
                 "First collision: %s"),
          arm, {
            dup <- duplicated(mangled) | duplicated(mangled, fromLast = TRUE)
            paste(utils::head(rownames(mat)[dup], 2), collapse = " and ")
          })
id_key <- stats::setNames(rownames(mat), mangled)

assigned <- data.frame(
  feature_id = unname(id_key[as.character(deg$df$genes)]),
  cluster = as.integer(deg$df$cluster),
  assign_method = "core",
  assign_cor = NA_real_,
  stringsAsFactors = FALSE
)
assert_no_na(assigned, "feature_id", "cluster assignment after un-mangling IDs")

# --- assign the held-out remainder by correlation --------------------------
centroids <- NULL
if (length(rest)) {
  cl <- split(assigned$feature_id, assigned$cluster)
  centroids <- t(vapply(cl, function(ids) colMeans(profiles[ids, , drop = FALSE]),
                        numeric(ncol(profiles))))
  r <- stats::cor(t(profiles[rest, , drop = FALSE]), t(centroids))
  hit <- max.col(r, ties.method = "first")
  conf <- r[cbind(seq_len(nrow(r)), hit)]
  min_cor <- cfg_req("cluster.assign_min_cor")
  ok <- conf >= min_cor
  assigned <- rbind(assigned, data.frame(
    feature_id = rest,
    cluster = ifelse(ok, as.integer(rownames(centroids))[hit], NA_integer_),
    assign_method = ifelse(ok, "correlation", "unassigned_low_cor"),
    assign_cor = conf,
    stringsAsFactors = FALSE
  ))
  message(sprintf(
    "[cluster] assigned %d of %d held-out by correlation (r >= %s); %d left unassigned",
    sum(ok), length(rest), format(min_cor), sum(!ok)))
}

# --- re-apply minc to the FINAL sizes --------------------------------------
# DEGreport's own semantics: a cluster smaller than minc does not survive, and
# its members are dropped rather than merged elsewhere. Marking them
# Unassigned keeps them in the output, which is what lets a user see that a
# feature was clustered and then discarded instead of never appearing.
sizes <- table(assigned$cluster[!is.na(assigned$cluster)])
too_small <- names(sizes)[sizes <= minc]
n_minc_dropped <- 0L
if (length(too_small)) {
  drop <- !is.na(assigned$cluster) & as.character(assigned$cluster) %in% too_small
  n_minc_dropped <- sum(drop)
  assigned$cluster[drop] <- NA_integer_
  assigned$assign_method[drop] <- "unassigned_minc"
  message(sprintf(
    "[cluster] %d feature(s) in %d cluster(s) below the final minc of %d marked Unassigned",
    n_minc_dropped, length(too_small), minc))
}

n_final_clusters <- length(unique(assigned$cluster[!is.na(assigned$cluster)]))
message(sprintf("[cluster] arm '%s': %d clusters, %d assigned, %d unassigned",
                arm, n_final_clusters, sum(!is.na(assigned$cluster)), sum(is.na(assigned$cluster))))

qc <- data.frame(
  metric = c("n_dynamic_in", "n_zero_variance_dropped", "n_clustered", "n_core", "n_remainder",
             "n_assigned_by_correlation", "n_unassigned_low_cor", "n_unassigned_minc",
             "n_clusters", "minc_requested", "minc_core", "minc_effective",
             "max_features", "assign_min_cor", "seed", "degpatterns_seconds"),
  value = c(length(dynamic_ids), n_zero_var, nrow(assigned), length(core), length(rest),
            sum(assigned$assign_method == "correlation"),
            sum(assigned$assign_method == "unassigned_low_cor"),
            n_minc_dropped, n_final_clusters,
            if (identical(minc_cfg, "auto")) NA_integer_ else as.integer(minc_cfg),
            minc_core, minc, max_features, cfg_req("cluster.assign_min_cor"), seed,
            round(elapsed, 1)),
  stringsAsFactors = FALSE
)
write_tsv_atomic(qc, snakemake@output$qc)
for (i in seq_len(nrow(qc))) cd_note(qc$metric[i], qc$value[i])

save_obj(list(
  deg = deg,                       # full DEGreport object, including $plot and $normalized
  id_key = id_key,                 # make.names() -> original. Never re-parse an ID.
  assignment = assigned,           # feature_id, cluster, assign_method, assign_cor
  profiles = profiles,             # z-scored timepoint means, ALL clustered features
  centroids = centroids,
  meta = meta,
  times = d$times,
  timepoint_means = d$timepoint_means,
  samples = samples,
  minc_effective = minc,
  seed = seed,
  arm = arm,
  degreport_version = as.character(utils::packageVersion("DEGreport"))
), snakemake@output$obj)
cd_finish()
