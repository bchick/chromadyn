# ---------------------------------------------------------------------------
# 07_cluster.R: group the dynamic features by the shape of their trajectory.
#
# cluster.method chooses the clusterer: DEGreport::degPatterns (the default,
# and the source analysis) or k-means. Both see the same representation,
# per-timepoint means z-scored across time, and return the same assignment
# table, so everything downstream is shared.
#
# The core of the method, and the step with the most ways to be silently
# wrong. Four of them are handled explicitly here, the first three specific
# to degPatterns:
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
#   Scale. degPatterns builds an all-pairs correlation matrix and its cost
#   grows with the cube of the feature count. Above cluster.max_features we
#   cluster a stratified subsample and assign the remainder by correlation to
#   cluster centroids, reporting both counts. k-means needs no subsample.
# ---------------------------------------------------------------------------

if (!exists("snakemake")) {
  source("workflow/scripts/lib/standalone.R")
  snakemake <- cd_standalone(
    rule       = "cluster",
    wildcards  = list(arm = "WT"),
    configfile = "config/demo.yaml",
    input      = list(dds = "obj_dds", transformed = "obj_transformed",
                      diff_obj = "obj_differential", results = "differential"),
    output     = list(obj = "obj_cluster", qc = "assignment_qc",
                      selection = "cluster_selection"))
}
source(file.path(snakemake@scriptdir, "lib", "common.R"))
source(file.path(snakemake@scriptdir, "lib", "namespace_fix.R"))
cd_init(snakemake)

method <- cfg_req("cluster.method")
suppressPackageStartupMessages({
  library(DESeq2)
  if (method == "degpatterns") library(DEGreport)
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

# --- profiles used for clustering, centroids and remainder assignment ------
# The same representation degPatterns clusters on: summarize by timepoint,
# then z-score each feature across timepoints. k-means clusters these
# directly.
profiles <- zscore_rows(timepoint_means(mat, times_num))

# Clusters numbered by size, largest first, so that k-means labels do not
# depend on the arbitrary order in which the algorithm found them.
relabel_by_size <- function(x) {
  sizes <- sort(table(x), decreasing = TRUE)
  match(x, as.integer(names(sizes)))
}

# --- guardrail 9.9: scale --------------------------------------------------
# degPatterns cannot finish on a large input, so above max_features it
# clusters a stratified subsample and the rest is assigned by correlation.
# k-means scales linearly and always uses every feature.
max_features <- cfg_req("cluster.max_features")
if (method == "degpatterns" && nrow(mat) > max_features) {
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
deg <- NULL
t0 <- Sys.time()

if (method == "degpatterns") {
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

  deg <- do.call(DEGreport::degPatterns, deg_args)
  k_found <- length(unique(deg$df$cluster))
  message(sprintf(
    "[cluster] degPatterns returned %d cluster(s) over %d features in %.1fs",
    k_found, nrow(deg$df), as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  # --- undo the make.names() mangling --------------------------------------
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

  selection <- data.frame(
    k = k_found,
    criterion = if (is.null(n_clusters)) "diana_cut" else "n_clusters",
    value = NA_real_, tot_withinss = NA_real_, selected = TRUE,
    stringsAsFactors = FALSE
  )
} else if (method == "kmeans") {
  kc <- cfg_req("cluster.kmeans")
  k_lo <- kc$k_range[[1]]
  k_hi <- kc$k_range[[2]]
  cd_assert(k_lo <= k_hi, "cluster.kmeans.k_range must be [low, high] with low <= high.")
  # The whole range is always fitted, so the silhouette sweep is reported even
  # when k is fixed. It is cheap next to degPatterns.
  ks <- sort(unique(c(seq(k_lo, k_hi), kc$k)))
  n_distinct <- nrow(unique(round(profiles, 10)))
  cd_assert(max(ks) < n_distinct,
            paste0("arm '%s': k-means cannot make %d clusters from %d distinct profiles. ",
                   "Lower cluster.kmeans.k or the top of k_range."),
            arm, max(ks), n_distinct)

  seed <- cd_set_seed("kmeans")
  message(sprintf("[cluster] k-means on %d features, k in %d..%d%s, nstart = %d, seed = %d",
                  nrow(profiles), min(ks), max(ks),
                  if (is.null(kc$k)) " (chosen by silhouette)" else sprintf(" (k = %d used)", kc$k),
                  kc$nstart, seed))
  fits <- lapply(ks, function(k) {
    # Reseeded per k, so each fit is the same whether or not the others ran.
    set.seed(seed)
    stats::kmeans(profiles, centers = k, nstart = kc$nstart, iter.max = kc$iter_max)
  })

  # Silhouette on a seeded sample: it needs every pairwise distance.
  cd_set_seed("kmeans_silhouette")
  idx <- if (nrow(profiles) > kc$silhouette_sample) {
    sort(sample(nrow(profiles), kc$silhouette_sample))
  } else {
    seq_len(nrow(profiles))
  }
  dsub <- stats::dist(profiles[idx, , drop = FALSE])
  sil <- vapply(fits, function(f) {
    lab <- f$cluster[idx]
    if (length(unique(lab)) < 2L) return(NA_real_)
    mean(cluster::silhouette(lab, dsub)[, "sil_width"])
  }, numeric(1))

  if (is.null(kc$k)) {
    cd_assert(any(is.finite(sil)),
              "arm '%s': no k in cluster.kmeans.k_range gave a usable silhouette.", arm)
    best <- which.max(sil)
    if (ks[best] == k_hi && k_lo < k_hi) {
      warning(sprintf(paste0("arm '%s': silhouette chose k = %d, the top of ",
                             "cluster.kmeans.k_range, so a larger k may fit better. Consider ",
                             "widening the range."), arm, ks[best]), call. = FALSE)
    }
  } else {
    best <- match(kc$k, ks)
  }
  fit <- fits[[best]]
  k_found <- ks[best]
  message(sprintf("[cluster] k-means used k = %d (mean silhouette %.3f on %d features) in %.1fs",
                  k_found, sil[best], length(idx),
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  assigned <- data.frame(
    feature_id = rownames(profiles),
    cluster = relabel_by_size(unname(fit$cluster)),
    assign_method = "core",
    assign_cor = NA_real_,
    stringsAsFactors = FALSE
  )
  selection <- data.frame(
    k = ks, criterion = "mean_silhouette", value = round(sil, 4),
    tot_withinss = round(vapply(fits, function(f) f$tot.withinss, numeric(1)), 1),
    selected = seq_along(ks) == best,
    stringsAsFactors = FALSE
  )
}
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# --- assign the held-out remainder by correlation --------------------------
# degPatterns on the subsample path only.
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
# DEGreport's own semantics, applied to every method: a cluster smaller than
# minc does not survive, and its members are dropped rather than merged
# elsewhere. Marking them Unassigned keeps them in the output, which is what
# lets a user see that a feature was clustered and then discarded instead of
# never appearing.
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
             "k_selected", "n_clusters", "minc_requested", "minc_core",
             "minc_effective", "max_features", "assign_min_cor", "seed", "cluster_seconds"),
  value = c(length(dynamic_ids), n_zero_var, nrow(assigned), length(core), length(rest),
            sum(assigned$assign_method == "correlation"),
            sum(assigned$assign_method == "unassigned_low_cor"),
            n_minc_dropped, k_found, n_final_clusters,
            if (identical(minc_cfg, "auto")) NA_integer_ else as.integer(minc_cfg),
            minc_core, minc, max_features, cfg_req("cluster.assign_min_cor"), seed,
            round(elapsed, 1)),
  stringsAsFactors = FALSE
)
write_tsv_atomic(qc, snakemake@output$qc)
for (i in seq_len(nrow(qc))) cd_note(qc$metric[i], qc$value[i])
cd_note("cluster_method", method)

write_tsv_atomic(cbind(method = method, selection), snakemake@output$selection)

save_obj(list(
  method = method,
  deg = deg,                       # full DEGreport object when method is degpatterns, else NULL
  id_key = if (method == "degpatterns") id_key else NULL,  # make.names() -> original
  assignment = assigned,           # feature_id, cluster, assign_method, assign_cor
  selection = selection,           # how the cluster count was reached
  profiles = profiles,             # z-scored timepoint means, ALL clustered features
  centroids = centroids,
  meta = meta,
  times = d$times,
  timepoint_means = d$timepoint_means,
  samples = samples,
  minc_effective = minc,
  seed = seed,
  arm = arm,
  degreport_version = if (method == "degpatterns") {
    as.character(utils::packageVersion("DEGreport"))
  } else {
    NA_character_
  }
), snakemake@output$obj)
cd_finish()
