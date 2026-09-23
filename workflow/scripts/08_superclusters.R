# ---------------------------------------------------------------------------
# 08_superclusters.R: collapse many clusters into a few named trajectory
# classes, and optionally split one of them.
#
# In the source analysis this was done by eye: a dendrogram was plotted, a cut
# was chosen by squinting at it, and fourteen cluster numbers were typed into
# a tribble() by hand. That map is the reason the original result could not
# be reproduced from the code alone. Here the cut is made by a stated rule,
# the dendrogram and a k-selection diagnostic are emitted so the rule can be
# checked, and an explicit override map is still available for reproducing a
# historical hand assignment exactly.
#
# Two stages, because one broad class usually contains several real shapes:
# an "Increasing" group is often a mixture of transient, sustained and late
# rises that only separate when you split it again.
# ---------------------------------------------------------------------------

if (!exists("snakemake")) {
  source("workflow/scripts/lib/standalone.R")
  snakemake <- cd_standalone(
    rule       = "superclusters",
    wildcards  = list(arm = "WT"),
    configfile = "config/demo.yaml",
    input      = list(obj = "obj_cluster", features = "demo/features.bed"),
    output     = list(clusters = "clusters", profiles = "cluster_profiles",
                      sizes = "supercluster_sizes", kdiag = "kdiag",
                      fit = "obj_supercluster"))
}
source(file.path(snakemake@scriptdir, "lib", "common.R"))
source(file.path(snakemake@scriptdir, "lib", "namespace_fix.R"))
source(file.path(snakemake@scriptdir, "lib", "classify.R"))
cd_init(snakemake)

arm <- snakemake@wildcards$arm
cl <- load_obj(snakemake@input$obj)
assignment <- cl$assignment
profiles <- cl$profiles
times <- as.numeric(colnames(profiles))

clustered <- assignment[!is.na(assignment$cluster), , drop = FALSE]
cd_assert(nrow(clustered) > 0, "arm '%s': no features survived clustering.", arm)
cluster_ids <- sort(unique(clustered$cluster))

# --- cluster mean profiles -------------------------------------------------
# Computed from the assignment rather than from deg$normalized, so that
# features assigned by correlation on the subsample path are included.
cluster_profiles <- t(vapply(
  cluster_ids,
  function(k) colMeans(profiles[clustered$feature_id[clustered$cluster == k], , drop = FALSE]),
  numeric(ncol(profiles))))
rownames(cluster_profiles) <- as.character(cluster_ids)
colnames(cluster_profiles) <- colnames(profiles)

# --- the dendrogram, and evidence for where to cut it ----------------------
# Emitted whatever the assignment method is. The point of the diagnostic is
# that you can see whether the cut you took was defensible, including when
# the cut was not what produced the assignment.
dd <- stats::dist(cluster_profiles)
hc <- stats::hclust(dd, method = cfg("superclusters.linkage", "complete"))

kdiag <- NULL
if (length(cluster_ids) >= 3L) {
  ks <- 2:min(length(cluster_ids) - 1L, 10L)
  kdiag <- do.call(rbind, lapply(ks, function(k) {
    memb <- stats::cutree(hc, k = k)
    sil <- if (requireNamespace("cluster", quietly = TRUE) && k < length(cluster_ids)) {
      mean(cluster::silhouette(memb, dd)[, "sil_width"])
    } else NA_real_
    wss <- sum(vapply(split(seq_len(nrow(cluster_profiles)), memb), function(i) {
      sub <- cluster_profiles[i, , drop = FALSE]
      sum((sub - rep(colMeans(sub), each = nrow(sub)))^2)
    }, numeric(1)))
    data.frame(k = k, silhouette = sil, wss = wss)
  }))
  best <- kdiag$k[which.max(kdiag$silhouette)]
  message(sprintf(
    "[superclusters] k diagnostics over %d clusters: best silhouette at k = %d (%.3f)",
    length(cluster_ids), best, max(kdiag$silhouette, na.rm = TRUE)))
} else {
  kdiag <- data.frame(k = integer(), silhouette = numeric(), wss = numeric())
  message("[superclusters] fewer than 3 clusters, k diagnostics skipped")
}

# --- stage 1 ---------------------------------------------------------------
method <- cfg_req("superclusters.method")
k <- cfg("superclusters.k")
h <- cfg("superclusters.h")
auto_label <- cfg("superclusters.auto_label", "shape")
split_cfg <- cfg("superclusters.split")

# Stage 1 and stage 2 have independent thresholds. Sharing them looks tidy
# and is wrong: a `cross` margin tuned to separate a late rise from a
# sustained one also stops stage 1 recognising a transient, because the same
# comparison means something different when the candidates are all rising.
shape_cfg <- cfg("superclusters.shape")
stage1_params <- list(
  mid_fraction = shape_cfg$mid_fraction,
  cross = shape_cfg$cross
)
split_params <- list(
  late_ratio_cutoff = split_cfg$late_ratio_cutoff,
  late_peak_fraction = split_cfg$late_peak_fraction,
  mid_time = split_cfg$mid_time,
  mid_fraction = split_cfg$mid_fraction,
  cross = split_cfg$cross
)

cluster_to_sc <- NULL

if (method == "overrides") {
  ov <- cfg_req("superclusters.overrides")
  # Guardrail 9.5. A join against an incomplete map silently drops the
  # unmapped clusters; the source analysis did exactly that.
  assert_map_covers(ov, cluster_ids, "superclusters.overrides")
  cluster_to_sc <- stats::setNames(unlist(ov[as.character(cluster_ids)]),
                                   as.character(cluster_ids))
  message(sprintf("[superclusters] overrides map: %d clusters onto %d classes",
                  length(cluster_ids), length(unique(cluster_to_sc))))

} else if (method == "shape") {
  # Label each cluster from its OWN centroid, then let clusters sharing a
  # label form a class. The reference classes this reproduces are a
  # per-feature shape classification, not a property of dendrogram structure,
  # and cutting the tree instead folds Transient into Decreasing whenever the
  # two branches happen to be close.
  lab <- classify_centroids(cluster_profiles, scheme = auto_label, vocab = "stage1",
                            params = stage1_params, times = times)
  cluster_to_sc <- stats::setNames(lab$label, rownames(cluster_profiles))
  message(sprintf("[superclusters] shape labelling: %s",
                  paste(sprintf("%s=%d", names(table(cluster_to_sc)), table(cluster_to_sc)),
                        collapse = ", ")))

} else if (method %in% c("hclust", "kmeans")) {
  memb <- if (method == "hclust") {
    if (!is.null(h)) stats::cutree(hc, h = h) else stats::cutree(hc, k = k)
  } else {
    cd_set_seed("supercluster_kmeans")
    stats::kmeans(cluster_profiles, centers = k, nstart = 25, iter.max = 100)$cluster
  }
  grp_cent <- t(vapply(split(seq_len(nrow(cluster_profiles)), memb), function(i) {
    colMeans(cluster_profiles[i, , drop = FALSE])
  }, numeric(ncol(cluster_profiles))))
  colnames(grp_cent) <- colnames(cluster_profiles)
  lab <- classify_centroids(grp_cent, scheme = auto_label, vocab = "stage1",
                            params = stage1_params, times = times)
  names_by_group <- stats::setNames(lab$label, rownames(grp_cent))
  cluster_to_sc <- stats::setNames(names_by_group[as.character(memb)],
                                   rownames(cluster_profiles))
  message(sprintf("[superclusters] %s cut into %d group(s): %s",
                  method, length(unique(memb)), paste(unique(lab$label), collapse = ", ")))
} else {
  stop(sprintf("unknown superclusters.method '%s'", method), call. = FALSE)
}

assignment$supercluster <- NA_character_
idx <- !is.na(assignment$cluster)
assignment$supercluster[idx] <- unname(cluster_to_sc[as.character(assignment$cluster[idx])])
assignment$supercluster[is.na(assignment$supercluster)] <- "Unassigned"
assignment$supercluster_label <- assignment$supercluster

# --- stage 2: split one class ----------------------------------------------
split_done <- FALSE
if (!is.null(split_cfg) && !is.null(split_cfg$target) && method != "overrides") {
  target <- split_cfg$target
  in_target <- which(assignment$supercluster == target)
  ksplit <- as.integer(split_cfg$k)
  if (length(in_target) == 0L) {
    message(sprintf("[superclusters] split target '%s' has no features, skipping stage 2", target))
  } else if (length(in_target) < ksplit * 2L) {
    message(sprintf(
      "[superclusters] split target '%s' has only %d features for k = %d, skipping",
      target, length(in_target), ksplit))
  } else {
    ids <- assignment$feature_id[in_target]
    z <- profiles[ids, , drop = FALSE]
    seed <- cd_set_seed("supercluster_split")
    km <- stats::kmeans(z, centers = ksplit, nstart = 25, iter.max = 100)
    cen <- km$centers
    colnames(cen) <- colnames(profiles)
    lab <- classify_centroids(cen, scheme = auto_label, vocab = "split",
                              params = split_params, prefix = target, times = times)
    assignment$supercluster_label[in_target] <- lab$label[km$cluster]
    split_done <- TRUE
    message(sprintf("[superclusters] split '%s' (%d features) into %d: %s (seed %d)",
                    target, length(ids), ksplit,
                    paste(sprintf("%s=%d", lab$label, as.integer(table(km$cluster))),
                          collapse = ", "), seed))
  }
}

# --- coordinates, carried as real columns ----------------------------------
features_path <- cfg("input.features")
if (!is.null(features_path)) {
  bed <- read_bed(features_path, feature_ids = assignment$feature_id)
  m <- match(assignment$feature_id, bed$name)
  assignment$chr <- bed$chr[m]
  assignment$start <- bed$start[m]
  assignment$end <- bed$end[m]
  # Guardrail 9.4: never write coordinates that were silently lost on a join.
  assert_coords(assignment[!is.na(assignment$chr), ], "cluster table coordinates")
  cd_assert(!anyNA(assignment$chr),
            "%d feature(s) have no BED record after the join; coordinates would be NA.",
            sum(is.na(assignment$chr)))
}

# --- write -----------------------------------------------------------------
cols <- c("feature_id", "cluster", "supercluster", "supercluster_label",
          "assign_method", "assign_cor",
          if (!is.null(features_path)) c("chr", "start", "end"))
out <- assignment[order(assignment$supercluster_label, assignment$cluster,
                        assignment$feature_id), cols]
write_tsv_atomic(out, snakemake@output$clusters)

prof_long <- do.call(rbind, lapply(seq_along(cluster_ids), function(i) {
  kk <- cluster_ids[i]
  data.frame(cluster = kk, time = times, mean_z = as.numeric(cluster_profiles[i, ]),
             n = sum(clustered$cluster == kk),
             supercluster = unname(cluster_to_sc[as.character(kk)]),
             stringsAsFactors = FALSE)
}))
write_tsv_atomic(prof_long, snakemake@output$profiles)

sizes <- as.data.frame(table(out$supercluster_label), stringsAsFactors = FALSE)
colnames(sizes) <- c("supercluster_label", "n")
sizes <- sizes[order(-sizes$n), ]
write_tsv_atomic(sizes, snakemake@output$sizes)
write_tsv_atomic(kdiag, snakemake@output$kdiag)

message(sprintf("[superclusters] arm '%s': %s",
                arm, paste(sprintf("%s=%d", sizes$supercluster_label, sizes$n), collapse = ", ")))
cd_note("method", method)
cd_note("n_clusters", length(cluster_ids))
cd_note("n_superclusters", nrow(sizes))
cd_note("split_applied", split_done)
cd_note("supercluster_sizes", stats::setNames(as.list(sizes$n), sizes$supercluster_label))

# Cached for the figures, so the dendrogram is drawn from the same tree the
# diagnostics were computed on rather than a second, possibly different one.
save_obj(list(hc = hc, dist = dd, cluster_profiles = cluster_profiles,
              cluster_to_sc = cluster_to_sc, kdiag = kdiag, times = times,
              assignment = assignment, arm = arm, split_applied = split_done),
         snakemake@output$fit)
cd_finish()
