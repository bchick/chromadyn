# ---------------------------------------------------------------------------
# plots.R: every figure timecourse-patterns draws.
#
# These exist as functions because in the source project each one was written
# inline, two or three times, in different notebooks, and the copies drifted.
# Nothing here reads the config or writes a file: a plot function takes data
# and returns a ggplot, and the calling script decides where it goes.
#
# Ported and generalized from the source notebooks. The generalizations are
# the point: the originals hard-coded one experiment's timepoint grid, one
# experiment's treatment colours, and an assumption of exactly two replicates.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
})

# plots.R depends on lib/theme.R for the palettes (cat_palette, cor_palette,
# col_accent) and for theme_publication(). It does not source it itself,
# because theme.R must be sourced exactly once per session and the scripts
# already do that. Fail loudly rather than at the first undefined symbol
# halfway through drawing a figure.
if (!exists("theme_publication", mode = "function")) {
  stop("lib/plots.R requires lib/theme.R. Source theme.R first.", call. = FALSE)
}

# --- QC ---------------------------------------------------------------------

#' PCA of samples, coloured by group and shaped by replicate.
#'
#' Ported from 00_preprocessing.Rmd:416-435. That version fixed the time factor
#' to c("0m","30m","60m","120m","240m") and the colour scale to this project's
#' three treatments, so it could not be reused. Here the grid comes from the data
#' and the palette from the group count.
#'
#' @param mat Transformed matrix, features x samples.
#' @param coldata data.frame with rownames == colnames(mat) and columns
#'   group, time, replicate.
#' @return list(df, var_exp, plot)
run_pca <- function(mat, coldata, label = NULL) {
  stopifnot(identical(rownames(coldata), colnames(mat)))
  pca <- stats::prcomp(t(mat), center = TRUE, scale. = FALSE)
  var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)

  df <- data.frame(
    sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    group = coldata$group,
    time = coldata$time,
    replicate = as.character(coldata$replicate),
    stringsAsFactors = FALSE
  )

  groups <- unique(df$group)
  pal <- stats::setNames(rep(cat_palette, length.out = length(groups)), groups)

  p <- ggplot(df, aes(x = .data$PC1, y = .data$PC2,
                      colour = .data$group, shape = .data$replicate)) +
    geom_point(size = 2.4) +
    geom_text(aes(label = .data$time), vjust = -0.9, size = 2.2, show.legend = FALSE) +
    scale_colour_manual(values = pal) +
    labs(x = sprintf("PC1 (%s%%)", var_exp[1]),
         y = sprintf("PC2 (%s%%)", var_exp[2]),
         colour = "Group", shape = "Replicate",
         title = label %||% sprintf("PCA, %s features", format(nrow(mat), big.mark = ","))) +
    theme_publication()

  list(df = df, var_exp = var_exp, plot = p)
}

#' Within-condition replicate correlations, for any number of replicates.
#'
#' Ported from 00_preprocessing.Rmd:464-474, which returned NULL unless a
#' condition had exactly two libraries. Any timepoint with three replicates was
#' silently dropped from the QC table, with no message. This version takes all
#' pairwise combinations, so a 2-2-2-3 design reports 1+1+1+3 correlations.
#'
#' @return data.frame(group, time, sample_a, sample_b, r)
compute_rep_cor <- function(mat, coldata, method = "pearson") {
  stopifnot(identical(rownames(coldata), colnames(mat)))
  key <- paste(coldata$group, coldata$time, sep = "_")
  out <- lapply(unique(key), function(k) {
    s <- rownames(coldata)[key == k]
    if (length(s) < 2L) return(NULL)
    pairs <- utils::combn(s, 2L)
    data.frame(
      group = coldata$group[key == k][1],
      time = coldata$time[key == k][1],
      sample_a = pairs[1, ],
      sample_b = pairs[2, ],
      r = apply(pairs, 2, function(p) stats::cor(mat[, p[1]], mat[, p[2]], method = method)),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  if (is.null(out)) {
    return(data.frame(group = character(), time = numeric(), sample_a = character(),
                      sample_b = character(), r = numeric()))
  }
  out[order(out$group, out$time), ]
}

#' Replicate correlation by timepoint, as points over a reference line.
plot_replicate_correlation <- function(rep_cor, floor_r = 0.9) {
  if (!nrow(rep_cor)) {
    return(ggplot() + annotate("text", 0, 0, label = "no replicated conditions") +
             theme_publication())
  }
  ggplot(rep_cor, aes(x = factor(.data$time), y = .data$r, colour = .data$group)) +
    # Reference lines are accent gold, never red. See lib/theme.R.
    geom_hline(yintercept = floor_r, linetype = "dashed", colour = col_accent, linewidth = 0.3) +
    geom_point(size = 2, position = position_dodge(width = 0.3)) +
    scale_colour_manual(values = stats::setNames(
      rep(cat_palette, length.out = length(unique(rep_cor$group))), unique(rep_cor$group))) +
    labs(x = "Time", y = "Replicate correlation (r)", colour = "Group",
         title = "Within-condition replicate correlation") +
    theme_publication()
}

#' Sample-to-sample correlation, as a tile plot ordered by time then replicate.
plot_sample_correlation <- function(mat, coldata, method = "pearson") {
  stopifnot(identical(rownames(coldata), colnames(mat)))
  ord <- order(coldata$group, coldata$time, coldata$replicate)
  m <- stats::cor(mat[, ord], method = method)
  lv <- colnames(m)
  df <- data.frame(
    a = factor(rep(lv, times = length(lv)), levels = lv),
    b = factor(rep(lv, each = length(lv)), levels = rev(lv)),
    r = as.vector(m)
  )
  ggplot(df, aes(.data$a, .data$b, fill = .data$r)) +
    geom_tile() +
    scale_fill_gradientn(colours = cor_palette(100), limits = c(min(m), 1)) +
    labs(x = NULL, y = NULL, fill = sprintf("%s r", method),
         title = "Sample correlation") +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}



# --- trajectories -----------------------------------------------------------

#' Long-format z-profiles joined to their trajectory class.
#'
#' The source plots read deg$normalized, which only contains the features
#' degPatterns itself clustered. On the subsample path that silently omits
#' every feature assigned by correlation, which can be most of them. This
#' builds the same long frame from the per-feature profiles and the final
#' assignment, so a figure always shows what the cluster table says.
supercluster_long <- function(profiles, assignment, include_unassigned = FALSE) {
  a <- assignment[!is.na(assignment$supercluster_label), , drop = FALSE]
  if (!include_unassigned) a <- a[a$supercluster_label != "Unassigned", , drop = FALSE]
  a <- a[a$feature_id %in% rownames(profiles), , drop = FALSE]
  times <- as.numeric(colnames(profiles))
  m <- profiles[a$feature_id, , drop = FALSE]
  data.frame(
    feature_id = rep(a$feature_id, times = length(times)),
    cluster = rep(a$cluster, times = length(times)),
    supercluster_label = rep(a$supercluster_label, times = length(times)),
    time = rep(times, each = nrow(a)),
    value = as.vector(m),
    stringsAsFactors = FALSE
  )
}

#' Facet labels carrying the class size, in canonical order.
.sc_facets <- function(long) {
  n <- tapply(long$feature_id, long$supercluster_label, function(x) length(unique(x)))
  lv <- levels(supercluster_factor(names(n)))
  lab <- stats::setNames(sprintf("%s (n = %s)", lv, format(n[lv], big.mark = ",", trim = TRUE)), lv)
  factor(lab[long$supercluster_label], levels = unname(lab))
}

.time_breaks <- function(times) sort(unique(as.numeric(times)))

#' Per-feature spaghetti, a boxplot per timepoint, and a loess trend.
#'
#' Ported from 01_degpatterns_clustering.Rmd:414-446. The original fixed the
#' class order to three MCF7 names, fixed the x breaks to 0/30/60/120/240 and
#' fixed the boxplot width to 8, which is meaningless on any other grid: on a
#' 0 to 8 day timecourse a width of 8 is the entire axis.
plot_supercluster_traces <- function(long, title = NULL, time_unit = NULL) {
  long$facet <- .sc_facets(long)
  brk <- .time_breaks(long$time)
  box_w <- max(diff(range(brk)) / 20, .Machine$double.eps)
  ggplot(long, aes(x = .data$time, y = .data$value)) +
    geom_line(aes(group = .data$feature_id), alpha = 0.05, linewidth = 0.3) +
    geom_boxplot(aes(group = .data$time), width = box_w, outlier.shape = NA,
                 fill = NA, linewidth = 0.4) +
    stat_smooth(aes(group = 1), method = "loess", formula = y ~ x, se = FALSE,
                colour = col_accent, linewidth = 1) +
    facet_wrap(~facet, nrow = 1) +
    scale_x_continuous(breaks = brk) +
    labs(x = if (is.null(time_unit)) "Time" else sprintf("Time (%s)", time_unit),
         y = "Z-score", title = title) +
    theme_publication()
}

#' The publication panel: median line over interquartile and 10-90 bands.
#'
#' Ported from 01_degpatterns_clustering.Rmd:468-514, with the same
#' generalizations, and taking its colours from lib/theme.R rather than from
#' a local three-class vector that disagreed with the shared palette.
plot_supercluster_ribbon <- function(long, title = NULL, time_unit = NULL) {
  long$sc <- supercluster_factor(long$supercluster_label)
  long$facet <- .sc_facets(long)
  brk <- .time_breaks(long$time)
  key <- interaction(long$sc, long$time, drop = TRUE)
  agg <- do.call(rbind, lapply(split(seq_len(nrow(long)), key), function(i) {
    v <- long$value[i]
    data.frame(sc = long$sc[i][1], facet = long$facet[i][1], time = long$time[i][1],
               median = stats::median(v, na.rm = TRUE),
               q25 = stats::quantile(v, 0.25, na.rm = TRUE),
               q75 = stats::quantile(v, 0.75, na.rm = TRUE),
               q10 = stats::quantile(v, 0.10, na.rm = TRUE),
               q90 = stats::quantile(v, 0.90, na.rm = TRUE))
  }))
  pal <- supercluster_palette(levels(long$sc))
  ggplot(agg, aes(x = .data$time)) +
    geom_ribbon(aes(ymin = .data$q10, ymax = .data$q90, fill = .data$sc), alpha = 0.15) +
    geom_ribbon(aes(ymin = .data$q25, ymax = .data$q75, fill = .data$sc), alpha = 0.30) +
    geom_line(aes(y = .data$median, colour = .data$sc), linewidth = 1) +
    geom_point(aes(y = .data$median, colour = .data$sc), size = 1.6) +
    facet_wrap(~facet, nrow = 1) +
    scale_colour_manual(values = pal) +
    scale_fill_manual(values = pal) +
    scale_x_continuous(breaks = brk) +
    labs(x = if (is.null(time_unit)) "Time" else sprintf("Time (%s)", time_unit),
         y = "Z-score", title = title) +
    theme_publication() +
    theme(legend.position = "none")
}

#' Dendrogram of cluster mean profiles, with the cut drawn.
#'
#' Ported from 01_degpatterns_clustering.Rmd:332-354, which returned a base
#' plot and drew no cut, because in that analysis the cut was made by eye
#' afterwards. Drawing the height that was actually used is the point: it is
#' what makes the assignment checkable.
plot_cluster_dendrogram <- function(hc, cluster_to_sc = NULL, cut_height = NULL,
                                    title = NULL) {
  dd <- stats::as.dendrogram(hc)
  ord <- stats::order.dendrogram(dd)
  labs <- hc$labels[ord]
  seg <- ggdendro_segments(hc)
  tip <- data.frame(x = seq_along(labs), label = labs,
                    sc = if (is.null(cluster_to_sc)) NA_character_
                         else unname(cluster_to_sc[labs]),
                    stringsAsFactors = FALSE)
  p <- ggplot() +
    geom_segment(data = seg, aes(x = .data$x, y = .data$y,
                                 xend = .data$xend, yend = .data$yend),
                 linewidth = 0.3) +
    scale_x_continuous(breaks = tip$x, labels = tip$label) +
    labs(x = "Cluster", y = "Height", title = title, colour = NULL) +
    theme_publication()
  if (!is.null(cut_height)) {
    p <- p + geom_hline(yintercept = cut_height, linetype = "dashed",
                        colour = col_accent, linewidth = 0.4)
  }
  if (!all(is.na(tip$sc))) {
    p <- p + geom_point(data = tip, aes(x = .data$x, y = 0, colour = .data$sc), size = 2.2) +
      scale_colour_manual(values = supercluster_palette(unique(tip$sc)))
  }
  p
}

#' Dendrogram segments as a data frame, so the tree can be drawn in ggplot
#' without taking a dependency on ggdendro for one function.
ggdendro_segments <- function(hc) {
  merge <- hc$merge
  height <- hc$height
  ord <- order(stats::order.dendrogram(stats::as.dendrogram(hc)))
  xpos <- numeric(nrow(merge))
  ypos <- numeric(nrow(merge))
  segs <- list()
  leaf_x <- function(i) ord[-i]
  for (k in seq_len(nrow(merge))) {
    kids <- merge[k, ]
    cx <- numeric(2); cy <- numeric(2)
    for (j in 1:2) {
      if (kids[j] < 0) { cx[j] <- leaf_x(kids[j]); cy[j] <- 0 }
      else { cx[j] <- xpos[kids[j]]; cy[j] <- ypos[kids[j]] }
    }
    xpos[k] <- mean(cx); ypos[k] <- height[k]
    segs[[length(segs) + 1L]] <- data.frame(
      x = c(cx[1], cx[1], cx[2]), y = c(cy[1], height[k], height[k]),
      xend = c(cx[1], cx[2], cx[2]), yend = c(height[k], height[k], cy[2]))
  }
  do.call(rbind, segs)
}

#' Dynamic against static features, as a labelled stacked bar.
plot_dynamic_vs_static <- function(res, arm, title = NULL) {
  df <- data.frame(
    arm = arm,
    status = factor(c("Dynamic", "Static"), levels = c("Static", "Dynamic")),
    n = c(sum(res$dynamic), sum(!res$dynamic))
  )
  df$label <- format(df$n, big.mark = ",", trim = TRUE)
  ggplot(df, aes(x = .data$arm, y = .data$n, fill = .data$status)) +
    geom_col(width = 0.55) +
    geom_text(aes(label = .data$label), position = position_stack(vjust = 0.5), size = 2.6) +
    scale_fill_manual(values = c(Static = "#BBBBBB", Dynamic = cat_palette[1])) +
    labs(x = NULL, y = "Features", fill = NULL, title = title) +
    theme_publication()
}

invisible(TRUE)
