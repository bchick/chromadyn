# ---------------------------------------------------------------------------
# plots.R: every figure chromadyn draws.
#
# These exist as functions because in the source project each one was written
# inline, two or three times, in different notebooks, and the copies drifted.
# Nothing here reads the config or writes a file: a plot function takes data
# and returns a ggplot, and the calling script decides where it goes.
#
# Ported and generalized from the source notebooks. The generalizations are
# the point: the originals hard-coded one experiment's timepoint grid, one
# experiment's ligand colours, and an assumption of exactly two replicates.
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
#' three ligands, so it could not be reused. Here the grid comes from the data
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

`%||%` <- function(a, b) if (is.null(a)) b else a

invisible(TRUE)
