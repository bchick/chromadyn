# ---------------------------------------------------------------------------
# classify.R: name a trajectory from the shape of its centroid.
#
# This is the half of the method that was done by eye in the source analysis.
# Fourteen degPatterns clusters were mapped onto three names by a hand-written
# tribble(), and the three were later split into five by a k-means whose
# labelling rules were hard-coded to one experiment's timepoint grid: the
# literal strings "60" and "240" appear in the conditions.
#
# One function does both stages. It computes the shape statistics once,
# evaluates an ordered list of predicates to a slot, then maps the slot to a
# name through a vocabulary. Stage 1 names broad classes; stage 2 names the
# pieces of a class that has been sub-split. The predicates are the source's,
# generalized off the literal timepoints.
# ---------------------------------------------------------------------------

#' Derive the reference timepoints of a grid.
#'
#' `mid` is the observed timepoint nearest a fixed FRACTION of the last
#' timepoint, not the middle post-baseline index. On the source grid
#' c(0, 30, 60, 120, 240) a fraction of 0.25 gives exactly 60, reproducing the
#' value hard-coded there.
#'
#' The index-based reading looks equivalent on that grid and is not. The
#' source grid is roughly log spaced (30, 60, 120, 240), so its middle index
#' sits a quarter of the way through in TIME; on a near-linear grid such as
#' c(0, 3, 5, 8) the middle index sits at 62% of the way through, by which
#' point a genuinely late-rising trajectory has already risen and the
#' late-versus-sustained test silently stops working.
grid_refs <- function(times, mid_time = NULL, mid_fraction = 0.25) {
  times <- sort(unique(as.numeric(times)))
  last <- max(times)
  post <- times[-1]
  mid <- if (!is.null(mid_time)) {
    as.numeric(mid_time)
  } else {
    post[which.min(abs(post - mid_fraction * last))]
  }
  list(first = times[1], last = last, mid = mid, all = times)
}

#' Shape statistics for a matrix of centroids (k rows, one column per timepoint).
centroid_stats <- function(centroids, times = as.numeric(colnames(centroids)),
                           mid_time = NULL, mid_fraction = 0.25) {
  refs <- grid_refs(times, mid_time, mid_fraction)
  peak_time <- times[apply(centroids, 1, which.max)]
  peak_val <- apply(centroids, 1, max)
  last_val <- centroids[, which.min(abs(times - refs$last)), drop = TRUE]
  mid_val <- centroids[, which.min(abs(times - refs$mid)), drop = TRUE]

  # late_ratio divides by the maximum, which is meaningless when the maximum
  # is not positive: on z-scored profiles a monotonically falling trajectory
  # can peak at a negative value, and the ratio then flips sign. Mark those
  # NA and let the caller fall back rather than emit a confident wrong label.
  late_ratio <- ifelse(peak_val > 0, last_val / peak_val, NA_real_)

  data.frame(
    row = rownames(centroids) %||% as.character(seq_len(nrow(centroids))),
    peak_time = peak_time, peak_value = peak_val,
    mid_value = as.numeric(mid_val), last_value = as.numeric(last_val),
    late_ratio = as.numeric(late_ratio),
    stringsAsFactors = FALSE
  )
}

# Slot to name, per stage. Keeping the vocabularies as data rather than as
# branches is what lets both stages share one predicate evaluator.
CD_VOCAB <- list(
  stage1 = c(down = "Decreasing", transient = "Transient", up = "Increasing"),
  split = c(transient_up = "Transient Increasing",
            late_up = "Late Increasing",
            sustained_up = "Sustained Increasing"),
  simple = c(down = "Decreasing", up = "Increasing")
)

#' Label centroids by shape.
#'
#' @param centroids Numeric matrix, k x T, colnames parseable as timepoints.
#' @param scheme "shape" for the three-way stage 1 and the three-way split,
#'   "shape_simple" for the two-class variant from the CUT&RUN analysis,
#'   "numeric" for Cluster-1..k, "none" for the bare cluster ids.
#' @param vocab "stage1" or "split".
#' @param params late_ratio_cutoff, late_peak_fraction, mid_time, cross,
#'   simple_cutoff.
#' @param prefix Base name for the numeric fallback.
#' @param unique_labels Whether two centroids sharing a label is a failure.
#'   TRUE for the split stage, where the k pieces of one class must each get a
#'   distinct name. FALSE for stage 1, where several clusters sharing a label
#'   is the whole point: that is how many clusters become a few classes.
#' @return data.frame with the stats plus `slot` and `label`.
classify_centroids <- function(centroids,
                               scheme = c("shape", "shape_simple", "numeric", "none"),
                               vocab = c("stage1", "split"),
                               params = list(),
                               prefix = "Cluster",
                               unique_labels = (match.arg(vocab) == "split"),
                               times = as.numeric(colnames(centroids))) {
  scheme <- match.arg(scheme)
  vocab <- match.arg(vocab)
  p <- utils::modifyList(
    list(late_ratio_cutoff = 0.5, late_peak_fraction = 0.5, mid_fraction = 0.25,
         mid_time = NULL, cross = 0, simple_cutoff = 0.7),
    params[!vapply(params, is.null, logical(1))]
  )
  refs <- grid_refs(times, p$mid_time, p$mid_fraction)
  st <- centroid_stats(centroids, times, p$mid_time, p$mid_fraction)

  st$slot <- NA_character_
  if (scheme == "shape") {
    if (vocab == "stage1") {
      # Ordered: first match wins.
      st$slot <- ifelse(
        st$peak_time <= refs$first, "down",
        ifelse(st$peak_time <= refs$mid & st$last_value < p$cross, "transient", "up"))
    } else {
      # An NA late_ratio means the centroid never rises above its own mean,
      # so none of the three increasing shapes describes it. Leave the slot
      # unset and let the fallback number it, rather than defaulting a
      # falling profile into "Sustained Increasing".
      st$slot <- ifelse(
        is.na(st$late_ratio), NA_character_,
        ifelse(st$late_ratio < p$late_ratio_cutoff, "transient_up",
               ifelse(st$peak_time >= p$late_peak_fraction * refs$last & st$mid_value < p$cross,
                      "late_up", "sustained_up")))
    }
  } else if (scheme == "shape_simple") {
    # The two-class variant: anything still high at the end is Increasing,
    # anything peaking at baseline is Decreasing, everything else Increasing.
    st$slot <- ifelse(
      !is.na(st$late_ratio) & st$late_ratio >= p$simple_cutoff, "up",
      ifelse(st$peak_time <= refs$first, "down", "up"))
    vocab <- "simple"
  }

  st$label <- if (scheme %in% c("numeric", "none")) {
    NA_character_
  } else {
    unname(CD_VOCAB[[vocab]][st$slot])
  }

  # Fallback, kept as the source had it: if the rules cannot name the
  # centroids, say so and number them rather than emit names that silently
  # merge two different shapes.
  #
  # What counts as a failure depends on the stage. In the split, two pieces of
  # one class sharing a name means the rules could not tell them apart. In
  # stage 1, several clusters sharing a name is the intended result, so only
  # an unlabelled centroid is a problem.
  bad <- is.na(st$label)
  if (unique_labels) {
    bad <- bad | duplicated(st$label) | duplicated(st$label, fromLast = TRUE)
  }
  if (scheme %in% c("numeric", "none") || any(bad)) {
    if (scheme %in% c("shape", "shape_simple")) {
      message(sprintf(paste0("[classify] shape labelling was ambiguous (%d of %d centroids ",
                             "share or lack a label); using numeric labels instead."),
                      sum(bad), nrow(st)))
      message(sprintf("[classify]   reference timepoints: first=%s mid=%s last=%s",
                      refs$first, refs$mid, refs$last))
      for (i in seq_len(nrow(st))) {
        message(sprintf(paste0("[classify]   %-14s peak_time=%-6s mid_value=%+.2f ",
                               "last_value=%+.2f late_ratio=%s"),
                        st$row[i], st$peak_time[i], st$mid_value[i], st$last_value[i],
                        if (is.na(st$late_ratio[i])) "NA" else sprintf("%+.2f", st$late_ratio[i])))
      }
      message(paste0("[classify]   To separate these, tune superclusters.split.cross ",
                     "(a margin below the feature mean, currently ", format(p$cross),
                     ") or late_ratio_cutoff. See docs/parameters.md."))
    }
    st$label <- if (scheme == "none") st$row else paste0(prefix, "-", seq_len(nrow(st)))
  }
  st
}

`%||%` <- function(a, b) if (is.null(a)) b else a

invisible(TRUE)
