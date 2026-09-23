# ---------------------------------------------------------------------------
# namespace_fix.R: undo Bioconductor's masking of tidyverse verbs.
#
# Loading DESeq2, DEGreport, AnnotationDbi, mclust or S4Vectors quietly
# replaces functions you thought you were calling. dplyr::select becomes
# AnnotationDbi::select, purrr::map becomes mclust::map, base::paste becomes
# BiocGenerics::paste, and so on. None of it errors; the code just starts
# doing something else. DEGreport pulls several of these in at once.
#
# Sourced by every timecourse-patterns script, immediately after lib/common.R.
#
# Assignment MUST be into .GlobalEnv. A plain `<-` inside this file's own
# evaluation frame, or sourcing this file from inside a function, does not
# propagate to the calling script.
#
# Ported from mcf7_project/scripts/rmd_setup.R, which already covered a
# superset of the functions listed in the timecourse-patterns specification.
# ---------------------------------------------------------------------------

local({
  # function -> the package it must resolve to. Ordered most-masked last, so
  # that if two entries ever collide the later one wins.
  pairs <- list(
    c("count",     "dplyr"),
    c("filter",    "dplyr"),   # masked by stats
    c("select",    "dplyr"),   # masked by AnnotationDbi, MASS
    c("rename",    "dplyr"),   # masked by S4Vectors
    c("map",       "purrr"),   # masked by mclust
    c("first",     "dplyr"),   # masked by S4Vectors
    c("slice",     "dplyr"),
    c("desc",      "dplyr"),
    c("paste",     "base"),    # masked by BiocGenerics
    c("setdiff",   "base"),
    c("union",     "base"),
    c("intersect", "base")
  )

  for (p in pairs) {
    fn <- p[[1]]
    pkg <- p[[2]]
    if (requireNamespace(pkg, quietly = TRUE)) {
      val <- tryCatch(getExportedValue(pkg, fn), error = function(e) NULL)
      if (!is.null(val)) assign(fn, val, envir = .GlobalEnv)
    }
  }
})

invisible(TRUE)
