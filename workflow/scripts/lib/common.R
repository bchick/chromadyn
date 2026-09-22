# ---------------------------------------------------------------------------
# common.R: config access, path resolution, IO, assertions and provenance.
#
# Sourced first by every chromadyn script. Everything a script needs to talk to
# the outside world lives here, so that the scripts themselves contain method
# and nothing else.
#
# Two rules this file exists to enforce:
#   1. No script contains a numeric threshold. Scripts call cfg(); the value
#      comes from config.yaml. tests/test_static.sh greps for violations.
#   2. No script builds an output path by pasting strings. Paths come from
#      workflow/paths.yaml via cd_path(), which the Snakefile reads too.
#
# House style follows mcf7_project/analyses/40_brg1_inaccessible_classes/R/common.R
# (path constants, one parameter list, small resolvers) with config accessors
# in place of here::here().
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(yaml)
  library(jsonlite)
})

# Session state. Populated by cd_init(); every accessor below reads from here.
.CD <- new.env(parent = emptyenv())
.CD$notes <- list()

# --- config ----------------------------------------------------------------

#' Read a config value by dotted path, for example cfg("cluster.degpatterns.minc").
#'
#' Returns `default` when the key is absent. A key that is present but YAML-null
#' returns NULL, which is meaningful in this config (features: null means gene
#' mode) and is deliberately distinct from absent.
cfg <- function(key, default = NULL) {
  stopifnot(is.character(key), length(key) == 1L)
  node <- .CD$config
  if (is.null(node)) stop("cfg(): no config loaded. Call cd_init(snakemake) first.", call. = FALSE)
  for (part in strsplit(key, ".", fixed = TRUE)[[1]]) {
    if (!is.list(node) || !(part %in% names(node))) return(default)
    node <- node[[part]]
  }
  node
}

#' As cfg(), but a missing or NULL value is a hard error naming the key.
cfg_req <- function(key) {
  val <- cfg(key)
  if (is.null(val)) {
    stop(sprintf("cfg_req(): required config key '%s' is missing or null.", key), call. = FALSE)
  }
  val
}

#' The seed for one stochastic step, derived from cluster.seed so that every
#' step is reproducible and no two steps share a stream. Recorded in the
#' manifest by cd_finish().
cd_seed <- function(step) {
  base <- as.integer(cfg("cluster.seed", 42L))
  # Stable per-step offset from the step name, so adding a step never shifts
  # the seed of an existing one.
  offset <- sum(utf8ToInt(step) * seq_along(utf8ToInt(step))) %% 10000L
  seed <- base + offset
  .CD$seeds[[step]] <- seed
  seed
}

#' Set the seed for a step and record it. Use immediately before the
#' stochastic call, never once at the top of a script.
cd_set_seed <- function(step) {
  seed <- cd_seed(step)
  set.seed(seed)
  invisible(seed)
}

# --- paths ------------------------------------------------------------------

#' Resolve an output path template from workflow/paths.yaml.
#'
#' cd_path("clusters", arm = "WT") -> "results/clusters/WT_clusters.tsv"
#' Unfilled placeholders are an error, so a typo cannot silently write a file
#' with a literal "{arm}" in its name.
cd_path <- function(key, ...) {
  tmpl <- .CD$paths[[key]]
  if (is.null(tmpl)) {
    stop(sprintf("cd_path(): '%s' is not in workflow/paths.yaml. Known keys: %s",
                 key, paste(sort(names(.CD$paths)), collapse = ", ")), call. = FALSE)
  }
  subs <- c(list(out = cfg("output.dir", "results")), list(...))
  for (nm in names(subs)) {
    tmpl <- gsub(paste0("{", nm, "}"), as.character(subs[[nm]]), tmpl, fixed = TRUE)
  }
  if (grepl("\\{[a-z_]+\\}", tmpl)) {
    stop(sprintf("cd_path('%s'): unfilled placeholder in '%s'. Supply it as a named argument.",
                 key, tmpl), call. = FALSE)
  }
  tmpl
}

#' Every path a figure key expands to, one per configured format.
cd_fig <- function(key, ..., formats = cfg("figures.formats")) {
  vapply(formats, function(f) cd_path(key, ..., fmt = f), character(1), USE.NAMES = FALSE)
}

#' Filesystem-safe version of a label that may contain spaces, for filenames.
#' "Late Increasing" -> "late_increasing"
cd_slug <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

cd_mkdir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

# --- assertions -------------------------------------------------------------

#' Fail with a formatted message. The workhorse; every other assert calls it.
cd_assert <- function(cond, fmt, ...) {
  if (!isTRUE(all(cond))) stop(sprintf(fmt, ...), call. = FALSE)
  invisible(TRUE)
}

assert_cols <- function(df, cols, what = "table") {
  missing <- setdiff(cols, colnames(df))
  cd_assert(length(missing) == 0L,
            "%s is missing required column(s): %s\nPresent: %s",
            what, paste(missing, collapse = ", "), paste(colnames(df), collapse = ", "))
}

#' Guardrail 9.2. Metadata rownames must equal matrix colnames in the SAME
#' ORDER, not merely as sets. A reordered merge() attached positionally is a
#' known way to produce a confident wrong answer, and degPatterns will not
#' notice.
assert_aligned <- function(metadata, mat, what = "degPatterns input") {
  rn <- rownames(metadata)
  cn <- colnames(mat)
  if (identical(rn, cn)) return(invisible(TRUE))
  if (setequal(rn, cn)) {
    stop(sprintf(paste0("%s: metadata rownames and matrix colnames contain the same ",
                        "samples but in DIFFERENT ORDER.\n  metadata: %s\n  matrix:   %s\n",
                        "Reorder the metadata to match the matrix; do not rely on ",
                        "positional attachment."),
                 what, paste(rn, collapse = ", "), paste(cn, collapse = ", ")), call. = FALSE)
  }
  stop(sprintf(paste0("%s: metadata rownames and matrix colnames differ.\n",
                      "  only in metadata: %s\n  only in matrix:   %s"),
               what,
               paste(setdiff(rn, cn), collapse = ", "),
               paste(setdiff(cn, rn), collapse = ", ")), call. = FALSE)
}

#' Guardrail 9.1. The single most destructive silent failure in this method.
#'
#' factor(c(0,15,30,60,120,240)) sorts numerically. factor(as.character(...))
#' sorts lexicographically to 0, 120, 15, 240, 30, 60, every trajectory becomes
#' noise, and nothing errors. Assert the levels are numeric and increasing.
assert_time_levels <- function(f, what = "time factor") {
  cd_assert(is.factor(f), "%s is not a factor (class: %s).", what, paste(class(f), collapse = "/"))
  lv <- levels(f)
  num <- suppressWarnings(as.numeric(lv))
  cd_assert(!any(is.na(num)),
            "%s has non-numeric levels: %s\nTimepoints must parse as numbers.",
            what, paste(lv[is.na(num)], collapse = ", "))
  cd_assert(!is.unsorted(num),
            paste0("%s levels are not in increasing numeric order: %s\n",
                   "This is the lexicographic-sort trap: every trajectory would be ",
                   "reordered into noise with no error raised. Build the factor with ",
                   "levels = sort(unique(as.numeric(time)))."),
            what, paste(lv, collapse = ", "))
  invisible(TRUE)
}

#' Guardrail 9.5. A supercluster map that misses a cluster silently drops those
#' features on the join. Name the unmapped clusters and fail.
assert_map_covers <- function(map, clusters, what = "superclusters.overrides") {
  unmapped <- setdiff(as.character(unique(clusters)), as.character(names(map)))
  cd_assert(length(unmapped) == 0L,
            paste0("%s does not cover every cluster.\n  unmapped: %s\n",
                   "  mapped:   %s\nAdd them to the map, or use a method that labels ",
                   "clusters automatically."),
            what, paste(sort(unmapped), collapse = ", "),
            paste(sort(names(map)), collapse = ", "))
}

assert_no_na <- function(df, cols, what = "table") {
  for (cl in cols) {
    n <- sum(is.na(df[[cl]]))
    cd_assert(n == 0L, "%s: column '%s' has %d NA value(s) where none are allowed.", what, cl, n)
  }
  invisible(TRUE)
}

#' Guardrail 9.4. Coordinates must be real, usable numbers before anything is
#' written to BED. The source analysis wrote thousands of rows of
#' "chr1.7401731.7402231  NA  NA" because a feature ID was parsed with a regex
#' that did not match after make.names() had replaced the separators.
assert_coords <- function(df, what = "coordinates") {
  assert_cols(df, c("chr", "start", "end"), what)
  assert_no_na(df, c("chr", "start", "end"), what)
  cd_assert(is.numeric(df$start) && is.numeric(df$end),
            "%s: start and end must be numeric (got %s and %s).",
            what, class(df$start)[1], class(df$end)[1])
  bad <- which(df$start >= df$end)
  cd_assert(length(bad) == 0L,
            "%s: %d row(s) have start >= end, first at row %d (%s:%s-%s).",
            what, length(bad), bad[1], df$chr[bad[1]], df$start[bad[1]], df$end[bad[1]])
  cd_assert(all(nzchar(as.character(df$chr))), "%s: empty chromosome name(s).", what)
  invisible(TRUE)
}

# --- IO ---------------------------------------------------------------------

#' Read a counts matrix from .tsv, .csv, .rds (matrix, data.frame,
#' SummarizedExperiment or DESeqDataSet). Returns a base matrix with feature
#' rownames and sample colnames.
read_counts <- function(path) {
  cd_assert(file.exists(path), "counts file not found: %s", path)
  ext <- tolower(tools::file_ext(path))

  obj <- switch(
    ext,
    rds = readRDS(path),
    csv = utils::read.csv(path, row.names = 1, check.names = FALSE),
    utils::read.delim(path, row.names = 1, check.names = FALSE)   # tsv, txt, anything else
  )

  if (inherits(obj, c("SummarizedExperiment", "DESeqDataSet", "RangedSummarizedExperiment"))) {
    cd_assert(requireNamespace("SummarizedExperiment", quietly = TRUE),
              "%s holds a SummarizedExperiment but the package is not installed.", path)
    an <- SummarizedExperiment::assayNames(obj)
    # Prefer raw counts explicitly. Taking assay(obj) blindly can hand back a
    # vst assay that happens to be stored alongside.
    which_assay <- if ("counts" %in% an) "counts" else an[1]
    if (which_assay != "counts") {
      message(sprintf("read_counts(): no 'counts' assay in %s, using '%s'. ",
                      basename(path), which_assay),
              "Set transform.method: none if this is already normalized.")
    }
    obj <- SummarizedExperiment::assay(obj, which_assay)
  }

  mat <- as.matrix(obj)
  cd_assert(!is.null(rownames(mat)), "counts matrix %s has no rownames (feature IDs).", path)
  cd_assert(!is.null(colnames(mat)), "counts matrix %s has no colnames (sample IDs).", path)
  mat
}

#' Read and minimally validate the samplesheet. Column semantics are checked
#' properly in 01_validate.R; this only guarantees the shape.
read_samplesheet <- function(path) {
  cd_assert(file.exists(path), "samplesheet not found: %s", path)
  ss <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  assert_cols(ss, c("sample", "group", "time", "replicate"), "samplesheet")
  ss$sample <- as.character(ss$sample)
  ss$group <- as.character(ss$group)
  ss$replicate <- as.character(ss$replicate)
  ss$time <- suppressWarnings(as.numeric(ss$time))
  bad <- which(is.na(ss$time))
  cd_assert(length(bad) == 0L,
            paste0("samplesheet: `time` must be numeric; row(s) %s are not.\n",
                   "Units go in config input.time_unit, never in the value."),
            paste(bad, collapse = ", "))
  ss
}

#' Read BED3 or BED6. For BED6 the name column maps to counts rownames; for
#' BED3 row order does. Returns a data.frame with chr/start/end/name.
read_bed <- function(path, feature_ids = NULL) {
  cd_assert(file.exists(path), "features BED not found: %s", path)
  bed <- utils::read.delim(path, header = FALSE, stringsAsFactors = FALSE, comment.char = "#")
  cd_assert(ncol(bed) >= 3L, "features BED %s has %d columns, need at least 3.", path, ncol(bed))
  out <- data.frame(chr = as.character(bed[[1]]),
                    start = as.numeric(bed[[2]]),
                    end = as.numeric(bed[[3]]),
                    stringsAsFactors = FALSE)
  if (ncol(bed) >= 4L) {
    out$name <- as.character(bed[[4]])
  } else {
    cd_assert(!is.null(feature_ids),
              paste0("features BED %s is BED3, so rows map to features by order, ",
                     "but no feature IDs were supplied."), path)
    cd_assert(nrow(out) == length(feature_ids),
              paste0("features BED %s is BED3 with %d rows but there are %d ",
                     "features; order mapping is ambiguous."),
              path, nrow(out), length(feature_ids))
    out$name <- feature_ids
  }
  assert_coords(out, sprintf("features BED %s", basename(path)))
  out
}

read_tsv_strict <- function(path, cols = NULL, what = basename(path)) {
  cd_assert(file.exists(path), "file not found: %s", path)
  df <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!is.null(cols)) assert_cols(df, cols, what)
  df
}

#' Write a TSV via a temporary file in the same directory, then rename.
#' An interrupted run leaves no half-written table for Snakemake to treat as
#' complete on the next invocation.
write_tsv_atomic <- function(df, path, ...) {
  cd_mkdir(dirname(path))
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  gz <- grepl("\\.gz$", path, ignore.case = TRUE)
  con <- if (gz) gzfile(tmp, "w") else file(tmp, "w")
  utils::write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, ...)
  close(con)
  ok <- file.rename(tmp, path)
  cd_assert(ok, "failed to move temporary file into place: %s", path)
  invisible(path)
}

save_obj <- function(obj, path) {
  cd_mkdir(dirname(path))
  saveRDS(obj, path, version = 3)
  invisible(path)
}

load_obj <- function(path) {
  cd_assert(file.exists(path), "cached object not found: %s", path)
  readRDS(path)
}

# --- provenance -------------------------------------------------------------

#' Record a fact about this run. Merged into run_manifest.json by 12_manifest.R.
cd_note <- function(key, value) {
  .CD$notes[[key]] <- value
  invisible(value)
}

#' Start a script: load config and paths off the snakemake object, open a log,
#' and begin a provenance fragment. Every script calls this before doing work.
cd_init <- function(snakemake) {
  .CD$config <- snakemake@config
  .CD$rule <- tryCatch(snakemake@rule, error = function(e) "unknown")
  .CD$wildcards <- tryCatch(as.list(snakemake@wildcards), error = function(e) list())
  .CD$seeds <- list()
  .CD$notes <- list()
  .CD$started <- Sys.time()

  paths_file <- file.path(snakemake@scriptdir, "..", "paths.yaml")
  cd_assert(file.exists(paths_file), "workflow/paths.yaml not found at %s", paths_file)
  .CD$paths <- yaml::read_yaml(paths_file)

  message(sprintf("[chromadyn] %s%s starting",
                  .CD$rule,
                  if (!is.null(.CD$wildcards$arm)) paste0(" (", .CD$wildcards$arm, ")") else ""))
  invisible(TRUE)
}

#' Close a script: write the provenance fragment 12_manifest.R will merge.
cd_finish <- function() {
  arm <- .CD$wildcards$arm
  frag <- cd_path("prov_fragment",
                  rule = .CD$rule,
                  arm_suffix = if (is.null(arm)) "" else paste0("_", arm))
  cd_mkdir(dirname(frag))
  payload <- list(
    rule = .CD$rule,
    arm = if (is.null(arm)) NA_character_ else arm,
    started = format(.CD$started, "%Y-%m-%dT%H:%M:%S%z"),
    finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    elapsed_sec = round(as.numeric(difftime(Sys.time(), .CD$started, units = "secs")), 2),
    seeds = .CD$seeds,
    notes = .CD$notes,
    r_version = R.version.string
  )
  writeLines(jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null", pretty = TRUE), frag)
  message(sprintf("[chromadyn] %s done in %.1fs", .CD$rule, payload$elapsed_sec))
  invisible(frag)
}

# --- small shared helpers ---------------------------------------------------

#' Per-timepoint means of a matrix, columns ordered by increasing numeric time.
#' Used by the range gate, the k-means sub-split and the heatmaps, so it is
#' defined once here rather than three times inline as in the source notebooks.
timepoint_means <- function(mat, times) {
  cd_assert(ncol(mat) == length(times),
            "timepoint_means(): %d columns but %d time values.", ncol(mat), length(times))
  idx <- split(seq_along(times), factor(times, levels = sort(unique(as.numeric(times)))))
  out <- vapply(idx, function(i) rowMeans(mat[, i, drop = FALSE]), numeric(nrow(mat)))
  rownames(out) <- rownames(mat)
  out
}

#' Row-wise z-score, returning a matrix. Rows with zero variance become 0
#' rather than NaN, and the count of such rows is returned as an attribute.
zscore_rows <- function(mat) {
  z <- t(scale(t(mat)))
  flat <- which(!is.finite(rowSums(z)))
  if (length(flat)) z[flat, ] <- 0
  attr(z, "n_zero_variance") <- length(flat)
  z
}

invisible(TRUE)
