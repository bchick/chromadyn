#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# build_demo.R: derive the bundled chromadyn demo from the published T cell
# ATAC timecourse.
#
# Source data: McDonald, Chick et al., Immunity 56(6):1303-1319.e5 (2023),
# doi:10.1016/j.immuni.2023.05.005. GEO SuperSeries GSE228381, ATAC sub-series
# GSE228171. Counts are the nf-core/atacseq 2.1.2 consensus-peak
# featureCounts matrix (129,314 peaks x 62 libraries, GRCm39).
#
# This script is committed so the provenance of demo/counts.tsv is auditable:
# anyone can see exactly which libraries were chosen, which rows were kept, and
# with what seed. The large source matrix is NOT committed (see .gitignore).
#
# Run from the repo root:
#   Rscript demo/build_demo.R [path/to/consensus_peaks.mRp.clN.featureCounts.txt]
#
# Outputs: demo/counts.tsv, demo/samplesheet.tsv, demo/features.bed
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(stats)
})

# --- parameters ------------------------------------------------------------
# Everything that makes this subsample what it is, in one place.
#
# The mix matters. An early version drew the static pool from the bottom half
# of the range distribution, which on ATAC data this deep is not static at all:
# 89% of the demo came back dynamic and there was nothing for the
# dynamic-versus-static figure to contrast. The static pool is now the bottom
# 15%, and a middle band is sampled too, so the demo is not purely bimodal.
PARAMS <- list(
  seed            = 42L,
  n_dynamic       = 2500L,   # clearly time-varying, shape-stratified
  n_flat          = 1500L,   # genuinely static, so the false-positive rate is measurable
  n_middle        = 500L,    # the ambiguous middle, so the demo is not purely bimodal
  n_lowcount      = 500L,    # below the pipeline's prefilter, so prefilter has work to do
  shape_k         = 6L,      # k-means groups used ONLY to spread the dynamic sample
  dynamic_top_q   = 0.80,    # dynamic pool: range above this quantile
  flat_max_q      = 0.15,    # flat pool: range below this quantile
  min_mean_counts = 10,      # matches config filter.min_mean_counts
  lowcount_min    = 1,       # ignore features that are essentially all zero
  std_chroms      = paste0("chr", c(1:19, "X", "Y"))
)

DEFAULT_MATRIX <- file.path(
  "/data/bchick/wproj/tcell_project/results/atac/bowtie2/merged_replicate",
  "macs2/narrow_peak/consensus/consensus_peaks.mRp.clN.featureCounts.txt")

# The nine libraries, and the design they encode.
#
# Naive / D3 / D5 are unsorted total CD8; each has exactly two libraries and
# there is nothing to choose. D8 is the three terminal-effector (TE) libraries
# from experiment 2: the deepest TE set, all within one experiment, taken as
# real libraries rather than pseudobulked across sorted subsets. That leaves
# D8 with three replicates against two elsewhere, which DESeq2 handles and
# which gives the demo a realistically unbalanced design.
#
# The 48h timepoint is deliberately absent: FRiP 0.115 and 0.314 against
# 0.38-0.65 for every other library. Every temporal analysis in the source
# project drops it too.
SAMPLES <- data.frame(
  sample    = c("Naive_WT_REP1", "Naive_WT_REP2",
                "D3_WT_REP1", "D3_WT_REP2",
                "D5_WT_REP1", "D5_WT_REP2",
                "D8_WT_TE_Exp2_REP1", "D8_WT_TE_Exp2_REP2", "D8_WT_TE_Exp2_REP3"),
  group     = "WT",
  time      = c(0, 0, 3, 3, 5, 5, 8, 8, 8),   # days post infection
  replicate = c("r1", "r2", "r1", "r2", "r1", "r2", "r1", "r2", "r3"),
  stringsAsFactors = FALSE
)
# Deliberately no `batch` column. The only batch variable available (Exp1 vs
# Exp2) is perfectly confounded with time here, because every D8 library is
# Exp2 and no earlier timepoint has an experiment label at all. Shipping it
# would invite a model that cannot be fit.

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  mpath <- if (length(args) >= 1) args[1] else DEFAULT_MATRIX
  outdir <- "demo"
  stopifnot(dir.exists(outdir))
  if (!file.exists(mpath)) {
    stop("Source matrix not found: ", mpath,
         "\nPass the path as the first argument. See demo/PROVENANCE.md.", call. = FALSE)
  }

  set.seed(PARAMS$seed)
  message("Reading ", basename(mpath), " ...")
  # Line 1 is the featureCounts command banner.
  fc <- utils::read.delim(mpath, skip = 1, check.names = FALSE, stringsAsFactors = FALSE)

  # nf-core leaves the aligner suffix on the column names even in the
  # merged-replicate matrix. Strip it so names match the samplesheet exactly.
  colnames(fc) <- sub("\\.mLb\\.clN\\.sorted\\.bam$", "", colnames(fc))

  missing <- setdiff(SAMPLES$sample, colnames(fc))
  if (length(missing)) stop("columns absent from the matrix: ", paste(missing, collapse = ", "))

  # featureCounts in SAF mode with -O emits semicolon-joined coordinates for
  # intervals that merged. Collapse to one span: first chromosome, min start,
  # max end. Never parse coordinates back out of the feature ID later.
  first_of  <- function(x) vapply(strsplit(x, ";", fixed = TRUE), function(p) p[1], character(1))
  min_of    <- function(x) vapply(strsplit(x, ";", fixed = TRUE), function(p) min(as.numeric(p)), numeric(1))
  max_of    <- function(x) vapply(strsplit(x, ";", fixed = TRUE), function(p) max(as.numeric(p)), numeric(1))

  coords <- data.frame(
    feature_id = as.character(fc$Geneid),
    chr        = first_of(as.character(fc$Chr)),
    start      = min_of(as.character(fc$Start)),
    end        = max_of(as.character(fc$End)),
    stringsAsFactors = FALSE
  )

  counts <- as.matrix(fc[, SAMPLES$sample, drop = FALSE])
  rownames(counts) <- coords$feature_id
  storage.mode(counts) <- "integer"
  message("  matrix: ", nrow(counts), " features x ", ncol(counts), " libraries")

  # --- restrict to standard chromosomes ------------------------------------
  # By chromosome, never by row order. The contig naming is mixed UCSC and
  # Ensembl, unplaced scaffolds sort first, and Interval_1 is on GL456210.1:
  # head() of this matrix is 5,000 rows of scaffold junk.
  keep_chr <- coords$chr %in% PARAMS$std_chroms
  message("  standard chromosomes: ", sum(keep_chr), " of ", nrow(counts),
          " (dropping ", sum(!keep_chr), " on unplaced scaffolds)")
  counts <- counts[keep_chr, , drop = FALSE]
  coords <- coords[keep_chr, , drop = FALSE]

  # --- signal used only to stratify the subsample --------------------------
  rowmean <- rowMeans(counts)
  cpm  <- t(t(counts) / (colSums(counts) / 1e6))
  lcpm <- log2(cpm + 1)
  tp_means <- vapply(split(seq_len(ncol(lcpm)), factor(SAMPLES$time, levels = sort(unique(SAMPLES$time)))),
                     function(i) rowMeans(lcpm[, i, drop = FALSE]), numeric(nrow(lcpm)))
  rng <- apply(tp_means, 1, function(x) max(x) - min(x))

  expressed <- rowmean >= PARAMS$min_mean_counts
  q <- quantile(rng[expressed], c(PARAMS$flat_max_q, PARAMS$dynamic_top_q), na.rm = TRUE)

  dyn_pool  <- which(expressed & rng >= q[2])
  flat_pool <- which(expressed & rng <= q[1])
  mid_pool  <- which(expressed & rng > q[1] & rng < q[2])
  low_pool  <- which(rowmean >= PARAMS$lowcount_min & rowmean < PARAMS$min_mean_counts)
  message("  pools: dynamic ", length(dyn_pool), ", flat ", length(flat_pool),
          ", middle ", length(mid_pool), ", low-count ", length(low_pool))

  # --- shape-stratify the dynamic pool -------------------------------------
  # Sampling the dynamic pool at random would over-represent whichever
  # trajectory shape happens to be most common and could drop a rare one
  # entirely. Cluster the z-scored profiles first, then sample within each
  # shape in proportion to its size, so every shape survives the subsample.
  zp <- t(scale(t(tp_means[dyn_pool, , drop = FALSE])))
  zp <- zp[is.finite(rowSums(zp)), , drop = FALSE]
  km <- kmeans(zp, centers = PARAMS$shape_k, nstart = 25, iter.max = 100)
  message("  dynamic-pool shapes: ", paste(km$size, collapse = ", "))

  take_prop <- function(ids, groups, n_total) {
    per <- split(ids, groups)
    want <- round(n_total * lengths(per) / length(ids))
    unlist(Map(function(g, k) sample(g, min(k, length(g))), per, want), use.names = FALSE)
  }
  sel_dyn  <- take_prop(rownames(zp), km$cluster, PARAMS$n_dynamic)
  sel_flat <- sample(rownames(counts)[flat_pool], min(PARAMS$n_flat, length(flat_pool)))
  sel_mid  <- sample(rownames(counts)[mid_pool],  min(PARAMS$n_middle, length(mid_pool)))
  sel_low  <- sample(rownames(counts)[low_pool],  min(PARAMS$n_lowcount, length(low_pool)))

  sel <- unique(c(sel_dyn, sel_flat, sel_mid, sel_low))
  # Keep genomic order: a demo BED that is sorted is easier to eyeball and to
  # load into a browser.
  sel <- rownames(counts)[rownames(counts) %in% sel]
  message("  selected ", length(sel), " features (",
          length(sel_dyn), " dynamic, ", length(sel_flat), " flat, ",
          length(sel_mid), " middle, ", length(sel_low), " low-count)")

  out_counts <- counts[sel, , drop = FALSE]
  out_coords <- coords[match(sel, coords$feature_id), ]
  stopifnot(identical(rownames(out_counts), out_coords$feature_id))

  # --- write ---------------------------------------------------------------
  cdf <- data.frame(feature_id = rownames(out_counts), out_counts,
                    check.names = FALSE, stringsAsFactors = FALSE)
  utils::write.table(cdf, file.path(outdir, "counts.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(SAMPLES, file.path(outdir, "samplesheet.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)

  # BED6. Coordinates carried as real columns from the source matrix, never
  # reconstructed from the feature ID (guardrail 9.4). featureCounts Start is
  # 1-based inclusive; BED start is 0-based half-open.
  bed <- data.frame(chr = out_coords$chr,
                    start = as.integer(out_coords$start - 1L),
                    end = as.integer(out_coords$end),
                    name = out_coords$feature_id,
                    score = 0L,
                    strand = ".",
                    stringsAsFactors = FALSE)
  stopifnot(!any(is.na(bed$start)), !any(is.na(bed$end)), all(bed$start < bed$end))
  utils::write.table(bed, file.path(outdir, "features.bed"),
                     sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  message("\nWrote:")
  for (f in c("counts.tsv", "samplesheet.tsv", "features.bed")) {
    p <- file.path(outdir, f)
    message(sprintf("  %-18s %8.1f KB", p, file.size(p) / 1024))
  }
}

main()
