#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# build_inputs.R: the full WT T cell ATAC timecourse as chromadyn inputs.
#
# Same source, libraries and conventions as demo/build_demo.R (see
# demo/PROVENANCE.md), but every peak on a standard chromosome rather than a
# 5,000-peak subsample, so that each trajectory class holds enough regions for
# motif enrichment to be meaningful.
#
#   Rscript examples/tcell_motifs/build_inputs.R <featureCounts matrix> <outdir>
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
mpath <- if (length(args) >= 1) args[1] else "demo/_source_featureCounts.txt"
outdir <- if (length(args) >= 2) args[2] else "examples/tcell_motifs/work/inputs"
if (!file.exists(mpath)) {
  stop("Source matrix not found: ", mpath,
       "\nPass the path as the first argument. See demo/PROVENANCE.md.", call. = FALSE)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# The demo's nine libraries and design, reused verbatim so the two runs are
# directly comparable.
ss <- utils::read.delim("demo/samplesheet.tsv", stringsAsFactors = FALSE)
std_chroms <- paste0("chr", c(1:19, "X", "Y"))

fc <- utils::read.delim(mpath, skip = 1, check.names = FALSE, stringsAsFactors = FALSE)
colnames(fc) <- sub("\\.mLb\\.clN\\.sorted\\.bam$", "", colnames(fc))
stopifnot(all(ss$sample %in% colnames(fc)))

part <- function(x, f) vapply(strsplit(as.character(x), ";", fixed = TRUE), f, numeric(1))
chr <- vapply(strsplit(as.character(fc$Chr), ";", fixed = TRUE), `[`, character(1), 1)
start <- part(fc$Start, function(p) min(as.numeric(p)))
end <- part(fc$End, function(p) max(as.numeric(p)))

keep <- chr %in% std_chroms            # by chromosome, never by row order
counts <- as.matrix(fc[keep, ss$sample])
storage.mode(counts) <- "integer"
ids <- as.character(fc$Geneid[keep])

utils::write.table(data.frame(feature_id = ids, counts, check.names = FALSE),
                   file.path(outdir, "counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
file.copy("demo/samplesheet.tsv", file.path(outdir, "samplesheet.tsv"), overwrite = TRUE)
bed <- data.frame(chr = chr[keep], start = as.integer(start[keep] - 1L), end = as.integer(end[keep]),
                  name = ids, score = 0L, strand = ".")
stopifnot(all(bed$start < bed$end), !anyNA(bed$start))
utils::write.table(bed, file.path(outdir, "features.bed"), sep = "\t", quote = FALSE,
                   row.names = FALSE, col.names = FALSE)
message(sprintf("wrote %d peaks x %d libraries to %s", nrow(counts), ncol(counts), outdir))
