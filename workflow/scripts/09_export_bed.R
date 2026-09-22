# ---------------------------------------------------------------------------
# 09_export_bed.R: one BED6 per trajectory class. Region mode only.
#
# Guardrail 9.4 in its most literal form. A greenlist variant of the source
# analysis stored peak_id as "chr1.7401731.7402231", because degPatterns had
# run make.names() over the rownames, then parsed coordinates back out with
# sub(":.*", "", peak_id). That regex does not match dot separators, so both
# exported BED files were thousands of rows of
# "chr1.7401731.7402231  NA  NA" and nobody noticed.
#
# chromadyn never parses an ID. Coordinates arrive as real columns from the
# features BED and are asserted valid immediately before writing.
# ---------------------------------------------------------------------------

if (!exists("snakemake")) {
  source("workflow/scripts/lib/standalone.R")
  snakemake <- cd_standalone(
    rule       = "export_bed",
    wildcards  = list(arm = "WT"),
    configfile = "config/demo.yaml",
    input      = list(clusters = "clusters"),
    output     = list(dir = "bed_dir", manifest = "bed_manifest"))
}
source(file.path(snakemake@scriptdir, "lib", "common.R"))
source(file.path(snakemake@scriptdir, "lib", "namespace_fix.R"))
cd_init(snakemake)

arm <- snakemake@wildcards$arm
cl <- read_tsv_strict(snakemake@input$clusters,
                      c("feature_id", "supercluster_label", "chr", "start", "end"))
outdir <- snakemake@output$dir
cd_mkdir(outdir)

cl <- cl[cl$supercluster_label != "Unassigned", , drop = FALSE]
assert_coords(cl, sprintf("arm '%s' cluster table", arm))

rows <- list()
for (sc in sort(unique(cl$supercluster_label))) {
  sub <- cl[cl$supercluster_label == sc, , drop = FALSE]
  sub <- sub[order(sub$chr, sub$start), ]
  bed <- data.frame(chr = sub$chr, start = as.integer(sub$start), end = as.integer(sub$end),
                    name = sub$feature_id, score = 0L, strand = ".",
                    stringsAsFactors = FALSE)
  # Assert again on exactly what is about to be written, not on what was read.
  cd_assert(!anyNA(bed$start) && !anyNA(bed$end),
            "arm '%s', class '%s': NA coordinates would be written to BED.", arm, sc)
  cd_assert(all(bed$start < bed$end),
            "arm '%s', class '%s': %d row(s) with start >= end.",
            arm, sc, sum(bed$start >= bed$end))
  path <- file.path(outdir, sprintf("%s_%s.bed", arm, cd_slug(sc)))
  utils::write.table(bed, path, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)
  rows[[length(rows) + 1L]] <- data.frame(
    arm = arm, supercluster_label = sc, file = basename(path), n = nrow(bed),
    stringsAsFactors = FALSE)
}

manifest <- do.call(rbind, rows)
write_tsv_atomic(manifest, snakemake@output$manifest)
message(sprintf("[export_bed] arm '%s': %d BED files, %d regions total",
                arm, nrow(manifest), sum(manifest$n)))
cd_note("n_bed_files", nrow(manifest))
cd_finish()
