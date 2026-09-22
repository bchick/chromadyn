# ---------------------------------------------------------------------------
# 09b_annotate.R: nearest gene and distance to TSS, per feature.
#
# Region mode with annotate.txdb set. Annotation packages are genome-specific
# and large, so they are not in the environment: install the one your genome
# needs and name it in the config. When txdb is null this rule is not part of
# the DAG at all.
# ---------------------------------------------------------------------------

if (!exists("snakemake")) {
  source("workflow/scripts/lib/standalone.R")
  snakemake <- cd_standalone(
    rule       = "annotate",
    wildcards  = list(arm = "WT"),
    configfile = "config/demo.yaml",
    input      = list(clusters = "clusters"),
    output     = list(annot = "annot"))
}
source(file.path(snakemake@scriptdir, "lib", "common.R"))
source(file.path(snakemake@scriptdir, "lib", "namespace_fix.R"))
cd_init(snakemake)

arm <- snakemake@wildcards$arm
txdb_name <- cfg_req("annotate.txdb")
orgdb_name <- cfg("annotate.orgdb")
window <- cfg_req("annotate.tss_window")

cd_assert(requireNamespace(txdb_name, quietly = TRUE),
          paste0("annotate.txdb is '%s' but that package is not installed. Annotation ",
                 "packages are genome-specific and are deliberately not in the chromadyn ",
                 "environment; install it, or set annotate.txdb: null to skip annotation."),
          txdb_name)
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicFeatures)
})

cl <- read_tsv_strict(snakemake@input$clusters,
                      c("feature_id", "supercluster_label", "chr", "start", "end"))
assert_coords(cl, sprintf("arm '%s' cluster table", arm))

txdb <- get(txdb_name, envir = asNamespace(txdb_name))
gr <- GenomicRanges::GRanges(cl$chr,
                             IRanges::IRanges(cl$start + 1L, cl$end),
                             feature_id = cl$feature_id)
tss <- GenomicFeatures::promoters(GenomicFeatures::genes(txdb), upstream = 0, downstream = 1)

hits <- GenomicRanges::distanceToNearest(gr, tss, ignore.strand = TRUE)
cl$nearest_gene_id <- NA_character_
cl$distance_to_tss <- NA_integer_
q <- S4Vectors::queryHits(hits)
sj <- S4Vectors::subjectHits(hits)
cl$nearest_gene_id[q] <- names(tss)[sj]
cl$distance_to_tss[q] <- S4Vectors::mcols(hits)$distance
cl$is_promoter <- !is.na(cl$distance_to_tss) & cl$distance_to_tss <= window

if (!is.null(orgdb_name) && requireNamespace(orgdb_name, quietly = TRUE)) {
  orgdb <- get(orgdb_name, envir = asNamespace(orgdb_name))
  sym <- tryCatch(
    AnnotationDbi::mapIds(orgdb, keys = unique(stats::na.omit(cl$nearest_gene_id)),
                          column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"),
    error = function(e) { message("[annotate] symbol lookup failed: ", conditionMessage(e)); NULL })
  if (!is.null(sym)) cl$nearest_gene_symbol <- unname(sym[cl$nearest_gene_id])
}

write_tsv_atomic(cl, snakemake@output$annot)
message(sprintf("[annotate] arm '%s': %d features, %d within %d bp of a TSS",
                arm, nrow(cl), sum(cl$is_promoter, na.rm = TRUE), window))
cd_note("n_promoter_proximal", sum(cl$is_promoter, na.rm = TRUE))
cd_note("txdb", txdb_name)
cd_finish()
