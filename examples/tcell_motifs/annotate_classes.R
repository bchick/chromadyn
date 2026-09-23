#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# annotate_classes.R: where in the genome does each trajectory class sit?
#
#   Rscript examples/tcell_motifs/annotate_classes.R [results_dir] [gencode.gtf]
#
# Annotates every peak with ChIPseeker against GENCODE vM35, collapses the
# categories to five that can be defended from sequence annotation alone, and
# writes a per-class composition table. Static peaks are included as the
# reference: a class is only notably distal or promoter-bound relative to
# accessible chromatin that does not change.
#
# No category is called "enhancer". That needs histone marks or perturbation
# data this analysis does not use; intronic and distal-intergenic peaks are
# reported as putative distal regulatory elements, which is what they are.
#
# Needs ChIPseeker and txdbmaker (Bioconductor); not part of the timecourse-patterns
# environment, since this is a demonstration rather than a pipeline rule.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ChIPseeker); library(GenomicRanges); library(txdbmaker)
})

args <- commandArgs(trailingOnly = TRUE)
res <- if (length(args) >= 1) args[1] else "examples/tcell_motifs/work/results"
gtf <- if (length(args) >= 2) args[2] else
  "examples/tcell_motifs/work/reference/gencode.vM35.primary_assembly.annotation.gtf"
work <- dirname(res)
PROMOTER_BP <- 1000       # promoter = within 1 kb of a TSS
FLANK_BP <- 3000          # promoter flank = 1 to 3 kb

txdb_cache <- file.path(work, "gencode_vM35.txdb.sqlite")
txdb <- if (file.exists(txdb_cache)) {
  AnnotationDbi::loadDb(txdb_cache)
} else {
  message("building TxDb from ", basename(gtf), " (once; cached)")
  t <- makeTxDbFromGFF(gtf, format = "gtf", organism = "Mus musculus")
  AnnotationDbi::saveDb(t, txdb_cache); t
}

cl <- utils::read.delim(file.path(res, "clusters", "WT_clusters.tsv"), stringsAsFactors = FALSE)
cl <- cl[cl$supercluster_label != "Unassigned", ]
d <- utils::read.delim(gzfile(file.path(res, "differential", "WT_results.tsv.gz")),
                       stringsAsFactors = FALSE)
bed <- utils::read.delim(file.path(work, "inputs", "features.bed"), header = FALSE,
                         col.names = c("chr", "start", "end", "name", "score", "strand"))
static <- bed[bed$name %in% d$feature_id[!as.logical(d$dynamic)], ]

regions <- rbind(
  data.frame(feature_id = cl$feature_id, chr = cl$chr, start = cl$start, end = cl$end,
             class = cl$supercluster_label),
  data.frame(feature_id = static$name, chr = static$chr, start = static$start,
             end = static$end, class = "Static")
)
gr <- GRanges(regions$chr, IRanges(regions$start + 1L, regions$end), feature_id = regions$feature_id)

ann <- annotatePeak(gr, TxDb = txdb, tssRegion = c(-FLANK_BP, FLANK_BP), verbose = FALSE)
a <- as.data.frame(ann)
stopifnot(identical(a$feature_id, regions$feature_id))

cat5 <- function(x, dist) {
  ifelse(grepl("^Promoter", x) & abs(dist) <= PROMOTER_BP, "Promoter (<=1 kb)",
  ifelse(grepl("^Promoter", x), "Promoter flank (1-3 kb)",
  ifelse(grepl("Exon|UTR", x), "Exon / UTR",
  ifelse(grepl("Intron", x), "Intronic",
         "Distal intergenic"))))
}
regions$category <- cat5(a$annotation, a$distanceToTSS)
regions$distance_to_tss <- a$distanceToTSS

tab <- as.data.frame(table(class = regions$class, category = regions$category),
                     stringsAsFactors = FALSE)
tab$fraction <- tab$Freq / ave(tab$Freq, tab$class, FUN = sum)
utils::write.table(tab, file.path("examples/tcell_motifs", "genomic_annotation.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(regions, file.path(work, "annotation_regions.rds"))

w <- stats::reshape(tab[, c("class", "category", "fraction")], idvar = "class",
                    timevar = "category", direction = "wide")
names(w) <- sub("^fraction\\.", "", names(w))
print(format(w, digits = 2), row.names = FALSE)
