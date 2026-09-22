# ---------------------------------------------------------------------------
# 02_prefilter.R: build the DESeqDataSet and drop features no assay could see.
#
# The design here is deliberately ~ 1. Each arm gets its own model, fitted in
# 05_differential.R against the shared baseline; this object exists to hold
# the counts, the colData and the prefilter decision in one place.
# ---------------------------------------------------------------------------

if (!exists("snakemake")) {
  source("workflow/scripts/lib/standalone.R")
  snakemake <- cd_standalone(
    rule       = "prefilter",
    configfile = "config/demo.yaml",
    input      = list(samplesheet = "demo/samplesheet.tsv", counts = "demo/counts.tsv"),
    output     = list(dds = "obj_dds"))
}
source(file.path(snakemake@scriptdir, "lib", "common.R"))
source(file.path(snakemake@scriptdir, "lib", "namespace_fix.R"))
cd_init(snakemake)

suppressPackageStartupMessages(library(DESeq2))

ss <- read_samplesheet(snakemake@input$samplesheet)
cts <- read_counts(snakemake@input$counts)

# Column order is the samplesheet's order, fixed once, here. Everything
# downstream indexes off this object rather than re-deriving an order.
cd_assert(setequal(ss$sample, colnames(cts)),
          "samplesheet and counts disagree on samples; 01_validate should have caught this.")
cts <- cts[, ss$sample, drop = FALSE]

# Guardrail 9.1 once more, at the point the factor is actually built. The
# validator checked the samplesheet; this checks the object the model sees.
coldata <- data.frame(
  sample = ss$sample,
  group = factor(ss$group),
  time = factor(ss$time, levels = sort(unique(ss$time))),
  replicate = factor(ss$replicate),
  row.names = ss$sample,
  stringsAsFactors = FALSE
)
for (extra in setdiff(colnames(ss), c("sample", "group", "time", "replicate"))) {
  coldata[[extra]] <- ss[[extra]]         # carried through for custom designs
}
assert_time_levels(coldata$time, "colData time")
assert_aligned(coldata, cts, "colData against counts")

storage.mode(cts) <- "integer"
dds <- DESeqDataSetFromMatrix(cts, colData = coldata, design = ~ 1)

min_mean <- cfg_req("filter.min_mean_counts")
keep <- rowMeans(counts(dds)) >= min_mean
n_in <- nrow(dds)
dds <- dds[keep, ]

message(sprintf("[prefilter] rowMeans >= %s keeps %d of %d features (dropped %d)",
                format(min_mean), sum(keep), n_in, n_in - sum(keep)))
cd_assert(nrow(dds) > 0,
          "filter.min_mean_counts = %s removed every feature.", format(min_mean))

cd_note("n_features_in", n_in)
cd_note("n_features_kept", as.integer(sum(keep)))
cd_note("n_features_dropped", as.integer(n_in - sum(keep)))
cd_note("min_mean_counts", min_mean)

save_obj(dds, snakemake@output$dds)
cd_finish()
