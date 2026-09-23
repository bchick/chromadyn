# ---------------------------------------------------------------------------
# 12_manifest.R: everything needed to explain how this result happened.
#
# The fully resolved config, every package version, every seed, sessionInfo()
# and SHA-256 of each input file. A result that cannot be traced back to its
# parameters is not a result, and the source analysis could not be: its
# clustering was reproducible only through a cached RDS, because degPatterns
# was never seeded even though the methods section said it was.
# ---------------------------------------------------------------------------

if (!exists("snakemake")) {
  source("workflow/scripts/lib/standalone.R")
  snakemake <- cd_standalone(
    rule       = "manifest",
    configfile = "config/demo.yaml",
    output     = list(manifest = "manifest"))
}
source(file.path(snakemake@scriptdir, "lib", "common.R"))
source(file.path(snakemake@scriptdir, "lib", "namespace_fix.R"))
cd_init(snakemake)

sha256 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  if (requireNamespace("digest", quietly = TRUE)) {
    digest::digest(path, algo = "sha256", file = TRUE)
  } else NA_character_
}

inputs <- Filter(Negate(is.null), list(
  samplesheet = cfg("input.samplesheet"),
  counts = cfg("input.counts"),
  features = cfg("input.features")
))
input_info <- lapply(names(inputs), function(k) {
  p <- inputs[[k]]
  list(role = k, path = p, exists = file.exists(p),
       bytes = if (file.exists(p)) as.numeric(file.size(p)) else NA_real_,
       sha256 = sha256(p))
})

pkgs <- c("DESeq2", "DEGreport", "ComplexHeatmap", "limma", "SummarizedExperiment",
          "GenomicRanges", "rtracklayer", "ggplot2", "cluster", "mclust",
          "digest", "yaml", "jsonlite")
versions <- stats::setNames(lapply(pkgs, function(p) {
  tryCatch(as.character(utils::packageVersion(p)), error = function(e) NULL)
}), pkgs)
versions <- versions[!vapply(versions, is.null, logical(1))]

# Merge the per-rule provenance fragments each script wrote as it finished.
frag_dir <- file.path(cfg("output.dir"), "logs", "prov")
frags <- list()
if (dir.exists(frag_dir)) {
  for (f in sort(list.files(frag_dir, pattern = "\\.json$", full.names = TRUE))) {
    frags[[sub("\\.json$", "", basename(f))]] <-
      tryCatch(jsonlite::fromJSON(f), error = function(e) list(error = conditionMessage(e)))
  }
}
seeds <- Filter(length, lapply(frags, function(x) x$seeds))

manifest <- list(
  tool = list(name = "timecourse-patterns", version = "0.1.0"),
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  config = snakemake@config,
  inputs = input_info,
  seeds = seeds,
  steps = frags,
  versions = c(list(R = R.version.string), versions),
  platform = list(os = R.version$os, arch = R.version$arch,
                  sysname = unname(Sys.info()[["sysname"]])),
  session_info = utils::capture.output(utils::sessionInfo())
)

cd_mkdir(dirname(snakemake@output$manifest))
writeLines(jsonlite::toJSON(manifest, auto_unbox = TRUE, null = "null",
                            pretty = TRUE, force = TRUE),
           snakemake@output$manifest)
message(sprintf("[manifest] %d steps, %d input file(s), %d package versions",
                length(frags), length(input_info), length(versions)))
cd_finish()
