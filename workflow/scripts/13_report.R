# ---------------------------------------------------------------------------
# 13_report.R: render the HTML report.
#
# The script resolves paths and hands them to the Rmd as parameters, so the
# Rmd itself contains no path logic and can be knitted by hand against any
# results directory.
# ---------------------------------------------------------------------------

if (!exists("snakemake")) {
  source("workflow/scripts/lib/standalone.R")
  snakemake <- cd_standalone(
    rule       = "report",
    configfile = "config/demo.yaml",
    output     = list(html = "report"),
    params     = list(arms = "WT"))
}
source(file.path(snakemake@scriptdir, "lib", "common.R"))
source(file.path(snakemake@scriptdir, "lib", "namespace_fix.R"))
cd_init(snakemake)

rmd <- normalizePath(file.path(snakemake@scriptdir, "..", "report", "report.Rmd"),
                     mustWork = TRUE)
out <- snakemake@output$html
cd_mkdir(dirname(out))

params <- list(
  title = cfg("report.title", "timecourse-patterns results"),
  results_dir = normalizePath(cfg("output.dir"), mustWork = TRUE),
  arms = unlist(snakemake@params$arms),
  time_unit = cfg("input.time_unit"),
  config = snakemake@config
)

# Knit in a temporary directory so that intermediates never land in results/
# and a failed render cannot leave a half-written report behind.
tmp <- tempfile("timecourse-patterns-report-")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

rmarkdown::render(
  rmd,
  output_file = basename(out),
  output_dir = normalizePath(dirname(out)),
  intermediates_dir = tmp,
  knit_root_dir = normalizePath(dirname(out)),
  params = params,
  envir = new.env(parent = globalenv()),
  quiet = TRUE
)

cd_assert(file.exists(out), "report rendering produced no file at %s", out)
message(sprintf("[report] wrote %s (%.0f KB)", out, file.size(out) / 1024))
cd_note("report_bytes", as.numeric(file.size(out)))
cd_finish()
