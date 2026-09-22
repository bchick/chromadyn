# ---------------------------------------------------------------------------
# standalone.R: run any chromadyn script outside Snakemake, for debugging.
#
# Snakemake injects an S4 `snakemake` object holding @input, @output, @params,
# @wildcards, @config and @scriptdir. When you want to step through a script in
# an R session instead, its header block builds an equivalent object with this:
#
#   if (!exists("snakemake")) {
#     source("workflow/scripts/lib/standalone.R")
#     snakemake <- cd_standalone(
#       rule      = "cluster",
#       wildcards = list(arm = "WT"),
#       configfile = "config/demo.yaml",
#       input     = list(transformed = "obj_transformed"),
#       output    = list(clusters = "clusters", obj = "obj_degpatterns"))
#   }
#
# `input` and `output` name keys from workflow/paths.yaml, resolved against the
# same config and wildcards the rule would use. A value containing "/" is taken
# as a literal path instead, which is how input files named in the config
# (counts, samplesheet, features) are passed.
#
# The point of naming paths.yaml keys rather than writing paths out is that the
# debug entry point cannot drift from the rule: if a path changes in
# paths.yaml, both move together.
#
# Run from the repository root.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(yaml))

if (!methods::isClass("Snakemake")) {
  methods::setClass(
    "Snakemake",
    representation(
      input = "list", output = "list", params = "list", wildcards = "list",
      threads = "numeric", log = "list", resources = "list", config = "list",
      rule = "character", bench_iteration = "numeric", scriptdir = "character",
      source = "function"
    )
  )
}

#' Deep-merge `over` onto `base`, as Snakemake merges a --configfile over the
#' workflow's own configfile.
.cd_merge <- function(base, over) {
  for (k in names(over)) {
    base[[k]] <- if (is.list(over[[k]]) && is.list(base[[k]])) {
      .cd_merge(base[[k]], over[[k]])
    } else {
      over[[k]]
    }
  }
  base
}

cd_standalone <- function(rule,
                          wildcards = list(),
                          configfile = NULL,
                          input = list(),
                          output = list(),
                          params = list(),
                          repo_root = ".",
                          threads = 1) {

  base_cfg <- file.path(repo_root, "config", "config.yaml")
  if (!file.exists(base_cfg)) {
    stop("cd_standalone(): ", base_cfg, " not found. Run from the repository root, ",
         "or pass repo_root=.", call. = FALSE)
  }
  config <- yaml::read_yaml(base_cfg)
  if (!is.null(configfile)) {
    cf <- if (file.exists(configfile)) configfile else file.path(repo_root, configfile)
    if (!file.exists(cf)) stop("cd_standalone(): configfile not found: ", configfile, call. = FALSE)
    config <- .cd_merge(config, yaml::read_yaml(cf))
  }

  paths <- yaml::read_yaml(file.path(repo_root, "workflow", "paths.yaml"))

  resolve <- function(spec) {
    if (grepl("/", spec, fixed = TRUE)) return(spec)          # literal path
    tmpl <- paths[[spec]]
    if (is.null(tmpl)) {
      stop("cd_standalone(): '", spec, "' is neither a path (no '/') nor a key in ",
           "workflow/paths.yaml.", call. = FALSE)
    }
    subs <- c(list(out = config$output$dir), wildcards)
    for (nm in names(subs)) {
      tmpl <- gsub(paste0("{", nm, "}"), as.character(subs[[nm]]), tmpl, fixed = TRUE)
    }
    if (grepl("\\{[a-z_]+\\}", tmpl)) {
      stop("cd_standalone(): unfilled placeholder resolving '", spec, "' -> '", tmpl,
           "'. Add it to `wildcards`.", call. = FALSE)
    }
    tmpl
  }

  methods::new(
    "Snakemake",
    input     = lapply(input,  resolve),
    output    = lapply(output, resolve),
    params    = params,
    wildcards = wildcards,
    threads   = as.numeric(threads),
    log       = list(),
    resources = list(),
    config    = config,
    rule      = rule,
    bench_iteration = NA_real_,
    scriptdir = normalizePath(file.path(repo_root, "workflow", "scripts"), mustWork = TRUE),
    source    = function(...) invisible(NULL)
  )
}

invisible(TRUE)
