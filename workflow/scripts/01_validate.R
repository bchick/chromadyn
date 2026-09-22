# ---------------------------------------------------------------------------
# 01_validate.R: check the inputs before anything expensive runs.
#
# Every check here corresponds to a real silent failure from the source
# analyses this method was extracted from. The rule is simple: an input problem
# must surface here, named, with the offending values printed, rather than
# three rules later as a confident wrong answer.
#
# Writes results/qc/validation_report.tsv with one row per check (PASS, WARN or
# FAIL and a message) and results/qc/library_sizes.tsv, then exits non-zero if
# anything FAILed.
# ---------------------------------------------------------------------------

if (!exists("snakemake")) {                              # standalone debug entry
  source("workflow/scripts/lib/standalone.R")
  snakemake <- cd_standalone(
    rule       = "validate",
    configfile = "config/demo.yaml",
    input      = list(samplesheet = "demo/samplesheet.tsv",
                      counts      = "demo/counts.tsv"),
    output     = list(report = "validation", libsizes = "libsizes"))
}
source(file.path(snakemake@scriptdir, "lib", "common.R"))
source(file.path(snakemake@scriptdir, "lib", "namespace_fix.R"))
cd_init(snakemake)

# --- check accumulator ------------------------------------------------------
CHECKS <- list()
record <- function(name, status, message) {
  CHECKS[[length(CHECKS) + 1L]] <<- data.frame(
    check = name, status = status, message = message, stringsAsFactors = FALSE)
  invisible(NULL)
}
pass <- function(name, fmt, ...) record(name, "PASS", sprintf(fmt, ...))
warn <- function(name, fmt, ...) record(name, "WARN", sprintf(fmt, ...))
fail <- function(name, fmt, ...) record(name, "FAIL", sprintf(fmt, ...))

# --- load -------------------------------------------------------------------
ss  <- read_samplesheet(snakemake@input$samplesheet)
cts <- read_counts(snakemake@input$counts)
pass("load", "samplesheet %d rows; counts %d features x %d samples",
     nrow(ss), nrow(cts), ncol(cts))

# --- sample identity: never silently drop a sample --------------------------
only_ss  <- setdiff(ss$sample, colnames(cts))
only_cts <- setdiff(colnames(cts), ss$sample)
if (length(only_ss) || length(only_cts)) {
  fail("sample_match",
       paste0("samplesheet and counts columns disagree. In samplesheet only: %s. ",
              "In counts only: %s. Sample IDs must match exactly; chromadyn will not ",
              "guess a mapping or quietly analyse a subset."),
       if (length(only_ss)) paste(only_ss, collapse = ", ") else "(none)",
       if (length(only_cts)) paste(only_cts, collapse = ", ") else "(none)")
} else {
  pass("sample_match", "all %d samplesheet IDs are present as counts columns", nrow(ss))
}

if (anyDuplicated(ss$sample)) {
  fail("sample_unique", "duplicated sample IDs in samplesheet: %s",
       paste(unique(ss$sample[duplicated(ss$sample)]), collapse = ", "))
} else {
  pass("sample_unique", "sample IDs are unique")
}

# --- guardrail 9.1: numeric time ordering -----------------------------------
# The destructive one. A character time column sorts lexicographically to
# 0, 120, 15, 240, 30, 60; every trajectory becomes noise and nothing errors.
times <- sort(unique(ss$time))
tf <- factor(ss$time, levels = times)
ok <- tryCatch({ assert_time_levels(tf, "samplesheet time"); TRUE },
               error = function(e) { fail("time_numeric", "%s", conditionMessage(e)); FALSE })
if (ok) {
  pass("time_numeric", "%d timepoints parse as numbers and are ordered: %s (%s)",
       length(times), paste(times, collapse = ", "), cfg("input.time_unit", "unitless"))
}
if (length(times) < 3L) {
  fail("time_count",
       "only %d distinct timepoint(s): %s. Trajectory clustering needs at least 3 to have a shape.",
       length(times), paste(times, collapse = ", "))
} else if (length(times) < 4L) {
  warn("time_count",
       "only 3 timepoints. Shapes will be coarse and the transient and late classes may not separate.")
} else {
  pass("time_count", "%d timepoints", length(times))
}

# --- arms and the shared baseline ------------------------------------------
baseline <- cfg_req("input.baseline_group")
arms <- setdiff(unique(ss$group), baseline)
n_baseline <- sum(ss$group == baseline)
if (!length(arms)) {
  fail("arms", "every row has group == baseline_group ('%s'), so there is nothing to analyse.",
       baseline)
} else {
  pass("arms", "%d arm(s): %s%s", length(arms), paste(arms, collapse = ", "),
       if (n_baseline) sprintf("; %d shared baseline libraries joined to each", n_baseline) else "")
}
bad_names <- arms[!grepl("^[A-Za-z0-9._+-]+$", arms)]
if (length(bad_names)) {
  fail("arm_names", "arm name(s) unusable in filenames: %s. Use letters, digits, dot, underscore, plus or hyphen.",
       paste(bad_names, collapse = ", "))
} else {
  pass("arm_names", "arm names are safe as filename components")
}

# --- per-arm replication, and guardrail 9.10 --------------------------------
# DESeq(test = "LRT") on ~ time needs replication within timepoint. A
# single-replicate timecourse is a common input shape and fails obscurely
# inside DESeq2, so catch it here and name the supported alternative.
for (arm in arms) {
  sub <- ss[ss$group %in% c(arm, baseline), ]
  tab <- table(sub$time)
  missing_tp <- setdiff(times, as.numeric(names(tab)))
  if (length(missing_tp)) {
    fail(sprintf("arm_%s_timepoints", arm),
         "arm '%s' has no libraries at timepoint(s) %s.",
         arm, paste(missing_tp, collapse = ", "))
  }
  if (max(tab) == 1L) {
    if (cfg("differential.method") == "lrt") {
      fail(sprintf("arm_%s_replication", arm),
           paste0("arm '%s' has one library per timepoint, so DESeq2's LRT over ~ time ",
                  "cannot estimate dispersion. Set differential.method: none to rank by ",
                  "variance and take differential.top_n instead. That path is supported ",
                  "and documented in docs/parameters.md."),
           arm)
    } else {
      warn(sprintf("arm_%s_replication", arm),
           "arm '%s' has one library per timepoint. Fine for differential.method: %s, but no dispersion is estimated.",
           arm, cfg("differential.method"))
    }
  } else if (min(tab) == 1L) {
    warn(sprintf("arm_%s_replication", arm),
         "arm '%s' has %d replicate(s) at some timepoints and %d at others (%s). DESeq2 handles this, but the sparse timepoints carry less weight.",
         arm, min(tab), max(tab), paste(sprintf("t=%s:n=%d", names(tab), as.integer(tab)), collapse = ", "))
  } else {
    pass(sprintf("arm_%s_replication", arm),
         "arm '%s': %d libraries, %s", arm, nrow(sub),
         paste(sprintf("t=%s:n=%d", names(tab), as.integer(tab)), collapse = ", "))
  }
}

# --- counts sanity ----------------------------------------------------------
if (anyNA(cts)) {
  fail("counts_na", "counts contain %d NA value(s). chromadyn will not impute them.", sum(is.na(cts)))
} else {
  pass("counts_na", "no NA values in counts")
}
if (cfg("transform.method") != "none") {
  nonint <- sum(abs(cts - round(cts)) > 1e-9)
  if (nonint > 0L) {
    fail("counts_integer",
         paste0("%d count value(s) are not whole numbers, but transform.method is '%s', ",
                "which expects raw counts. Set transform.method: none if this matrix is ",
                "already normalized."),
         nonint, cfg("transform.method"))
  } else {
    pass("counts_integer", "counts are whole numbers, as %s requires", cfg("transform.method"))
  }
  if (min(cts) < 0) {
    fail("counts_negative", "counts contain negative values (min %s).", format(min(cts)))
  } else {
    pass("counts_negative", "no negative counts")
  }
}
if (anyDuplicated(rownames(cts))) {
  fail("feature_unique", "duplicated feature IDs in counts: %d",
       sum(duplicated(rownames(cts))))
} else {
  pass("feature_unique", "%d unique feature IDs", nrow(cts))
}

# --- guardrail 9.8: minc against input size ---------------------------------
# minc = 50 is meaningless on 2,000 features. The check is against the feature
# count that will actually reach clustering, approximated here by the count
# surviving the prefilter, since select_dynamic has not run yet.
n_after_prefilter <- sum(rowMeans(cts) >= cfg("filter.min_mean_counts"))
minc <- cfg("cluster.degpatterns.minc")
if (identical(minc, "auto")) {
  pass("minc", "minc: auto resolves to %d for %d prefiltered features",
       max(15L, round(n_after_prefilter / 100)), n_after_prefilter)
} else {
  minc <- as.integer(minc)
  ceiling_minc <- n_after_prefilter / 20
  if (minc > ceiling_minc) {
    fail("minc",
         paste0("cluster.degpatterns.minc is %d but only about %d features survive the ",
                "prefilter, so a cluster would need %.0f%% of the data to be kept. Try ",
                "minc: %d, or minc: auto."),
         minc, n_after_prefilter, 100 * minc / n_after_prefilter,
         max(15L, round(n_after_prefilter / 100)))
  } else {
    pass("minc", "minc %d is workable against roughly %d prefiltered features (cap %.0f)",
         minc, n_after_prefilter, ceiling_minc)
  }
}

# --- prefilter would not empty the matrix -----------------------------------
if (n_after_prefilter == 0L) {
  fail("prefilter", "filter.min_mean_counts = %s removes every feature.",
       format(cfg("filter.min_mean_counts")))
} else if (n_after_prefilter < nrow(cts) * cfg("filter.warn_below_fraction")) {
  warn("prefilter", "filter.min_mean_counts = %s keeps only %d of %d features (%.1f%%).",
       format(cfg("filter.min_mean_counts")), n_after_prefilter, nrow(cts),
       100 * n_after_prefilter / nrow(cts))
} else {
  pass("prefilter", "prefilter keeps %d of %d features (%.1f%%)",
       n_after_prefilter, nrow(cts), 100 * n_after_prefilter / nrow(cts))
}

# --- region mode ------------------------------------------------------------
features <- cfg("input.features")
if (is.null(features)) {
  pass("region_mode", "gene mode: no coordinates supplied, BED export and annotation are skipped")
} else {
  ok <- tryCatch({
    bed <- read_bed(features, feature_ids = rownames(cts))
    unmatched <- setdiff(rownames(cts), bed$name)
    if (length(unmatched)) {
      fail("region_mode",
           "features BED is missing %d of %d feature IDs, first: %s",
           length(unmatched), nrow(cts), paste(utils::head(unmatched, 5), collapse = ", "))
    } else {
      pass("region_mode", "region mode: %d BED records cover every feature (%s)",
           nrow(bed), paste(utils::head(sort(unique(bed$chr)), 3), collapse = ", "))
    }
    TRUE
  }, error = function(e) { fail("region_mode", "%s", conditionMessage(e)); FALSE })
}

# --- batch column, only if it will be used ----------------------------------
if (cfg("batch_correct.method") != "none") {
  bcol <- cfg_req("batch_correct.column")
  if (!bcol %in% colnames(ss)) {
    fail("batch", "batch_correct.method is '%s' but the samplesheet has no '%s' column.",
         cfg("batch_correct.method"), bcol)
  } else if (length(unique(ss[[bcol]])) < 2L) {
    fail("batch", "batch column '%s' has a single value, so there is nothing to correct.", bcol)
  } else {
    # A batch that cannot be separated from time must not be "corrected": the
    # correction would delete the trajectory instead of the artefact. The
    # general test is whether ~ batch + time is full rank, which is exactly
    # what limma and DESeq2 would choke on. Counting timepoints per batch is
    # not enough: a batch spanning t=0,3,5 against another holding only t=8 is
    # still perfectly confounded, because that one indicator column equals the
    # t=8 indicator.
    mm <- stats::model.matrix(~ b + t,
                              data = data.frame(b = factor(ss[[bcol]]),
                                                t = factor(ss$time, levels = times)))
    if (qr(mm)$rank < ncol(mm)) {
      fail("batch",
           paste0("batch column '%s' is confounded with time: the model matrix for ",
                  "~ batch + time is rank %d with %d columns, so the two cannot be ",
                  "separated. Correcting this batch would remove the trajectory along ",
                  "with it. Batch levels per timepoint: %s"),
           bcol, qr(mm)$rank, ncol(mm),
           paste(vapply(split(as.character(ss[[bcol]]), ss$time),
                        function(x) paste0(unique(x), collapse = "+"), character(1)),
                 collapse = ", "))
    } else {
      pass("batch", "batch column '%s' has %d levels and ~ batch + time is full rank",
           bcol, length(unique(ss[[bcol]])))
    }
  }
} else {
  pass("batch", "batch correction off")
}

# --- write ------------------------------------------------------------------
report <- do.call(rbind, CHECKS)
write_tsv_atomic(report, snakemake@output$report)

libsizes <- data.frame(
  sample = colnames(cts),
  group = ss$group[match(colnames(cts), ss$sample)],
  time = ss$time[match(colnames(cts), ss$sample)],
  replicate = ss$replicate[match(colnames(cts), ss$sample)],
  total_counts = colSums(cts),
  detected_features = colSums(cts > 0),
  stringsAsFactors = FALSE
)
write_tsv_atomic(libsizes, snakemake@output$libsizes)

n_fail <- sum(report$status == "FAIL")
n_warn <- sum(report$status == "WARN")
cd_note("n_checks", nrow(report))
cd_note("n_fail", n_fail)
cd_note("n_warn", n_warn)
cd_note("arms", arms)
cd_note("timepoints", times)

message(sprintf("[validate] %d checks: %d PASS, %d WARN, %d FAIL",
                nrow(report), sum(report$status == "PASS"), n_warn, n_fail))
if (n_warn) for (i in which(report$status == "WARN")) {
  message(sprintf("  WARN %s: %s", report$check[i], report$message[i]))
}
if (n_fail) {
  for (i in which(report$status == "FAIL")) {
    message(sprintf("  FAIL %s: %s", report$check[i], report$message[i]))
  }
  cd_finish()
  stop(sprintf("validation failed: %d check(s). See %s", n_fail, snakemake@output$report),
       call. = FALSE)
}
cd_finish()
