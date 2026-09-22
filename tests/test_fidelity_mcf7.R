#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Fidelity: does chromadyn still produce the MCF7 trajectory classes?
#
#   Rscript tests/test_fidelity_mcf7.R
#
# OPT-IN. Needs the source project the method was extracted from, and SKIPS
# cleanly when it is absent, which is always the case in CI and on a fresh
# clone. Nothing from that project is committed here.
#
#   CHROMADYN_MCF7=/path/to/mcf7_project   override the default location
#   KEEP_SANDBOX=1                         keep the working directory
#
# What this checks, and what it does not.
#
# The published assignment maps 14 degPatterns clusters onto classes. Cluster
# IDs are meaningful only within one run over one feature set, so they cannot
# be compared across runs and the published override map cannot be applied to
# a subsample. What IS comparable is the per-feature class label, which is a
# statement about the shape of that feature's trajectory. So the test asks:
# starting from the same counts and the same parameters, does chromadyn put
# each peak in the same trajectory class the published analysis put it in?
#
# The comparison runs on a stratified subsample, because degPatterns cost is
# cubic and the full 12,703 peaks would take over an hour.
# ---------------------------------------------------------------------------

SRC <- Sys.getenv("CHROMADYN_MCF7", "/data/bchick/wproj/mcf7_project")
N_SUBSAMPLE <- 8000L
SEED <- 42L

DDS <- file.path(SRC, "analyses/05_atac_temporal_clustering/output/dds_consensus.rds")
PUB <- file.path(SRC, paste0("analyses/05_atac_temporal_clustering/output/tables/",
                             "atac_hrg_superclusters_v2.tsv"))

skip <- function(why) {
  cat(sprintf("SKIP  MCF7 fidelity: %s\n", why))
  cat("      This test is opt-in and needs the source project. Set CHROMADYN_MCF7\n")
  cat("      to its path if it lives elsewhere.\n")
  quit(status = 0)
}
if (!dir.exists(SRC)) skip(sprintf("source project not found at %s", SRC))
for (f in c(DDS, PUB)) if (!file.exists(f)) skip(sprintf("missing %s", f))

suppressPackageStartupMessages({
  library(SummarizedExperiment)
})

P <- 0L; F <- 0L
ok <- function(m) { cat(sprintf("  PASS %s\n", m)); P <<- P + 1L }
no <- function(m) { cat(sprintf("  FAIL %s\n", m)); F <<- F + 1L }
chk <- function(c, m) if (isTRUE(c)) ok(m) else no(m)

cat("chromadyn MCF7 fidelity\n")
cat(sprintf("  source: %s\n", SRC))

pub <- read.delim(PUB, stringsAsFactors = FALSE)
stopifnot(all(c("peak_id", "supercluster", "supercluster_v2") %in% colnames(pub)))
cat(sprintf("  published: %d peaks, %d v2 classes (%s)\n", nrow(pub),
            length(unique(pub$supercluster_v2)),
            paste(names(sort(table(pub$supercluster_v2), decreasing = TRUE)), collapse = ", ")))

dds <- readRDS(DDS)
cd <- as.data.frame(colData(dds))

# The HRG arm: the two shared t=0 libraries plus every hrg library. This is
# the split the source made at 01_degpatterns_clustering.Rmd:98.
keep <- cd$treatment %in% c("shared", "hrg")
cd <- cd[keep, ]
cts <- counts(dds)[, rownames(cd), drop = FALSE]

# Sample from ALL peaks, not only the published dynamic ones.
#
# An earlier version drew the subsample exclusively from the 12,703 peaks the
# published analysis called dynamic. That looks like the efficient choice and
# it invalidates the comparison: DESeq2 fits its dispersion trend across the
# features it is given, and a set composed entirely of strongly time-varying
# peaks produces an inflated trend and much less power. Forty per cent of
# already-published-dynamic peaks then failed to reach padj < 0.01, which
# says nothing about chromadyn and everything about the fixture.
#
# Sampling from the whole matrix keeps the dynamic fraction near what the
# original run saw. Peaks that are in the published table are the ones the
# comparison is made on; the rest are there to make the model honest.
set.seed(SEED)
pub <- pub[pub$peak_id %in% rownames(cts), ]
sel <- sample(rownames(cts), min(N_SUBSAMPLE, nrow(cts)))
# Guarantee every published class is represented even if the draw is unlucky.
per <- split(pub$peak_id, pub$supercluster_v2)
floor_n <- 150L
boost <- unlist(lapply(per, function(g) {
  have <- intersect(g, sel)
  if (length(have) >= floor_n) character(0) else sample(setdiff(g, have), floor_n - length(have))
}), use.names = FALSE)
sel <- unique(c(sel, boost))
cts <- cts[sel, , drop = FALSE]
truth <- setNames(pub$supercluster_v2[match(sel, pub$peak_id)], sel)
truth <- truth[!is.na(truth)]
cat(sprintf("  subsample: %d peaks x %d libraries, of which %d carry a published label (%.0f%%)\n",
            nrow(cts), ncol(cts), length(truth), 100 * length(truth) / nrow(cts)))

sandbox <- file.path(dirname(tempdir()), sprintf("chromadyn-fidelity-mcf7-%d", Sys.getpid()))
dir.create(file.path(sandbox, "data"), recursive = TRUE, showWarnings = FALSE)
finish <- function() {
  if (F == 0L && !identical(Sys.getenv("KEEP_SANDBOX"), "1")) {
    unlink(sandbox, recursive = TRUE)
  } else cat(sprintf("\nSandbox kept: %s\n", sandbox))
  cat(sprintf("\n----------------------------------------\n%d passed, %d failed\n", P, F))
  quit(status = as.integer(F > 0))
}

ss <- data.frame(
  sample = rownames(cd),
  group = ifelse(cd$treatment == "shared", "shared", "HRG"),
  time = as.numeric(as.character(cd$time_min)),
  replicate = as.character(cd$rep),
  stringsAsFactors = FALSE
)
ss <- ss[order(ss$time, ss$replicate), ]
cts <- cts[, ss$sample, drop = FALSE]
write.table(data.frame(feature_id = rownames(cts), cts, check.names = FALSE),
            file.path(sandbox, "data", "counts.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(ss, file.path(sandbox, "data", "samplesheet.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

cfg <- file.path(sandbox, "config.yaml")
preset <- readLines("config/presets/mcf7_as28.yaml")
preset <- sub("^  samplesheet: .*", sprintf("  samplesheet: %s",
              file.path(sandbox, "data", "samplesheet.tsv")), preset)
preset <- sub("^  counts: .*", sprintf("  counts: %s",
              file.path(sandbox, "data", "counts.tsv")), preset)
writeLines(c(preset, "output:", sprintf("  dir: %s", file.path(sandbox, "results"))), cfg)

cat("\nRunning chromadyn with config/presets/mcf7_as28.yaml\n")
runner <- if (nzchar(Sys.getenv("CONDA_PREFIX"))) character(0) else c("pixi", "run", "--frozen")
logf <- file.path(sandbox, "run.log")
st <- system2(c(runner, "snakemake")[1], c(tail(c(runner, "snakemake"), -1),
              "--configfile", cfg, "-j", "4"), stdout = logf, stderr = logf)
if (st != 0) {
  cat(paste(tail(readLines(logf), 25), collapse = "\n"), "\n")
  no("pipeline completed"); finish()
}
ok("pipeline completed")

cl <- read.delim(file.path(sandbox, "results", "clusters", "HRG_clusters.tsv"),
                 stringsAsFactors = FALSE)
cl <- cl[cl$supercluster_label != "Unassigned", ]
cl$truth <- truth[cl$feature_id]
cl <- cl[!is.na(cl$truth), ]

adj_rand <- function(a, b) {
  t <- table(a, b); n <- sum(t); cn <- function(x) sum(choose(x, 2))
  i <- cn(as.vector(t)); ea <- cn(rowSums(t)); eb <- cn(colSums(t))
  e <- ea * eb / choose(n, 2); (i - e) / ((ea + eb) / 2 - e)
}

tab <- table(cl$truth, cl$supercluster_label)
cat("\n  published class (rows) against chromadyn class (columns)\n\n")
print(tab)

named <- intersect(rownames(tab), colnames(tab))
agree <- if (length(named)) sum(diag(tab[named, named, drop = FALSE])) / nrow(cl) else 0
ari <- adj_rand(cl$truth, cl$supercluster_label)
cat(sprintf("\n  %d peaks compared | exact label agreement %.3f | ARI %.3f\n",
            nrow(cl), agree, ari))

# Floors are set below the observed values so noise cannot make this flap,
# while a real regression still trips it.
#
# Observed at calibration, 8,000-peak subsample, seed 42:
#   exact per-peak label agreement  0.930
#   adjusted Rand index             0.811
#   published names reproduced      5 of 5
#   per-class: Decreasing 100%, Transient 100%, Transient Increasing 100%,
#              Late Increasing 92%, Sustained Increasing 81%
#
# This is not a claim of bit-exact reproduction and cannot be: the run is on
# a subsample, so degPatterns sees a different feature set and cuts its tree
# differently. It is the stronger claim that matters, which is that starting
# from raw counts and the documented parameters, chromadyn puts nine peaks in
# ten into the same trajectory class the published analysis did.
chk(agree >= 0.85,
    sprintf("exact per-peak label agreement is %.3f (>= 0.85)", agree))
chk(ari >= 0.70, sprintf("ARI against the published assignment is %.3f (>= 0.70)", ari))
chk(length(named) == 5,
    sprintf("all 5 published class names are reproduced (%d): %s",
            length(named), paste(named, collapse = ", ")))

# Direction fidelity is the part that must hold regardless of granularity.
rising <- cl[grepl("Increasing", cl$truth, fixed = TRUE), ]
falling <- cl[cl$truth == "Decreasing", ]
bad <- sum(falling$supercluster_label != "Decreasing" &
             grepl("Increasing", falling$supercluster_label, fixed = TRUE))
chk(bad / max(nrow(falling), 1) <= 0.05,
    sprintf("published Decreasing peaks are not called Increasing (%d of %d)",
            bad, nrow(falling)))
bad2 <- sum(rising$supercluster_label == "Decreasing")
chk(bad2 / max(nrow(rising), 1) <= 0.05,
    sprintf("published Increasing peaks are not called Decreasing (%d of %d)",
            bad2, nrow(rising)))

for (k in rownames(tab)) {
  best <- colnames(tab)[which.max(tab[k, ])]
  cat(sprintf("  published '%s' -> mostly '%s' (%.0f%%)\n",
              k, best, 100 * max(tab[k, ]) / sum(tab[k, ])))
}
finish()
