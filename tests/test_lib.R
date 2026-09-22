# ---------------------------------------------------------------------------
# tests/test_lib.R: unit tests for workflow/scripts/lib/common.R.
#
# These cover the guardrails from specification section 9, each of which is a
# real silent failure from the source analyses. Run from the repo root:
#   Rscript tests/test_lib.R
# Called by tests/test_static.sh (tier 1); needs no data and no pipeline run.
# ---------------------------------------------------------------------------

source("workflow/scripts/lib/common.R")
.CD$config <- yaml::read_yaml("config/config.yaml")
.CD$paths  <- yaml::read_yaml("workflow/paths.yaml")

P <- 0; F <- 0
chk <- function(name, expr, want_error = TRUE) {
  e <- tryCatch({expr; NULL}, error = function(e) conditionMessage(e))
  hit <- if (want_error) !is.null(e) else is.null(e)
  if (hit) { cat("  pass ", name, "\n"); P <<- P+1 } else {
    cat("  FAIL ", name, if (!is.null(e)) paste0(" -> ", e) else " (no error raised)", "\n"); F <<- F+1 }
}

cat("== cfg ==\n")
chk("cfg dotted access", stopifnot(cfg("cluster.degpatterns.minc") == 50), FALSE)
chk("cfg null key returns NULL", stopifnot(is.null(cfg("input.features"))), FALSE)
chk("cfg missing -> default", stopifnot(cfg("no.such.key", "dflt") == "dflt"), FALSE)
chk("cfg_req on null key errors", cfg_req("input.features"))

cat("== paths ==\n")
chk("cd_path fills arm", stopifnot(cd_path("clusters", arm="WT") == "results/clusters/WT_clusters.tsv"), FALSE)
chk("cd_path unfilled placeholder errors", cd_path("clusters"))
chk("cd_path unknown key errors", cd_path("nope"))
chk("cd_fig expands formats", stopifnot(length(cd_fig("pca")) == 3), FALSE)
chk("cd_slug", stopifnot(cd_slug("Late Increasing") == "late_increasing"), FALSE)

cat("== guardrail 9.1 time ordering ==\n")
chk("numeric-ordered levels OK", assert_time_levels(factor(c(0,30,120), levels=c(0,30,120))), FALSE)
chk("lexicographic levels caught", assert_time_levels(factor(as.character(c(0,30,120,240)))))
chk("non-numeric levels caught", assert_time_levels(factor(c("naive","d3"))))

cat("== guardrail 9.2 alignment ==\n")
m <- matrix(1:6, nrow=2, dimnames=list(c("f1","f2"), c("s1","s2","s3")))
md_ok  <- data.frame(t=1:3, row.names=c("s1","s2","s3"))
md_ord <- data.frame(t=1:3, row.names=c("s2","s1","s3"))
md_dif <- data.frame(t=1:3, row.names=c("s1","s2","sX"))
chk("aligned passes", assert_aligned(md_ok, m), FALSE)
chk("reordered caught", assert_aligned(md_ord, m))
chk("different samples caught", assert_aligned(md_dif, m))

cat("== guardrail 9.4 coordinates ==\n")
good <- data.frame(chr="chr1", start=c(100,200), end=c(150,250))
chk("good coords pass", assert_coords(good), FALSE)
chk("NA coords caught", assert_coords(data.frame(chr="chr1", start=c(100,NA), end=c(150,250))))
chk("start>=end caught", assert_coords(data.frame(chr="chr1", start=300, end=250)))
chk("missing cols caught", assert_coords(data.frame(chr="chr1", start=1)))

cat("== guardrail 9.5 map coverage ==\n")
chk("full map passes", assert_map_covers(list(`1`="A",`2`="B"), c(1,2,1)), FALSE)
chk("missing cluster caught", assert_map_covers(list(`1`="A"), c(1,2,3)))

cat("== helpers ==\n")
mm <- matrix(c(1,2, 3,4, 10,20, 30,40), nrow=2,
             dimnames=list(c("f1","f2"), c("a","b","c","d")))
tmres <- timepoint_means(mm, c(0,0,5,5))
chk("timepoint_means shape", stopifnot(identical(dim(tmres), c(2L,2L)), colnames(tmres)==c("0","5")), FALSE)
z <- zscore_rows(rbind(a=c(1,2,3), b=c(5,5,5)))
chk("zscore flat row -> 0", stopifnot(all(z["b",]==0), attr(z,"n_zero_variance")==1), FALSE)

cat(sprintf("\n%d passed, %d failed\n", P, F))
quit(status = as.integer(F > 0))
