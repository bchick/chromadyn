# ---------------------------------------------------------------------------
# 10_figures.R: every per-arm figure.
#
# The plotting itself lives in lib/plots.R. This script decides what gets
# drawn and where it goes, and nothing else, so that a figure can be redrawn
# in an R session without running a rule.
#
# Writes a flag file rather than naming every format of every figure as a
# Snakemake output, so that changing figures.formats does not multiply the
# DAG. The flag lists what was actually written.
# ---------------------------------------------------------------------------

if (!exists("snakemake")) {
  source("workflow/scripts/lib/standalone.R")
  snakemake <- cd_standalone(
    rule       = "figures",
    wildcards  = list(arm = "WT"),
    configfile = "config/demo.yaml",
    input      = list(cluster_obj = "obj_degpatterns", fit = "obj_supercluster",
                      results = "differential", transformed = "obj_transformed",
                      dds = "obj_dds"),
    output     = list(flag = "fig_flag"))
}
source(file.path(snakemake@scriptdir, "lib", "common.R"))
source(file.path(snakemake@scriptdir, "lib", "namespace_fix.R"))
source(file.path(snakemake@scriptdir, "lib", "theme.R"))
source(file.path(snakemake@scriptdir, "lib", "plots.R"))
cd_init(snakemake)

suppressPackageStartupMessages(library(ggplot2))

arm <- snakemake@wildcards$arm
cl <- load_obj(snakemake@input$cluster_obj)
fit <- load_obj(snakemake@input$fit)
res <- read_tsv_strict(snakemake@input$results, c("feature_id", "dynamic"))
res$dynamic <- as.logical(res$dynamic)
vst_all <- load_obj(snakemake@input$transformed)
dds <- load_obj(snakemake@input$dds)

formats <- cfg("figures.formats")
dpi <- cfg("figures.dpi")
unit <- cfg("input.time_unit")
written <- character(0)
emit <- function(paths) written <<- c(written, paths)

long <- supercluster_long(cl$profiles, fit$assignment)
n_sc <- length(unique(long$supercluster_label))

# --- dynamic against static -------------------------------------------------
emit(save_figure(
  plot_dynamic_vs_static(res, arm, sprintf("%s: dynamic features", arm)),
  cd_path("fig_dynamic", arm = arm, fmt = formats[1]),
  width = FIG_PANEL, height = FIG_PANEL, formats = formats, dpi = dpi))

# --- the native degPatterns grid -------------------------------------------
# deg$plot is a ggplot built by DEGreport. Kept as-is: it is the figure a
# reader of the original method will recognise.
if (!is.null(cl$deg$plot) && inherits(cl$deg$plot, "ggplot")) {
  emit(save_figure(
    cl$deg$plot + labs(title = sprintf("%s: degPatterns clusters", arm)),
    cd_path("fig_degpatterns", arm = arm, fmt = formats[1]),
    width = FIG_W_2COL, height = FIG_W_1COL, formats = formats, dpi = dpi))
} else {
  message("[figures] degPatterns returned no plot object, skipping that panel")
}

# --- dendrogram, with the cut drawn ----------------------------------------
cut_h <- NULL
if (cfg("superclusters.method") == "hclust") {
  kk <- cfg("superclusters.k")
  hh <- cfg("superclusters.h")
  cut_h <- if (!is.null(hh)) hh else {
    hts <- sort(fit$hc$height, decreasing = TRUE)
    if (!is.null(kk) && kk >= 2 && kk <= length(hts) + 1) mean(hts[c(kk - 1, kk)]) else NULL
  }
}
emit(save_figure(
  plot_cluster_dendrogram(fit$hc, fit$cluster_to_sc, cut_h,
                          sprintf("%s: cluster dendrogram", arm)),
  cd_path("fig_dendrogram", arm = arm, fmt = formats[1]),
  width = FIG_W_1_5COL, height = FIG_PANEL, formats = formats, dpi = dpi))

# --- trajectories -----------------------------------------------------------
panel_w <- max(FIG_W_1COL, FIG_PANEL * n_sc)
emit(save_figure(
  plot_supercluster_traces(long, sprintf("%s: superclusters", arm), unit),
  cd_path("fig_traces", arm = arm, fmt = formats[1]),
  width = panel_w, height = FIG_PANEL, formats = formats, dpi = dpi))

emit(save_figure(
  plot_supercluster_ribbon(long, sprintf("%s: temporal superclusters", arm), unit),
  cd_path("fig_ribbon", arm = arm, fmt = formats[1]),
  width = panel_w, height = FIG_PANEL, formats = formats, dpi = dpi))

# --- heatmap ----------------------------------------------------------------
# Z-scored for display, outside degPatterns, exactly as the source did:
# columns ordered by time then replicate, rows split by trajectory class.
if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  a <- fit$assignment[fit$assignment$supercluster_label != "Unassigned", , drop = FALSE]
  a <- a[order(supercluster_factor(a$supercluster_label), a$cluster), ]
  coldata <- as.data.frame(SummarizedExperiment::colData(dds))[cl$samples, , drop = FALSE]
  ord <- order(as.numeric(as.character(coldata$time)), coldata$replicate)
  z <- zscore_rows(vst_all[a$feature_id, cl$samples[ord], drop = FALSE])
  rs <- supercluster_factor(a$supercluster_label)
  cs <- factor(as.numeric(as.character(coldata$time))[ord],
               levels = sort(unique(as.numeric(as.character(coldata$time)))))
  hm <- ComplexHeatmap::Heatmap(
    z, name = "z",
    col = circlize_ramp(z),
    row_split = rs, column_split = cs,
    cluster_columns = FALSE, cluster_row_slices = FALSE, cluster_rows = FALSE,
    show_row_names = FALSE, show_column_names = TRUE,
    column_names_gp = grid::gpar(fontsize = 5),
    row_title_gp = grid::gpar(fontsize = 6),
    column_title_gp = grid::gpar(fontsize = 7),
    use_raster = TRUE, raster_quality = 2,
    heatmap_legend_param = list(labels_gp = grid::gpar(fontsize = 5),
                                title_gp = grid::gpar(fontsize = 6)))
  emit(save_grid_figure(function() ComplexHeatmap::draw(hm),
                        cd_path("fig_heatmap", arm = arm, fmt = formats[1]),
                        width = FIG_W_1_5COL, height = FIG_W_1_5COL,
                        formats = formats, dpi = dpi))
} else {
  message("[figures] ComplexHeatmap not installed, skipping the heatmap")
}

flag <- snakemake@output$flag
cd_mkdir(dirname(flag))
writeLines(sort(written), flag)
message(sprintf("[figures] arm '%s': wrote %d files across %d panels",
                arm, length(written), length(written) / length(formats)))
cd_note("n_figure_files", length(written))
cd_finish()
