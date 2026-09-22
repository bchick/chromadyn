#===============================================================================
# theme.R: the one place chromadyn defines how a figure looks.
#===============================================================================
# Ported from mcf7_project/scripts/project_colors.R, with the dataset-specific
# palettes removed (ligand_colors, timepoint_colors, br2_condition_colors,
# experiment_colors, treatment_colors). What is kept is generic: the PRGn
# diverging system, the cnsplots geometry, and supercluster_colors.
#
# ONE DEFINITION PER PALETTE. The project this came from redefined a palette in
# a notebook after sourcing the shared file, with different hex values, so two
# panels of the same figure disagreed about what a condition's colour was.
# Nothing outside this file may define a palette. If you need a new one, add it
# here.
#
# THRESHOLD AND REFERENCE LINES use col_accent (#e6ab02), never red. Red reads
# as "error" and as "downregulated" at the same time, and it is the worst
# choice for the most common form of colour blindness.
#
# FIGURE SIZES ARE IN PIXELS AT 72 PX PER INCH, not inches. This is inherited
# from cnsplots and it is the single easiest thing to get wrong here:
# save_figure(p, path, width = 6.4) does not give you a 6.4 inch figure, it
# gives you a 25 x 13 PIXEL image. Use the FIG_W_* constants, or multiples of
# FIG_PANEL. Note also that `bg` defaults to transparent.
#
# COLOUR (from the source project, unchanged):
#   Diverging: PRGn (purple-green) for quantitative diverging data, z-scores
#   Categorical: teal / indigo / purple / green for discrete groups
#   Accent: harvest gold for significance and highlights
#
# STRUCTURE (ported from cnsplots, github.com/faridrashidi/cnsplots):
#   Journal-submission geometry and typography. cnsplots is a Python package,
#   so nothing is imported; its style spec is reimplemented here for ggplot2.
#
# Author: Brent Chick
#===============================================================================

# === PRGn Diverging Scale (full 11-class ColorBrewer) ===
prgn <- list(
  purple_dark = "#40004b",
  purple      = "#762a83",
  purple_mid  = "#9970ab",
  purple_lt   = "#c2a5cf",
  purple_pale = "#e7d4e8",
  neutral     = "#f7f7f7",
  green_pale  = "#d9f0d3",
  green_lt    = "#a6dba0",
  green_mid   = "#5aae61",
  green       = "#1b7837",
  green_dark  = "#00441b"
)

# === Structural & Accent Colors ===
col_bg        <- "#fdfdfd"   # Snow white (opaque fill; see theme bg = "paper")
col_text      <- "#333333"   # Charcoal (text & axes)
col_grid      <- "#F8F9FA"   # Cool grey (gridlines, when explicitly enabled)
col_accent    <- "#e6ab02"   # Harvest gold (significance markers, highlights)

# === Categorical Palette (for discrete groups/clusters) ===
cat_palette <- c(
  "#1b9e77",   # Deep Teal
  "#7570b3",   # Slate Indigo
  "#762a83",   # Deep Purple (matches PRGn low)
  "#1b7837"    # Deep Green (matches PRGn high)
)

# === DE direction colors ===
direction_colors <- c(
  Up   = prgn$green,     # #1b7837 (sustained / gained)
  Down = prgn$purple,    # #762a83 (transient / lost)
  NS   = "#CCCCCC"
)

# === Diverging heatmap palette function ===
# Full PRGn spectrum for smooth gradients (LFC, Z-scores)
div_palette <- colorRampPalette(c(
  prgn$purple, prgn$purple_mid, prgn$purple_lt, prgn$purple_pale,
  prgn$neutral,
  prgn$green_pale, prgn$green_lt, prgn$green_mid, prgn$green
))

# === Correlation / sequential heatmap palette ===
# Purple sequential, for non-diverging quantitative data (correlation, enrichment)
cor_palette <- colorRampPalette(c("#fcfbfd", "#dadaeb", "#9e9ac8", "#6a51a3", "#3f007d"))

#===============================================================================
# TYPOGRAPHY: cnsplots font stack
#===============================================================================
# cnsplots preference order is Helvetica > Helvetica Neue > Arial > DejaVu Sans.
# Resolved against installed families at source() time; Nimbus Sans is inserted
# because it is metrically identical to Helvetica and is what most Linux hosts
# actually ship. Falls back to the device default ("sans") if none are found.

col_font <- local({
  preferred <- c("Helvetica", "Helvetica Neue", "Arial",
                 "Nimbus Sans", "DejaVu Sans")
  installed <- tryCatch(
    if (requireNamespace("systemfonts", quietly = TRUE)) {
      unique(systemfonts::system_fonts()$family)
    } else character(0),
    error = function(e) character(0)
  )
  hit <- preferred[preferred %in% installed]
  if (length(hit)) hit[1] else "sans"
})

#===============================================================================
# FIGURE GEOMETRY: cnsplots pixel convention
#===============================================================================
# cnsplots specifies figure sizes in PIXELS at 72 px/inch, so 1 px == 1 pt and
# a pixel width maps directly onto a journal column width in mm. Use fig_px()
# to convert for ggsave()/knitr fig.width, or use save_figure() which does it.

PX_PER_INCH  <- 72

FIG_W_1COL   <- 255   # ~90 mm  single column (Nature 89 mm, Cell 85 mm)
FIG_W_1_5COL <- 390   # ~137 mm, one-and-a-half column
FIG_W_2COL   <- 540   # ~190 mm, full width (cnsplots multipanel_max_width)
FIG_PANEL    <- 150   # ~53 mm  cnsplots default single panel (150 x 150 px)

#' Convert cnsplots pixels to inches (for ggsave / knitr fig.width)
#'
#' @param px Numeric size in pixels at 72 px/inch.
#' @return Size in inches.
fig_px <- function(px) px / PX_PER_INCH

#===============================================================================
# PUBLICATION THEME: PRGn color on cnsplots structure
#===============================================================================
#' Publication ggplot2 theme
#'
#' Structural spec ported from cnsplots: no gridlines, transparent background,
#' left/bottom spines only at 0.5 pt, 2 pt ticks at 0.6 pt width, tight label
#' padding, frameless compact legend, bold centered title.
#'
#' Font sizes now DERIVE from `base_size` (they were hard-coded before, so
#' `base_size = 10` had no effect on axis text). cnsplots ratio: titles and axis
#' labels at `base_size`, tick labels and legend text at 0.875 * `base_size`.
#'
#' @param base_size Title / axis-label size in pt. Default 8 (cnsplots
#'   `title_fontsize`). Tick and legend text derive at 0.875 * base_size.
#' @param base_family Font family. Defaults to the resolved `col_font`.
#' @param grid Gridlines: "none" (cnsplots default), "y", "x", or "both".
#'   Restores the pre-2026-07-27 cool-grey grid for a single panel that needs it.
#' @param bg Background: "transparent" (cnsplots default) or "paper" (opaque
#'   `col_bg` snow white, the pre-2026-07-27 behavior).
#' @return A ggplot2 theme object.
theme_publication <- function(base_size   = 8,
                              base_family = col_font,
                              grid        = c("none", "y", "x", "both"),
                              bg          = c("transparent", "paper")) {

  grid <- match.arg(grid)
  bg   <- match.arg(bg)
  `%+replace%` <- ggplot2::`%+replace%`

  # cnsplots sizes: title_fontsize -> titles + axis labels; legend_fontsize
  # (7 of 8) -> tick labels, legend text, colorbar ticks.
  size_small <- base_size * 0.875

  # matplotlib line widths are pt; ggplot2 element_line(linewidth) is mm.
  pt_to_mm      <- function(pt) pt / ggplot2::.pt
  lw_spine      <- pt_to_mm(0.5)   # cnsplots axes_linewidth
  lw_tick       <- pt_to_mm(0.6)   # cnsplots {x,y}tick_major_width
  lw_grid       <- pt_to_mm(0.3)

  bg_fill  <- if (bg == "paper") col_bg else NA
  grid_maj <- ggplot2::element_line(color = col_grid, linewidth = lw_grid)

  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    ggplot2::theme(
      # --- Background (cnsplots savefig_transparent = TRUE) ---
      plot.background   = ggplot2::element_rect(fill = bg_fill, color = NA),
      panel.background  = ggplot2::element_rect(fill = bg_fill, color = NA),

      # --- Grid (cnsplots axes_grid = FALSE) ---
      panel.grid.major.y = if (grid %in% c("y", "both")) grid_maj else ggplot2::element_blank(),
      panel.grid.major.x = if (grid %in% c("x", "both")) grid_maj else ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),

      # --- Axes: left + bottom spines only (top/right off) ---
      axis.line         = ggplot2::element_line(color = col_text, linewidth = lw_spine),
      axis.ticks        = ggplot2::element_line(color = col_text, linewidth = lw_tick),
      axis.ticks.length = grid::unit(2, "pt"),          # cnsplots tick_major_size

      # tick_major_pad = 1 pt; axes_labelpad = 2 pt
      axis.text.x  = ggplot2::element_text(color = col_text, size = size_small,
                                           margin = ggplot2::margin(t = 1)),
      axis.text.y  = ggplot2::element_text(color = col_text, size = size_small,
                                           margin = ggplot2::margin(r = 1), hjust = 1),
      axis.title.x = ggplot2::element_text(color = col_text, size = base_size,
                                           margin = ggplot2::margin(t = 2)),
      axis.title.y = ggplot2::element_text(color = col_text, size = base_size,
                                           margin = ggplot2::margin(r = 2), angle = 90),

      # --- Strips (facets). No matplotlib analog; kept left-aligned bold. ---
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(color = col_text, size = size_small,
                                               face = "bold", hjust = 0,
                                               margin = ggplot2::margin(b = 2)),

      # --- Legend: frameless, compact (cnsplots legend_* settings) ---
      legend.background   = ggplot2::element_rect(fill = NA, color = NA),
      legend.key          = ggplot2::element_rect(fill = NA, color = NA),
      legend.key.size     = grid::unit(base_size * 1.1, "pt"),
      legend.key.spacing  = grid::unit(2, "pt"),        # handletextpad 0.3 * size
      legend.text         = ggplot2::element_text(color = col_text, size = size_small),
      legend.title        = ggplot2::element_text(color = col_text, size = base_size,
                                                  hjust = 0),
      legend.margin       = ggplot2::margin(0, 0, 0, 0),

      # --- Titles (cnsplots axes_titlelocation = "center", titleweight bold) ---
      plot.title    = ggplot2::element_text(color = col_text, size = base_size,
                                            face = "bold", hjust = 0.5,
                                            margin = ggplot2::margin(b = 4)),
      plot.subtitle = ggplot2::element_text(color = col_text, size = size_small,
                                            hjust = 0.5,
                                            margin = ggplot2::margin(b = 4)),
      plot.caption  = ggplot2::element_text(color = col_text, size = size_small,
                                            hjust = 1,
                                            margin = ggplot2::margin(t = 4)),

      # --- Panel labels (patchwork tags): bold, top-left ---
      plot.tag          = ggplot2::element_text(color = col_text,
                                                size = base_size + 2,
                                                face = "bold", hjust = 0),
      plot.tag.position = "topleft",

      # cnsplots savefig_pad_inches = 0.01 in ~ 0.7 pt
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
}

#===============================================================================
# EXPORT: cnsplots savefig conventions
#===============================================================================
#' Save a figure at journal size in Illustrator-editable formats
#'
#' Implements the cnsplots export spec: sizes in pixels at 72 px/inch,
#' transparent background, tight margins, 288 dpi rasters, and vector text left
#' as TEXT (svglite for SVG, cairo_pdf for PDF) so Illustrator can edit it
#' rather than receiving outlined paths.
#'
#' @param plot A ggplot (or any object ggsave accepts).
#' @param path Output path. Any extension is stripped; one file per `formats`.
#' @param width,height Size in PIXELS at 72 px/inch (see FIG_W_* constants).
#' @param formats Extensions to write. Default PDF + SVG + PNG.
#' @param dpi Raster resolution. Default 288 (cnsplots savefig_dpi = 72 * 4).
#' @param bg Background. Default "transparent" (cnsplots). Pass "white" for a
#'   PNG destined for Slack or a dark slide, where transparent renders as
#'   invisible dark-on-dark text.
#' @param fit If TRUE, ignore `width`/`height` and size the canvas to exactly
#'   fit the assembled figure. Only meaningful with `panel_px()`, since panels
#'   sized in relative units have no intrinsic size to measure. Falls back to
#'   `width`/`height` with a warning if measurement fails.
#' @return Character vector of written paths (invisibly).
save_figure <- function(plot, path,
                        width   = FIG_PANEL,
                        height  = FIG_PANEL,
                        formats = c("pdf", "svg", "png"),
                        dpi     = 288,
                        bg      = "transparent",
                        fit     = FALSE) {

  stem <- sub("\\.(pdf|svg|png|jpg|jpeg|tiff|eps)$", "", path, ignore.case = TRUE)
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)

  if (isTRUE(fit)) {
    got <- fig_size_px(plot)
    # A figure built from relative panels measures near zero: its panels are
    # null units that only resolve at draw time. Refuse to write a 4 px canvas.
    if (all(is.finite(got)) && got[["width"]] >= 20 && got[["height"]] >= 20) {
      width  <- got[["width"]]
      height <- got[["height"]]
    } else {
      warning("fit = TRUE could not measure the figure (are panels sized with ",
              "panel_px()?); falling back to width/height", call. = FALSE)
    }
  }

  w_in <- fig_px(width)
  h_in <- fig_px(height)

  written <- vapply(formats, function(fmt) {
    out <- paste0(stem, ".", fmt)
    dev <- switch(
      fmt,
      # svglite writes <text> elements with a font-family attribute, the
      # analog of cnsplots svg.fonttype = "none" (editable text, not paths).
      svg = if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite else NULL,
      # cairo_pdf subsets and embeds the font (~ pdf.fonttype = 42) and keeps
      # text as text; the base pdf() device would need font metrics we lack.
      pdf = grDevices::cairo_pdf,
      png = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else NULL,
      NULL
    )
    args <- list(filename = out, plot = plot, width = w_in, height = h_in,
                 units = "in", dpi = dpi, bg = bg)
    if (!is.null(dev)) args$device <- dev
    do.call(ggplot2::ggsave, args)
    out
  }, character(1))

  invisible(unname(written))
}

#===============================================================================
# MULTI-PANEL LAYOUT: cnsplots multipanel conventions
#===============================================================================
#' Lock panel boxes to an exact pixel size
#'
#' cnsplots `multipanel.panel(label, width, height)` sizes each panel in
#' absolute pixels. patchwork sizes panels relatively by default, as a share of
#' whatever is left after axis titles and legends are drawn, so a panel with a
#' long legend silently ends up smaller than its neighbor. This passes absolute
#' units through to `plot_layout()` so the PLOT BOX is exactly the size asked
#' for, independent of decorations around it.
#'
#' IMPORTANT: use a flat `design` grid, not nested `/` and `|` composition.
#' In a nested patchwork the absolute size applies to the outer cell and the
#' inner panels divide it up, which is not what the numbers say.
#'
#' @param width,height Panel size in PIXELS at 72 px/inch. Recycled across the
#'   grid, or pass a vector to size columns/rows individually.
#' @param design Optional patchwork design string, e.g. "AABB\\nCCDD". Strongly
#'   preferred over nested composition (see above).
#' @param ... Further arguments passed to `patchwork::plot_layout()`.
#' @return A patchwork layout spec to add to an assembled figure.
panel_px <- function(width = FIG_PANEL, height = FIG_PANEL, design = NULL, ...) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("panel_px() needs the patchwork package")
  }
  patchwork::plot_layout(
    widths  = grid::unit(width,  "pt"),   # 1 px == 1 pt at 72 px/inch
    heights = grid::unit(height, "pt"),
    design  = design,
    ...
  )
}

#' Measure an assembled figure in cnsplots pixels
#'
#' Returns the total size the figure wants, including axes, legends and panel
#' labels. Only meaningful once panels are sized with `panel_px()`; relative
#' panels are null units and measure as zero.
#'
#' @param plot A ggplot or patchwork figure.
#' @return Named numeric `c(width, height)` in pixels at 72 px/inch.
fig_size_px <- function(plot) {
  # Measuring needs an open device. Without one, R opens the default, which on
  # a headless Rscript run is pdf(): it does not know systemfonts families (so
  # it warns "font family not found in PostScript font database") and it leaves
  # an Rplots.pdf in the working directory. Measure on an in-memory ragg device
  # instead, and hand the active device back untouched.
  prev <- grDevices::dev.cur()
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_capture(width = 2000, height = 2000, units = "px", res = 72)
  } else {
    grDevices::pdf(NULL)
  }
  on.exit({
    grDevices::dev.off()
    if (prev > 1) grDevices::dev.set(prev)
  }, add = TRUE)

  g <- tryCatch(
    if (inherits(plot, "patchwork")) patchwork::patchworkGrob(plot)
    else ggplot2::ggplotGrob(plot),
    error = function(e) NULL
  )
  na <- c(width = NA_real_, height = NA_real_)
  if (is.null(g)) return(na)

  # Only the PANELS are null units under relative sizing; the surrounding axis
  # and label furniture always measures nonzero. So checking the total is not
  # enough to tell "absolutely sized" from "relative": measure the panel cells
  # and bail out if they resolve to nothing.
  cells <- g$layout[grepl("^panel", g$layout$name), , drop = FALSE]
  if (!nrow(cells)) return(na)
  panel_extent <- suppressWarnings(c(
    grid::convertWidth(sum(g$widths[min(cells$l):max(cells$r)]), "pt", valueOnly = TRUE),
    grid::convertHeight(sum(g$heights[min(cells$t):max(cells$b)]), "pt", valueOnly = TRUE)
  ))
  if (!all(is.finite(panel_extent)) || any(panel_extent < 1)) return(na)

  suppressWarnings(c(
    width  = grid::convertWidth(sum(g$widths),   "pt", valueOnly = TRUE),
    height = grid::convertHeight(sum(g$heights), "pt", valueOnly = TRUE)
  ))
}

#' Add A/B/C panel labels to a patchwork multi-panel figure
#'
#' cnsplots labels panels in bold Helvetica at the top-left of each panel.
#'
#' @param levels Tag sequence passed to patchwork. Default uppercase letters.
#' @return A patchwork annotation to add to an assembled plot.
panel_labels <- function(levels = "A") {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("panel_labels() needs the patchwork package")
  }
  patchwork::plot_annotation(tag_levels = levels)
}

#===============================================================================
# SIGNIFICANCE ANNOTATION: cnsplots pvalue_format
#===============================================================================
#' Significance stars (cnsplots pvalue_format = "star")
#'
#' @param p Numeric vector of p-values.
#' @return Character vector: ****, ***, **, * or ns.
p_stars <- function(p) {
  ifelse(is.na(p), NA_character_,
  ifelse(p <= 1e-4, "****",
  ifelse(p <= 1e-3, "***",
  ifelse(p <= 1e-2, "**",
  ifelse(p <= 0.05, "*", "ns")))))
}

#' Threshold p-value labels (cnsplots pvalue_format = "threshold")
#'
#' @param p Numeric vector of p-values.
#' @return Character vector: "P < 0.0001" ... "P > 0.05".
p_label <- function(p) {
  ifelse(is.na(p), NA_character_,
  ifelse(p <= 1e-4, "P < 0.0001",
  ifelse(p <= 1e-3, "P < 0.001",
  ifelse(p <= 1e-2, "P < 0.01",
  ifelse(p <= 0.05, "P < 0.05", "P > 0.05")))))
}

# === ATAC temporal supercluster colors (LOCKED, reuse in Figs 2–5) ===
# Semantics: purple = accessibility LOST; amber = transient spike;
# green gradient (light→dark) = de novo gain, increasing in how sustained/late.
# Sustained Increasing + Late Increasing are the duration-dependent de novo
# enhancer sets carried into Figs 3–5, keep these two greens fixed.
supercluster_colors <- c(
  "Decreasing"           = "#CC79A7",  # Wong pink   (lost)
  "Transient"            = "#E69F00",  # Wong amber  (transient spike)
  "Transient Increasing" = "#56B4E9",  # Wong sky    (gains, not held)
  "Sustained Increasing" = "#009E73",  # Wong green  (de novo, sustained; HERO)
  "Late Increasing"      = "#0072B2"   # Wong blue   (de novo, late; HERO)
)


#===============================================================================
# chromadyn additions
#===============================================================================

# === Canonical supercluster display order ===
# Every consumer orders classes by this, so a legend, a heatmap's row split and
# a ribbon panel in the same figure cannot disagree. Classes present in the
# data but absent here are appended alphabetically after these.
supercluster_order <- c(
  "Decreasing",
  "Transient",
  "Transient Increasing",
  "Sustained Increasing",
  "Late Increasing"
)

#' Order a vector of supercluster labels canonically.
#'
#' @param x Character or factor of supercluster labels.
#' @return A factor with levels in canonical order, unknown labels appended
#'   alphabetically and "Unassigned" always last.
supercluster_factor <- function(x) {
  x <- as.character(x)
  known <- supercluster_order[supercluster_order %in% x]
  rest <- sort(setdiff(unique(x), c(known, "Unassigned")))
  lv <- c(known, rest, if ("Unassigned" %in% x) "Unassigned")
  factor(x, levels = lv)
}

#' Colours for an arbitrary set of supercluster labels.
#'
#' Named classes keep their locked colour; anything else is filled from
#' cat_palette, and "Unassigned" is always neutral grey. This is what lets the
#' palette survive a run that produces classes the canonical five do not cover,
#' without silently recolouring the ones that matter.
supercluster_palette <- function(labels) {
  labels <- levels(supercluster_factor(labels))
  out <- stats::setNames(rep(NA_character_, length(labels)), labels)
  known <- intersect(labels, names(supercluster_colors))
  out[known] <- supercluster_colors[known]
  if ("Unassigned" %in% labels) out[["Unassigned"]] <- "#BBBBBB"
  gaps <- names(out)[is.na(out)]
  if (length(gaps)) {
    out[gaps] <- rep(cat_palette, length.out = length(gaps))
  }
  out
}

#' A sequential colour per timepoint, for any timepoint grid.
#'
#' Replaces the source project's hard-coded 0m/5m/.../240m vector, which only
#' worked for one experiment. Walks the PRGn ramp so that early is purple and
#' late is green, matching div_palette.
#'
#' @param times Numeric vector of timepoints; sorted and de-duplicated here.
timepoint_palette <- function(times) {
  times <- sort(unique(as.numeric(times)))
  cols <- div_palette(max(length(times), 2L))
  stats::setNames(cols[seq_along(times)], as.character(times))
}

#' Symmetric diverging colour ramp for a z-scored matrix.
#'
#' Clipped at the 99th percentile of |z| so that a handful of extreme features
#' do not flatten the whole heatmap to one colour.
circlize_ramp <- function(z, q = 0.99) {
  lim <- stats::quantile(abs(z[is.finite(z)]), q, na.rm = TRUE)
  if (!is.finite(lim) || lim <= 0) lim <- 1
  if (requireNamespace("circlize", quietly = TRUE)) {
    circlize::colorRamp2(c(-lim, 0, lim), c(prgn$purple, prgn$neutral, prgn$green))
  } else {
    grDevices::colorRampPalette(c(prgn$purple, prgn$neutral, prgn$green))(64)
  }
}

#' save_figure() for grid graphics, such as a ComplexHeatmap.
#'
#' ggsave() cannot take a Heatmap object, so this opens the same devices with
#' the same pixels-at-72-dpi geometry and draws into them. `draw_fn` is called
#' with no arguments and must do the drawing.
save_grid_figure <- function(draw_fn, path,
                             width = FIG_PANEL, height = FIG_PANEL,
                             formats = c("pdf", "svg", "png"),
                             dpi = 288, bg = "white") {
  stem <- sub("\\.(pdf|svg|png)$", "", path, ignore.case = TRUE)
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  w_in <- fig_px(width)
  h_in <- fig_px(height)

  supported <- c("pdf", "svg", "png")
  bad <- setdiff(formats, supported)
  if (length(bad)) {
    stop("save_grid_figure(): unsupported format(s): ", paste(bad, collapse = ", "),
         ". Supported: ", paste(supported, collapse = ", "), call. = FALSE)
  }

  open_device <- function(fmt, out) {
    switch(
      fmt,
      pdf = grDevices::cairo_pdf(out, width = w_in, height = h_in, bg = bg),
      svg = if (requireNamespace("svglite", quietly = TRUE)) {
        svglite::svglite(out, width = w_in, height = h_in, bg = bg)
      } else {
        grDevices::svg(out, width = w_in, height = h_in, bg = bg)
      },
      png = if (requireNamespace("ragg", quietly = TRUE)) {
        ragg::agg_png(out, width = w_in, height = h_in, units = "in",
                      res = dpi, background = bg)
      } else {
        grDevices::png(out, width = w_in, height = h_in, units = "in", res = dpi, bg = bg)
      }
    )
  }

  written <- character(0)
  for (fmt in formats) {
    out <- paste0(stem, ".", fmt)
    open_device(fmt, out)
    # Close THIS device whatever happens, without accumulating handlers across
    # iterations: an on.exit(add = TRUE) inside the loop would fire once per
    # format at function exit and close devices this call never opened.
    ok <- tryCatch({ draw_fn(); TRUE },
                   error = function(e) { grDevices::dev.off(); stop(e) })
    if (isTRUE(ok)) grDevices::dev.off()
    written <- c(written, out)
  }
  invisible(written)
}

invisible(TRUE)
