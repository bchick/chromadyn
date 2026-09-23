#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# make_figure.R: one figure, three independent views of each trajectory class.
#
#   A  the trajectory itself (chromadyn's ribbon panel)
#   B  where in the genome the class sits, against static peaks
#   C  which known motifs are enriched in it, against static peaks
#
# B and C use information the clustering never saw: it was given counts and
# timepoints only. If classes defined purely by the shape of their
# accessibility over time also differ in genomic context and in the
# transcription factors whose motifs they carry, the classes are biology and
# not an artefact of the clustering.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
source("workflow/scripts/lib/theme.R")
source("workflow/scripts/lib/plots.R")

work <- "examples/tcell_motifs/work"
odir <- "examples/tcell_motifs"

# Split pieces renamed by their actual shape; see README.md.
RENAME <- c("Increasing-1" = "Transient Increasing", "Increasing-2" = "Late Increasing",
            "Increasing-3" = "Dip and recovery")
ORDER <- c("Decreasing", "Transient", "Transient Increasing", "Late Increasing",
           "Dip and recovery")
relabel <- function(x) ifelse(x %in% names(RENAME), RENAME[x], x)

fit <- readRDS(file.path(work, "results/objects/WT_supercluster_fit.rds"))
cl <- readRDS(file.path(work, "results/objects/WT_degpatterns.rds"))
a <- fit$assignment
a$supercluster_label <- relabel(a$supercluster_label)
pal <- supercluster_palette(ORDER)

# --- A: trajectories ---------------------------------------------------------
long <- supercluster_long(cl$profiles, a)
long$supercluster_label <- factor(long$supercluster_label, levels = ORDER)
pA <- plot_supercluster_ribbon(long, NULL, "days") +
  scale_colour_manual(values = pal) + scale_fill_manual(values = pal) +
  facet_wrap(~facet, nrow = 1) +
  labs(tag = "A", title = "Trajectory classes, full WT timecourse (54,282 dynamic peaks)")

# --- B: genomic context --------------------------------------------------------
ann <- utils::read.delim(file.path(odir, "genomic_annotation.tsv"), stringsAsFactors = FALSE)
ann$class <- relabel(ann$class)
ann$class <- factor(ann$class, levels = rev(c(ORDER, "Static")))
cats <- c("Promoter (<=1 kb)", "Promoter flank (1-3 kb)", "Exon / UTR", "Intronic",
          "Distal intergenic")
ann$category <- factor(ann$category, levels = rev(cats))
cat_cols <- stats::setNames(c(prgn$purple, prgn$purple_lt, "#BBBBBB", prgn$green_lt, prgn$green),
                            cats)
static_prom <- sum(ann$fraction[ann$class == "Static" &
                                  grepl("^Promoter \\(", ann$category)])
pB <- ggplot(ann, aes(.data$fraction, .data$class, fill = .data$category)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = static_prom, linetype = "dashed", colour = col_accent, linewidth = 0.4) +
  scale_fill_manual(values = cat_cols, breaks = cats, name = NULL) +
  scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"), expand = c(0, 0)) +
  labs(x = "Fraction of peaks", y = NULL, tag = "B",
       title = "Genomic context",
       subtitle = "Dashed line: promoter fraction of static peaks") +
  theme_publication() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2, reverse = FALSE))

# --- C: motifs -----------------------------------------------------------------
mot <- utils::read.delim(file.path(odir, "motif_enrichment.tsv"), stringsAsFactors = FALSE)
Q_MAX <- 1e-5; TOP <- 6
sig <- mot[mot$q < Q_MAX, ]
top <- unique(unlist(lapply(split(sig, sig$class), function(d) utils::head(d$motif[order(d$q)], TOP))))
g <- expand.grid(motif = top, class = ORDER, stringsAsFactors = FALSE)
m <- merge(g, mot[, c("class", "motif", "enr", "q")], all.x = TRUE)
m$l2 <- log2(pmax(ifelse(is.na(m$enr), 1, m$enr), 1e-3))
m$sig <- !is.na(m$q) & m$q < Q_MAX
best <- tapply(seq_len(nrow(m)), m$motif, function(i) {
  j <- i[which.max(m$l2[i])]; match(m$class[j], ORDER) * 100 - m$l2[j] })
m$motif <- factor(m$motif, levels = rev(names(sort(best))))
m$class <- factor(m$class, levels = ORDER)
lim <- max(abs(m$l2))
pC <- ggplot(m, aes(.data$class, .data$motif, fill = .data$l2)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_point(data = m[m$sig, ], size = 0.6, colour = col_text) +
  scale_fill_gradientn(colours = div_palette(100), limits = c(-lim, lim),
                       name = "log2 enrichment\nvs static") +
  labs(x = NULL, y = NULL, tag = "C", title = "Known motif enrichment",
       subtitle = sprintf("SEA, JASPAR2024; dot: q < %g", Q_MAX)) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        axis.text.y = element_text(size = 6))

fig <- pA / (pB | pC) + plot_layout(heights = c(1, 2.1))
save_figure(fig, file.path(odir, "tcell_classes_overview"),
            width = FIG_W_2COL * 1.35, height = FIG_W_2COL * 1.15,
            formats = c("png", "pdf"), dpi = 300, bg = "white")
dir.create("docs/img", showWarnings = FALSE)
file.copy(file.path(odir, "tcell_classes_overview.png"), "docs/img/tcell_classes_overview.png",
          overwrite = TRUE)
message("wrote ", file.path(odir, "tcell_classes_overview.png"))
