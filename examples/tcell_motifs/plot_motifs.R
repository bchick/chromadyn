#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# plot_motifs.R: summarize SEA results as one motif-by-class figure.
#
#   Rscript examples/tcell_motifs/plot_motifs.R [motifs_dir] [out_dir]
#
# For each class, the top motifs by SEA q-value (enriched over static peaks)
# are pooled; the figure shows the enrichment ratio of every pooled motif in
# every class, so a motif specific to one temporal program stands out as a
# single bright cell in its column.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(ggplot2))
source("workflow/scripts/lib/theme.R")

args <- commandArgs(trailingOnly = TRUE)
mdir <- if (length(args) >= 1) args[1] else "examples/tcell_motifs/work/motifs"
odir <- if (length(args) >= 2) args[2] else "examples/tcell_motifs"
TOP_PER_CLASS <- 6
Q_MAX <- 1e-5

# the workflow's labels for this run, and the shape each actually has. The split
# named its three pieces numerically because two share a late rise and one is
# not an increasing shape at all; see README.md for the profiles.
LABELS <- c(
  "decreasing"   = "Decreasing",
  "transient"    = "Transient",
  "increasing-1" = "Transient Increasing",
  "increasing-2" = "Late Increasing",
  "increasing-3" = "Dip and recovery"
)
ORDER <- c("Decreasing", "Transient", "Transient Increasing", "Late Increasing",
           "Dip and recovery")

read_sea <- function(slug) {
  f <- file.path(mdir, paste0("sea_", slug), "sea.tsv")
  if (!file.exists(f)) return(NULL)
  x <- utils::read.delim(f, comment.char = "#", stringsAsFactors = FALSE)
  if (!nrow(x)) return(NULL)
  data.frame(class = LABELS[[slug]], motif = x$ALT_ID, id = x$ID,
             enr = x$ENR_RATIO, q = x$QVALUE, tp = x$TP, fp = x$FP,
             stringsAsFactors = FALSE)
}
all <- do.call(rbind, lapply(names(LABELS), read_sea))
stopifnot(!is.null(all))
utils::write.table(all[order(all$class, all$q), ], file.path(odir, "motif_enrichment.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

sig <- all[all$q < Q_MAX, ]
top <- unique(unlist(lapply(split(sig, sig$class), function(d)
  utils::head(d$motif[order(d$q)], TOP_PER_CLASS))))
grid <- expand.grid(motif = top, class = ORDER, stringsAsFactors = FALSE)
m <- merge(grid, all[, c("class", "motif", "enr", "q")], all.x = TRUE)
m$enr[is.na(m$enr)] <- 1
m$log2enr <- log2(pmax(m$enr, 1e-3))
m$sig <- !is.na(m$q) & m$q < Q_MAX
# Order motifs by the class in which they are most enriched, so the figure
# reads as a staircase down the diagonal.
best <- tapply(seq_len(nrow(m)), m$motif, function(i) {
  j <- i[which.max(m$log2enr[i])]; match(m$class[j], ORDER) * 100 - m$log2enr[j]
})
m$motif <- factor(m$motif, levels = rev(names(sort(best))))
m$class <- factor(m$class, levels = ORDER)

lim <- max(abs(m$log2enr))
p <- ggplot(m, aes(.data$class, .data$motif, fill = .data$log2enr)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_point(data = m[m$sig, ], size = 0.7, colour = col_text) +
  scale_fill_gradientn(colours = div_palette(100), limits = c(-lim, lim),
                       name = "log2 enrichment\nvs static peaks") +
  labs(x = NULL, y = NULL,
       title = "Known motifs enriched per trajectory class",
       subtitle = sprintf("SEA, JASPAR2024 vertebrates; dot = q < %g", Q_MAX)) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_figure(p, file.path(odir, "motif_enrichment"),
            width = FIG_W_1_5COL, height = max(FIG_W_1COL, 11 * length(top) + 90),
            formats = c("png", "pdf"), bg = "white")
message(sprintf("%d motifs across %d classes -> %s", length(top), length(ORDER),
                file.path(odir, "motif_enrichment.png")))
