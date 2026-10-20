## Shared pie-chart helpers for R/plot_pies.R (step 19) and extra/compare_pies.R (one-off comparisons).
## Needs `cfg` (args_config()) and utils.R. Defines: CT, G (host nuclear chromosome sizes), order_all, viral, excluded, nuc,
## pal_for(), fmt(), draw_pie().
CT <- rd(file.path(cfg$work, "annot", "contigs.tsv"))
G <- rd(file.path(cfg$work, "annot", "genome.sizes"), header = FALSE, col.names = c("contig", "length"))
if (is.null(CT) || is.null(G)) stop("work/annot/contigs.tsv missing: run step 1 first")

## genome order: chr1..chr22, chrX, chrY, other host, excluded host contigs (chrM), viral genomes
chrom_rank <- function(x) {
  n <- suppressWarnings(as.numeric(sub("^chr", "", x)))
  ifelse(!is.na(n), n, ifelse(x == "chrX", 100, ifelse(x == "chrY", 101, 200)))
}
CT[, rank := chrom_rank(contig) + fifelse(group == "excluded", 1000, 0)]
CT[group == "viral", rank := 2000 + seq_len(.N)]                  # viral genomes in config order
CT <- CT[order(rank, contig)]
order_all <- CT$contig
viral <- CT[group == "viral", contig]; excluded <- CT[group == "excluded", contig]
nuc <- CT[group == "nuclear", contig]
VCOL <- c("#d62728", "#ff7f0e", "#e377c2", "#bcbd22", "#17becf")
pal_for <- function(contigs) {
  hosts <- setdiff(contigs, c(viral, excluded))
  ## alternate light and dark tones so neighbouring chromosomes stay distinguishable
  hc <- grDevices::hcl(h = seq(15, 375, length.out = length(hosts) + 1)[-1], c = 45, l = rep(c(62, 78), length.out = length(hosts)))
  cols <- c(setNames(hc, hosts), setNames(rep("#4d4d4d", length(intersect(contigs, excluded))), intersect(contigs, excluded)),
            setNames(VCOL[seq_along(intersect(contigs, viral))], intersect(contigs, viral)))
  cols[contigs]
}
fmt <- function(x) ifelse(x >= 1, sprintf("%.1f%%", x), ifelse(x >= 0.01, sprintf("%.2f%%", x), ifelse(x > 0, sprintf("%.1e%%", x), "0%")))

## d: panel, contig, value (slice size), length; one pie per panel. ring = draw the length-share ring inside
draw_pie <- function(d, contigs, title, subtitle, file, ring, value_label = "share", legend_extra = NULL) {
  d <- copy(d)[contig %in% contigs]
  d[, share := 100 * value / sum(value), by = panel]
  d[, len_share := 100 * length / sum(length), by = panel]
  d[, contig := factor(contig, contigs)]
  setorder(d, panel, contig)
  ## legend label: contig and its share in every panel
  lg <- dcast(d, contig ~ panel, value.var = "share")
  pn <- setdiff(names(lg), "contig")
  lg[, lab := paste0(contig, ":  ", do.call(paste, c(lapply(pn, function(p) paste(if (length(pn) > 1) p else "", fmt(get(p)))), sep = " | ")))]
  if (!is.null(legend_extra)) lg[, lab := paste0(lab, legend_extra[as.character(contig)])]
  labs_map <- setNames(trimws(gsub("\\s+", " ", lg$lab)), as.character(lg$contig))
  cols <- pal_for(contigs)
  ## slices as explicit rectangles (angle = cumulative share), clockwise from 12 o'clock in genome order, so the labels
  ## sit on the right slices. Outer ring = values; inner ring = length shares.
  d[, `:=`(y1 = cumsum(share) - share, y2 = cumsum(share)), by = panel]
  d[, `:=`(l1 = cumsum(len_share) - len_share, l2 = cumsum(len_share)), by = panel]
  ox <- if (ring) c(1.05, 2) else c(0, 2)
  p <- ggplot(d) + geom_rect(aes(xmin = ox[1], xmax = ox[2], ymin = y1, ymax = y2, fill = contig), colour = "white", linewidth = 0.2)
  if (ring) p <- p + geom_rect(aes(xmin = 0.35, xmax = 1, ymin = l1, ymax = l2, fill = contig), colour = "white", linewidth = 0.2, alpha = 0.55)
  ## names on slices of at least 2%
  big <- d[share >= 2]
  if (nrow(big)) p <- p + geom_text(data = big, aes(x = if (ring) 1.55 else 1.35, y = (y1 + y2) / 2, label = paste0(contig, "\n", fmt(share))),
                                    size = 2.3, lineheight = 0.85, colour = "grey10")
  n_pan <- uniqueN(d$panel); n_col <- min(n_pan, 4); n_row <- ceiling(n_pan / n_col)
  p <- p + coord_polar(theta = "y", direction = 1) + facet_wrap(~panel, ncol = n_col) +
    scale_x_continuous(limits = c(0, 2.05), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 100 + 1e-6), expand = c(0, 0)) +
    scale_fill_manual(values = cols, labels = labs_map[names(cols)], name = NULL, drop = FALSE) +
    labs(title = title, subtitle = paste(c(subtitle, if (ring) "inner ring: share of contig length (an even spread)"), collapse = "; ")) +
    theme_void(10) + theme(legend.text = element_text(size = 7), legend.key.size = unit(0.35, "cm"), strip.text = element_text(face = "bold", size = 10),
                           plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 8))
  n_leg <- length(contigs)
  ggsave(file, p, width = 5.2 * n_col + 4, height = max(5.5 * n_row, 0.17 * n_leg + 1.5), dpi = 150, bg = "white", limitsize = FALSE)
  fwrite(d[, .(panel, contig, value, share_pct = signif(share, 5), length, length_share_pct = signif(len_share, 5))],
         sub("\\.png$", ".tsv", file), sep = "\t")
}

