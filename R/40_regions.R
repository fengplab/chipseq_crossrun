#!/usr/bin/env Rscript
## Q4: Which regions recur across samples, and where do samples overlap most?
## Tables:  table05_recurrent_regions (top N), table05_all_regions, table05b_jaccard_matrix, table05c_hotspot_bins
## Figures: fig10_region_recurrence_hist, fig11_jaccard_heatmap, fig12_upset_<group>, fig13_hotspot_genome
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "utils.R"))
suppressPackageStartupMessages(library(pheatmap))
cfg <- args_config()
S <- samples_dt(cfg)
D <- file.path(cfg$work, "regions")
mr <- rd(file.path(D, "merged_regions.tsv"))
if (is.null(mr)) stop("run step 08 first")
top_n <- cfg$regions$top_n %||% 50

## per-region metadata breadth
long <- mr[, .(sample_id = strsplit(samples, ",", fixed = TRUE)[[1]]), by = region_id]
cv <- intersect(unlist(cfg$condition_vars), names(S))
long <- merge(long, S[, c("sample_id", cv), with = FALSE], by = "sample_id", all.x = TRUE)
if (length(cv)) long[, condition := do.call(paste, c(.SD, sep = "|")), .SDcols = cv] else long[, condition := sample_id]
br <- long[, c(list(n_conditions = uniqueN(condition)), lapply(.SD, uniqueN)), by = region_id, .SDcols = cv]
setnames(br, cv, paste0("n_", cv))
reg <- merge(mr, br, by = "region_id")
reg[, width := end - start]

rl <- rd(file.path(D, "region_labels.tsv"))
if (!is.null(rl) && nrow(rl)) {
  lab <- rl[, .(labels = paste(unique(label), collapse = ",")), by = .(region_id, feature_set)]
  reg <- merge(reg, dcast(lab, region_id ~ feature_set, value.var = "labels"), by = "region_id", all.x = TRUE)
}
ng <- rd(file.path(D, "region_nearest_gene.tsv"))
if (!is.null(ng)) reg <- merge(reg, ng[, .(region_id, nearest_gene = gene, distance_to_gene = distance)], by = "region_id", all.x = TRUE)
io <- rd(file.path(D, "region_input_overlap.tsv"))
reg <- merge(reg, if (!is.null(io) && nrow(io)) io else data.table(region_id = character(), n_input_peak_sets = integer()),
             by = "region_id", all.x = TRUE)
reg[is.na(n_input_peak_sets), n_input_peak_sets := 0L]
reg[, likely_artifact := n_input_peak_sets > 0]
setorder(reg, -n_samples, -n_conditions, chrom, start)
wt(reg, cfg, "table05_all_regions")
wt(head(reg[, !"samples"], top_n), cfg, "table05_recurrent_regions")

n_s <- uniqueN(long$sample_id)
h <- reg[, .N, by = n_samples]
p10 <- ggplot(h, aes(n_samples, N)) + geom_col(fill = "#2a78d6", width = 0.8) + scale_y_log10() +
  scale_x_continuous(breaks = seq_len(n_s)) +
  labs(x = "Number of samples sharing the region", y = "Merged regions (log scale)", title = "Region recurrence") + theme_x()
save_plot(p10, cfg, "fig10_region_recurrence_hist", max(7, 3 + 0.2 * n_s), 4.5)
fig_note(cfg, "fig10_region_recurrence_hist", sprintf("%d merged regions across %d samples.", nrow(reg), n_s))

## Jaccard heatmap
jc <- rd(file.path(D, "jaccard.tsv"))
if (!is.null(jc) && nrow(jc)) {
  jj <- rbind(jc, jc[sample_a != sample_b, .(sample_a = sample_b, sample_b = sample_a, jaccard)])
  w <- dcast(jj, sample_a ~ sample_b, value.var = "jaccard")
  wt(w, cfg, "table05b_jaccard_matrix")
  m <- as.matrix(w[, -1]); rownames(m) <- w$sample_a
  if (nrow(m) > 1) {
    diag(m) <- NA
    save_pheatmap(function() pheatmap(m, annotation_col = sample_annotation(cfg, S, colnames(m)),
                                      clustering_distance_rows = as.dist(1 - replace(m, is.na(m), 1)),
                                      clustering_distance_cols = as.dist(1 - replace(m, is.na(m), 1)),
                                      color = colorRampPalette(c("white", "#2a78d6", "#0c447c"))(50), fontsize = 8,
                                      na_col = "grey90", main = "Pairwise peak Jaccard similarity"),
                  cfg, "fig11_jaccard_heatmap", width = max(9, 5 + 0.3 * ncol(m)), height = max(7, 4 + 0.3 * nrow(m)))
    fig_note(cfg, "fig11_jaccard_heatmap", sprintf("%d samples; diagonal (self-comparison) left blank.", nrow(m)))
  }
}

## UpSet by group
gv <- cfg$upset_group_var
if (!is.null(gv) && gv %in% names(long) && requireNamespace("UpSetR", quietly = TRUE)) {
  gl <- unique(long[!is.na(get(gv)), .(region_id, grp = as.character(get(gv)))])
  sets <- split(gl$region_id, gl$grp)
  if (length(sets) >= 2) {
    f <- function() print(UpSetR::upset(UpSetR::fromList(sets), nsets = length(sets), nintersects = 30,
                                        order.by = "freq", mainbar.y.label = "Shared regions",
                                        sets.x.label = paste("Regions per", gv), main.bar.color = "#2a78d6"))
    save_pheatmap(f, cfg, paste0("fig12_upset_", gv), width = 11, height = 6.5)
    fig_note(cfg, paste0("fig12_upset_", gv), sprintf("Samples grouped by '%s' (%d groups); up to 30 largest intersections shown.", gv, length(sets)))
  }
}

## Hotspot bins
bn <- rd(file.path(D, "bins.tsv"))
if (!is.null(bn) && nrow(bn)) {
  hb <- bn[, .(n_samples = uniqueN(sample_id), total_peaks = sum(n_peaks)), by = .(chrom, start, end)]
  setorder(hb, -n_samples, -total_peaks)
  wt(head(hb, 200), cfg, "table05c_hotspot_bins")
  ch <- unique(hb$chrom)
  ch <- ch[order(suppressWarnings(as.numeric(sub("^chr", "", ch))), ch)]
  ch <- ch[!grepl("_|^chrUn|random|alt|EBV|^chrM$", ch)]
  hb <- hb[chrom %in% ch]
  hb[, chrom := factor(chrom, levels = ch)]
  thr <- max(2, ceiling(0.5 * n_s))
  p13 <- ggplot(hb, aes(start / 1e6, n_samples)) +
    geom_point(aes(colour = n_samples >= thr), size = 0.7, alpha = 0.7) +
    facet_wrap(~chrom, scales = "free_x", ncol = 6) +
    scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "#e34948"), guide = "none") +
    labs(x = "Position (Mb)", y = "Samples with a peak in bin",
         title = sprintf("Recurrent %s-kb bins (red: >= %d samples)", cfg$regions$bin_size / 1000, thr)) +
    theme_x(9) + theme(panel.spacing = unit(0.4, "lines"))
  save_plot(p13, cfg, "fig13_hotspot_genome", 14, 1.5 + 2 * ceiling(length(ch) / 6))
  fig_note(cfg, "fig13_hotspot_genome", sprintf("Bins with at least one peak; red = bins with peaks in >= %d of %d samples.", thr, n_s))
}
