## Peak overview: number of peaks per sample and per chromosome, peak width and fold enrichment.
## Inputs (step 01): work/peaks/<sample>.info.tsv, work/peaks/peak_stats.tsv, work/annot/genome.sizes
## Tables:  table00_peak_summary, table00b_peaks_by_chromosome
## Figures: pk01_peaks_per_sample, pk02_peaks_by_chromosome, pk03_peak_width, pk04_fold_enrichment
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "utils.R"))
suppressPackageStartupMessages(library(patchwork))
cfg <- args_config()
clear_figures(cfg, "^pk0[0-9]_")

S <- samples_dt(cfg)
G <- rd(file.path(cfg$work, "annot", "genome.sizes"), header = FALSE, col.names = c("chrom", "length"))
if (is.null(G)) stop("work/annot/genome.sizes missing: run step 1 first")
P <- rbindlist(lapply(S$sample_id, function(s) {
  d <- rd(file.path(cfg$work, "peaks", paste0(s, ".info.tsv")))
  if (is.null(d)) return(NULL)
  d[, sample_id := s][]
}), fill = TRUE)
missing <- setdiff(S$sample_id, unique(P$sample_id))
if (length(missing) == nrow(S)) stop("no work/peaks/<sample>.info.tsv files (written by steps/09_peak_info.sh at the start of step 9)")
if (length(missing)) msg("no peak details for: ", paste(missing, collapse = ", "), "; they are shown with 0 peaks")
for (v in c("width", "signal_value", "neglog10_p", "neglog10_q")) set(P, j = v, value = suppressWarnings(as.numeric(P[[v]])))

## sample order: grouped by the categorical covariates (config order), then by name
covs <- intersect(unlist(cfg$categorical_covariates), names(S))
tv <- if (!is.null(cfg$time_var) && cfg$time_var %in% names(S)) cfg$time_var else NULL
ord <- do.call(order, c(lapply(c(covs, tv), function(v) S[[v]]), list(S$sample_id)))
s_levels <- S$sample_id[ord]
gv <- cfg$upset_group_var %||% covs[1] %||% NULL
grp <- if (!is.null(gv) && gv %in% names(S)) setNames(ifelse(is.na(S[[gv]]) | S[[gv]] == "", "NA", S[[gv]]), S$sample_id) else NULL
genome_bp <- sum(G$length)
chrom_levels <- G$chrom

## ---------------- table00: per-sample summary ----------------
q6 <- function(x, pre) {
  x <- x[!is.na(x)]
  v <- if (length(x)) c(min(x), quantile(x, c(0.25, 0.5), names = FALSE), mean(x), quantile(x, 0.75, names = FALSE), max(x)) else rep(NA_real_, 6)
  setNames(as.list(signif(v, 4)), paste0(pre, c("_min", "_q25", "_median", "_mean", "_q75", "_max")))
}
summ <- P[, c(list(n_peaks = .N, bp_covered = sum(as.numeric(width))),
              q6(width, "width"), q6(signal_value, "fold_enrichment"),
              list(n_fold_enrichment_ge_2 = if (all(is.na(signal_value))) NA_integer_ else sum(signal_value >= 2, na.rm = TRUE),
                   n_q_below_0.05 = if (all(is.na(neglog10_q))) NA_integer_ else sum(neglog10_q > -log10(0.05), na.rm = TRUE),
                   n_q_below_0.01 = if (all(is.na(neglog10_q))) NA_integer_ else sum(neglog10_q > 2, na.rm = TRUE))),
          by = sample_id]
top <- P[, .N, by = .(sample_id, chrom)][order(-N)][, .SD[1], by = sample_id]
summ <- merge(summ, top[, .(sample_id, top_chrom = chrom, n_top_chrom = N)], by = "sample_id", all.x = TRUE)
summ <- merge(S[, c("sample_id", covs, tv, intersect("peak_type_resolved", names(S))), with = FALSE], summ, by = "sample_id", all.x = TRUE)
summ[is.na(n_peaks), `:=`(n_peaks = 0L, bp_covered = 0)]
summ[, `:=`(pct_genome_covered = signif(100 * bp_covered / genome_bp, 3),
            pct_top_chrom = signif(100 * n_top_chrom / pmax(n_peaks, 1), 3))]
summ[, n_top_chrom := NULL]
front <- setdiff(names(summ), c("pct_genome_covered", "top_chrom", "pct_top_chrom"))
i <- match("bp_covered", front)
setcolorder(summ, c(front[seq_len(i)], "pct_genome_covered", front[-seq_len(i)], "top_chrom", "pct_top_chrom"))
summ <- summ[match(s_levels, sample_id)]
wt(summ, cfg, "table00_peak_summary")

## ---------------- table00b: peaks per chromosome ----------------
cnt <- P[, .(n_peaks = .N), by = .(sample_id, chrom)]
grid <- CJ(sample_id = S$sample_id, chrom = G$chrom)
pc <- merge(grid, cnt, by = c("sample_id", "chrom"), all.x = TRUE)
pc[is.na(n_peaks), n_peaks := 0L]
pc <- merge(pc, G, by = "chrom")
pc[, tot := sum(n_peaks), by = sample_id]
pc[, `:=`(pct_of_sample_peaks = ifelse(tot > 0, 100 * n_peaks / tot, NA_real_),
          peaks_per_mb = n_peaks / (length / 1e6))]
## density relative to an even spread over the genome: share of peaks / share of genome
pc[, log2_density_ratio := ifelse(n_peaks > 0 & tot > 0, log2((n_peaks / tot) / (length / genome_bp)), NA_real_)]
pc[, tot := NULL]
pc[, `:=`(chrom = factor(chrom, chrom_levels), sample_id = factor(sample_id, s_levels))]
setorder(pc, sample_id, chrom)
out_pc <- copy(pc)[, `:=`(chrom = as.character(chrom), sample_id = as.character(sample_id),
                          pct_of_sample_peaks = signif(pct_of_sample_peaks, 4), peaks_per_mb = signif(peaks_per_mb, 4),
                          log2_density_ratio = signif(log2_density_ratio, 4))]
setcolorder(out_pc, c("sample_id", "chrom", "length", "n_peaks", "pct_of_sample_peaks", "peaks_per_mb", "log2_density_ratio"))
setnames(out_pc, "length", "chrom_length")
wt(out_pc, cfg, "table00b_peaks_by_chromosome")

## ---------------- figures ----------------
ns <- length(s_levels)
lab_s <- function(x) short_id(x)
fill_scale <- if (!is.null(grp)) list(scale_fill_brewer(palette = "Set2", name = gv)) else list()
fill_aes <- if (!is.null(grp)) aes(fill = group) else aes()

## pk01: peaks per sample (count and base pairs covered)
d1 <- summ[, .(sample_id = factor(sample_id, rev(s_levels)), n_peaks, mb = bp_covered / 1e6,
               group = if (!is.null(grp)) grp[sample_id] else NA_character_)]
bar <- function(y, xl, fmt) {
  ggplot(d1, aes(x = .data[[y]], y = sample_id)) +
    (if (!is.null(grp)) geom_col(aes(fill = group), width = 0.75) else geom_col(fill = "#2a78d6", width = 0.75)) +
    geom_text(aes(label = fmt(.data[[y]])), hjust = -0.1, size = 3) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) + scale_y_discrete(labels = lab_s) +
    fill_scale + labs(x = xl, y = NULL) + theme_x(10) + theme(panel.grid.major.y = element_blank())
}
p1 <- bar("n_peaks", "Number of peaks", function(v) format(v, big.mark = ",")) +
  bar("mb", "Genome covered by peaks (Mb)", function(v) sprintf("%.2f", v)) + theme(axis.text.y = element_blank()) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
save_plot(p1, cfg, "pk01_peaks_per_sample", 11, 1.6 + 0.28 * ns)
fig_note(cfg, "pk01_peaks_per_sample", sprintf("%d samples; %s peaks in total (median %s per sample).", ns,
         format(sum(summ$n_peaks), big.mark = ","), format(median(summ$n_peaks), big.mark = ",")))

## pk02: peaks per chromosome (count in tile, colour = density relative to chromosome length)
lim <- max(2, ceiling(max(abs(pc$log2_density_ratio), na.rm = TRUE) * 2) / 2)
lim <- min(lim, 6)
p2 <- ggplot(pc, aes(chrom, sample_id)) +
  geom_tile(aes(fill = pmax(pmin(log2_density_ratio, lim), -lim)), colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0, limits = c(-lim, lim),
                       na.value = "grey85", name = "log2(share of peaks /\nshare of genome)") +
  scale_y_discrete(limits = rev(s_levels), labels = lab_s) +
  labs(x = NULL, y = NULL) + theme_x(9) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.grid = element_blank(),
        legend.position = "right")
if (length(chrom_levels) * ns <= 1500)
  p2 <- p2 + geom_text(aes(label = ifelse(n_peaks > 0, n_peaks, "")), size = if (ns > 30) 1.8 else 2.3)
save_plot(p2, cfg, "pk02_peaks_by_chromosome", max(9, 3 + 0.36 * length(chrom_levels)), 1.6 + 0.3 * ns)
zero <- pc[, .(z = sum(n_peaks == 0)), by = chrom][z == ns, as.character(chrom)]
fig_note(cfg, "pk02_peaks_by_chromosome", sprintf("All %d chromosomes of the genome are shown.%s", length(chrom_levels),
         if (length(zero)) sprintf(" No sample has peaks on: %s.", paste(zero, collapse = ", ")) else ""))

## shared layout for the two distribution figures: distribution panel + a text column with the median and n
dist_plot <- function(d, col, xl, log_base, ref = NULL) {
  d <- d[!is.na(get(col)) & get(col) > 0]
  d[, sample_id := factor(sample_id, rev(s_levels))]
  if (!is.null(grp)) d[, group := grp[as.character(sample_id)]]
  med <- d[, .(m = as.numeric(median(get(col))), n = .N), by = sample_id]
  p <- ggplot(d, aes(x = .data[[col]], y = sample_id)) +
    (if (!is.null(grp)) geom_violin(aes(fill = group), scale = "width", colour = NA, alpha = 0.55)
     else geom_violin(fill = "#9ec5ee", scale = "width", colour = NA)) +
    geom_boxplot(width = 0.22, outlier.shape = NA, fill = "white", linewidth = 0.3) +
    scale_x_continuous(trans = if (log_base == 2) "log2" else "log10",
                       breaks = scales::breaks_log(8, base = log_base),
                       labels = function(b) formatC(b, format = "fg", big.mark = ",")) +
    scale_y_discrete(labels = lab_s, drop = FALSE) + fill_scale +
    labs(x = xl, y = NULL) + theme_x(10) + theme(panel.grid.major.y = element_blank())
  if (!is.null(ref)) p <- p + geom_vline(xintercept = ref, linetype = "dashed", colour = "grey40")
  txt <- ggplot(med, aes(x = 0, y = sample_id, label = sprintf("%s  (n = %s)", formatC(m, format = "fg", digits = 3, big.mark = ","),
                                                                format(n, big.mark = ",")))) +
    geom_text(hjust = 0, size = 3) + scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_discrete(limits = levels(d$sample_id)) + labs(title = "Median  (peaks)") + theme_void(10) +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0))
  p + txt + plot_layout(widths = c(5, 1.1), guides = "collect") & theme(legend.position = "bottom")
}
p3 <- dist_plot(P, "width", "Peak width (bp, log scale)", 10)
save_plot(p3, cfg, "pk03_peak_width", 10, 1.8 + 0.3 * ns)
fig_note(cfg, "pk03_peak_width", sprintf("Median peak width over all samples: %s bp; range of sample medians %s-%s bp.",
         format(median(P$width), big.mark = ","), format(min(summ$width_median, na.rm = TRUE), big.mark = ","),
         format(max(summ$width_median, na.rm = TRUE), big.mark = ",")))

if (any(!is.na(P$signal_value))) {
  p4 <- dist_plot(P, "signal_value", "Fold enrichment over control (log2 scale)", 2, ref = 2)
  save_plot(p4, cfg, "pk04_fold_enrichment", 10, 1.8 + 0.3 * ns)
  nofe <- summ[is.na(fold_enrichment_median), sample_id]
  fig_note(cfg, "pk04_fold_enrichment", paste0("Dashed line: 2-fold enrichment.",
           if (length(nofe)) sprintf(" No fold enrichment in the peak file of: %s.", paste(nofe, collapse = ", ")) else ""))
} else {
  fig_note(cfg, "pk04_fold_enrichment__summary", "Not shown: no peak file has a fold-enrichment column (column 7 of narrowPeak/broadPeak).")
}
msg("peak overview done")
