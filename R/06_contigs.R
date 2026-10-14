## Reads on the host nuclear genome, excluded host contigs (chrM) and viral genomes, per IP and control library.
## Inputs: work/contigs/bams.tsv and the per-BAM count tables (steps/09_peak_info.sh), work/annot/contigs.tsv
## Tables:  table00c_contig_reads (sample x library x contig), qc07_viral_mapq_summary (MAPQ before the viral filter)
## Figure:  pk05_viral_reads (Figure P5)
## Also writes work/contigs/contig_reads.tsv, which the pie charts (step 19) reuse.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "utils.R"))
suppressPackageStartupMessages(library(patchwork))
cfg <- args_config()
clear_figures(cfg, "^pk05_")
fig_note(cfg, "pk05__summary", "")
vmq <- as.integer((cfg$contigs %||% list())$viral_mapq %||% 10)

B <- rd(file.path(cfg$work, "contigs", "bams.tsv"))
CT <- rd(file.path(cfg$work, "annot", "contigs.tsv"))
if (is.null(B) || !nrow(B) || is.null(CT)) {
  fig_note(cfg, "pk05__summary", "Not available: reads per contig were not counted (samtools missing or no BAM files in the manifest).")
  msg("no per-contig read counts; Figure P5 skipped"); quit(save = "no")
}
S <- samples_dt(cfg)
H <- rbindlist(lapply(unique(B$counts_file), function(f) { x <- rd(f, colClasses = list(character = c("contig", "bam_contig"))); x[, counts_file := f][] }))
H <- merge(H, CT[, .(contig, length, group)], by = "contig", all.x = TRUE)
H[is.na(group), group := "other"]
## MAPQ filter applies to the viral genomes only
H[, pass := group != "viral" | mapq >= vmq]

## ---------------- per library x contig ----------------
per <- H[, .(reads_all = sum(reads), reads = sum(reads[pass])), by = .(counts_file, contig, group, length)]
grid <- CJ(counts_file = unique(B$counts_file), contig = CT$contig)
grid <- merge(grid, CT[, .(contig, group, length)], by = "contig")
per <- rbind(per, grid[!per, on = c("counts_file", "contig")][, `:=`(reads_all = 0, reads = 0)], use.names = TRUE)
L <- merge(B[, .(sample_id, library, counts_file)], per, by = "counts_file", allow.cartesian = TRUE)
L[, total := sum(reads), by = .(sample_id, library)]
nuc <- L[group == "nuclear", .(nuc_reads = sum(reads), nuc_len = sum(as.numeric(length))), by = .(sample_id, library)]
L <- merge(L, nuc, by = c("sample_id", "library"), all.x = TRUE)
L[, `:=`(pct_reads = 100 * reads / pmax(total, 1),
         reads_per_mb = reads / (length / 1e6),
         relative_coverage = (reads / length) / (nuc_reads / nuc_len))]   # 1 = same depth as the host nuclear genome
fwrite(L[, .(sample_id, library, contig, group, length, reads_all, reads, total)], file.path(cfg$work, "contigs", "contig_reads.tsv"), sep = "\t")

## display groups for the figure: host nuclear (all chromosomes together), each excluded contig, each virus, other
L[, unit := fifelse(group == "nuclear", "Host nuclear", fifelse(group == "other", "Other contigs", contig))]
U <- L[, .(reads = sum(reads), reads_all = sum(reads_all), length = sum(as.numeric(length)), total = total[1],
           nuc_reads = nuc_reads[1], nuc_len = nuc_len[1]), by = .(sample_id, library, unit, group)]
U[, `:=`(pct_reads = 100 * reads / pmax(total, 1), relative_coverage = (reads / length) / (nuc_reads / nuc_len))]
ipc <- dcast(U, sample_id + unit ~ library, value.var = "pct_reads")
if (!"control" %in% names(ipc)) ipc[, control := NA_real_]
if (!"IP" %in% names(ipc)) ipc[, IP := NA_real_]
ipc[, log2_ip_vs_control := log2((IP + 1e-6) / (control + 1e-6))]
ipc[is.na(IP) | is.na(control) | (IP == 0 & control == 0), log2_ip_vs_control := NA_real_]

## table00c: one row per sample x library x contig, with the IP/control ratio for contig groups
T0c <- L[, .(sample_id, library, contig, group, length, reads_all, reads, pct_reads = signif(pct_reads, 4),
             reads_per_mb = signif(reads_per_mb, 4), relative_coverage = signif(relative_coverage, 4))]
T0c <- merge(T0c, ipc[, .(sample_id, contig = unit, log2_ip_vs_control = signif(log2_ip_vs_control, 4))], by = c("sample_id", "contig"), all.x = TRUE)
T0c[group == "nuclear", log2_ip_vs_control := NA_real_]
setorder(T0c, sample_id, library, contig)
wt(T0c, cfg, "table00c_contig_reads")

## ---------------- qc07: MAPQ before the viral filter ----------------
q_stats <- function(d) {
  d <- d[, .(n = sum(reads)), by = mapq][order(mapq)]
  tot <- sum(d$n)
  med <- if (tot) d$mapq[which(cumsum(d$n) >= tot / 2)[1]] else NA_real_
  list(reads_before_filter = tot, mapq_0 = d[mapq == 0, sum(n)], mapq_1_9 = d[mapq >= 1 & mapq <= 9, sum(n)],
       mapq_10_29 = d[mapq >= 10 & mapq <= 29, sum(n)], mapq_30_plus = d[mapq >= 30, sum(n)],
       reads_passing = d[mapq >= vmq, sum(n)], pct_passing = if (tot) round(100 * d[mapq >= vmq, sum(n)] / tot, 2) else NA_real_,
       median_mapq = med, mean_mapq = if (tot) round(sum(d$mapq * d$n) / tot, 2) else NA_real_)
}
HB <- merge(B[, .(sample_id, library, counts_file)], H, by = "counts_file", allow.cartesian = TRUE)
HB[, unit := fifelse(group == "nuclear", "Host nuclear (reference)", fifelse(group == "viral", contig, paste0(contig, " (reference)")))]
Q7 <- HB[group != "other", q_stats(.SD), by = .(sample_id, library, contig = unit, group)]
## every viral genome for every library, also when it has no reads at all
## (no viral genomes configured -> nothing to fill in)
vg <- unique(B[, .(sample_id, library)])[, .(contig = CT[group == "viral", contig]), by = .(sample_id, library)]
miss7 <- if (nrow(vg)) vg[!Q7, on = c("sample_id", "library", "contig")] else vg
if (nrow(miss7)) Q7 <- rbind(Q7, miss7[, .(sample_id, library, contig, group = "viral", reads_before_filter = 0, mapq_0 = 0, mapq_1_9 = 0,
                                            mapq_10_29 = 0, mapq_30_plus = 0, reads_passing = 0, pct_passing = NA_real_,
                                            median_mapq = NA_real_, mean_mapq = NA_real_)], use.names = TRUE)
Q7[, viral_mapq_cutoff := vmq]
setcolorder(Q7, c("sample_id", "library", "contig", "group", "viral_mapq_cutoff"))
Q7 <- Q7[order(sample_id, library, match(group, c("viral", "excluded", "nuclear")), contig)]
wt(Q7, cfg, "qc07_viral_mapq_summary")

## ---------------- Figure P5 ----------------
covs <- intersect(unlist(cfg$categorical_covariates), names(S))
tv <- if (!is.null(cfg$time_var) && cfg$time_var %in% names(S)) cfg$time_var else NULL
ord <- do.call(order, c(lapply(c(covs, tv), function(v) S[[v]]), list(S$sample_id)))
s_levels <- intersect(S$sample_id[ord], unique(U$sample_id))
viral_names <- CT[group == "viral", contig]
unit_levels <- c("Host nuclear", CT[group == "excluded", contig], viral_names, "Other contigs")
unit_levels <- intersect(unit_levels, unique(U$unit))
U[, `:=`(sample_id = factor(sample_id, rev(s_levels)), unit = factor(unit, unit_levels),
         library = factor(ifelse(library == "IP", "IP", "Control"), c("IP", "Control")))]
lab_s <- function(x) short_id(x)
fmt_pct <- function(x) ifelse(x == 0, "0", ifelse(x >= 10, sprintf("%.0f", x), ifelse(x >= 0.01, formatC(x, digits = 2, format = "fg"), formatC(x, digits = 1, format = "e"))))
unit_cols <- setNames(c("#7f7f7f", "#b07aa1", "#e15759", "#4e79a7", "#59a14f", "#f28e2b", "#76b7b2", "#edc948")[seq_along(unit_levels)], unit_levels)
unit_cols["Host nuclear"] <- "#7f7f7f"

## A: % of reads (tiles), IP and control side by side; a strip with each sample's virus metadata
gv <- cfg$upset_group_var %||% covs[1]
strip <- NULL
if (!is.null(gv) && gv %in% names(S)) {
  sd <- data.table(sample_id = factor(s_levels, rev(s_levels)), v = S[[gv]][match(s_levels, S$sample_id)])
  sd[is.na(v) | v == "", v := "NA"]
  strip <- ggplot(sd, aes(x = gv, y = sample_id, fill = v)) + geom_tile(colour = "white") +
    scale_fill_brewer(palette = "Set2", name = gv) + scale_y_discrete(labels = lab_s) +
    labs(x = NULL, y = NULL) + theme_x(9) + theme(panel.grid = element_blank(), legend.position = "bottom",
                                                  axis.text.x = element_text(angle = 40, hjust = 1))
}
pA <- ggplot(U, aes(unit, sample_id)) +
  geom_tile(aes(fill = log10(pmax(pct_reads, 1e-4))), colour = "white") +
  geom_text(aes(label = fmt_pct(pct_reads)), size = 2.4) +
  scale_fill_gradient(low = "#f4f1fa", high = "#6a3d9a", limits = c(-4, 2), breaks = c(-4, -2, 0, 2),
                      labels = c("0.0001", "0.01", "1", "100"), name = "% of reads", oob = scales::squish) +
  facet_grid(. ~ library) + scale_y_discrete(labels = if (is.null(strip)) lab_s else NULL) +
  labs(x = NULL, y = NULL, title = "A. Share of reads") + theme_x(9) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "bottom",
        legend.key.width = grid::unit(1, "cm"))
if (!is.null(strip)) pA <- pA + theme(axis.text.y = element_blank())

## a panel saying why it is empty (log axes cannot be drawn without data)
empty_panel <- function(title, why) ggplot() + annotate("text", x = 0, y = 0, label = why, size = 3, colour = "grey40") +
  labs(title = title, x = NULL, y = NULL) + theme_void(9) + theme(plot.title = element_text(face = "bold", size = 9))
no_units <- if (!length(viral_names) && !nrow(CT[group == "excluded"]))
  "No viral genomes or excluded contigs in the genome\n(contigs.viral, contigs.host_exclude, chrom_sizes)" else "No reads on viral genomes\nor excluded contigs"

## B: relative coverage (depth / host nuclear depth), viral genomes and excluded contigs
UB <- U[!unit %in% c("Host nuclear", "Other contigs") & reads > 0]
pB <- if (!nrow(UB)) empty_panel("B. Relative coverage", no_units) else ggplot(UB, aes(relative_coverage, sample_id, colour = unit, shape = library)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") + geom_point(size = 2.2, stroke = 0.9) +
  scale_x_log10(labels = function(b) formatC(b, format = "fg", big.mark = ",")) +
  scale_shape_manual(values = c(IP = 16, Control = 1), name = NULL) + scale_colour_manual(values = unit_cols, name = NULL, drop = TRUE) +
  scale_y_discrete(drop = FALSE, labels = NULL) +
  labs(x = "Depth relative to host nuclear genome (log)", y = NULL, title = "B. Relative coverage") + theme_x(9) +
  theme(legend.position = "bottom")

## C: IP enrichment over control per viral genome / excluded contig
UC <- ipc[unit %in% setdiff(unit_levels, c("Host nuclear", "Other contigs")) & !is.na(log2_ip_vs_control)]
UC[, `:=`(sample_id = factor(sample_id, rev(s_levels)), unit = factor(unit, unit_levels))]
pC <- if (!nrow(UC)) empty_panel("C. IP vs control", no_units) else ggplot(UC, aes(log2_ip_vs_control, sample_id, colour = unit)) +
  geom_vline(xintercept = 0, colour = "grey50") + geom_point(size = 2.2) +
  scale_colour_manual(values = unit_cols, name = NULL, drop = TRUE) + scale_y_discrete(drop = FALSE, labels = NULL) +
  labs(x = "log2(IP share / control share)", y = NULL, title = "C. IP vs control") + theme_x(9) + theme(legend.position = "bottom")

## D: pooled over all samples
P <- U[, .(reads = sum(reads), length = length[1]), by = .(library, unit)]
tot <- unique(U[, .(sample_id, library, total)])[, .(total = sum(total)), by = library]
P <- merge(P, tot, by = "library")
nucP <- P[unit == "Host nuclear", .(library, nuc_reads = reads, nuc_len = length)]
P <- merge(P, nucP, by = "library")
P[, `:=`(pct_reads = 100 * reads / total, relative_coverage = (reads / length) / (nuc_reads / nuc_len))]
## points, not bars: on a log axis a bar has no natural baseline. Zero values are drawn at the bottom as open symbols.
floor_pct <- 10^floor(log10(max(1e-6, min(P[pct_reads > 0, pct_reads], 1)) / 3))
P[, `:=`(y1 = fifelse(pct_reads > 0, pct_reads, floor_pct), zero = pct_reads == 0)]
dodge <- position_dodge(0.5)
pD1 <- ggplot(P, aes(unit, y1, colour = library, shape = zero)) + geom_point(size = 2.8, position = dodge) +
  geom_text(aes(label = fmt_pct(pct_reads)), position = dodge, vjust = -1.1, size = 2.4, show.legend = FALSE) +
  scale_y_log10(labels = function(b) formatC(b, format = "fg"), expand = expansion(mult = c(0.05, 0.15))) +
  scale_colour_manual(values = c(IP = "#6a3d9a", Control = "#8c8c8c"), name = NULL) + scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
  labs(x = NULL, y = "% of reads (log)", title = "D. All samples pooled: share of reads", caption = "Open symbols at the bottom: no reads") + theme_x(9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
P2 <- P[!unit %in% c("Host nuclear", "Other contigs") & reads > 0]
pD2 <- if (!nrow(P2)) empty_panel("All samples pooled: relative coverage", no_units) else ggplot(P2, aes(unit, relative_coverage, colour = library)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") + geom_point(size = 2.8, position = dodge) +
  geom_text(aes(label = formatC(relative_coverage, format = "fg", digits = 2)), position = dodge, vjust = -1.1, size = 2.4, show.legend = FALSE) +
  scale_y_log10(labels = function(b) formatC(b, format = "fg", big.mark = ","), expand = expansion(mult = c(0.08, 0.15))) +
  scale_colour_manual(values = c(IP = "#6a3d9a", Control = "#8c8c8c"), name = NULL) +
  labs(x = NULL, y = "Depth relative to host nuclear (log)", title = "All samples pooled: relative coverage", caption = "Dashed: same depth as the host nuclear genome") +
  theme_x(9) + theme(axis.text.x = element_text(angle = 30, hjust = 1))

top <- if (!is.null(strip)) (strip | pA | pB | pC) + plot_layout(widths = c(0.35, 2.2, 1.3, 1.3)) else (pA | pB | pC) + plot_layout(widths = c(2.2, 1.3, 1.3))
fig <- (top / (pD1 | pD2)) + plot_layout(heights = c(max(2.5, 0.28 * length(s_levels)), 2.2))
save_plot(fig, cfg, "pk05_viral_reads", 17, 4.5 + 0.28 * length(s_levels) + 2.6)
vir_tot <- U[unit %in% viral_names & library == "IP", .(r = sum(reads)), by = unit]
if (!length(viral_names)) msg("WARNING: no viral genomes configured (contigs.viral is empty), so Figure P5 shows host contigs only")
fig_note(cfg, "pk05_viral_reads", sprintf("%sViral reads counted at MAPQ >= %d (contigs.viral_mapq; see qc07_viral_mapq_summary); host contigs are not MAPQ-filtered. %s",
         if (!length(viral_names)) "No viral genomes were configured (contigs.viral), so only host contigs are shown. " else "", vmq, if (nrow(vir_tot)) paste0("Viral reads in IP libraries: ", paste(sprintf("%s %s", vir_tot$unit, format(vir_tot$r, big.mark = ",")), collapse = ", "), ".") else ""))
msg("viral / excluded contig read summary done")
