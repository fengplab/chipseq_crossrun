## Reads on the host nuclear genome, excluded host contigs (chrM) and viral genomes, per IP and control library.
## Inputs: work/contigs/tables.tsv (the pipeline's per-library contig_reads.tsv and its viral BAMs) and work/contigs/bams.tsv
##         (per-BAM counts of the filtered BAMs), both from steps/09_peak_info.sh; work/annot/contigs.tsv
## chrM and viral reads come from the pipeline's contig_reads.tsv when it exists (counted before host filtering, duplicates
## removed); shares are of all reads (total_reads_dedup, unaligned included). Libraries without it fall back to the filtered BAM.
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
TB <- rd(file.path(cfg$work, "contigs", "tables.tsv"), colClasses = "character")
CT <- rd(file.path(cfg$work, "annot", "contigs.tsv"))
if ((is.null(B) || !nrow(B)) && (is.null(TB) || !nrow(TB)) || is.null(CT)) {
  fig_note(cfg, "pk05__summary", "Not available: reads per contig were not counted (samtools missing or no BAM files in the manifest).")
  msg("no per-contig read counts; Figure P5 skipped"); quit(save = "no")
}
if (is.null(B)) B <- data.table(sample_id = character(), library = character(), bam = character(), counts_file = character(), paired = integer())
if (is.null(TB)) TB <- data.table(sample_id = character(), library = character(), contig_table = character(), viral_bam = character(),
                                  viral_counts_file = character())
S <- samples_dt(cfg)
AL <- rd(file.path(cfg$work, "annot", "chrom_alias.tsv"), header = FALSE, col.names = c("from", "to"), colClasses = "character")
canon <- function(x) { i <- match(x, AL$from); ifelse(is.na(i), x, AL$to[i]) }
read_counts <- function(files) {
  files <- unique(files[!is.na(files) & nzchar(files)])
  if (!length(files)) return(data.table(contig = character(), bam_contig = character(), mapq = integer(), reads = numeric(),
                                        fragments = numeric(), counts_file = character()))
  rbindlist(lapply(files, function(f) { x <- rd(f, colClasses = list(character = c("contig", "bam_contig"))); x[, counts_file := f][] }))
}
H <- read_counts(B$counts_file)
H <- merge(H, CT[, .(contig, length, group)], by = "contig", all.x = TRUE)
H[is.na(group), group := "other"]
## MAPQ filter applies to the viral genomes only
H[, pass := group != "viral" | mapq >= vmq]

## ---------------- libraries with the pipeline's contig_reads.tsv ----------------
## The pipeline counts chrM and the viral genomes before host filtering (duplicates removed: reads_dedup), and gives the
## whole-sample totals, so every library splits into parts that add up to total_reads_dedup:
##   host nuclear = mapped_reads - duplicates - sum(reads_dedup of the listed contigs); unaligned = total_reads - mapped_reads.
## Host nuclear reads are split over the chromosomes in proportion to the filtered BAM's reads per chromosome.
need_cols <- c("contig", "contig_type", "length_bp", "total_reads", "duplicates", "total_reads_dedup", "mapped_reads", "reads",
               "reads_dedup", "mean_depth_dedup", "pct_genome_covered")
VC <- read_counts(TB$viral_counts_file)
if (nrow(VC)) VC <- merge(VC, CT[, .(contig, group)], by = "contig", all.x = TRUE)
lib_key <- function(s, l) paste(s, l, sep = "\r")
P <- rbindlist(lapply(seq_len(nrow(TB)), function(i) {
  s <- TB$sample_id[i]; lib <- TB$library[i]
  x <- fread(TB$contig_table[i], sep = "\t", quote = "", colClasses = list(character = intersect(c("sample", "contig", "contig_type", "description"), names(fread(TB$contig_table[i], nrows = 0)))))
  miss <- setdiff(need_cols, names(x))
  if (length(miss)) stop(sprintf("%s: column(s) %s missing", TB$contig_table[i], paste(miss, collapse = ", ")))
  if (!"description" %in% names(x)) x[, description := NA_character_]
  tot <- unique(x[, .(total_reads, duplicates, total_reads_dedup, mapped_reads)])
  if (nrow(tot) != 1) stop(sprintf("%s: whole-sample columns differ between rows", TB$contig_table[i]))
  x <- x[, .(sample_id = s, library = lib, contig = canon(contig), table_contig = contig, group = contig_type, length = as.numeric(length_bp),
             reads_raw = as.numeric(reads), reads_all = as.numeric(reads_dedup), reads = as.numeric(reads_dedup),
             mean_depth = as.numeric(mean_depth_dedup), pct_genome_covered = as.numeric(pct_genome_covered), description)]
  x[!group %in% c("viral", "excluded"), group := fifelse(contig %in% CT[group == "viral", contig], "viral", "excluded")]
  ## optional MAPQ cutoff on the viral genomes, from the pipeline's duplicate-marked viral BAM (primary, non-duplicate reads)
  x[, mapq_source := fifelse(group == "viral", if (vmq > 0) "unfiltered (no viral BAM)" else "unfiltered", "unfiltered")]
  vcf <- TB$viral_counts_file[i]
  if (!is.na(vcf) && nzchar(vcf) && nrow(VC)) {
    v <- VC[counts_file == vcf, .(all = sum(reads), hi = sum(reads[mapq >= vmq])), by = contig]
    for (j in which(x$group == "viral")) {
      k <- match(x$contig[j], v$contig)
      a <- if (is.na(k)) 0 else v$all[k]
      if (abs(a - x$reads_all[j]) > max(5, 0.01 * x$reads_all[j]))
        msg(sprintf("NOTE: %s %s %s: %s reads in the viral BAM vs %s reads_dedup in contig_reads.tsv", s, lib, x$contig[j],
                    format(a, big.mark = ","), format(x$reads_all[j], big.mark = ",")))
      if (vmq > 0) { x$reads[j] <- if (is.na(k)) 0 else v$hi[k]; x$mapq_source[j] <- sprintf("viral BAM, MAPQ >= %d", vmq) }
    }
  }
  host <- tot$mapped_reads - tot$duplicates - sum(x$reads_all)
  unal <- tot$total_reads - tot$mapped_reads
  if (host < 0) stop(sprintf("%s: mapped_reads - duplicates is smaller than the reads on the listed contigs", TB$contig_table[i]))
  ## host nuclear, split over chromosomes like the filtered BAM of the same library
  hb <- B[sample_id == s & library == lib, counts_file][1]
  hc <- if (!is.na(hb) && nrow(H)) H[counts_file == hb & group == "nuclear", .(n = sum(reads)), by = .(contig, length)] else NULL
  nuc_ct <- CT[group == "nuclear"]
  if (!is.null(hc) && nrow(hc) && sum(hc$n) > 0) {
    hc <- merge(nuc_ct[, .(contig, length)], hc[, .(contig, n)], by = "contig", all.x = TRUE)[is.na(n), n := 0]
    hn <- hc[, .(sample_id = s, library = lib, contig, table_contig = NA_character_, group = "nuclear", length = as.numeric(length),
                 reads_raw = NA_real_, reads_all = host * n / sum(n), reads = host * n / sum(n), mean_depth = NA_real_,
                 pct_genome_covered = NA_real_, description = NA_character_, mapq_source = "split by filtered BAM")]
  } else {
    hn <- data.table(sample_id = s, library = lib, contig = "Host nuclear", table_contig = NA_character_, group = "nuclear",
                     length = sum(as.numeric(nuc_ct$length)), reads_raw = NA_real_, reads_all = host, reads = host, mean_depth = NA_real_,
                     pct_genome_covered = NA_real_, description = NA_character_, mapq_source = "no filtered BAM to split by")
  }
  un <- data.table(sample_id = s, library = lib, contig = "Unaligned", table_contig = NA_character_, group = "unaligned", length = NA_real_,
                   reads_raw = NA_real_, reads_all = unal, reads = unal, mean_depth = NA_real_, pct_genome_covered = NA_real_,
                   description = "not aligned, or dropped by Allo", mapq_source = NA_character_)
  out <- rbind(x, hn, un, use.names = TRUE)
  out[, `:=`(total = tot$total_reads_dedup, source = "pipeline contig_reads.tsv")][]
}), use.names = TRUE)

## ---------------- libraries without it: counts from the filtered BAMs ----------------
## (chrM and viral reads are zero there when the pipeline removed them before writing those BAMs)
per <- H[, .(reads_all = sum(reads), reads = sum(reads[pass])), by = .(counts_file, contig, group, length)]
grid <- CJ(counts_file = unique(B$counts_file), contig = CT$contig)
grid <- merge(grid, CT[, .(contig, group, length)], by = "contig")
per <- rbind(per, grid[!per, on = c("counts_file", "contig")][, `:=`(reads_all = 0, reads = 0)], use.names = TRUE)
LB <- merge(B[, .(sample_id, library, counts_file)], per, by = "counts_file", allow.cartesian = TRUE)
if (nrow(P)) LB <- LB[!lib_key(sample_id, library) %in% lib_key(P$sample_id, P$library)]
LB[, `:=`(total = sum(reads), source = "filtered BAM", table_contig = NA_character_, reads_raw = NA_real_, mean_depth = NA_real_,
          pct_genome_covered = NA_real_, description = NA_character_,
          mapq_source = fifelse(group == "viral", sprintf("filtered BAM, MAPQ >= %d", vmq), "filtered BAM")), by = .(sample_id, library)]
LB[, counts_file := NULL]
L <- rbind(P, LB, use.names = TRUE, fill = TRUE)
if (nrow(P) && nrow(LB)) msg(sprintf("WARNING: %d of %d libraries have no contig_reads.tsv; their chrM/viral reads come from the filtered BAMs",
                                     uniqueN(LB[, .(sample_id, library)]), uniqueN(L[, .(sample_id, library)])))
if (nrow(P)) msg(sprintf("chrM and viral reads from the pipeline's contig_reads.tsv for %d libraries", uniqueN(P[, .(sample_id, library)])))

nuc <- L[group == "nuclear", .(nuc_reads = sum(reads), nuc_len = sum(as.numeric(length))), by = .(sample_id, library)]
L <- merge(L, nuc, by = c("sample_id", "library"), all.x = TRUE)
L[, `:=`(pct_reads = 100 * reads / pmax(total, 1),
         reads_per_mb = reads / (length / 1e6),
         relative_coverage = (reads / length) / (nuc_reads / nuc_len))]   # 1 = same depth as the host nuclear genome
fwrite(L[, .(sample_id, library, contig, group, length, reads_all, reads, total, source)], file.path(cfg$work, "contigs", "contig_reads.tsv"), sep = "\t")

## display groups for the figure: host nuclear (all chromosomes together), each excluded contig, each virus, other, unaligned
L[, unit := fifelse(group == "nuclear", "Host nuclear", fifelse(group == "other", "Other contigs", fifelse(group == "unaligned", "Unaligned", contig)))]
U <- L[, .(reads = sum(reads), reads_all = sum(reads_all), length = sum(as.numeric(length)), total = total[1],
           nuc_reads = nuc_reads[1], nuc_len = nuc_len[1], pct_genome_covered = pct_genome_covered[1]), by = .(sample_id, library, unit, group)]
U[, `:=`(pct_reads = 100 * reads / pmax(total, 1), relative_coverage = (reads / length) / (nuc_reads / nuc_len))]
ipc <- dcast(U, sample_id + unit ~ library, value.var = "pct_reads")
if (!"control" %in% names(ipc)) ipc[, control := NA_real_]
if (!"IP" %in% names(ipc)) ipc[, IP := NA_real_]
ipc[, log2_ip_vs_control := log2((IP + 1e-6) / (control + 1e-6))]
ipc[is.na(IP) | is.na(control) | (IP == 0 & control == 0), log2_ip_vs_control := NA_real_]

## table00c: one row per sample x library x contig, with the IP/control ratio for contig groups
T0c <- L[, .(sample_id, library, contig, group, length, reads_all, reads = round(reads, 1), pct_reads = signif(pct_reads, 4),
             reads_per_mb = signif(reads_per_mb, 4), relative_coverage = signif(relative_coverage, 4), mean_depth, pct_genome_covered,
             count_source = source, viral_mapq = mapq_source, description)]
T0c[, reads_all := round(reads_all, 1)]
T0c <- merge(T0c, ipc[, .(sample_id, contig = unit, log2_ip_vs_control = signif(log2_ip_vs_control, 4))], by = c("sample_id", "contig"), all.x = TRUE)
T0c[group %in% c("nuclear", "unaligned"), log2_ip_vs_control := NA_real_]
setorder(T0c, sample_id, library, contig)
wt(T0c, cfg, "table00c_contig_reads")

## ---------------- qc07: MAPQ of the viral reads before the viral filter ----------------
q_stats <- function(d) {
  d <- d[, .(n = sum(reads)), by = mapq][order(mapq)]
  tot <- sum(d$n)
  med <- if (tot) d$mapq[which(cumsum(d$n) >= tot / 2)[1]] else NA_real_
  list(reads_before_filter = tot, mapq_0 = d[mapq == 0, sum(n)], mapq_1_9 = d[mapq >= 1 & mapq <= 9, sum(n)],
       mapq_10_29 = d[mapq >= 10 & mapq <= 29, sum(n)], mapq_30_plus = d[mapq >= 30, sum(n)],
       reads_passing = d[mapq >= vmq, sum(n)], pct_passing = if (tot) round(100 * d[mapq >= vmq, sum(n)] / tot, 2) else NA_real_,
       median_mapq = med, mean_mapq = if (tot) round(sum(d$mapq * d$n) / tot, 2) else NA_real_)
}
## libraries with a pipeline table: MAPQ from the duplicate-marked viral BAM (viral genomes only)
QP <- if (nrow(VC)) {
  vb <- merge(TB[!is.na(viral_counts_file) & nzchar(viral_counts_file), .(sample_id, library, counts_file = viral_counts_file)],
              VC[group %in% "viral"], by = "counts_file", allow.cartesian = TRUE)
  if (nrow(vb)) vb[, q_stats(.SD), by = .(sample_id, library, contig, group)][, source := "pipeline viral BAM"][] else NULL
} else NULL
## libraries without one: the filtered BAM (also gives chrM and the host reference rows)
QB <- if (nrow(H) && nrow(LB)) {
  HB <- merge(unique(LB[, .(sample_id, library)])[B, on = c("sample_id", "library"), nomatch = 0][, .(sample_id, library, counts_file)],
              H, by = "counts_file", allow.cartesian = TRUE)
  HB[, unit := fifelse(group == "nuclear", "Host nuclear (reference)", fifelse(group == "viral", contig, paste0(contig, " (reference)")))]
  if (nrow(HB)) HB[group != "other", q_stats(.SD), by = .(sample_id, library, contig = unit, group)][, source := "filtered BAM"][] else NULL
} else NULL
Q7 <- rbind(QP, QB, use.names = TRUE, fill = TRUE)
if (is.null(Q7) || !nrow(Q7)) Q7 <- data.table(sample_id = character(), library = character(), contig = character(), group = character(),
                                              reads_before_filter = numeric(), source = character())
## every viral genome for every library, also when it has no reads at all (or no viral BAM: NA, see source)
libs <- unique(rbind(B[, .(sample_id, library)], TB[, .(sample_id, library)]))
vg <- libs[, .(contig = CT[group == "viral", contig]), by = .(sample_id, library)]
miss7 <- if (nrow(vg)) vg[!Q7, on = c("sample_id", "library", "contig")] else vg
if (nrow(miss7)) {
  has_vb <- lib_key(miss7$sample_id, miss7$library) %in% lib_key(TB$sample_id, TB$library)[!is.na(TB$viral_counts_file) & nzchar(TB$viral_counts_file)]
  in_tb <- lib_key(miss7$sample_id, miss7$library) %in% lib_key(TB$sample_id, TB$library)
  z <- fifelse(in_tb & !has_vb, NA_real_, 0)
  Q7 <- rbind(Q7, miss7[, .(sample_id, library, contig, group = "viral", reads_before_filter = z, mapq_0 = z, mapq_1_9 = z,
                            mapq_10_29 = z, mapq_30_plus = z, reads_passing = z, pct_passing = NA_real_,
                            median_mapq = NA_real_, mean_mapq = NA_real_,
                            source = fifelse(in_tb & !has_vb, "no viral BAM", fifelse(in_tb, "pipeline viral BAM", "filtered BAM")))],
              use.names = TRUE, fill = TRUE)
}
Q7[, viral_mapq_cutoff := vmq]
setcolorder(Q7, intersect(c("sample_id", "library", "contig", "group", "viral_mapq_cutoff"), names(Q7)))
Q7 <- Q7[order(sample_id, library, match(group, c("viral", "excluded", "nuclear")), contig)]
wt(Q7, cfg, "qc07_viral_mapq_summary")

## ---------------- Figure P5 ----------------
covs <- intersect(unlist(cfg$categorical_covariates), names(S))
tv <- if (!is.null(cfg$time_var) && cfg$time_var %in% names(S)) cfg$time_var else NULL
ord <- do.call(order, c(lapply(c(covs, tv), function(v) S[[v]]), list(S$sample_id)))
s_levels <- intersect(S$sample_id[ord], unique(U$sample_id))
viral_names <- CT[group == "viral", contig]
unit_levels <- c("Host nuclear", CT[group == "excluded", contig], viral_names,
                 setdiff(unique(L[group %in% c("excluded", "viral"), contig]), c(CT$contig)), "Other contigs", "Unaligned")
unit_levels <- intersect(unit_levels, unique(U$unit))
U[, `:=`(sample_id = factor(sample_id, rev(s_levels)), unit = factor(unit, unit_levels),
         library = factor(ifelse(library == "IP", "IP", "Control"), c("IP", "Control")))]
lab_s <- function(x) short_id(x)
fmt_pct <- function(x) ifelse(x == 0, "0", ifelse(x >= 10, sprintf("%.0f", x), ifelse(x >= 0.01, formatC(x, digits = 2, format = "fg"), formatC(x, digits = 1, format = "e"))))
unit_cols <- setNames(c("#7f7f7f", "#b07aa1", "#e15759", "#4e79a7", "#59a14f", "#f28e2b", "#76b7b2", "#edc948")[seq_along(unit_levels)], unit_levels)
unit_cols["Host nuclear"] <- "#7f7f7f"
if ("Unaligned" %in% names(unit_cols)) unit_cols["Unaligned"] <- "#c7c7c7"
if (anyNA(unit_cols)) unit_cols[is.na(unit_cols)] <- grDevices::hcl.colors(sum(is.na(unit_cols)), "Dark 3")

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
NONCTG <- c("Host nuclear", "Other contigs", "Unaligned")
UB <- U[!unit %in% NONCTG & reads > 0]
pB <- if (!nrow(UB)) empty_panel("B. Relative coverage", no_units) else ggplot(UB, aes(relative_coverage, sample_id, colour = unit, shape = library)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") + geom_point(size = 2.2, stroke = 0.9) +
  scale_x_log10(labels = function(b) formatC(b, format = "fg", big.mark = ",")) +
  scale_shape_manual(values = c(IP = 16, Control = 1), name = NULL) + scale_colour_manual(values = unit_cols, name = NULL, drop = TRUE) +
  scale_y_discrete(drop = FALSE, labels = NULL) +
  labs(x = "Depth relative to host nuclear genome (log)", y = NULL, title = "B. Relative coverage") + theme_x(9) +
  theme(legend.position = "bottom")

## C: IP enrichment over control per viral genome / excluded contig
UC <- ipc[unit %in% setdiff(unit_levels, NONCTG) & !is.na(log2_ip_vs_control)]
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
P2 <- P[!unit %in% NONCTG & reads > 0]
pD2 <- if (!nrow(P2)) empty_panel("All samples pooled: relative coverage", no_units) else ggplot(P2, aes(unit, relative_coverage, colour = library)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") + geom_point(size = 2.8, position = dodge) +
  geom_text(aes(label = formatC(relative_coverage, format = "fg", digits = 2)), position = dodge, vjust = -1.1, size = 2.4, show.legend = FALSE) +
  scale_y_log10(labels = function(b) formatC(b, format = "fg", big.mark = ","), expand = expansion(mult = c(0.08, 0.15))) +
  scale_colour_manual(values = c(IP = "#6a3d9a", Control = "#8c8c8c"), name = NULL) +
  labs(x = NULL, y = "Depth relative to host nuclear (log)", title = "All samples pooled: relative coverage", caption = "Dashed: same depth as the host nuclear genome") +
  theme_x(9) + theme(axis.text.x = element_text(angle = 30, hjust = 1))

## E: share of each viral genome / excluded contig covered by at least one read (pipeline contig_reads.tsv only)
UE <- U[!unit %in% NONCTG & !is.na(pct_genome_covered)]
pE <- if (!nrow(UE)) NULL else ggplot(UE, aes(pct_genome_covered, sample_id, colour = unit, shape = library)) +
  geom_point(size = 2.2, stroke = 0.9) + scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  scale_shape_manual(values = c(IP = 16, Control = 1), name = NULL) + scale_colour_manual(values = unit_cols, name = NULL, drop = TRUE) +
  scale_y_discrete(drop = FALSE, labels = lab_s) +
  labs(x = "% of genome covered by >= 1 read", y = NULL, title = "E. Genome covered") + theme_x(9) + theme(legend.position = "bottom")

top <- if (!is.null(strip)) (strip | pA | pB | pC) + plot_layout(widths = c(0.35, 2.2, 1.3, 1.3)) else (pA | pB | pC) + plot_layout(widths = c(2.2, 1.3, 1.3))
bottom <- if (is.null(pE)) (pD1 | pD2) else (pD1 | pD2 | pE) + plot_layout(widths = c(1.4, 1.1, 1.1))
fig <- (top / bottom) + plot_layout(heights = c(max(2.5, 0.28 * length(s_levels)), max(2.2, if (is.null(pE)) 0 else 0.22 * length(s_levels))))
save_plot(fig, cfg, "pk05_viral_reads", 17, 4.5 + 0.28 * length(s_levels) + 2.6)
vir_tot <- U[group == "viral" & library == "IP", .(r = round(sum(reads))), by = unit]
src_note <- if (nrow(P)) sprintf("chrM and viral reads are the pipeline's counts before host filtering (contig_reads.tsv, duplicates removed); shares are of all reads, unaligned included (total_reads_dedup); host nuclear = mapped - duplicates - listed contigs. %s",
                                 if (nrow(LB)) "Libraries without contig_reads.tsv use the filtered BAM. " else "") else ""
if (!length(viral_names)) msg("WARNING: no viral genomes configured (contigs.viral is empty), so Figure P5 shows host contigs only")
fig_note(cfg, "pk05_viral_reads", sprintf("%s%sViral reads %s (contigs.viral_mapq; see qc07_viral_mapq_summary); host contigs are not MAPQ-filtered. %s",
         if (!length(viral_names)) "No viral genomes were configured (contigs.viral), so only host contigs are shown. " else "", src_note,
         if (vmq > 0) sprintf("counted at MAPQ >= %d", vmq) else "not MAPQ-filtered", if (nrow(vir_tot)) paste0("Viral reads in IP libraries: ", paste(sprintf("%s %s", vir_tot$unit, format(vir_tot$r, big.mark = ",")), collapse = ", "), ".") else ""))
msg("viral / excluded contig read summary done")
