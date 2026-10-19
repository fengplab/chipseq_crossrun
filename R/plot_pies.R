## Pie charts of reads and peaks per contig (step 19; written to <outdir>/extra/pies/, not part of the report).
##   per sample:  coverage_<sample>_{with_viruses,host}_{ring,noring}.png   reads per contig, IP and control side by side
##                peaks_<sample>_{ring,noring}.png                          peaks per host chromosome
##   pooled:      pooled_coverage_rpkm_{with_viruses,host}.png               all samples' reads pooled, as RPKM shares (no ring)
##                pooled_peaks_{ring,noring}.png                             all samples' peaks pooled
## "with_viruses" = host chromosomes + excluded host contigs (chrM) + viral genomes; "host" = host chromosomes + chrM.
## Counts come from work/contigs/contig_reads.tsv (step 9). With the pipeline's contig_reads.tsv, chrM and viral reads are its
## counts before host filtering (duplicates removed) and the host nuclear total is split over the chromosomes like the filtered
## BAM; otherwise everything is counted from the filtered BAM. Viral reads at MAPQ >= contigs.viral_mapq when that is > 0.
## Inner ring: each contig's share of the total length of the contigs in the pie (what an even spread would give).
## Every figure has a TSV of its numbers next to it.
## --peaks-only: draw only the peak pies, plus peaks_per_chromosome.tsv (needs steps 0-1 only).
## --group-by var1,var2 (or extra_plots.pie_group_by in the config): also pool samples by each level of a manifest column
##   (e.g. virus, cell_type), one figure per column with one pie per level:
##                peaks_by_<var>_{ring,noring}.png                          peaks per host chromosome, summed over the level's samples
##                coverage_rpkm_by_<var>_{with_viruses,host}.png             IP reads pooled per level, as RPKM shares (no ring)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "utils.R"))
cfg <- args_config()
ca <- commandArgs(TRUE)
peaks_only <- "--peaks-only" %in% ca
gi <- match("--group-by", ca)
if (!is.na(gi) && (gi == length(ca) || startsWith(ca[gi + 1], "--"))) stop("--group-by needs a comma-separated list of manifest columns")
group_by <- unique(c(unlist((cfg$extra_plots %||% list())$pie_group_by), if (!is.na(gi)) trimws(strsplit(ca[gi + 1], ",")[[1]])))
group_by <- group_by[nzchar(group_by)]
out <- file.path(cfg$outdir, "extra", "pies"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
CR <- rd(file.path(cfg$work, "contigs", "contig_reads.tsv"))
CT <- rd(file.path(cfg$work, "annot", "contigs.tsv"))
G <- rd(file.path(cfg$work, "annot", "genome.sizes"), header = FALSE, col.names = c("contig", "length"))
if (is.null(CT) || is.null(G)) stop("work/annot/contigs.tsv missing: run step 1 first")
S <- samples_dt(cfg)
bad_g <- setdiff(group_by, names(S))
if (length(bad_g)) stop("--group-by / extra_plots.pie_group_by: not a manifest column: ", paste(bad_g, collapse = ", "),
                        ". Columns: ", paste(setdiff(names(S), "sample_id"), collapse = ", "))
## level of each sample for a grouping column; blank or missing -> "NA" (logged)
level_of <- function(v) {
  x <- as.character(S[[v]]); x[is.na(x) | !nzchar(trimws(x))] <- "NA"
  if (any(x == "NA")) msg(sprintf("pie grouping by %s: %d sample(s) without a value, pooled as 'NA': %s", v, sum(x == "NA"),
                                  paste(S$sample_id[x == "NA"], collapse = ", ")))
  setNames(x, S$sample_id)
}
safe <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
group_panel <- function(lv, ids) {                 # "<level> (n = number of samples)", ordered by level; ids may repeat
  n <- table(lv[unique(ids)]); lab <- setNames(sprintf("%s (n = %d)", names(n), as.integer(n)), names(n))
  factor(lab[lv[ids]], lab[order(names(n))])
}

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

vmq <- as.integer((cfg$contigs %||% list())$viral_mapq %||% 10)
sets <- list(with_viruses = order_all, host = c(nuc, excluded))
n_files <- 0

## ---------------- coverage pies per sample ----------------
pipe_src <- !is.null(CR) && "source" %in% names(CR) && any(CR$source == "pipeline contig_reads.tsv")
src_txt <- if (pipe_src) "chrM and viral reads: pipeline counts before host filtering; chromosomes: host total split as in the filtered BAM" else NULL
mq_txt <- if (vmq > 0) sprintf("viral reads at MAPQ >= %d", vmq) else "viral reads not MAPQ-filtered"
if (peaks_only) {
  msg("--peaks-only: coverage pies skipped")
} else if (!is.null(CR) && nrow(CR)) {
  CR <- CR[contig %in% order_all]
  CR[, panel := ifelse(library == "IP", "IP", "Control")]
  for (s in intersect(S$sample_id, unique(CR$sample_id))) {
    d <- CR[sample_id == s, .(panel, contig, value = reads, length)]
    d[, panel := factor(panel, intersect(c("IP", "Control"), unique(panel)))]
    for (nm in names(sets)) for (ring in c(TRUE, FALSE)) {
      cs <- intersect(sets[[nm]], d$contig)
      f <- file.path(out, sprintf("coverage_%s_%s_%s.png", s, nm, if (ring) "ring" else "noring"))
      draw_pie(d, cs, sprintf("Reads per contig: %s", s),
               paste(c(sprintf("%s; %s", if (nm == "host") "host chromosomes and chrM" else "host chromosomes, chrM and viral genomes", mq_txt), src_txt), collapse = "\n"),
               f, ring)
      n_files <- n_files + 1
    }
  }
  ## pooled: reads summed over all samples, as RPKM (reads per kb of contig per million pooled reads); no ring
  pooled <- CR[, .(reads = sum(reads)), by = .(panel, contig, length)]
  for (nm in names(sets)) {
    pd <- pooled[contig %in% sets[[nm]]]
    pd[, rpkm := reads / (length / 1e3) / (sum(reads) / 1e6), by = panel]
    pd[, panel := factor(panel, intersect(c("IP", "Control"), unique(panel)))]
    extra <- pd[panel == "IP", setNames(sprintf("  (IP RPKM %s)", formatC(rpkm, format = "fg", digits = 3, big.mark = ",")), contig)]
    f <- file.path(out, sprintf("pooled_coverage_rpkm_%s.png", nm))
    draw_pie(pd[, .(panel, contig, value = rpkm, length)], intersect(sets[[nm]], pd$contig),
             sprintf("All %d samples pooled: coverage per contig (RPKM)", uniqueN(CR$sample_id)),
             paste(c(sprintf("slice = share of RPKM (pooled reads per kb of contig per million pooled reads), i.e. relative depth; %s; %s",
                     if (nm == "host") "host chromosomes and chrM" else "with viral genomes", mq_txt), src_txt), collapse = "\n"),
             f, ring = FALSE, legend_extra = extra)
    fwrite(pd[, .(panel, contig, length, reads = round(reads, 1), rpkm = signif(rpkm, 5))], sub("\\.png$", "_rpkm.tsv", f), sep = "\t")
    n_files <- n_files + 1
  }
  ## grouped: IP reads pooled over each level's samples, as RPKM within the level (controls are shared between samples)
  for (v in group_by) {
    lv <- level_of(v)
    ip <- CR[panel == "IP" & sample_id %in% names(lv)]
    ip[, grp := group_panel(lv, sample_id)]
    for (nm in names(sets)) {
      pd <- ip[contig %in% sets[[nm]], .(reads = sum(reads)), by = .(panel = grp, contig, length)]
      pd[, rpkm := reads / (length / 1e3) / (sum(reads) / 1e6), by = panel]
      f <- file.path(out, sprintf("coverage_rpkm_by_%s_%s.png", safe(v), nm))
      draw_pie(pd[, .(panel, contig, value = rpkm, length)], intersect(sets[[nm]], pd$contig),
               sprintf("IP coverage per contig (RPKM), samples pooled by %s", v),
               paste(c(sprintf("slice = share of RPKM within each group (IP reads pooled over the group's samples); %s; %s",
                       if (nm == "host") "host chromosomes and chrM" else "with viral genomes", mq_txt), src_txt), collapse = "\n"),
               f, ring = FALSE)
      fwrite(pd[, .(group = panel, contig, length, reads = round(reads, 1), rpkm = signif(rpkm, 5))], sub("\\.png$", "_rpkm.tsv", f), sep = "\t")
      n_files <- n_files + 1
    }
  }
} else msg("no per-contig read counts (work/contigs/contig_reads.tsv, step 9): coverage pies skipped")

## ---------------- peak pies (host chromosomes; the analyses have no peaks on chrM or viral genomes) ----------------
P <- rbindlist(lapply(S$sample_id, function(s) {
  b <- rd(file.path(cfg$work, "peaks", paste0(s, ".bed5")), header = FALSE)
  if (is.null(b) || !nrow(b)) return(NULL)
  b[, .(n = .N), by = .(contig = V1)][, sample_id := s][]
}))
if (!is.null(P) && nrow(P)) {
  hosts <- intersect(order_all, G$contig)
  ## every peak counted must be on a host nuclear chromosome (step 1 drops chrM and viral peaks and logs them),
  ## and the per-chromosome counts must add up to step 1's per-sample totals
  off <- P[!contig %in% hosts]
  if (nrow(off)) stop("peaks on contigs outside the host nuclear chromosomes: ", paste(unique(off$contig), collapse = ", "))
  PS <- rd(file.path(cfg$work, "peaks", "peak_stats.tsv"))
  if (!is.null(PS)) {
    chk <- merge(P[, .(n = sum(n)), by = sample_id], PS[, .(sample_id, n_peaks)], by = "sample_id")
    if (chk[n != n_peaks, .N]) stop("peak counts per chromosome do not add up to work/peaks/peak_stats.tsv for: ", paste(chk[n != n_peaks, sample_id], collapse = ", "))
  }
  no_peaks <- setdiff(S$sample_id, P$sample_id)
  if (length(no_peaks)) msg("samples with no peaks on host chromosomes (no peak pie): ", paste(no_peaks, collapse = ", "))
  wide <- dcast(P, contig ~ sample_id, value.var = "n", fill = 0)
  wide <- merge(G[contig %in% hosts], wide, by = "contig", all.x = TRUE)
  for (j in setdiff(names(wide), c("contig", "length"))) set(wide, which(is.na(wide[[j]])), j, 0L)
  wide <- wide[order(match(contig, hosts))][, all_samples := rowSums(.SD), .SDcols = setdiff(names(wide), c("contig", "length"))]
  fwrite(wide, file.path(out, "peaks_per_chromosome.tsv"), sep = "\t")
  msg(sprintf("peaks per chromosome: %s peaks in %d samples on %d host chromosomes -> %s", format(sum(P$n), big.mark = ","),
              uniqueN(P$sample_id), length(hosts), file.path(out, "peaks_per_chromosome.tsv")))
  for (s in unique(P$sample_id)) for (ring in c(TRUE, FALSE)) {
    d <- merge(G[contig %in% hosts], P[sample_id == s, .(contig, value = n)], by = "contig", all.x = TRUE)
    d[is.na(value), value := 0][, panel := "Peaks"]
    draw_pie(d, hosts, sprintf("Peaks per chromosome: %s", s), sprintf("%s peaks", format(sum(d$value), big.mark = ",")),
             file.path(out, sprintf("peaks_%s_%s.png", s, if (ring) "ring" else "noring")), ring)
    n_files <- n_files + 1
  }
  for (ring in c(TRUE, FALSE)) {
    d <- merge(G[contig %in% hosts], P[, .(value = sum(n)), by = contig], by = "contig", all.x = TRUE)
    d[is.na(value), value := 0][, panel := "All samples"]
    draw_pie(d, hosts, sprintf("All %d samples pooled: peaks per chromosome", uniqueN(P$sample_id)), sprintf("%s peaks", format(sum(d$value), big.mark = ",")),
             file.path(out, sprintf("pooled_peaks_%s.png", if (ring) "ring" else "noring")), ring)
    n_files <- n_files + 1
  }
}
if (!is.null(P) && nrow(P)) for (v in group_by) {
  lv <- level_of(v)
  pp <- P[sample_id %in% names(lv)][, grp := group_panel(lv, sample_id)]
  agg <- pp[, .(value = sum(n)), by = .(panel = grp, contig)]
  d <- CJ(panel = levels(pp$grp)[levels(pp$grp) %in% pp$grp], contig = hosts)
  d <- merge(merge(d, G, by = "contig"), agg, by = c("panel", "contig"), all.x = TRUE)[is.na(value), value := 0]
  d[, panel := factor(panel, levels(pp$grp))]
  tot <- d[, .(n = sum(value)), by = panel]
  for (ring in c(TRUE, FALSE)) {
    draw_pie(d, hosts, sprintf("Peaks per chromosome, samples pooled by %s", v),
             sprintf("peaks summed over each group's samples (%s)", paste(sprintf("%s: %s", sub(" \\(n = .*", "", tot$panel), format(tot$n, big.mark = ",", trim = TRUE)), collapse = "; ")),
             file.path(out, sprintf("peaks_by_%s_%s.png", safe(v), if (ring) "ring" else "noring")), ring)
    n_files <- n_files + 1
  }
  if (sum(d$value) != sum(pp$n)) stop("grouped peak counts do not add up for ", v)
  msg(sprintf("peak pies by %s: %d groups, %s peaks", v, nrow(tot), format(sum(tot$n), big.mark = ",")))
}
msg(sprintf("%d pie charts written to %s", n_files, out))
