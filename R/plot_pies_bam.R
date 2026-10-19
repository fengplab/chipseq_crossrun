#!/usr/bin/env Rscript
## Pie charts of reads per contig for ONE sample, straight from its filtered BAM (no crossrun run needed).
##
## Usage:
##   Rscript R/plot_pies_bam.R --bam S1.bam [--control S1_input.bam] [--peaks S1_peaks.narrowPeak] [options]
##
## Options:
##   --bam FILE            indexed BAM of the sample (host chromosomes + chrM + viral genomes)
##   --control FILE        indexed control BAM: drawn as a second pie beside the IP
##   --peaks FILE          narrowPeak/broadPeak/BED of the sample: adds peaks-per-chromosome pies
##   --label NAME          sample name for titles and file names (default: BAM file name)
##   --viral LIST          viral genomes as DISPLAY=CONTIG, comma-separated
##                         (default HSV-1=NC_001806.2,KSHV=NC_009333.1,ADV=AC_000008.1)
##   --exclude LIST        host contigs kept out of the "host chromosomes" set but shown in the pies (default chrM)
##   --viral-mapq N        count viral reads only at MAPQ >= N (default 0; host contigs are never filtered)
##   --out DIR             output folder (default pies_<label>)
##
## Reads counted: primary, non-duplicate, QC-passed alignments (samtools view -F 0xF04), each mate separately,
## in one streaming pass. Host chromosomes = chr1-22, chrX, chrY (also 1-22/X/Y or T2T RefSeq/GenBank accessions).
## Contigs that are neither host chromosomes, excluded contigs nor listed viruses (unplaced scaffolds) are left out.
##
## Output (each with a TSV of its numbers):
##   coverage_<label>_with_viruses_{ring,noring}.png   host chromosomes, chrM and viral genomes
##   coverage_<label>_host_{ring,noring}.png           host chromosomes and chrM
##   coverage_rpkm_<label>_{with_viruses,host}.png     share of RPKM (relative depth), no ring
##   peaks_<label>_{ring,noring}.png                   peaks per host chromosome (with --peaks)
## Inner ring: each contig's share of the total length of the contigs in the pie (what an even spread would give).
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

a <- commandArgs(trailingOnly = TRUE)
self <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
if (!length(a) || any(a %in% c("-h", "--help"))) {
  cat(sub("^## ?", "", grep("^##", readLines(self), value = TRUE)), sep = "\n"); quit(status = if (length(a)) 0 else 1)
}
known <- c("--bam", "--control", "--peaks", "--label", "--viral", "--exclude", "--viral-mapq", "--out")
bad <- setdiff(grep("^--", a, value = TRUE), known)
die <- function(...) { message("ERROR: ", ...); quit(status = 1) }
if (length(bad)) die("unknown option(s): ", paste(bad, collapse = ", "), " (see --help)")
opt <- function(n, d = NULL) { i <- match(n, a); if (is.na(i)) d else a[i + 1] }
say <- function(...) message("[pies] ", ...)
if (!nzchar(Sys.which("samtools"))) die("samtools not found on PATH")

`%||%` <- function(x, y) if (is.null(x) || !length(x) || !nzchar(x)) y else x
bam <- opt("--bam") %||% die("--bam is required")
ctl <- opt("--control"); pk <- opt("--peaks")
for (f in c(bam, ctl, pk)) if (!file.exists(f)) die("file not found: ", f)
label <- opt("--label") %||% sub("\\.mLb\\.clN(\\.sorted)?$", "", sub("\\.bam$", "", basename(bam)))
vir <- strsplit(strsplit(opt("--viral", "HSV-1=NC_001806.2,KSHV=NC_009333.1,ADV=AC_000008.1"), ",")[[1]], "=")
if (any(lengths(vir) != 2)) die("--viral must look like HSV-1=NC_001806.2,KSHV=NC_009333.1")
vir <- setNames(sapply(vir, `[`, 1), sapply(vir, `[`, 2))          # contig name in the BAM -> display name
excl <- strsplit(opt("--exclude", "chrM"), ",")[[1]]
vmq <- as.integer(opt("--viral-mapq", 0))
out <- opt("--out", paste0("pies_", label)); dir.create(out, showWarnings = FALSE, recursive = TRUE)

## ---------------- contig names: chr-style for host chromosomes ----------------
t2t <- c(setNames(paste0("chr", 1:22), sprintf("NC_%06d.1", 60924 + 1:22)), setNames(paste0("chr", 1:22), sprintf("CP%06d.2", 68278 - 1:22)),
         NC_060947.1 = "chrX", CP068255.2 = "chrX", NC_060948.1 = "chrY", CP086569.2 = "chrY", NC_012920.1 = "chrM", CP068254.1 = "chrM")
canon <- function(x) {
  y <- ifelse(x %in% names(vir), vir[x], ifelse(x %in% names(t2t), t2t[x], x))
  y <- ifelse(grepl("^([0-9]{1,2}|X|Y)$", y), paste0("chr", y), y)
  ifelse(y %in% c("MT", "M", "chrMT"), "chrM", y)
}

## ---------------- reads per contig (one pass per BAM) ----------------
count_bam <- function(f) {
  h <- fread(cmd = sprintf("samtools idxstats '%s'", f), header = FALSE, col.names = c("bam_contig", "length", "mapped", "unmapped"))[length > 0]
  x <- fread(cmd = sprintf("samtools view -F 0xF04 '%s' | cut -f3,5 | awk 'BEGIN{FS=OFS=\"\\t\"} {n[$1 FS $2]++} END {for (k in n) print k, n[k]}'", f),
             header = FALSE, col.names = c("bam_contig", "mapq", "n"))
  if (!nrow(x)) die(basename(f), ": no primary, non-duplicate alignments")
  x[, contig := canon(bam_contig)]
  x[, keep := !(contig %in% vir) | mapq >= vmq]
  r <- x[, .(reads_all = sum(n), reads = sum(n[keep])), by = .(bam_contig, contig)]
  r <- merge(h[, .(bam_contig, length)], r, by = "bam_contig", all.x = TRUE)
  r[is.na(reads), `:=`(reads = 0, reads_all = 0)][, contig := canon(bam_contig)][]
}
libs <- list(IP = bam); if (!is.null(ctl)) libs$Control <- ctl
CR <- rbindlist(lapply(names(libs), function(l) { say("counting ", basename(libs[[l]])); count_bam(libs[[l]])[, panel := l] }))
host_re <- "^chr([0-9]{1,2}|X|Y)$"
CR[, group := fifelse(grepl(host_re, contig), "nuclear", fifelse(contig %in% excl, "excluded", fifelse(contig %in% vir, "viral", "other")))]
miss_v <- setdiff(vir, CR$contig); if (length(miss_v)) say("note: not in the BAM header: ", paste(miss_v, collapse = ", "))
n_other <- CR[group == "other" & panel == "IP", sum(reads)]
if (n_other > 0) say(sprintf("note: %s reads on %d other contigs (unplaced scaffolds etc.) are left out", format(n_other, big.mark = ","),
                             CR[group == "other" & panel == "IP", .N]))
CR <- CR[group != "other"]

## genome order: chr1..chr22, chrX, chrY, excluded contigs, viral genomes in --viral order
rank_of <- function(x, g) { n <- suppressWarnings(as.numeric(sub("^chr", "", x)))
  ifelse(g == "viral", 2000 + match(x, vir), ifelse(g == "excluded", 1000, ifelse(!is.na(n), n, ifelse(x == "chrX", 100, 101)))) }
CT <- unique(CR[, .(contig, group, length)])[, rank := rank_of(contig, group)][order(rank, contig)]
order_all <- CT$contig
nuc <- CT[group == "nuclear", contig]; exc <- CT[group == "excluded", contig]; vrl <- CT[group == "viral", contig]

## ---------------- drawing (same style as the crossrun pies) ----------------
VCOL <- c("#d62728", "#ff7f0e", "#e377c2", "#bcbd22", "#17becf")
pal_for <- function(contigs) {
  hosts <- setdiff(contigs, c(vrl, exc))
  hc <- grDevices::hcl(h = seq(15, 375, length.out = length(hosts) + 1)[-1], c = 45, l = rep(c(62, 78), length.out = length(hosts)))
  cols <- c(setNames(hc, hosts), setNames(rep("#4d4d4d", length(intersect(contigs, exc))), intersect(contigs, exc)),
            setNames(VCOL[seq_along(intersect(contigs, vrl))], intersect(contigs, vrl)))
  cols[contigs]
}
fmt <- function(x) ifelse(x >= 1, sprintf("%.1f%%", x), ifelse(x >= 0.01, sprintf("%.2f%%", x), ifelse(x > 0, sprintf("%.1e%%", x), "0%")))
draw_pie <- function(d, contigs, title, subtitle, file, ring) {
  d <- copy(d)[contig %in% contigs]
  d[, share := 100 * value / sum(value), by = panel]
  d[, len_share := 100 * length / sum(length), by = panel]
  d[, contig := factor(contig, contigs)][, panel := factor(panel, intersect(c("IP", "Control", "Peaks"), unique(panel)))]
  setorder(d, panel, contig)
  lg <- dcast(d, contig ~ panel, value.var = "share"); pn <- setdiff(names(lg), "contig")
  lg[, lab := paste0(contig, ":  ", do.call(paste, c(lapply(pn, function(p) paste(if (length(pn) > 1) p else "", fmt(get(p)))), sep = " | ")))]
  labs_map <- setNames(trimws(gsub("\\s+", " ", lg$lab)), as.character(lg$contig))
  cols <- pal_for(contigs)
  d[, `:=`(y1 = cumsum(share) - share, y2 = cumsum(share)), by = panel]
  d[, `:=`(l1 = cumsum(len_share) - len_share, l2 = cumsum(len_share)), by = panel]
  ox <- if (ring) c(1.05, 2) else c(0, 2)
  p <- ggplot(d) + geom_rect(aes(xmin = ox[1], xmax = ox[2], ymin = y1, ymax = y2, fill = contig), colour = "white", linewidth = 0.2)
  if (ring) p <- p + geom_rect(aes(xmin = 0.35, xmax = 1, ymin = l1, ymax = l2, fill = contig), colour = "white", linewidth = 0.2, alpha = 0.55)
  big <- d[share >= 2]
  if (nrow(big)) p <- p + geom_text(data = big, aes(x = if (ring) 1.55 else 1.35, y = (y1 + y2) / 2, label = paste0(contig, "\n", fmt(share))),
                                    size = 2.3, lineheight = 0.85, colour = "grey10")
  p <- p + coord_polar(theta = "y", direction = 1) + facet_wrap(~panel) +
    scale_x_continuous(limits = c(0, 2.05), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 100 + 1e-6), expand = c(0, 0)) +
    scale_fill_manual(values = cols, labels = labs_map[names(cols)], name = NULL, drop = FALSE) +
    labs(title = title, subtitle = paste(c(subtitle, if (ring) "inner ring: share of contig length (an even spread)"), collapse = "; ")) +
    theme_void(10) + theme(legend.text = element_text(size = 7), legend.key.size = unit(0.35, "cm"), strip.text = element_text(face = "bold", size = 10),
                           plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 8))
  ggsave(file, p, width = 5.2 * uniqueN(d$panel) + 4, height = max(5.5, 0.17 * length(contigs) + 1.5), dpi = 150, bg = "white")
  fwrite(d[, .(panel, contig, value, share_pct = signif(share, 5), length, length_share_pct = signif(len_share, 5))], sub("\\.png$", ".tsv", file), sep = "\t")
  say("wrote ", file)
}

mq_txt <- if (vmq > 0) sprintf("viral reads at MAPQ >= %d", vmq) else "viral reads not MAPQ-filtered"
sets <- list(with_viruses = order_all, host = c(nuc, exc))
set_txt <- c(with_viruses = "host chromosomes, chrM and viral genomes", host = "host chromosomes and chrM")
for (nm in names(sets)) {
  cs <- intersect(sets[[nm]], CR$contig)
  for (ring in c(TRUE, FALSE))
    draw_pie(CR[, .(panel, contig, value = reads, length)], cs, sprintf("Reads per contig: %s", label),
             sprintf("%s; %s", set_txt[[nm]], mq_txt), file.path(out, sprintf("coverage_%s_%s_%s.png", label, nm, if (ring) "ring" else "noring")), ring)
  ## RPKM: reads per kb of contig per million reads in the pie, i.e. relative depth (no ring)
  pd <- CR[contig %in% cs][, rpkm := reads / (length / 1e3) / (sum(reads) / 1e6), by = panel]
  draw_pie(pd[, .(panel, contig, value = rpkm, length)], cs, sprintf("Coverage per contig (RPKM): %s", label),
           sprintf("slice = share of RPKM (relative depth); %s; %s", set_txt[[nm]], mq_txt),
           file.path(out, sprintf("coverage_rpkm_%s_%s.png", label, nm)), ring = FALSE)
  fwrite(pd[, .(panel, contig, length, reads, rpkm = signif(rpkm, 5))], file.path(out, sprintf("coverage_rpkm_%s_%s_values.tsv", label, nm)), sep = "\t")
}

## ---------------- peaks per host chromosome ----------------
if (!is.null(pk)) {
  P <- fread(pk, header = FALSE, select = 1, col.names = "chrom")[, chrom := canon(chrom)][, .(value = .N), by = .(contig = chrom)]
  off <- P[!contig %in% nuc, sum(value)]
  if (off > 0) say(sprintf("note: %d peaks not on host chromosomes are left out of the peak pies", off))
  d <- merge(CT[group == "nuclear", .(contig, length)], P, by = "contig", all.x = TRUE)[is.na(value), value := 0][, panel := "Peaks"]
  for (ring in c(TRUE, FALSE))
    draw_pie(d, intersect(order_all, nuc), sprintf("Peaks per chromosome: %s", label), sprintf("%s peaks", format(sum(d$value), big.mark = ",")),
             file.path(out, sprintf("peaks_%s_%s.png", label, if (ring) "ring" else "noring")), ring)
}
say("done: ", out)
