#!/usr/bin/env Rscript
## Spacing histogram (as in Figure 3) for chosen motif families, from an existing crossrun output folder.
##
## Usage:
##   Rscript R/plot_spacing.R <crossrun_outdir> <family_id[,family_id...]> [options]
## Options:
##   --samples a,b,c   only these samples (default: all)
##   --by-sample       one panel per sample instead of pooling samples
##   --max-bp N        x-axis limit in bp (default: motifs.spacing_max used in step 05, usually 400)
##   --min-bp N        smallest spacing shown (default 2, as in Figure 3)
##   --period          dashed lines at the most common spacing and its multiples
##   --out FILE        output file; .png, .pdf or .svg (default: <outdir>/figures/spacing_<families>.png)
##   --width W --height H   size in inches (default: chosen from the number of panels)
##
## Example:
##   Rscript R/plot_spacing.R crossrun_out_human_broad F001 --period --max-bp 60
##   Rscript R/plot_spacing.R crossrun_out_human_broad F001,F015 --by-sample --samples THP1_HSV1_1h_IP-cGAS,THP1_HSV1_8h_IP-cGAS
## Reads work/scan/<sample>.peaks.spacing.tsv (step 05) and work/motifs/family_rep_features.tsv (step 04).
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 2 || a[1] %in% c("-h", "--help")) {
  cat(sub("^## ?", "", grep("^##", readLines(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), value = TRUE)), sep = "\n")
  quit(status = if (length(a) < 2) 1 else 0)
}
outdir <- a[1]; fams <- strsplit(a[2], ",")[[1]]; a <- a[-(1:2)]
opt <- function(name, default = NULL) { i <- match(name, a); if (is.na(i)) default else a[i + 1] }
flag <- function(name) name %in% a
samples <- opt("--samples"); if (!is.null(samples)) samples <- strsplit(samples, ",")[[1]]
min_bp <- as.numeric(opt("--min-bp", 2)); max_bp <- opt("--max-bp"); if (!is.null(max_bp)) max_bp <- as.numeric(max_bp)
by_sample <- flag("--by-sample")

scan_dir <- file.path(outdir, "work", "scan")
files <- list.files(scan_dir, "\\.peaks\\.spacing\\.tsv$", full.names = TRUE)
if (!length(files)) stop("no spacing files in ", scan_dir, " (run step 05 first, or check the output folder)")
sp <- rbindlist(lapply(files, fread), fill = TRUE)
if (!nrow(sp)) stop("the spacing files are empty")
miss <- setdiff(fams, unique(sp$family_id))
if (length(miss)) stop("no spacing data for: ", paste(miss, collapse = ", "),
                       " (families need at least two same-strand sites in one peak; check the ID in table01)")
if (!is.null(samples)) {
  bad <- setdiff(samples, unique(sp$sample_id))
  if (length(bad)) stop("unknown sample(s): ", paste(bad, collapse = ", "))
  sp <- sp[sample_id %in% samples]
}
sp <- sp[family_id %in% fams & spacing >= min_bp]
if (!is.null(max_bp)) sp <- sp[spacing <= max_bp]
if (!nrow(sp)) stop("nothing left to plot after the filters")

feat_f <- file.path(outdir, "work", "motifs", "family_rep_features.tsv")
cons <- if (file.exists(feat_f)) { f <- fread(feat_f); setNames(f$pwm_consensus, f$family_id) } else character()
lab <- function(id) ifelse(is.na(cons[id]), id, paste(id, cons[id]))

keys <- c("family_id", if (by_sample) "sample_id")
d <- sp[, .(n = sum(n)), by = c(keys, "spacing")]
d[, family := factor(lab(family_id), levels = lab(fams))]
empty <- setdiff(fams, unique(d$family_id))
if (length(empty)) message("no site pairs in the selected samples/range for: ", paste(empty, collapse = ", "), " (shown as empty panels)")
xm <- if (is.null(max_bp)) max(20, max(d$spacing)) else max_bp
p <- ggplot(d, aes(spacing, n)) + geom_col(width = 1, fill = "#2a78d6") +
  labs(x = "Distance between adjacent same-strand sites (bp)", y = "Site pairs",
       title = "Motif site spacing",
       subtitle = sprintf("%s; %s", if (by_sample) "per sample" else sprintf("summed over %d samples", uniqueN(sp$sample_id)),
                          "peaks only")) +
  coord_cartesian(xlim = c(max(0, min_bp - 1), xm + 0.5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_bw(11) + theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 9))
if (flag("--period")) {
  md <- d[, .(mode = spacing[which.max(n)]), by = keys]
  per <- md[mode >= 3, .(spacing = seq(mode, xm, by = mode)), by = c(keys, "mode")]
  per[, family := factor(lab(family_id), levels = lab(fams))]
  p <- p + geom_vline(data = per, aes(xintercept = spacing), linetype = "dashed", colour = "#d1495b", linewidth = 0.4)
  print(md)
}
npan <- if (by_sample) length(fams) * uniqueN(d$sample_id) else length(fams)
if (by_sample) p <- p + facet_grid(sample_id ~ family, scales = "free_y", drop = FALSE) else if (length(fams) > 1) p <- p + facet_wrap(~family, scales = "free_y", drop = FALSE)
ncol <- if (by_sample) length(fams) else min(length(fams), 3)
nrow <- ceiling(npan / ncol)
w <- as.numeric(opt("--width", min(16, 3 + 4 * ncol))); h <- as.numeric(opt("--height", min(40, 1.5 + 2.4 * nrow)))
out <- opt("--out", file.path(outdir, "figures", sprintf("spacing_%s%s.png", paste(fams, collapse = "_"), if (by_sample) "_by_sample" else "")))
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
ggsave(out, p, width = w, height = h, dpi = 200, limitsize = FALSE)
cat("wrote", out, "\n")
