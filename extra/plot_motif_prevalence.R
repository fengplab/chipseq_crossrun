## Motif family prevalence for chosen samples, drawn for a figure panel at its printed size. Same data as the report's
## fig01_motif_prevalence_heatmap (tables/table02_motif_prevalence.tsv): the share of each sample's peaks that contain at
## least one site of a motif family, plus whether the family is detected in that sample.
## Usage: Rscript extra/plot_motif_prevalence.R <crossrun_out>/tables/table02_motif_prevalence.tsv [options]
##   --labels-file TSV       columns sample_id, label (tab-separated, with that header). Chooses the samples: only those
##                           listed are drawn, left to right in the file's order (same file as extra/plot_chrom_heatmap.R)
##   --samples a,b,...       which samples and in what order (overrides the --labels-file selection; its labels still apply)
##   --labels a=Label,...    column labels on the command line, added to or replacing --labels-file labels
##   --families F001,F007    draw exactly these families, top to bottom in this order
##   --top N                 otherwise: the N families most enriched over the genomic background in each chosen sample
##                           (fold enrichment, then % of peaks; as in the report), combined (default 5)
##   --max-families M        at most M rows: the families with the highest mean % of peaks over the chosen samples (default 15)
##   --style dots|heatmap    dots: size = % of peaks, filled = detected (default); heatmap: colour = % of peaks
##   --row-order prevalence|cluster|given   rows by mean % of peaks (default), by clustering, or as in --families
##   --highlight a,b         shade these samples' columns (sample_ids)
##   --width W --height H    inches (default 3.6 wide; height fitted to the rows)
##   --font PT               text size (default 9); --family NAME (default Arial); --dpi N (default 600)
##   --out FILE              .png, .pdf or .svg (default motif_prevalence.pdf); the plotted values go to <out>.tsv
## "Detected" is the table's own call (column detected). Families never hit in a chosen sample still show, as 0%.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
a <- commandArgs(TRUE)
if (!length(a) || a[1] %in% c("-h", "--help")) {
  f <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)); h <- readLines(f)
  cat(sub("^## ?", "", h[seq_len(which(!grepl("^##", h))[1] - 1)]), sep = "\n"); quit(status = 0)
}
opt <- function(n, d = NULL) { i <- match(n, a); if (is.na(i)) d else a[i + 1] }
die <- function(...) { message("ERROR: ", ...); quit(status = 1) }
say <- function(...) message("[prevalence] ", ...)
known <- c("--labels-file", "--samples", "--labels", "--families", "--top", "--max-families", "--style", "--row-order",
           "--highlight", "--width", "--height", "--font", "--family", "--dpi", "--out")
bad <- setdiff(grep("^--", a, value = TRUE), known); if (length(bad)) die("unknown option(s): ", paste(bad, collapse = ", "))
if (!file.exists(a[1])) die("not found: ", a[1])
P <- fread(a[1])
need <- c("family_id", "pwm_consensus", "sample_id", "n_hit", "pct_peaks", "enrichment", "detected")
miss <- setdiff(need, names(P)); if (length(miss)) die(a[1], " lacks columns: ", paste(miss, collapse = ", "))
P[, `:=`(sample_id = as.character(sample_id), family_id = as.character(family_id),
         detected = toupper(as.character(detected)) == "TRUE")]
say(sprintf("read %d rows: %d families x %d samples", nrow(P), uniqueN(P$family_id), uniqueN(P$sample_id)))

## ---- samples: --labels-file selects and orders (as in plot_chrom_heatmap.R), --samples overrides
L <- NULL
if (!is.null(opt("--labels-file"))) {
  if (!file.exists(opt("--labels-file"))) die("not found: ", opt("--labels-file"))
  L <- fread(opt("--labels-file"), sep = "\t", colClasses = "character")
  if (!all(c("sample_id", "label") %in% names(L))) die("--labels-file needs a header with columns sample_id and label (tab-separated)")
  L[, `:=`(sample_id = trimws(sample_id), label = trimws(label))]; L <- L[nzchar(sample_id)]
  if (!nrow(L)) die("--labels-file lists no samples")
  dup <- unique(L$sample_id[duplicated(L$sample_id)]); if (length(dup)) die("--labels-file repeats: ", paste(dup, collapse = ", "))
  emp <- L$sample_id[is.na(L$label) | !nzchar(L$label)]; if (length(emp)) die("--labels-file has no label for: ", paste(emp, collapse = ", "))
  miss <- setdiff(L$sample_id, P$sample_id)
  if (length(miss)) die("--labels-file samples not in the table (check spelling and suffixes such as _REP1): ", paste(miss, collapse = ", "))
}
sm <- if (!is.null(opt("--samples"))) strsplit(opt("--samples"), ",")[[1]] else if (!is.null(L)) L$sample_id else unique(P$sample_id)
miss <- setdiff(sm, P$sample_id); if (length(miss)) die("samples not in the table: ", paste(miss, collapse = ", "))
if (anyDuplicated(sm)) die("--samples repeats: ", paste(unique(sm[duplicated(sm)]), collapse = ", "))
left <- setdiff(unique(P$sample_id), sm)
say(sprintf("drawing %d of %d samples, left to right: %s", length(sm), uniqueN(P$sample_id), paste(sm, collapse = ", ")))
if (length(left)) say("left out: ", paste(left, collapse = ", "))
lab <- setNames(sm, sm)
if (!is.null(L)) { hit <- intersect(L$sample_id, sm); lab[hit] <- L$label[match(hit, L$sample_id)] }
if (!is.null(opt("--labels"))) { kv <- strsplit(strsplit(opt("--labels"), ",")[[1]], "=")
  k <- sapply(kv, `[`, 1); miss <- setdiff(k, sm); if (length(miss)) die("--labels for samples not shown: ", paste(miss, collapse = ", "))
  lab[k] <- sapply(kv, `[`, 2) }
unl <- names(lab)[lab == names(lab)]; if (length(unl)) say("no label given for (shown as the sample ID): ", paste(unl, collapse = ", "))
hl <- if (!is.null(opt("--highlight"))) strsplit(opt("--highlight"), ",")[[1]] else character()
miss <- setdiff(hl, sm); if (length(miss)) die("--highlight samples not drawn: ", paste(miss, collapse = ", "))

## ---- families: given, or the report's rule applied to the chosen samples only
D <- P[sample_id %in% sm]
cons <- unique(P[, .(family_id, pwm_consensus)])[, .SD[1], by = family_id]
if (!is.null(opt("--families"))) {
  fam <- strsplit(opt("--families"), ",")[[1]]
  miss <- setdiff(fam, P$family_id); if (length(miss)) die("families not in the table: ", paste(miss, collapse = ", "))
  rule <- "as given"
} else {
  topn <- as.integer(opt("--top", 5)); if (is.na(topn) || topn < 1) die("--top must be a positive whole number")
  fam <- unique(D[n_hit > 0][order(-enrichment, -pct_peaks), head(family_id, topn), by = sample_id]$V1)
  if (!length(fam)) die("no family has a site in the chosen samples' peaks")
  rule <- sprintf("the %d most enriched over the genomic background in each chosen sample, combined", topn)
}
G <- CJ(family_id = fam, sample_id = sm)                                      # every family x sample, missing = no site
G <- merge(G, D[, .(family_id, sample_id, pct_peaks, detected, enrichment, n_hit)], by = c("family_id", "sample_id"), all.x = TRUE)
G[is.na(pct_peaks), `:=`(pct_peaks = 0, detected = FALSE, n_hit = 0L)]
mx <- as.integer(opt("--max-families", 15))
if (is.null(opt("--families")) && length(fam) > mx) {
  keep <- G[, .(m = mean(pct_peaks)), by = family_id][order(-m)][seq_len(mx), family_id]
  say(sprintf("%d families qualify; keeping the %d with the highest mean %% of peaks (--max-families)", length(fam), mx))
  G <- G[family_id %in% keep]; fam <- fam[fam %in% keep]; rule <- paste0(rule, sprintf("; top %d by mean %% of peaks", mx))
}
ro <- opt("--row-order", if (!is.null(opt("--families"))) "given" else "prevalence")
if (!ro %in% c("prevalence", "cluster", "given")) die("--row-order must be prevalence, cluster or given")
fam_ord <- switch(ro, given = fam,
  prevalence = G[, .(m = mean(pct_peaks)), by = family_id][order(-m), family_id],
  cluster = { M <- as.matrix(dcast(G, family_id ~ sample_id, value.var = "pct_peaks")[, -1]); rownames(M) <- sort(unique(G$family_id))
              if (nrow(M) > 2) rownames(M)[hclust(dist(M), "average")$order] else rownames(M) })
G[, cons := cons$pwm_consensus[match(family_id, cons$family_id)]]
rowlab <- setNames(paste(fam_ord, cons$pwm_consensus[match(fam_ord, cons$family_id)]), fam_ord)
G[, `:=`(family = family_id, sample = sample_id)]

## ---- draw
fs <- as.numeric(opt("--font", 9)); fam_font <- opt("--family", "Arial"); style <- opt("--style", "dots")
if (!style %in% c("dots", "heatmap")) die("--style must be dots or heatmap")
W <- as.numeric(opt("--width", 3.6)); lh <- fs * 1.25 / 72
H <- as.numeric(opt("--height", length(fam_ord) * max(0.22, lh * 1.5) + max(nchar(lab)) * fs * 0.55 / 72 + 0.4 + 5 * lh))   # rows + column labels + two key rows below
out <- opt("--out", "motif_prevalence.pdf"); ext <- tolower(tools::file_ext(out))
if (!ext %in% c("png", "pdf", "svg")) die("--out must end in .png, .pdf or .svg")
th <- theme_minimal(base_size = fs, base_family = fam_font) +
  theme(text = element_text(size = fs, colour = "black"), axis.text = element_text(size = fs, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid = element_blank(),
        legend.title = element_text(size = fs), legend.text = element_text(size = fs), plot.margin = margin(2, 2, 2, 2),
        legend.position = "bottom", legend.box = "vertical", legend.box.just = "left", legend.justification = "left", legend.margin = margin(0, 0, 0, 0), legend.spacing.x = unit(4, "pt"))
p <- ggplot(G, aes(sample, family))
if (length(hl)) p <- p + geom_tile(data = G[sample_id %in% hl], width = 1, height = 1, fill = "#F3E3D6")
if (style == "dots") {
  p <- p + geom_point(aes(size = pct_peaks, fill = detected), shape = 21, colour = "#14213D", stroke = 0.3) +
    scale_size_area(max_size = 5, name = "% of peaks\nwith a site", breaks = function(l) pretty(c(0, l[2]), 3)[-1]) +
    scale_fill_manual(values = c(`TRUE` = "#2E5A88", `FALSE` = "#FFFFFF"), breaks = c(TRUE, FALSE), labels = c("yes", "no"), name = "Detected") +
    guides(fill = guide_legend(order = 1, override.aes = list(size = 3)), size = guide_legend(order = 2, title.position = "left")) +
    th + theme(panel.grid.major = element_line(colour = "#ECE8DF", linewidth = 0.2))
} else {
  p <- p + geom_tile(aes(fill = pct_peaks), colour = "white", linewidth = 0.3, width = if (length(hl)) 0.9 else 1, height = if (length(hl)) 0.9 else 1) +
    scale_fill_gradient(low = "#FFFFFF", high = "#0C447C", name = "% of peaks\nwith a site") + th
}
p <- p + scale_x_discrete(limits = sm, labels = lab) + scale_y_discrete(limits = rev(fam_ord), labels = rowlab) + labs(x = NULL, y = NULL)   # limits fix the order (a highlight layer would otherwise reorder)
ggsave(out, p, width = W, height = H, dpi = as.numeric(opt("--dpi", 600)), bg = "white",
       device = if (ext == "pdf" && capabilities("cairo")) grDevices::cairo_pdf else NULL)
tsv <- sub("\\.[A-Za-z]+$", ".tsv", out)
fwrite(G[order(match(family_id, fam_ord), match(sample_id, sm)),
         .(family_id, pwm_consensus = cons, sample_id, label = unname(lab[sample_id]), pct_peaks = signif(pct_peaks, 4), detected, n_hit,
           enrichment = signif(enrichment, 4))], tsv, sep = "\t")
say(sprintf("%d families (%s) x %d samples -> %s (%.2f x %.2f in, %s %g pt, %s); values -> %s",
            length(fam_ord), rule, length(sm), out, W, H, fam_font, fs, style, tsv))
