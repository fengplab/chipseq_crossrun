## Supplementary figure: peaks of FLAG-cGAS versus endogenous cGAS. Makes the three panels from a crossrun output folder,
## each at its printed size (text at --font points), and puts them together:
##   A  pies: share of peaks per chromosome for the FLAG-cGAS sample and the endogenous sample(s), pooled into one pie (extra/redraw_pies.R)
##   B  pairwise similarity of peak sets (Jaccard, from tables/table05b_jaccard_matrix.tsv), the matched pair outlined
##   C  motif family prevalence in the chosen samples (extra/plot_motif_prevalence.R), the FLAG-cGAS column shaded
##   supp_figure  A over B | C, panel letters added
## Usage:
##   Rscript extra/make_supp_flag_endogenous.R --crossrun <outdir> --flag <sample_id> --endo <sample_id> \
##       --labels-file samples_B.tsv [--motif-labels-file samples_C.tsv] [options]
##   --crossrun DIR          crossrun output folder (reads tables/table00b_peaks_by_chromosome.tsv, table05b_jaccard_matrix.tsv,
##                           table02_motif_prevalence.tsv)
##   --flag ID               the FLAG-cGAS sample (sample_id)
##   --endo ID[,ID,...]      the endogenous-cGAS sample, or several to pool into one pie (e.g. all THP1 HSV-1 time points)
##   --endo-label TEXT       title of the endogenous pie when several are pooled (default "Endogenous cGAS"); " (n = k)" is added
##   --pool sum|mean         pooling: sum = peak counts added over the pooled samples, so samples with more peaks weigh more
##                           (default); mean = each sample's per-chromosome share averaged, every sample weighing the same
##   --b-pair ID             the endogenous sample outlined with --flag in panel B (default: the first --endo sample)
##   --labels-file TSV       sample_id, label (tab-separated, with header): the samples of panel B, in this order, and the
##                           labels used in every panel (--flag and --b-pair must be listed)
##   --motif-labels-file TSV the samples of panel C, left to right (default: the --labels-file samples). Panel C is
##                           meant for one cell line, e.g. the THP1 samples
##   --b-order cluster|file  panel B order: clustered by 1 - Jaccard (default) or as in --labels-file
##   --b-values              print the Jaccard index in each panel-B tile (only readable for about 8 samples or fewer)
##   --outdir DIR            default supp_flag_endogenous
##   --width W               figure width, inches (default 7.5); panel A is this wide, B and C share it (--c-width, default 3.55)
##   --pie-height H          height of panel A (default 2.9); --bc-height H: height of B and C (default 4.5)
##   --font PT               text size (default 9); --family NAME (default Arial); --dpi N (default 600)
## Pie values are peak counts per host nuclear chromosome from table00b (viral genomes and chrM are not in it). Each panel
## is written as .png and .pdf with a .tsv of the values drawn; supp_figure.png/.pdf places the panel PNGs 1:1 (needs the
## R package png; without it the panels are still written and the assembly is skipped).
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
a <- commandArgs(TRUE)
if (!length(a) || a[1] %in% c("-h", "--help")) {
  f <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)); h <- readLines(f)
  cat(sub("^## ?", "", h[seq_len(which(!grepl("^##", h))[1] - 1)]), sep = "\n"); quit(status = 0)
}
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
opt <- function(n, d = NULL) { i <- match(n, a); if (is.na(i)) d else a[i + 1] }
has <- function(n) n %in% a
die <- function(...) { message("ERROR: ", ...); quit(status = 1) }
say <- function(...) message("[supp] ", ...)
known <- c("--crossrun", "--flag", "--endo", "--endo-label", "--pool", "--b-pair", "--labels-file", "--motif-labels-file", "--b-order", "--b-values", "--outdir", "--width",
           "--c-width", "--pie-height", "--bc-height", "--font", "--family", "--dpi")
bad <- setdiff(grep("^--", a, value = TRUE), known); if (length(bad)) die("unknown option(s): ", paste(bad, collapse = ", "))
for (o in c("--crossrun", "--flag", "--endo", "--labels-file")) if (is.null(opt(o))) die(o, " is required (see --help)")
cr <- opt("--crossrun"); fl <- opt("--flag"); en_all <- trimws(strsplit(opt("--endo"), ",")[[1]]); en_all <- en_all[nzchar(en_all)]
if (anyDuplicated(en_all) || fl %in% en_all) die("--endo must list distinct samples, none of them --flag")
en <- opt("--b-pair", en_all[1]); pool <- opt("--pool", "sum"); if (!pool %in% c("sum", "mean")) die("--pool must be sum or mean")
tab <- function(n) { f <- file.path(cr, "tables", paste0(n, ".tsv")); if (!file.exists(f)) die("not found: ", f, " (run the crossrun steps that write it)"); f }
f00b <- tab("table00b_peaks_by_chromosome"); f05b <- tab("table05b_jaccard_matrix"); f02 <- tab("table02_motif_prevalence")
od <- opt("--outdir", "supp_flag_endogenous"); dir.create(od, showWarnings = FALSE, recursive = TRUE)
fs <- as.numeric(opt("--font", 9)); fam <- opt("--family", "Arial"); dpi <- as.integer(opt("--dpi", 600))
W <- as.numeric(opt("--width", 7.5)); Wc <- as.numeric(opt("--c-width", 3.55)); Wb <- W - Wc
Ha <- as.numeric(opt("--pie-height", 2.9)); Hbc <- as.numeric(opt("--bc-height", 4.5))
if (Wb < 2) die("--c-width leaves less than 2 in for panel B")
rs <- function(script, args) {                                     # run a sibling script; stop on failure
  st <- system2("Rscript", shQuote(c(file.path(here, script), as.character(args)))); if (st != 0) die(script, " failed (exit ", st, ")") }

read_labels <- function(f) {
  if (!file.exists(f)) die("not found: ", f)
  L <- fread(f, sep = "\t", colClasses = "character")
  if (!all(c("sample_id", "label") %in% names(L))) die(f, " needs a header with columns sample_id and label (tab-separated)")
  L[, `:=`(sample_id = trimws(sample_id), label = trimws(label))][nzchar(sample_id)]
}
LB <- read_labels(opt("--labels-file"))
if (anyDuplicated(LB$sample_id)) die("--labels-file repeats a sample_id")
miss <- setdiff(c(fl, en), LB$sample_id); if (length(miss)) die("--labels-file must list --flag and --b-pair; missing: ", paste(miss, collapse = ", "))
lab <- setNames(LB$label, LB$sample_id)
fmc <- opt("--motif-labels-file", opt("--labels-file"))

## ---------------- A: pies from peak counts per chromosome
P <- fread(f00b)[sample_id %in% c(fl, en_all)]
miss <- setdiff(c(fl, en_all), P$sample_id); if (length(miss)) die("not in ", f00b, ": ", paste(miss, collapse = ", "))
k <- length(en_all)
en_lab <- if (k == 1) { if (en_all %in% names(lab)) unname(lab[en_all]) else en_all } else
  sprintf("%s (n = %d)", opt("--endo-label", "Endogenous cGAS"), k)
fl_lab <- unname(lab[fl])
PF <- P[sample_id == fl, .(panel = fl_lab, contig = chrom, value = as.numeric(n_peaks))]
PE <- P[sample_id %in% en_all]
per <- PE[, .(n = sum(n_peaks)), by = sample_id]
if (any(per$n == 0)) die("pooled samples without peaks: ", paste(per[n == 0, sample_id], collapse = ", "))
PEp <- if (pool == "sum" || k == 1) PE[, .(value = as.numeric(sum(n_peaks))), by = .(contig = chrom)] else
  PE[, sh := n_peaks / sum(n_peaks), by = sample_id][, .(value = mean(sh) * 1e6), by = .(contig = chrom)]   # mean share, x 1e6
PEp[, panel := en_lab]
PA <- rbind(PF[, sample_id := fl], PEp[, .(panel, contig, value, sample_id = paste(en_all, collapse = "+"))])   # sample_id: the pooled ids joined by +
tot <- c(sum(P[sample_id == fl, n_peaks]), sum(per$n)); names(tot) <- c(fl_lab, en_lab)
PA[, peaks_total := tot[panel]]
say(sprintf("A: %s %s peaks; %s %s peaks%s (host nuclear chromosomes, from table00b)", fl_lab, format(tot[1], big.mark = ","), en_lab,
            format(tot[2], big.mark = ","), if (k > 1) sprintf(", pooled by %s over %s (%s peaks)", pool, paste(per$sample_id, collapse = ", "),
                                                               paste(format(per$n, big.mark = ","), collapse = ", ")) else ""))
fA <- file.path(od, "A_pies.tsv"); fwrite(PA, fA, sep = "\t")
for (ext in c("png", "pdf"))
  rs("redraw_pies.R", c(fA, "--ncol", "2", "--width", W, "--height", Ha, "--font", fs, "--radius", "0.14", "--inside-min", "12",
                        "--dpi", dpi, "--panels", paste(c(fl_lab, en_lab), collapse = ","), "--out", file.path(od, paste0("A_pies.", ext))))

## ---------------- B: Jaccard similarity of peak sets
J <- fread(f05b); sid <- J$sample_a; M <- as.matrix(J[, -1]); rownames(M) <- sid
miss <- setdiff(LB$sample_id, intersect(rownames(M), colnames(M)))
if (length(miss)) die("not in ", f05b, ": ", paste(miss, collapse = ", "))
sm <- LB$sample_id; M <- M[sm, sm, drop = FALSE]
if (any(is.na(M[upper.tri(M)]))) die("the Jaccard matrix has missing pairs among the chosen samples")
ord <- if (opt("--b-order", "cluster") == "file" || length(sm) < 3) sm else sm[hclust(as.dist(1 - M), method = "average")$order]
DB <- as.data.table(as.table(M)); setnames(DB, c("a", "b", "jaccard")); DB[, `:=`(a = as.character(a), b = as.character(b))]
DB[a == b, jaccard := NA]
pairD <- DB[(a == fl & b == en) | (a == en & b == fl)]
say(sprintf("B: %d samples, order %s; Jaccard %s vs %s = %.3f (other pairs %.3f to %.3f)", length(sm), if (identical(ord, sm)) "as given" else "clustered",
            lab[fl], lab[en], M[fl, en], min(DB[!is.na(jaccard) & !(a %in% c(fl, en) & b %in% c(fl, en)), jaccard]),
            max(DB[!is.na(jaccard) & !(a %in% c(fl, en) & b %in% c(fl, en)), jaccard])))
th <- theme_minimal(base_size = fs, base_family = fam) +
  theme(text = element_text(size = fs, colour = "black"), axis.text = element_text(size = fs, colour = "black"),
        panel.grid = element_blank(), legend.title = element_text(size = fs), legend.text = element_text(size = fs), plot.margin = margin(2, 2, 2, 2))
jmax <- max(0.1, ceiling(max(DB$jaccard, na.rm = TRUE) * 10) / 10)
pB <- ggplot(DB, aes(a, b, fill = jaccard)) + geom_tile(colour = "white", linewidth = 0.3) +
  geom_tile(data = pairD, fill = NA, colour = "black", linewidth = 0.8) +
  scale_fill_gradient(low = "#F7F5F0", high = "#2E5A88", limits = c(0, jmax), na.value = "grey85", name = "Jaccard",
                      breaks = c(0, jmax / 2, jmax)) +
  scale_x_discrete(limits = ord, labels = lab, position = "top") + scale_y_discrete(limits = rev(ord), labels = lab) +
  labs(x = NULL, y = NULL) + coord_fixed() + th +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5), legend.position = "bottom",
        legend.key.height = unit(6, "pt"), legend.key.width = unit(18, "pt"), legend.title = element_text(vjust = 0.9)) +
  guides(fill = guide_colourbar(title.position = "left"))
if (has("--b-values")) pB <- pB + geom_text(aes(label = ifelse(is.na(jaccard), "", sprintf("%.2f", jaccard))), size = (fs - 3) / .pt, family = fam)
for (ext in c("png", "pdf")) ggsave(file.path(od, paste0("B_similarity.", ext)), pB, width = Wb, height = Hbc, dpi = dpi, bg = "white",
                                    device = if (ext == "pdf" && capabilities("cairo")) grDevices::cairo_pdf else NULL)
fwrite(DB[order(match(a, ord), match(b, ord))][, `:=`(label_a = unname(lab[a]), label_b = unname(lab[b]))], file.path(od, "B_similarity.tsv"), sep = "\t")

## ---------------- C: motif family prevalence in the chosen samples
for (ext in c("png", "pdf"))
  rs("plot_motif_prevalence.R", c(f02, "--labels-file", fmc, if (fl %in% read_labels(fmc)$sample_id) c("--highlight", fl),
                                  "--width", Wc, "--height", Hbc, "--font", fs, "--family", fam, "--dpi", dpi,
                                  "--out", file.path(od, paste0("C_motif_prevalence.", ext))))

## ---------------- assembly: panel PNGs placed 1:1 (no rescaling, so text keeps its size)
if (!requireNamespace("png", quietly = TRUE)) {
  say("R package png not installed: panels written, assembly skipped (place A_pies, B_similarity, C_motif_prevalence by hand)")
} else {
  im <- lapply(c(A = "A_pies.png", B = "B_similarity.png", C = "C_motif_prevalence.png"), function(f) png::readPNG(file.path(od, f)))
  dims <- sapply(im, function(x) dim(x)[2:1] / dpi)                 # width, height in inches at the panels' resolution
  top <- 0.15; gap <- 0.2; y2 <- top + dims[2, "A"] + gap   # gap: a row above B and C for their letters
  H <- y2 + max(dims[2, "B"], dims[2, "C"])
  draw <- function() {
    grid::grid.newpage()
    put <- function(x, l, t, w, h) grid::grid.raster(x, x = grid::unit(l, "in"), y = grid::unit(H - t, "in"), width = grid::unit(w, "in"),
                                                     height = grid::unit(h, "in"), just = c("left", "top"), interpolate = FALSE)
    put(im$A, 0, top, dims[1, "A"], dims[2, "A"])
    put(im$B, 0, y2, dims[1, "B"], dims[2, "B"])
    put(im$C, W - dims[1, "C"], y2, dims[1, "C"], dims[2, "C"])
    tag <- function(s, l, t) grid::grid.text(s, x = grid::unit(l, "in"), y = grid::unit(H - t, "in"), just = c("left", "top"),
                                             gp = grid::gpar(fontsize = fs + 3, fontface = "bold", fontfamily = fam))
    tag("A", 0.02, 0.02); tag("B", 0.02, top + dims[2, "A"]); tag("C", W - dims[1, "C"], top + dims[2, "A"])
  }
  png(file.path(od, "supp_figure.png"), width = W, height = H, units = "in", res = dpi, bg = "white"); draw(); invisible(dev.off())
  if (capabilities("cairo")) { grDevices::cairo_pdf(file.path(od, "supp_figure.pdf"), width = W, height = H, bg = "white"); draw(); invisible(dev.off()) }
  say(sprintf("supp_figure: %.2f x %.2f in (panels placed at their own size; the PDF holds them as %d dpi images: use the panel PDFs for vector art)", W, H, dpi))
}
say("done -> ", od)
