## Bin-size tuning (standalone; run by tune/run_tuning.sh, not by run_all.sh).
## Input: base-resolution bin tables of the tuning chromosomes (annotation, per-library signal, optional reference motif sites).
## For each candidate size w and each IP/control pair, separately for satellite (CenSat) and other bins:
##   R(w)   signal retained = covariance of the two read-name halves' binding values across bins (noise is independent
##          between halves, so this estimates real signal variance); reported relative to the floor (smallest size)
##   rho(w) split-half correlation (reproducibility)
##   beta(w) per reference motif: within-stratum slope of binding on site count, relative to the floor
## Chosen size: the largest w such that every size up to w keeps median R(w)/R(floor) >= 1 - tolerance in both groups and,
## with reference motifs, median beta(w)/beta(floor) >= 1 - tolerance.
## Args: <tuning config.yaml> <tuning outdir>
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "utils.R"))
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "bin_io.R"))
a <- commandArgs(TRUE); cfg <- args_config(); TO <- a[2]
tu <- cfg$tuning %||% list()
sizes <- sort(as.integer(unlist(tu$sizes %||% c(60, 100, 150, 200, 250, 300, 400, 500, 750, 1000, 2000))))
base <- as.integer(tu$base %||% 10); tol <- as.numeric(tu$tolerance %||% 0.05)
floor_w <- min(sizes); dchr <- tu$drift_chromosome %||% tail(unlist(tu$chromosomes), 1)
BD <- file.path(TO, "bins", paste0("w", base))
annot <- load_annot(BD)
libs <- fread(file.path(BD, "libraries.tsv"))
P <- lib_pairs(libs)
if (!nrow(P)) stop("no IP sample with a control library in the tuning set")
msg(sprintf("tuning: %d IP/control pairs, %d chromosomes (%s), %s bins at base %d bp; sizes %s; tolerance %.0f%%",
            nrow(P), uniqueN(annot$chrom), paste(unique(annot$chrom), collapse = ", "), format(nrow(annot), big.mark = ","),
            base, paste(sizes, collapse = ", "), 100 * tol))

## ---------------- signal retained and reproducibility ----------------
sf <- file.path(BD, "ref_sites.tsv.gz")
keep_cache <- file.exists(sf)                         # signal is kept in memory only when the motif criterion needs it
cache <- list()
curves <- rbindlist(lapply(seq_len(nrow(P)), function(i) {
  ip <- load_signal(P$ip_dir[i], cols = c("all", "h1", "h2")); ct <- load_signal(P$ctl_dir[i], cols = c("all", "h1", "h2"))
  if (keep_cache) cache[[P$sample_id[i]]] <<- list(ip = ip, ct = ct)
  r <- retention_curves(annot, ip, ct, sizes, base)
  r[, ratio := cov / cov[size == floor_w], by = group][, sample_id := P$sample_id[i]][]
}))
summ <- curves[, .(median_ratio = median(ratio, na.rm = TRUE), q25_ratio = quantile(ratio, 0.25, na.rm = TRUE),
                   q75_ratio = quantile(ratio, 0.75, na.rm = TRUE), median_cor = median(cor, na.rm = TRUE),
                   median_bins = median(n_bins)), by = .(size, group)]

## ---------------- reference motif effects (optional) ----------------
ME <- NULL
if (keep_cache) {
  S <- fread(cmd = paste("gzip -dc", shQuote(sf)), sep = "\t")
  if (nrow(S)) {
    fams <- unique(S$family)
    strat_cols <- intersect(c("class", "gc", "mappability"), names(annot))
    ME <- rbindlist(lapply(names(cache), function(s) {
      ip <- cache[[s]]$ip; ct <- cache[[s]]$ct; same_bins(annot, ip)
      base_dt <- cbind(annot[, c("chrom", "bin", "bin_len", strat_cols, if ("n_frac" %in% names(annot)) "n_frac"), with = FALSE],
                       ip = ip$all, inp = ct$all)
      if (!"n_frac" %in% names(base_dt)) base_dt[, n_frac := 0]
      for (f in fams) {                                   # site counts per base bin for each family
        x <- S[family == f, .(chrom, bin, n)]
        base_dt[x, (paste0("s_", f)) := as.numeric(i.n), on = .(chrom, bin)]
        set(base_dt, which(is.na(base_dt[[paste0("s_", f)]])), paste0("s_", f), 0)
      }
      Li <- attr(ip, "library"); Lc <- attr(ct, "library")
      rbindlist(lapply(sizes, function(w) {
        k <- w %/% base; scol <- paste0("s_", fams)
        g <- if (k == 1) copy(base_dt) else {
          m <- aggregate_bins(base_dt, k, c(intersect(c("gc", "mappability"), strat_cols), "n_frac", "ip", "inp"))
          ss <- base_dt[, lapply(.SD, sum), by = .(chrom, bin = bin %/% k), .SDcols = scol]   # sites are counts: summed
          cl <- base_dt[, .(class = class[which.max(bin_len)]), by = .(chrom, bin = bin %/% k)]  # class of the first full base bin
          m <- cbind(m, ss[, ..scol]); if ("class" %in% strat_cols) m[, class := cl$class]; m
        }
        g <- g[n_frac < 0.5 & inp > 0]
        g[, y := binding(ip, inp, Li$fragments, Lc$fragments)]
        sk <- "class"
        for (cc in intersect(c("gc", "mappability"), strat_cols)) {
          g[, (paste0(cc, "_dec")) := cut(get(cc), unique(quantile(get(cc), 0:10 / 10, na.rm = TRUE)), include.lowest = TRUE, labels = FALSE)]
          sk <- c(sk, paste0(cc, "_dec"))
        }
        g[, stratum := do.call(paste, c(.SD, sep = "|")), .SDcols = intersect(sk, names(g))]
        g[, yd := y - mean(y), by = stratum]
        rbindlist(lapply(fams, function(f) {
          x <- g[[paste0("s_", f)]]; xd <- x - ave(x, g$stratum)
          den <- sum(xd^2)
          data.table(size = w, family = f, sample_id = s, beta = if (den > 0) sum(xd * g$yd) / den else NA_real_,
                     n_site_bins = sum(x > 0))
        }))
      }))
    }))
    ME[, ratio := beta / beta[size == floor_w], by = .(family, sample_id)]
  }
}

## ---------------- choose the size ----------------
crit <- summ[, .(ok_signal = all(median_ratio >= 1 - tol, na.rm = TRUE) && !all(is.na(median_ratio))), by = size]
if (!is.null(ME) && nrow(ME)) {
  mm <- ME[is.finite(ratio), .(median_motif_ratio = median(ratio)), by = size]
  crit <- merge(crit, mm, by = "size", all.x = TRUE)[, ok_motif := is.na(median_motif_ratio) | median_motif_ratio >= 1 - tol]
} else crit[, `:=`(median_motif_ratio = NA_real_, ok_motif = TRUE)]
setorder(crit, size)
crit[, ok := ok_signal & ok_motif][, ok_all_below := cumprod(ok) == 1]
chosen <- if (any(crit$ok_all_below)) max(crit[ok_all_below == TRUE, size]) else floor_w
if (!any(crit$ok_all_below)) msg("WARNING: even the floor fails the criterion (no signal at any size?); using the floor")
rat <- summ[size == chosen, setNames(median_ratio, group)]

## ---------------- outputs ----------------
fwrite(curves[, .(sample_id, group, size, split_half_cov = signif(cov, 5), retained_vs_floor = signif(ratio, 4),
                  split_half_cor = signif(cor, 4), n_bins)], file.path(TO, "tuning_curves.tsv"), sep = "\t")
fwrite(merge(summ, crit, by = "size"), file.path(TO, "tuning_summary.tsv"), sep = "\t")
if (!is.null(ME)) fwrite(ME, file.path(TO, "tuning_motif_effects.tsv"), sep = "\t")
rec <- data.table(name = c("chosen_size", "floor_size", "tolerance", "base_size", "sizes", "chromosomes", "drift_chromosome",
                          "n_samples", "samples", "ratio_satellite", "ratio_other", "criterion", "date"),
                  value = c(chosen, floor_w, tol, base, paste(sizes, collapse = ","), paste(unique(annot$chrom), collapse = ","), dchr,
                            nrow(P), paste(P$sample_id, collapse = ","), signif(rat["satellite"], 4), signif(rat["other"], 4),
                            if (!is.null(ME)) "signal retained + reference motif effects" else "signal retained", format(Sys.Date())))
fwrite(rec, file.path(TO, "bin_size.tsv"), sep = "\t", col.names = FALSE)

pd <- copy(curves)[, group := factor(group, c("satellite", "other"), c("CenSat bins", "Other bins"))]
sm <- copy(summ)[, group := factor(group, c("satellite", "other"), c("CenSat bins", "Other bins"))]
p1 <- ggplot(pd, aes(size, ratio, group = sample_id)) + geom_hline(yintercept = 1 - tol, linetype = "dashed", colour = "#d6402a") +
  geom_line(colour = "grey70", linewidth = 0.4) + geom_line(data = sm, aes(size, median_ratio, group = NULL), colour = "black", linewidth = 0.9) +
  geom_vline(xintercept = chosen, colour = "#2a78d6") + scale_x_log10(breaks = sizes) + facet_wrap(~group) +
  labs(x = "Bin size (bp, log scale)", y = "Signal retained vs the floor\n(split-half covariance ratio)",
       title = sprintf("Chosen bin size: %d bp (largest within %.0f%% of %d bp)", chosen, 100 * tol, floor_w)) + theme_x(9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 <- ggplot(pd, aes(size, cor, group = sample_id)) + geom_line(colour = "grey70", linewidth = 0.4) +
  geom_line(data = sm, aes(size, median_cor, group = NULL), colour = "black", linewidth = 0.9) + geom_vline(xintercept = chosen, colour = "#2a78d6") +
  scale_x_log10(breaks = sizes) + facet_wrap(~group) + labs(x = "Bin size (bp, log scale)", y = "Split-half correlation") + theme_x(9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
suppressPackageStartupMessages(library(patchwork))
ggsave(file.path(TO, "tuning_curves.png"), p1 / p2, width = 10, height = 8, dpi = 150, bg = "white")
msg(sprintf("chosen bin size: %d bp (signal retained at that size: satellite %.3f, other %.3f of the %d bp floor) -> %s",
            chosen, rat["satellite"], rat["other"], floor_w, file.path(TO, "bin_size.tsv")))
