## Drift check (step 16): is the bin size chosen by the tuning run still within tolerance for this batch?
## Recomputes signal retained (split-half covariance of binding) at the tuning floor and at the chosen size on the drift
## chromosome, for every IP/control pair, and compares R(chosen) / R(floor) with the tuning tolerance.
## Args: config.yaml  <drift bins dir at the base size>  <libraries.tsv of the main bin table>  <tuning record>  <chosen size>
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "utils.R"))
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "bin_io.R"))
a <- commandArgs(TRUE)
cfg <- args_config()
base_dir <- a[2]; libs <- fread(a[3]); rec <- fread(a[4], header = FALSE, col.names = c("key", "value"), sep = "\t")
W <- as.integer(a[5]); getr <- function(k) rec[key == k, value][1]
floor_w <- as.integer(getr("floor_size")); tol <- as.numeric(getr("tolerance")); base <- bin_manifest(base_dir)$bin_size[1]
annot <- load_annot(base_dir, cols = c("class"))
P <- lib_pairs(libs[, .(sample_id, library, name, dir = file.path(base_dir, name))])
if (!nrow(P)) { msg("drift check: no IP/control pairs"); quit(status = 0) }
res <- rbindlist(lapply(seq_len(nrow(P)), function(i) {
  ip <- load_signal(P$ip_dir[i], cols = c("h1", "h2")); ct <- load_signal(P$ctl_dir[i], cols = c("h1", "h2"))
  r <- retention_curves(annot, ip, ct, unique(c(floor_w, W)), base)
  r[, ratio := cov / cov[size == floor_w], by = group][size == W][, sample_id := P$sample_id[i]][]
}))
summ <- res[, .(median_ratio = median(ratio, na.rm = TRUE), min_ratio = min(ratio, na.rm = TRUE), n_samples = sum(!is.na(ratio))), by = group]
summ[, `:=`(tuning_ratio = suppressWarnings(as.numeric(getr(paste0("ratio_", group))))), by = group]
summ[, `:=`(chosen_size = W, floor_size = floor_w, tolerance = tol,
            status = fifelse(is.na(median_ratio), "not computable", fifelse(median_ratio >= 1 - tol, "ok", "DRIFT: retune")))]
fwrite(rbind(res[, .(sample_id, group, size, ratio = signif(ratio, 4), cor = signif(cor, 4), n_bins)], fill = TRUE),
       file.path(cfg$work, "bins", "drift_check_samples.tsv"), sep = "\t")
fwrite(summ, file.path(cfg$work, "bins", "drift_check.tsv"), sep = "\t")
for (i in seq_len(nrow(summ)))
  msg(sprintf("%sdrift check, %s bins: R(%d bp) / R(%d bp) = %.3f (median of %d samples; tuning: %s; tolerance %.0f%%)",
              if (summ$status[i] == "DRIFT: retune") "WARNING: " else "", summ$group[i], W, floor_w, summ$median_ratio[i],
              summ$n_samples[i], format(summ$tuning_ratio[i], digits = 3), 100 * tol))
