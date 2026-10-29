## Binding models, Level 1: one sample, one read version.
##   binding (log2 IP/input per bin) ~ motif counts of every model family | stratum, standard errors clustered by
##   satellite array (blocks of binding_models.nonarray_block bp outside arrays).
## Fits written (work/models/level1/<sample>__<arm>.coef.tsv):
##   per_class          one fit per repeat class (CenSat class, RepeatMasker class or 'none'), strata GC x mappability
##   pooled_adjusted    all classes in one fit, strata repeat class x GC x mappability (bins compared within class)
##   pooled_unadjusted  all classes in one fit, strata GC x mappability only (the before-and-after comparison)
## Held-out arrays (arm in binding_models.validation_arms), <sample>__<arm>.valid.tsv: in each CenSat class with at least
##   validation_folds arrays, arrays are split into folds; the per-class model fitted on the other folds predicts the
##   held-out arrays, with all families and with each family left out. delta_mse = MSE without the family - MSE with it
##   (positive: the family improves prediction on arrays the model never saw). No threshold is applied.
## Usage: Rscript R/70_level1.R config.yaml <sample_id> <arm>      (arm: all | uniq)
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "utils.R")); source(file.path(here, "bin_io.R")); source(file.path(here, "models_lib.R"))
a <- commandArgs(TRUE); cfg <- args_config()
if (length(a) < 3) stop("usage: Rscript R/70_level1.R config.yaml <sample_id> <arm>")
sid <- a[2]; arm <- a[3]
s <- bm_settings(cfg); need_pkgs("fixest")
fixest::setFixest_nthreads(s$threads); fixest::setFixest_notes(FALSE)
B <- bin_dir(cfg); OD <- file.path(cfg$work, "models", "level1"); dir.create(OD, showWarnings = FALSE, recursive = TRUE)
pref <- file.path(OD, paste0(sid, "__", arm))
P <- model_libraries(cfg, B)[sample_id == sid]
if (!nrow(P)) stop(sid, ": not a modelled sample (needs IP and control libraries in ", B, " and a row in the design matrix)")
set.seed(s$seed)

A <- load_model_annot(cfg, B, s)
fam <- model_families(cfg, s, B)
b <- sample_binding(P, arm, s)
info <- data.table(sample_id = sid, arm = arm, bin_size = bin_manifest(B)$bin_size[1], bins_total = nrow(A),
                   bins_seq_ok = sum(A$keep_seq), n_families = length(fam), count_transform = s$count_transform)
if (is.null(b)) {
  msg(sid, " ", arm, ": no reads in this read version (library size 0); skipped")
  fwrite(info[, `:=`(status = "skipped: no reads in this read version", bins_used = 0L)], paste0(pref, ".info.tsv"), sep = "\t")
  for (x in c("coef", "valid")) fwrite(data.table(), paste0(pref, ".", x, ".tsv"), sep = "\t")
  quit(status = 0)
}
keep <- A$keep_seq & b$input_ok & is.finite(b$y)
info[, `:=`(bins_input_ok = sum(A$keep_seq & b$input_ok), bins_used = sum(keep), input_threshold_rpm = signif(b$input_threshold, 4),
            library_ip = b$n_ip, library_input = b$n_input)]
msg(sprintf("%s %s: %d bins of %d bp; %d pass the sequence filter, %d also pass the input filter (>= %.3g per million), %d used; %d families",
            sid, arm, nrow(A), info$bin_size, info$bins_seq_ok, info$bins_input_ok, b$input_threshold, info$bins_used, length(fam)))
if (sum(keep) < s$min_bins_per_class) stop(sid, " ", arm, ": only ", sum(keep), " bins pass the filters")

X <- motif_matrix(B, A, fam, s)[keep, , drop = FALSE]
vn <- make.names(fam, unique = TRUE); colnames(X) <- vn; v2f <- setNames(fam, vn)
d <- data.table(y = b$y[keep], class = A$class_name[keep], censat = A$is_censat[keep], cluster = A$cluster[keep],
                array = A$array_id[keep])
d[, stratum := strata_cols(A$gc[keep], A$mappability[keep], s)]
d <- cbind(d, as.data.table(X)); rm(X); invisible(gc())

## one fit; families absent from the fit (no variation within strata) are reported as not estimable
fit_one <- function(dd, fe, model, cls) {
  present <- vn[vapply(vn, function(v) any(dd[[v]] != 0), logical(1))]
  out <- data.table(sample_id = sid, arm = arm, model = model, class = cls, family_id = v2f[vn], estimate = NA_real_,
                    se = NA_real_, t = NA_real_, p = NA_real_, n_bins = nrow(dd), n_clusters = uniqueN(dd$cluster),
                    se_type = NA_character_, status = "not estimable in this class: no sites")
  if (!length(present)) return(out)
  ncl <- uniqueN(dd$cluster)
  vc <- if (ncl >= s$min_clusters) ~cluster else "hetero"   # too few clusters for a clustered SE: heteroskedasticity-robust
  f <- stats::as.formula(paste("y ~", paste(present, collapse = " + "), "|", fe))
  fit <- tryCatch(fixest::feols(f, data = dd, vcov = vc, lean = TRUE, warn = FALSE),
                  error = function(e) { msg("WARNING: ", sid, " ", arm, " ", model, " ", cls, ": ", conditionMessage(e)); NULL })
  if (is.null(fit)) { out[family_id %in% v2f[present], status := "fit failed"]; return(out) }
  ct <- fixest::coeftable(fit)
  got <- intersect(rownames(ct), present)
  i <- match(v2f[got], out$family_id)
  out[i, `:=`(estimate = ct[got, 1], se = ct[got, 2], t = ct[got, 3], p = ct[got, 4],
              se_type = if (is.character(vc)) paste0("heteroskedasticity-robust (fewer than ", s$min_clusters, " clusters)") else "clustered by array / block",
              status = "estimated")]
  out[family_id %in% v2f[setdiff(present, got)], status := "not estimable in this class: collinear with strata or other families"]
  out
}
cls_n <- d[, .N, by = class][N >= s$min_bins_per_class]
res <- list()
for (cl in cls_n$class) res[[length(res) + 1]] <- fit_one(d[class == cl], "stratum", "per_class", cl)
small <- d[, .N, by = class][N < s$min_bins_per_class]
if (nrow(small)) msg(sid, " ", arm, ": classes with fewer than ", s$min_bins_per_class, " bins get no per-class fit: ",
                     paste(sprintf("%s (%d)", small$class, small$N), collapse = ", "))
d[, cls_stratum := paste(class, stratum)]
res[[length(res) + 1]] <- fit_one(d, "cls_stratum", "pooled_adjusted", "ALL")
res[[length(res) + 1]] <- fit_one(d, "stratum", "pooled_unadjusted", "ALL")
coef <- rbindlist(res)
fwrite(coef, paste0(pref, ".coef.tsv"), sep = "\t")

## ---------------------------------------------------------------- held-out arrays
val <- data.table()
if (arm %in% s$validation_arms) {
  folds <- s$validation_folds
  for (cl in cls_n$class) {
    dc <- d[class == cl & censat == TRUE & array >= 0]
    arrs <- unique(dc$array)
    if (length(arrs) < folds) next
    present <- vn[vapply(vn, function(v) any(dc[[v]] != 0), logical(1))]
    if (!length(present)) next
    fid <- setNames(sample(rep_len(seq_len(folds), length(arrs))), arrs)
    dc[, fold := fid[as.character(array)]]
    mse <- function(tr, te, vars) {
      if (!length(vars)) return(NA_real_)
      f <- stats::as.formula(paste("y ~", paste(vars, collapse = " + "), "| stratum"))
      fit <- tryCatch(fixest::feols(f, data = tr, warn = FALSE, notes = FALSE), error = function(e) NULL)
      if (is.null(fit)) return(NA_real_)
      pr <- predict(fit, newdata = te)
      ok <- is.finite(pr); if (sum(ok) < 10) return(NA_real_)
      mean((te$y[ok] - pr[ok])^2)
    }
    per <- rbindlist(lapply(seq_len(folds), function(k) {
      tr <- dc[fold != k]; te <- dc[fold == k]
      full <- mse(tr, te, present)
      data.table(fold = k, family = present, mse_full = full,
                 mse_drop = vapply(present, function(v) mse(tr, te, setdiff(present, v)), numeric(1)))
    }))
    per[, delta := mse_drop - mse_full]
    val <- rbind(val, per[, .(sample_id = sid, arm = arm, class = cl, family_id = v2f[family], n_arrays = length(arrs),
                              folds_used = sum(is.finite(delta)), mse_full = mean(mse_full, na.rm = TRUE),
                              delta_mse = mean(delta, na.rm = TRUE), delta_mse_sd = sd(delta, na.rm = TRUE)), by = family])
  }
  if (nrow(val)) val[, family := NULL]
  msg(sid, " ", arm, ": held-out arrays in ", uniqueN(val$class), " CenSat class(es) with >= ", folds, " arrays")
}
fwrite(val, paste0(pref, ".valid.tsv"), sep = "\t")
info[, status := "done"]
fwrite(info, paste0(pref, ".info.tsv"), sep = "\t")
msg(sprintf("%s %s: %d per-class fits; %d family x class effects estimated", sid, arm, nrow(cls_n), coef[model == "per_class" & status == "estimated", .N]))
