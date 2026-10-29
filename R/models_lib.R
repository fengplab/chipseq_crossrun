## Shared code for the binding models (steps 18-20): settings, the bin table, the bin filter, binding values,
## motif-count matrices, strata, clusters and the model terms. Source after utils.R and bin_io.R.
suppressPackageStartupMessages(library(data.table))

## ---------------------------------------------------------------- settings (binding_models block, with defaults)
bm_settings <- function(cfg) {
  bm <- cfg$binding_models %||% list()
  s <- list(
    enabled = !isFALSE(bm$enabled),
    terms = bm$terms %||% "auto",                      # auto | list of covariates (e.g. [virus, time_bin, cgas_construct])
    arms = unlist(bm$arms %||% c("all", "uniq")),      # read versions: all reads, uniquely placed reads
    families = bm$families %||% "detected",            # detected | all | list of family ids
    max_families = as.integer(bm$max_families %||% 50),
    count_transform = bm$count_transform %||% "log2p1",  # log2p1 | raw | presence
    gc_bins = as.integer(bm$gc_bins %||% 10),          # GC strata: quantile bins over the filtered bins
    mappability_breaks = as.numeric(unlist(bm$mappability_breaks %||% c(0.25, 0.5, 0.9))),
    min_bins_per_class = as.integer(bm$min_bins_per_class %||% 1000),
    max_n_frac = as.numeric(bm$max_n_frac %||% 0.1),
    input_min_quantile = as.numeric(bm$input_min_quantile %||% 0.01),
    nonarray_block = as.numeric(bm$nonarray_block %||% 1e6),
    min_clusters = as.integer(bm$min_clusters %||% 10),    # fewer clusters in a fit: heteroskedasticity-robust SE instead
    validation_folds = as.integer(bm$validation_folds %||% 5),
    validation_arms = unlist(bm$validation_arms %||% "all"),
    similarity_arm = bm$similarity_arm %||% "all",
    similarity_max_bins = as.numeric(bm$similarity_max_bins %||% 2e6),
    permutations = as.integer(bm$permutations %||% 999),
    level1_mode = bm$level1_mode %||% "slurm",         # slurm (models/level1.sbatch array) | inline (inside step 18)
    threads = as.integer(bm$level1_threads %||% cfg$threads %||% 4),
    seed = as.integer(bm$seed %||% 7))
  bad <- setdiff(s$arms, c("all", "uniq")); if (length(bad)) stop("binding_models.arms: unknown read version(s) ", paste(bad, collapse = ", "), " (use all, uniq)")
  if (!s$count_transform %in% c("log2p1", "raw", "presence")) stop("binding_models.count_transform must be log2p1, raw or presence")
  if (!s$level1_mode %in% c("slurm", "inline")) stop("binding_models.level1_mode must be slurm or inline")
  if (s$validation_folds < 2) stop("binding_models.validation_folds must be at least 2")
  s
}

## the bin table at the configured size (work/config.env BINS_SIZE, written by step 0)
bin_dir <- function(cfg) {
  env <- readLines(file.path(cfg$work, "config.env"), warn = FALSE)
  w <- sub("^export BINS_SIZE=", "", grep("^export BINS_SIZE=", env, value = TRUE))
  w <- gsub("[\"']", "", w)
  if (!length(w) || !nzchar(w)) stop("no bin size: set bins.size or bins.tuning_record in the config and rerun step 0")
  d <- file.path(cfg$work, "bins", paste0("w", w))
  if (!file.exists(file.path(d, "manifest.tsv"))) stop("bin table ", d, " not found: run step 16 at this bin size first")
  d
}

## ---------------------------------------------------------------- design: samples and model terms
## Terms: 'auto' = every level step 17 marks estimable; a list of covariates = the estimable levels of those covariates.
## A covariate with no estimable level (all aliased or absent) is refused with step 17's reason.
model_terms <- function(cfg, s) {
  MD <- file.path(cfg$work, "models")
  et <- fread(file.path(MD, "estimable_terms.tsv"))
  est <- fread(file.path(cfg$tables, "design_estimability.tsv"))
  if (identical(s$terms, "auto") || identical(unlist(s$terms), "auto")) return(et)
  want <- unlist(s$terms)
  unknown <- setdiff(want, est$covariate)
  if (length(unknown)) stop("binding_models.terms: not a covariate of the design: ", paste(unknown, collapse = ", "),
                            " (covariates: ", paste(unique(est$covariate), collapse = ", "), ")")
  none <- setdiff(want, et$covariate)
  if (length(none)) {
    why <- est[covariate %in% none & status != "reference level", paste0(covariate, " ", level, ": ", status,
                 ifelse(is.na(aliased_with), "", paste0(" (", aliased_with, ")")))]
    stop("binding_models.terms: no estimable level for ", paste(none, collapse = ", "), ". Step 17 says: ", paste(why, collapse = "; "))
  }
  et[covariate %in% want]
}
design_matrix <- function(cfg) fread(file.path(cfg$work, "models", "design_matrix.tsv"))

## IP/control library directories of each modelled sample
model_libraries <- function(cfg, B) {
  libs <- fread(file.path(B, "libraries.tsv"))
  P <- lib_pairs(libs[, .(sample_id, library, name, dir = file.path(B, name))])
  dm <- design_matrix(cfg)
  P[sample_id %in% dm$sample_id]
}

## ---------------------------------------------------------------- bins, filter, binding
## annotation + filter flags shared by every sample (blacklist / include, N fraction)
load_model_annot <- function(cfg, B, s) {
  A <- load_annot(B, cols = c("class", "array_id", "edge_dist", "mappability", "gc", "n_frac"))
  m <- bin_manifest(B); bs <- m$bin_size[1]
  A[, `:=`(start = bin * bs, end = pmin(bin * bs + bin_len, bin * bs + bs))]
  A[, keep_seq := !is.na(n_frac) & n_frac <= s$max_n_frac & !is.na(gc) & !is.na(mappability)]
  ## include / exclude regions (BED, 0-based half-open; a bin is in a region when its midpoint is: mid = start + len/2)
  inreg <- function(bed) {
    r <- fread(bed, header = FALSE, select = 1:3, col.names = c("chrom", "s", "e"), colClasses = list(character = 1))
    mid <- A[, .(chrom, x = start + bin_len / 2)][, i := .I]
    setkey(r, chrom, s, e)
    hit <- foverlaps(mid[, .(chrom, s = x, e = x, i)], r, type = "within", nomatch = 0L)
    seq_len(nrow(A)) %in% hit$i
  }
  ann <- file.path(cfg$work, "annot")
  if (file.exists(f <- file.path(ann, "include.bed")) && file.size(f) > 0) A[, keep_seq := keep_seq & inreg(f)]
  if (file.exists(f <- file.path(ann, "exclude.bed")) && file.size(f) > 0) A[, keep_seq := keep_seq & !inreg(f)]
  ## clusters: a satellite array, or a fixed-size block of the chromosome outside arrays
  A[, cluster := fifelse(!is.na(array_id) & array_id >= 0, paste0("a", array_id), paste0(chrom, ":", start %/% s$nonarray_block))]
  A[, class_name := fifelse(is.na(class_name), "none", class_name)]
  A[, is_censat := !is.na(class_source) & class_source == "censat"]
  A
}

## binding of one sample in one read version, with the input filter; NULL when the version has no reads
sample_binding <- function(P_i, arm, s) {
  ip <- load_signal(P_i$ip_dir, cols = arm); ct <- load_signal(P_i$ctl_dir, cols = arm)
  Li <- attr(ip, "library"); Lc <- attr(ct, "library")
  key <- if (arm == "all") "fragments" else "unique"
  n_ip <- as.numeric(Li[[key]]); n_ct <- as.numeric(Lc[[key]])
  if (!length(n_ip) || !length(n_ct) || !is.finite(n_ip) || !is.finite(n_ct) || n_ip <= 0 || n_ct <= 0) return(NULL)
  a <- ip[[arm]]; b <- ct[[arm]]
  y <- binding(a, b, n_ip, n_ct)
  bn <- b / (n_ct / 1e6)                                 # input per million fragments
  thr <- as.numeric(quantile(bn[bn > 0], s$input_min_quantile, na.rm = TRUE))
  list(y = y, input_ok = !is.na(bn) & bn > 0 & bn >= thr, n_ip = n_ip, n_input = n_ct, input_threshold = thr)
}

## ---------------------------------------------------------------- motif counts
## families entering the models: detected in >= 1 sample (crossrun table02), capped at max_families by samples detected,
## then median % of peaks with a site
model_families <- function(cfg, s, B) {
  sites_fam <- unique(fread(cmd = sprintf("gzip -dc '%s'", file.path(B, "sites.tsv.gz")), select = "family")$family)
  if (is.character(s$families) && length(s$families) > 1 || !s$families %in% c("detected", "all")) {
    fam <- unlist(s$families); miss <- setdiff(fam, sites_fam)
    if (length(miss)) stop("binding_models.families: no genome-wide sites in the bin table for ", paste(miss, collapse = ", "))
    return(fam)
  }
  if (s$families == "all") return(sort(sites_fam))
  pv <- fread(file.path(cfg$tables, "table02_motif_prevalence.tsv"))
  pv[, det := toupper(as.character(detected)) == "TRUE"]
  st <- pv[, .(n_det = sum(det), med = median(pct_peaks, na.rm = TRUE)), by = family_id][n_det > 0 & family_id %in% sites_fam]
  setorder(st, -n_det, -med, family_id)
  if (!nrow(st)) stop("no motif family is detected in any sample (table02_motif_prevalence) and has genome-wide sites")
  head(st$family_id, s$max_families)
}

## bins x families matrix of transformed site counts, rows in annotation order
motif_matrix <- function(B, A, fam, s) {
  st <- fread(cmd = sprintf("gzip -dc '%s'", file.path(B, "sites.tsv.gz")), colClasses = list(character = c("family", "chrom")))
  st <- st[family %in% fam]
  idx <- A[, .(chrom, bin, row = .I)]
  st <- merge(st, idx, by = c("chrom", "bin"))
  X <- matrix(0, nrow(A), length(fam), dimnames = list(NULL, fam))
  X[cbind(st$row, match(st$family, fam))] <- switch(s$count_transform, log2p1 = log2(1 + st$n), raw = st$n, presence = 1)
  X
}

## strata: GC quantile bin x mappability bin (computed over the bins entering a fit)
strata_cols <- function(gc, mp, s) {
  qb <- unique(quantile(gc, seq(0, 1, length.out = s$gc_bins + 1), na.rm = TRUE))
  gcb <- if (length(qb) > 2) cut(gc, qb, include.lowest = TRUE, labels = FALSE) else rep(1L, length(gc))
  mb <- findInterval(mp, s$mappability_breaks, rightmost.closed = TRUE)
  paste0("g", gcb, "m", mb)
}

need_pkgs <- function(p) {
  miss <- p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("R package(s) not installed: ", paste(miss, collapse = ", "), " (see environment.yml)")
}
