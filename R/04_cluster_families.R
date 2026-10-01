#!/usr/bin/env Rscript
## Step 04b: cluster every de novo motif (all samples) into cross-run motif families.
## Distance between two motifs: d = min(1, E / family_evalue) using the better of the two Tomtom
## directions; pairs Tomtom did not report get d = 1. Average-linkage hclust, cut at family_cut_height.
## Output: work/motifs/families.tsv (motif_uid, sample_id, family_id, is_representative, ...)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "utils.R"))
cfg <- args_config()
thr <- cfg$motifs$family_evalue %||% 0.05
cut_h <- cfg$motifs$family_cut_height %||% 0.5

ids <- sub("^MOTIF\\s+(\\S+).*", "\\1", grep("^MOTIF", readLines(file.path(cfg$work, "motifs", "all_discovered.meme")), value = TRUE))
if (length(ids) == 0) stop("no de novo motifs found in any XSTREME output")
tt <- read_tomtom(file.path(cfg$work, "motifs", "tomtom_all_vs_all.tsv"))

n <- length(ids)
D <- matrix(1, n, n, dimnames = list(ids, ids))
diag(D) <- 0
if (!is.null(tt) && nrow(tt)) {
  tt <- tt[query %in% ids & target %in% ids & query != target]
  tt[, d := pmin(1, e_value / thr)]
  tt <- tt[, .(d = min(d)), by = .(a = pmin(query, target), b = pmax(query, target))]
  D[cbind(tt$a, tt$b)] <- tt$d
  D[cbind(tt$b, tt$a)] <- tt$d
}
cl <- if (n > 1) {
  hc <- hclust(as.dist(D), method = "average")
  hc$height <- cummax(round(hc$height, 12))  # guard against floating-point non-monotonicity with tied distances
  cutree(hc, h = cut_h)
} else setNames(1L, ids)

xm <- rd(file.path(cfg$tables, "xstreme_motifs.tsv"))
fam <- data.table(motif_uid = ids, cluster = cl[ids])
fam[, sample_id := sub("__.*$", "", motif_uid)]
fam[, motif_id := sub("^.*?__", "", motif_uid)]
if (!is.null(xm)) {
  xm <- unique(xm[, .(motif_uid, consensus, width, source, sea_pvalue = suppressWarnings(as.numeric(sea_pvalue)),
                      evalue = suppressWarnings(as.numeric(evalue)))], by = "motif_uid")
  fam <- merge(fam, xm, by = "motif_uid", all.x = TRUE)
} else {
  fam[, `:=`(consensus = NA_character_, width = NA_integer_, source = NA_character_, sea_pvalue = NA_real_, evalue = NA_real_)]
}

## representative: best SEA p-value (then E-value); family order: samples covered, then best p
fam[, rank_p := fifelse(is.na(sea_pvalue), fifelse(is.na(evalue), 1, evalue), sea_pvalue)]
setorder(fam, cluster, rank_p, -width)
fam[, is_representative := seq_len(.N) == 1L, by = cluster]
ord <- fam[, .(n_samples = uniqueN(sample_id), best = min(rank_p)), by = cluster][order(-n_samples, best)]
ord[, family_id := sprintf("F%03d", seq_len(.N))]
fam <- merge(fam, ord[, .(cluster, family_id, family_n_samples = n_samples)], by = "cluster")
fam[, family_size := .N, by = family_id]
setorder(fam, family_id, -is_representative, rank_p)
fam[, c("cluster", "rank_p") := NULL]
fwrite(fam, file.path(cfg$work, "motifs", "families.tsv"), sep = "\t", na = "NA", quote = FALSE)
msg(sprintf("%d motifs -> %d families", n, uniqueN(fam$family_id)))
