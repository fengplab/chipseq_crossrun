#!/usr/bin/env Rscript
## Q3: Which known DNA-binding-protein motifs do the ChIP'd sequences most resemble?
## Two complementary lines of evidence, both against the same database(s) for every sample:
##   (1) Tomtom similarity of de novo motifs (all discovered motifs and family representatives)
##   (2) SEA enrichment of every known motif in every sample's peaks
## Tables:  table04_known_motif_ranking, table04b_known_by_tf_class, table04c_family_known_matches,
##          table04d_sea_known_long
## Figures: fig08_known_motif_dotplot, fig09_tf_class_bar
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "utils.R"))
cfg <- args_config()
mc <- cfg$motifs
ke <- mc$known_evalue %||% 0.05
sq <- mc$sea_qvalue %||% 0.05
S <- samples_dt(cfg)

kidx <- rd(file.path(cfg$work, "known", "known_motif_index.tsv"))
if (is.null(kidx)) { msg("no known motif databases configured; skipping Q3"); quit(save = "no") }
if (!is.null(cfg$tf_class_table) && file.exists(cfg$tf_class_table)) {
  tc <- fread(cfg$tf_class_table)
  kidx <- merge(kidx, tc[, .(motif_id, tf_class, tf_family)], by = "motif_id", all.x = TRUE)
} else kidx[, `:=`(tf_class = NA_character_, tf_family = NA_character_)]
kidx[, tf_name := toupper(sub("^.*\\.", "", sub("\\.H1[0-9].*$", "", alt_id)))]  # HOCOMOCO 'SP1.H12CORE.0.P.B' -> SP1

## (1) Tomtom
td <- read_tomtom(file.path(cfg$work, "known", "tomtom_discovered_vs_known.tsv"))
tr <- read_tomtom(file.path(cfg$work, "known", "tomtom_reps_vs_known.tsv"))
tom <- if (!is.null(td)) td[e_value <= ke, .(sample_id = sub("__.*$", "", query), motif_uid = query, motif_id = target, e_value)] else
  data.table(sample_id = character(), motif_uid = character(), motif_id = character(), e_value = numeric())
tom_s <- if (nrow(tom)) tom[, .(n_samples_tomtom = uniqueN(sample_id), best_tomtom_evalue = min(e_value)), by = motif_id] else
  data.table(motif_id = character(), n_samples_tomtom = integer(), best_tomtom_evalue = numeric())
fm <- if (!is.null(tr)) tr[e_value <= ke][order(e_value)] else NULL
fam_match <- if (!is.null(fm) && nrow(fm)) fm[, .(matched_families = paste(unique(query), collapse = ",")), by = .(motif_id = target)] else NULL
if (!is.null(fm) && nrow(fm)) {
  fk <- merge(fm[, .(family_id = query, motif_id = target, e_value, orientation)], kidx[, .(motif_id, alt_id, tf_class)], by = "motif_id")
  fk <- fk[order(family_id, e_value)][, rank := seq_len(.N), by = family_id][rank <= 5]
  wt(fk, cfg, "table04c_family_known_matches")
}

## (2) SEA
sea_files <- list.files(file.path(cfg$work, "known", "sea"), "\\.sea\\.tsv$", full.names = TRUE)
if (!length(sea_files))
  msg("WARNING: no SEA results in work/known/sea; table04b, table04d and figures 8-9 cannot be made, ",
      "and table04 will only contain Tomtom evidence. Rerun step 5.")
sea <- rbindlist(lapply(sea_files, function(f) {
  l <- readLines(f); l <- l[!grepl("^#", l) & nzchar(l)]
  if (length(l) < 2) return(NULL)
  d <- fread(text = paste(l, collapse = "\n"))
  setnames(d, tolower(gsub("%", "_pct", gsub("[^A-Za-z0-9%]+", "_", names(d)))))
  d[, sample_id := sub("\\.sea\\.tsv$", "", basename(f))]
  d[, .(sample_id, motif_id = id, tp, tp_pct, fp_pct, enr_ratio = as.numeric(enr_ratio),
        pvalue = as.numeric(pvalue), qvalue = as.numeric(qvalue))]
}), fill = TRUE)
if (nrow(sea)) {
  sea <- merge(sea, kidx[, .(motif_id, alt_id, tf_name, tf_class)], by = "motif_id", all.x = TRUE)
  wt(sea, cfg, "table04d_sea_known_long")
  sea_s <- sea[qvalue <= sq, .(n_samples_sea = uniqueN(sample_id), median_enr_ratio = median(enr_ratio),
                               min_sea_q = min(qvalue)), by = motif_id]
} else sea_s <- data.table(motif_id = character(), n_samples_sea = integer(), median_enr_ratio = numeric(), min_sea_q = numeric())

n_s <- uniqueN(S[peaks_file != "", sample_id])
t4 <- merge(kidx, tom_s, by = "motif_id", all.x = TRUE)
t4 <- merge(t4, sea_s, by = "motif_id", all.x = TRUE)
if (!is.null(fam_match)) t4 <- merge(t4, fam_match, by = "motif_id", all.x = TRUE) else t4[, matched_families := NA_character_]
t4[is.na(n_samples_tomtom), n_samples_tomtom := 0L]
t4[is.na(n_samples_sea), n_samples_sea := 0L]
t4 <- t4[n_samples_tomtom > 0 | n_samples_sea > 0]
t4[, `:=`(frac_samples_sea = n_samples_sea / n_s, frac_samples_tomtom = n_samples_tomtom / n_s)]
setorder(t4, -n_samples_sea, -n_samples_tomtom, best_tomtom_evalue, min_sea_q)
t4[, rank := seq_len(.N)]
setcolorder(t4, c("rank", "motif_id", "alt_id", "tf_name", "tf_class", "tf_family", "n_samples_sea", "median_enr_ratio",
                  "min_sea_q", "n_samples_tomtom", "best_tomtom_evalue", "matched_families"))
wt(t4, cfg, "table04_known_motif_ranking")

## collapse redundant motifs: by TF class if available, else by TF name
t4[, grp := fifelse(is.na(tf_class) | tf_class == "NA", tf_name, tf_class)]
if (nrow(sea)) {
  sea[, grp := t4$grp[match(motif_id, t4$motif_id)]]
  bycls <- sea[qvalue <= sq & !is.na(grp), .(n_samples_sea = uniqueN(sample_id), n_motifs = uniqueN(motif_id),
                                             median_enr_ratio = median(enr_ratio)), by = grp]
  tomc <- merge(tom, t4[, .(motif_id, grp)], by = "motif_id")[, .(n_samples_tomtom = uniqueN(sample_id)), by = grp]
  bycls <- merge(bycls, tomc, by = "grp", all = TRUE)
  ## groups seen by only one method get 0 for the other (previously NA, which drew a label with no bar)
  for (v in c("n_samples_sea", "n_motifs", "n_samples_tomtom")) set(bycls, which(is.na(bycls[[v]])), v, 0L)
  bycls[, grouping := fifelse(grp %in% t4[!is.na(tf_class) & tf_class != "NA", tf_class], "TF class", "TF name (no class available)")]
  setnames(bycls, "grp", "tf_class_or_name")
  setorder(bycls, -n_samples_sea, -n_samples_tomtom)
  wt(bycls, cfg, "table04b_known_by_tf_class")
  top_c <- head(bycls, 20)
  pl <- melt(top_c[, .(tf_class_or_name, `SEA enrichment (q < threshold)` = n_samples_sea,
                       `De novo motif matches (Tomtom)` = n_samples_tomtom)],
             id.vars = "tf_class_or_name", variable.name = "evidence", value.name = "n_samples")
  pl[, tf_class_or_name := factor(tf_class_or_name, levels = rev(top_c$tf_class_or_name))]
  pb <- ggplot(pl, aes(tf_class_or_name, n_samples, fill = evidence)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) + coord_flip() +
    scale_fill_manual(values = c("#2a78d6", "#eb6834"), name = NULL) +
    labs(x = NULL, y = "Number of samples", title = "Known-motif evidence by TF class (top 20)") + theme_x()
  save_plot(pb, cfg, "fig09_tf_class_bar", 9, 2 + 0.35 * nrow(top_c))
  fig_note(cfg, "fig09_tf_class_bar", sprintf("Top %d of %d TF classes/names, ranked by samples with SEA enrichment, then Tomtom matches. %s",
           nrow(top_c), nrow(bycls), if (all(bycls$grouping != "TF class")) "No tf_class_table was supplied, so motifs are grouped by TF name." else ""))

  ## dot plot: top TFs x samples
  top <- head(t4[n_samples_sea > 0], 20)$motif_id
  dd <- sea[motif_id %in% top]
  dd[, label := paste0(tf_name, " (", motif_id, ")")]
  dd[, label := factor(label, levels = rev(unique(label[order(match(motif_id, top))])))]
  dd[, denovo_match := paste(sample_id, motif_id) %in% paste(tom$sample_id, tom$motif_id)]
  ordv <- intersect(c(unlist(cfg$categorical_covariates), cfg$time_var), names(S))
  so <- S[sample_id %in% dd$sample_id]
  if (length(ordv)) setorderv(so, ordv)
  dd[, sample := factor(short_id(sample_id), levels = unique(short_id(so$sample_id)))]
  pd <- ggplot(dd, aes(sample, label)) +
    geom_point(aes(size = pmin(-log10(qvalue), 50), colour = log2(enr_ratio))) +
    geom_point(data = dd[denovo_match == TRUE], aes(size = pmin(-log10(qvalue), 50)), shape = 21, colour = "black", stroke = 0.8, fill = NA) +
    scale_colour_gradient(low = "#f7c1c1", high = "#a32d2d", name = "log2 enrichment") +
    scale_size_continuous(name = "-log10 SEA q", range = c(1, 7)) +
    labs(x = NULL, y = NULL, title = "Known motifs enriched in peaks", subtitle = "black ring = a de novo motif in that sample matches (Tomtom)") +
    theme_x(9) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  save_plot(pd, cfg, "fig08_known_motif_dotplot", max(9, 4 + 0.35 * nlevels(dd$sample)), max(5, 2.5 + 0.3 * length(top)))
  fig_note(cfg, "fig08_known_motif_dotplot", sprintf("Top %d of %d known motifs, ranked by the number of samples with SEA enrichment.", length(top), nrow(t4)))
}
