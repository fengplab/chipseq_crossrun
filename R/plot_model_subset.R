#!/usr/bin/env Rscript
## Custom versions of Figure 18 (peak-level models) and Figure 19 (sample-level models) for chosen covariates,
## comparisons, motif families or responses, drawn from the tables a finished run wrote (step 17). Nothing is refitted.
##
## Usage:
##   Rscript R/plot_model_subset.R <crossrun_outdir> --figure 18|19 [filters] [options]
##
## Filters (combine freely):
##   --covariate virus[,time_bin]   only the comparisons of these covariates (name or label as in table09c/table11,
##                                  e.g. virus, cell_type, time_bin or "hpi"). For Figure 18 this switches to the
##                                  covariate view: how each family's peak enrichment changes with the covariate
##                                  (the "covariate x peak" terms), instead of the repeat-class panels
##   --levels "VACV vs HSV-1"       only comparisons whose label contains one of these texts (comma-separated)
##   --families F001,F015           only these motif families (Figure 19: responses about these families)
##   --types density,reads_class    Figure 19 only: response types (density, density_class, tandem, peak_share,
##                                  reads_class, reads_peaks); default all
##   --responses REGEX              Figure 19 only: responses whose name matches the regular expression
##   --with-classes                 Figure 18 with --covariate: also show the sequence and binding panels
##
## Which rows are shown:
##   default        for every column, all rows significant there, or its 5 smallest p-values if fewer are
##                  (the rule of the report figures); with --families, every requested family
##   --sig-only     only rows significant in at least one shown column
##   --top N        at most N rows (smallest p-value first)
##   --alpha A      significance threshold (default: models.alpha from the run's config, else 0.05)
##
## Output:
##   --out FILE     .png/.pdf/.svg (default <outdir>/figures/custom/fig18_<filters>.png); a TSV with the plotted
##                  values is written next to it. Long figures are split into parts of --max-rows rows (default 60)
##   --title TEXT, --width W, --height H (inches)
##
## Stars follow the report: Figure 18 covariate view = the covariate's joint test for the family (table09c) has
## q < alpha and this comparison has p < 0.05; Figure 19 = the covariate's F test has q < alpha and this level
## differs from the reference (p < 0.05).
##
## Examples:
##   Rscript R/plot_model_subset.R crossrun_out_human_broad --figure 19 --covariate virus --sig-only
##   Rscript R/plot_model_subset.R crossrun_out_human_broad --figure 18 --covariate virus
##   Rscript R/plot_model_subset.R crossrun_out_human_broad --figure 18 --families F001,F015,F027
##   Rscript R/plot_model_subset.R crossrun_out_human_broad --figure 19 --types density --levels "HSV-1 vs mock"
self <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(dirname(self), "utils.R"))
suppressPackageStartupMessages(library(patchwork))

a <- commandArgs(trailingOnly = TRUE)
if (!length(a) || any(a %in% c("-h", "--help"))) {
  cat(sub("^## ?", "", grep("^##", readLines(self), value = TRUE)), sep = "\n"); quit(status = if (length(a)) 0 else 1)
}
die <- function(...) { message("ERROR: ", ...); quit(status = 1) }
say <- function(...) message("[subset] ", ...)
outdir <- a[1]; a <- a[-1]
known <- c("--figure", "--covariate", "--levels", "--families", "--types", "--responses", "--with-classes", "--sig-only",
           "--top", "--alpha", "--out", "--title", "--width", "--height", "--max-rows")
bad <- setdiff(grep("^--", a, value = TRUE), known)
if (length(bad)) die("unknown option(s): ", paste(bad, collapse = ", "), " (see --help)")
opt <- function(name, default = NULL) { i <- match(name, a); if (is.na(i)) default else a[i + 1] }
has <- function(name) name %in% a
csv <- function(x) if (is.null(x)) character() else trimws(strsplit(x, ",")[[1]])

fig <- opt("--figure") %||% die("--figure 18 or --figure 19 is required")
if (!fig %in% c("18", "19")) die("--figure must be 18 or 19")
tabs <- file.path(outdir, "tables")
if (!dir.exists(tabs)) die("no tables folder in ", outdir)
alpha <- as.numeric(opt("--alpha", NA))
if (is.na(alpha)) {                     # the run's own setting: config path recorded in work/config.env
  env <- file.path(outdir, "work", "config.env")
  cf <- if (file.exists(env)) sub("^export CONFIG_YAML=", "", grep("^export CONFIG_YAML=", readLines(env), value = TRUE)[1]) else NA
  cy <- if (!is.na(cf) && file.exists(cf)) tryCatch(yaml::read_yaml(cf), error = function(e) NULL) else NULL
  alpha <- as.numeric(cy$models$alpha %||% 0.05)
}
covs_req <- csv(opt("--covariate")); lev_req <- csv(opt("--levels")); fam_req <- csv(opt("--families"))
top_n <- as.integer(opt("--top", NA)); max_rows <- as.integer(opt("--max-rows", 60))
need_tab <- function(name) {
  f <- file.path(tabs, paste0(name, ".tsv"))
  x <- rd(f); if (is.null(x) || !nrow(x)) die(name, ".tsv is missing or empty in ", tabs, " (run step 17 first)")
  x
}
## covariate names or labels -> covariate names, from a table with covariate and covariate_label columns
resolve_cov <- function(req, tab) {
  if (!length(req)) return(character())
  key <- unique(tab[, .(covariate, covariate_label)])
  out <- vapply(req, function(r) {
    hit <- key[tolower(covariate) == tolower(r) | tolower(covariate_label) == tolower(r)]
    if (!nrow(hit)) hit <- key[startsWith(tolower(covariate_label), tolower(r)) | startsWith(tolower(covariate), tolower(r))]
    if (!nrow(hit)) die("covariate '", r, "' not found; available: ", paste(sprintf("%s (%s)", key$covariate, key$covariate_label), collapse = ", "))
    hit$covariate[1]
  }, "")
  unique(out)
}
match_levels <- function(lbl) if (!length(lev_req)) rep(TRUE, length(lbl)) else Reduce(`|`, lapply(lev_req, function(l) grepl(l, lbl, fixed = TRUE)))
## rows to show: every requested family, or the report rule, then --sig-only and --top
pick_rows <- function(d, row) {
  rows <- if (length(fam_req) && fig == "18") intersect(fam_req, d[[row]]) else select_heatmap_rows(d, row, "col_key", "p_value", "is_sig", 5)
  if (has("--sig-only")) rows <- intersect(rows, unique(d[is_sig == TRUE][[row]]))
  if (!is.na(top_n) && length(rows) > top_n) {
    best <- d[get(row) %in% rows, .(p = min(p_value, na.rm = TRUE)), by = row][order(p)]
    rows <- head(best[[row]], top_n)
  }
  rows
}
slug <- function() {
  parts <- c(if (length(covs_req)) paste(covs_req, collapse = "-"), if (length(lev_req)) "levels",
             if (length(fam_req)) paste(head(fam_req, 4), collapse = "-"), if (fig == "19" && !is.null(opt("--types"))) gsub(",", "-", opt("--types")),
             if (has("--sig-only")) "sig")
  s <- if (length(parts)) paste(parts, collapse = "_") else "all"
  gsub("[^A-Za-z0-9_.-]+", "-", s)
}
out <- opt("--out", file.path(outdir, "figures", "custom", sprintf("fig%s_%s.png", fig, slug())))
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
part_name <- function(i, n) if (n > 1) sub("(\\.[A-Za-z]+)$", sprintf("_p%0*d\\1", nchar(n), i), out) else out
heat <- function(dp, xcol, ycol, legend, lim, title, subtitle, drop_y = FALSE) {
  ggplot(dp, aes(.data[[xcol]], .data[[ycol]])) +
    geom_tile(aes(fill = pmax(pmin(value, lim), -lim)), colour = "grey85") +
    geom_point(data = dp[is_sig == TRUE], shape = 8, size = 1.6) +
    scale_fill_gradient2(low = "#2a78d6", mid = "white", high = "#d6402a", midpoint = 0, limits = c(-lim, lim), name = legend) +
    scale_y_discrete(drop = drop_y) +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) + theme_x(9) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), panel.grid = element_blank(), legend.key.width = grid::unit(1.2, "cm"),
          strip.text.y = element_text(angle = 0, hjust = 0))
}
save_parts <- function(rows, make_plot, w, h_per_row, base_h) {
  pages <- split(rows, ceiling(seq_along(rows) / max_rows))
  for (i in seq_along(pages)) {
    p <- make_plot(pages[[i]], if (length(pages) > 1) sprintf("Part %d of %d", i, length(pages)) else NULL)
    f <- part_name(i, length(pages))
    ggsave(f, p, width = as.numeric(opt("--width", w)), height = as.numeric(opt("--height", base_h + h_per_row * length(pages[[i]]))),
           dpi = 200, bg = "white", limitsize = FALSE)
    say("wrote ", f)
  }
}

## =============================================================== Figure 18
if (fig == "18") {
  CO <- need_tab("table09b_peak_model_coefficients"); T9 <- need_tab("table09_motif_classification")
  fam_lab <- function(f) {
    cons <- T9$pwm_consensus[match(f, T9$family_id)]; cls <- T9$classification[match(f, T9$family_id)]
    paste0(f, " ", ifelse(is.na(cons), "", cons), "  [", cls, "]")
  }
  if (length(fam_req)) {
    miss <- setdiff(fam_req, T9$family_id)
    if (length(miss)) die("families not in table09: ", paste(miss, collapse = ", "))
    nofit <- setdiff(fam_req, CO$family_id)
    if (length(nofit)) say("no model coefficients (", paste(sprintf("%s: %s", nofit, T9$classification[match(nofit, T9$family_id)]), collapse = "; "),
                           ") - see table09 for the reason")
  }
  pieces <- list()
  covariate_view <- length(covs_req) > 0 || length(lev_req) > 0
  if (covariate_view) {
    CV <- need_tab("table09c_peak_model_covariates")
    cc <- CO[term_group == "covariate x peak"]
    if (!"covariate" %in% names(cc)) {          # tables from before the covariate column existed: recover it from the label
      cl <- unique(CV[, .(covariate, covariate_label)])
      cc[, covariate := NA_character_]
      for (i in seq_len(nrow(cl))) cc[startsWith(term_label, paste0("Peak enrichment, ", cl$covariate_label[i], ":")), covariate := cl$covariate[i]]
    }
    covs <- if (length(covs_req)) resolve_cov(covs_req, CV) else unique(CV$covariate)
    cc <- cc[covariate %in% covs]
    cc[, col := sub("^Peak enrichment, ", "", term_label)]
    cc <- cc[match_levels(col)]
    if (!nrow(cc)) die("no covariate x peak terms for ", paste(covs, collapse = ", "), if (length(lev_req)) " with those levels")
    cc <- merge(cc, CV[, .(family_id, covariate, cov_q = q_value)], by = c("family_id", "covariate"), all.x = TRUE)
    cc[, `:=`(is_sig = !is.na(cov_q) & cov_q < alpha & !is.na(p_value) & p_value < 0.05,
              panel = "Change in peak enrichment\n(vs reference level)", value = estimate)]
    pieces$cov <- cc
  }
  if (!covariate_view || has("--with-classes")) {
    pd <- CO[term_group %in% c("peak", "repeat class", "repeat class x peak")]
    pd <- merge(pd, T9[, .(family_id, repeat_q, interaction_q)], by = "family_id", all.x = TRUE)
    pd[, is_sig := !is.na(q_value) & q_value < alpha &
         (term_group == "peak" | (term_group == "repeat class" & !is.na(repeat_q) & repeat_q < alpha) |
          (term_group == "repeat class x peak" & !is.na(interaction_q) & interaction_q < alpha))]
    pd[, col := fifelse(term_group == "peak", "In peaks (outside repeats)", stringr_wrap(repeat_class_name %||% repeat_class, 22))]
    pd[, panel := fifelse(term_group == "repeat class", "Sequence effect:\nmotif frequency in the class", "Binding effect:\npeak vs background")]
    pd[, value := estimate]
    pieces$cls <- pd
  }
  D <- rbindlist(lapply(pieces, function(x) x[, .(family_id, col, panel, value, p_value, is_sig)]))
  if (length(fam_req)) D <- D[family_id %in% fam_req]
  if (!nrow(D)) die("nothing to plot after the filters")
  D[, col_key := paste(panel, col)]
  rows <- pick_rows(D, "family_id")
  if (!length(rows)) die("no family passes the filters", if (has("--sig-only")) sprintf(" (none significant at alpha = %g)", alpha))
  cls <- T9$classification[match(rows, T9$family_id)]
  if (!length(fam_req)) rows <- rows[order(cls, match(rows, rows))]
  D <- D[family_id %in% rows]
  panels <- c("Sequence effect:\nmotif frequency in the class", "Binding effect:\npeak vs background", "Change in peak enrichment\n(vs reference level)")
  D[, panel := factor(panel, intersect(panels, unique(panel)))]
  cols_order <- unique(c("In peaks (outside repeats)", D[order(panel), col]))
  D[, col := factor(col, intersect(cols_order, unique(col)))]
  lim <- max(2, min(8, max(abs(D$value), na.rm = TRUE)))
  cov_title <- if (covariate_view) paste(unique(CV$covariate_label[match(unique(pieces$cov$covariate), CV$covariate)]), collapse = ", ") else ""
  title <- opt("--title", if (covariate_view) sprintf("Peak-level models: change in peak enrichment with %s", cov_title)
                          else "Peak-level mixed logistic models: what explains motif presence")
  fwrite(D[, .(family_id, label = fam_lab(family_id), panel = gsub("\n", " ", panel), column = gsub("\n", " ", col), log_odds_ratio = value,
               odds_ratio = exp(value), p_value, significant = is_sig)], sub("\\.[A-Za-z]+$", ".tsv", out), sep = "\t", quote = FALSE)
  save_parts(rows, function(pr, subt) {
    dp <- D[family_id %in% pr]; dp[, fam := factor(fam_lab(family_id), levels = rev(fam_lab(pr)))]
    heat(dp, "col", "fam", "log odds ratio", lim, title,
         paste(c(subt, sprintf("star: %s", if (covariate_view) sprintf("covariate joint test q < %g and comparison p < 0.05", alpha)
                                            else sprintf("q < %g (and the family's joint test for repeat-class terms)", alpha))), collapse = "; ")) +
      facet_grid(. ~ panel, scales = "free_x", space = "free_x")
  }, max(7, 4.5 + 0.6 * uniqueN(D$col)), 0.3, 2.8)
}

## =============================================================== Figure 19
if (fig == "19") {
  C11 <- need_tab("table11b_covariate_coefficients"); T11 <- need_tab("table11_covariate_tests")
  covs <- if (length(covs_req)) resolve_cov(covs_req, T11) else unique(C11$covariate)
  cc <- C11[covariate %in% covs & !is.na(estimate)]
  ty <- csv(opt("--types"))
  if (length(ty)) { miss <- setdiff(ty, unique(C11$feature_type)); if (length(miss)) die("unknown --types: ", paste(miss, collapse = ", "),
                                                                                         "; available: ", paste(unique(C11$feature_type), collapse = ", "))
                    cc <- cc[feature_type %in% ty] }
  if (!is.null(opt("--responses"))) cc <- cc[grepl(opt("--responses"), response, perl = TRUE)]
  if (length(fam_req)) cc <- cc[response_id %in% fam_req | grepl(paste0("\\b(", paste(fam_req, collapse = "|"), ")\\b"), response)]
  cc <- cc[match_levels(term_label)]
  if (!nrow(cc)) die("nothing to plot after the filters")
  cc[, `:=`(is_sig = !is.na(covariate_q) & covariate_q < alpha & p_value < 0.05, value = estimate / log(2),
            col_key = term_label, row = paste(feature_type, response_id, sep = "::"))]
  type_lab <- unique(T11[, .(feature_type, feature_type_label)])
  cc[, type := factor(type_lab$feature_type_label[match(feature_type, type_lab$feature_type)])]
  rows <- pick_rows(cc, "row")
  if (!length(rows)) die("no response passes the filters", if (has("--sig-only")) sprintf(" (none significant at alpha = %g)", alpha))
  cc <- cc[row %in% rows]
  rows <- cc[, .(p = min(p_value)), by = .(row, type)][order(type, p), row]       # group by response type, strongest first
  cc[, tl := factor(stringr_wrap(term_label, 24), levels = unique(stringr_wrap(cc[order(match(covariate, covs)), term_label], 24)))]
  lim <- max(1, min(6, max(abs(cc$value), na.rm = TRUE)))
  models <- unique(cc$model)
  legend <- if (identical(models, "quasibinomial")) "log2 odds ratio" else if (identical(models, "quasipoisson")) "log2 rate ratio" else "log2 odds or\nrate ratio"
  title <- opt("--title", sprintf("Sample-level models: effects of %s",
                                  paste(unique(T11$covariate_label[match(unique(cc$covariate), T11$covariate)]), collapse = ", ")))
  ot <- cc[order(match(row, rows), tl), .(feature_type, response_id, response, covariate, comparison = term_label, log2_ratio = value,
                                          ratio = exp(estimate), p_value, covariate_q, significant = is_sig)]
  fwrite(ot, sub("\\.[A-Za-z]+$", ".tsv", out), sep = "\t", quote = FALSE)
  save_parts(rows, function(pr, subt) {
    dp <- cc[row %in% pr]
    lv <- rev(unique(stringr_wrap(dp$response[match(pr, dp$row)], 60)))
    dp[, resp := factor(stringr_wrap(response, 60), levels = lv)]
    heat(dp, "tl", "resp", legend, lim, title,
         paste(c(subt, sprintf("star: covariate F test q < %g and level differs from the reference (p < 0.05)", alpha)), collapse = "; "),
         drop_y = TRUE) +
      scale_x_discrete(drop = FALSE) + facet_grid(type ~ ., scales = "free_y", space = "free_y", labeller = label_wrap_gen(18))
  }, max(9, 7 + 0.75 * nlevels(cc$tl)), 0.5, 3)
}
