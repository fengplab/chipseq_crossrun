## One-off pie comparisons between hand-defined groups of samples (not part of run_all.sh).
## Usage: Rscript extra/compare_pies.R config.yaml comparisons.tsv [--peaks-only] [--equal-weight]
## comparisons.tsv (tab-separated, header): comparison  group  filter
##   filter = manifest column=value conditions: ';' = AND, '|' = OR between values, exact match after trimming,
##            e.g.  cell_type=THP1;virus=IAV|HSV-1     or  sample_id=S1|S2|S3
##   rows with the same comparison are drawn side by side, one pie per group, in file order. A sample may be in only one
##   group of a comparison; a group that matches no sample stops the script.
## Output: <outdir>/extra/pies/compare/
##   <comparison>_members.tsv                         which samples went into which group
##   <comparison>_peaks_{ring,noring}.{png,tsv}       peaks per host chromosome (chrM and viral genomes are never peak-calling contigs)
##   <comparison>_coverage_rpkm_{host,with_viruses}.{png,tsv}   IP reads per contig as RPKM shares (needs step 9; not with --peaks-only)
## Pooling: by default a group's peaks (reads) are summed over its samples, so samples with more peaks weigh more.
## --equal-weight: each sample's shares are computed first and averaged, so every sample counts the same
##   (file names get "_equal": <comparison>_peaks_equal_ring.png; the TSV 'value' is then the mean share x 1e6).
## Needs steps 0-1 (peaks) and step 9 (coverage) of the main run.
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "..", "R", "utils.R"))
cfg <- args_config()
ca <- commandArgs(TRUE)
if (length(ca) < 2 || startsWith(ca[2], "--")) stop("usage: Rscript extra/compare_pies.R config.yaml comparisons.tsv [--peaks-only] [--equal-weight]")
peaks_only <- "--peaks-only" %in% ca; equal_w <- "--equal-weight" %in% ca
source(file.path(here, "..", "R", "pie_lib.R"))
out <- file.path(cfg$outdir, "extra", "pies", "compare"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
S <- samples_dt(cfg)
CMP <- fread(ca[2], sep = "\t", colClasses = "character")
if (!all(c("comparison", "group", "filter") %in% names(CMP))) stop(ca[2], ": needs columns comparison, group, filter")
safe <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

## samples matching one filter string
match_filter <- function(f) {
  keep <- rep(TRUE, nrow(S))
  for (cond in trimws(strsplit(f, ";")[[1]])) {
    if (!nzchar(cond)) next
    kv <- strsplit(cond, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) stop("bad condition '", cond, "' (use column=value or column=v1|v2)")
    col <- trimws(kv[1]); vals <- trimws(strsplit(kv[2], "|", fixed = TRUE)[[1]])
    if (!col %in% names(S)) stop("'", col, "' is not a manifest column. Columns: ", paste(names(S), collapse = ", "))
    x <- trimws(as.character(S[[col]])); x[is.na(x)] <- ""
    keep <- keep & x %in% vals
  }
  S$sample_id[keep]
}

peaks_of <- function(s) {
  b <- rd(file.path(cfg$work, "peaks", paste0(s, ".bed5")), header = FALSE)
  if (is.null(b) || !nrow(b)) return(data.table(sample_id = s, contig = character(), n = integer()))
  b[, .(n = .N), by = .(contig = V1)][, sample_id := s][]
}
hosts <- intersect(order_all, G$contig)
CR <- if (peaks_only) NULL else rd(file.path(cfg$work, "contigs", "contig_reads.tsv"))
vmq <- as.integer((cfg$contigs %||% list())$viral_mapq %||% 10)
sets <- list(with_viruses = order_all, host = c(nuc, excluded))
sfx <- if (equal_w) "_equal" else ""                 # equal-weight figures get their own names
wt_txt <- if (equal_w) "mean of per-sample shares (each sample weighs the same)" else "summed over the group's samples"
n_files <- 0

for (cmp in unique(CMP$comparison)) {
  rows <- CMP[comparison == cmp]
  mem <- rbindlist(lapply(seq_len(nrow(rows)), function(i) {
    ids <- match_filter(rows$filter[i])
    if (!length(ids)) stop(sprintf("comparison '%s', group '%s': filter '%s' matches no sample", cmp, rows$group[i], rows$filter[i]))
    data.table(group = rows$group[i], sample_id = ids)
  }))
  dup <- mem[duplicated(sample_id), unique(sample_id)]
  if (length(dup)) stop(sprintf("comparison '%s': samples in more than one group: %s", cmp, paste(dup, collapse = ", ")))
  lab <- mem[, .(n = .N), by = group][, panel := sprintf("%s (n = %d)", group, n)]
  mem[, panel := factor(lab$panel[match(group, lab$group)], lab$panel)]
  fwrite(mem[, .(group, sample_id)], file.path(out, paste0(safe(cmp), "_members.tsv")), sep = "\t")
  msg(sprintf("%s: %s", cmp, paste(sprintf("%s = %s", lab$group, mem[, paste(sample_id, collapse = ", "), by = group]$V1), collapse = "  |  ")))

  ## ---- peaks
  P <- rbindlist(lapply(mem$sample_id, peaks_of))
  off <- P[!contig %in% hosts]
  if (nrow(off)) stop("peaks on contigs outside the host nuclear chromosomes: ", paste(unique(off$contig), collapse = ", "))
  empty <- setdiff(mem$sample_id, P$sample_id)
  if (length(empty)) msg(sprintf("%s: samples with no peaks (count as zero): %s", cmp, paste(empty, collapse = ", ")))
  P <- merge(P, mem, by = "sample_id")
  full <- merge(CJ(sample_id = mem$sample_id, contig = hosts), P[, .(sample_id, contig, n)], by = c("sample_id", "contig"), all.x = TRUE)
  full[is.na(n), n := 0L]; full <- merge(full, mem, by = "sample_id")
  d <- if (equal_w) {
    full[, sh := if (sum(n) > 0) n / sum(n) else 0, by = sample_id][, .(value = mean(sh) * 1e6), by = .(panel, contig)]
  } else full[, .(value = as.numeric(sum(n))), by = .(panel, contig)]
  d <- merge(d, G, by = "contig")
  tot <- full[, .(peaks = sum(n)), by = panel][order(panel)]
  if (sum(tot$peaks) != sum(P$n)) stop("peak counts do not add up in ", cmp)
  for (ring in c(TRUE, FALSE)) {
    f <- file.path(out, sprintf("%s_peaks%s_%s.png", safe(cmp), sfx, if (ring) "ring" else "noring"))
    draw_pie(d, hosts, sprintf("Peaks per chromosome: %s", cmp),
             sprintf("peaks %s; %s", wt_txt, paste(sprintf("%s: %s peaks", sub(" \\(n = .*", "", tot$panel), format(tot$peaks, big.mark = ",", trim = TRUE)), collapse = "; ")),
             f, ring)
    n_files <- n_files + 1
  }

  ## ---- IP coverage (RPKM shares)
  if (!is.null(CR) && nrow(CR)) {
    ip <- merge(CR[library == "IP" & contig %in% order_all], mem, by = "sample_id")
    miss <- setdiff(mem$sample_id, ip$sample_id)
    if (length(miss)) msg(sprintf("%s: no IP read counts for %s (left out of the coverage pies)", cmp, paste(miss, collapse = ", ")))
    for (nm in names(sets)) {
      x <- ip[contig %in% sets[[nm]]]
      pd <- if (equal_w) {
        x[, r := reads / (length / 1e3) / (sum(reads) / 1e6), by = sample_id][, r := r / sum(r), by = sample_id][
          , .(value = mean(r) * 1e6, reads = sum(reads)), by = .(panel, contig, length)]
      } else {
        x[, .(reads = sum(reads)), by = .(panel, contig, length)][, value := reads / (length / 1e3) / (sum(reads) / 1e6), by = panel]
      }
      f <- file.path(out, sprintf("%s_coverage_rpkm%s_%s.png", safe(cmp), sfx, nm))
      draw_pie(pd[, .(panel, contig, value, length)], intersect(sets[[nm]], pd$contig), sprintf("IP coverage per contig (RPKM): %s", cmp),
               sprintf("slice = share of RPKM, %s; %s; viral reads at MAPQ >= %d", wt_txt,
                       if (nm == "host") "host chromosomes and chrM" else "with viral genomes", vmq), f, ring = FALSE)
      n_files <- n_files + 1
    }
  } else if (!peaks_only) msg("no work/contigs/contig_reads.tsv (step 9): coverage pies skipped")
}
msg(sprintf("%d comparison pie charts written to %s", n_files, out))
