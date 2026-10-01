# chipseq_crossrun

Cross-run summaries of motif, repeat, region and metadata results from many runs of the
[fengplab/chipseq](https://github.com/fengplab/chipseq) pipeline (`--run_xstreme`, optionally `--annotate_features`).
It runs outside the pipeline: point it at the pipeline output directories through a manifest, and it
writes tables (TSV) and figures (PDF + PNG). Heavy lifting uses MEME Suite and bedtools; all tables
and figures are built in R.

## Install

```bash
conda env create -f environment.yml
conda activate chipseq_crossrun
```

If you also use R outside conda, packages in your personal library (`~/R/<platform>-library/<version>`) can shadow the env's copies and crash conda's R. `run_all.sh` hides that library whenever `Rscript` comes from the active conda env, and checks that every R package loads before any step runs. Set `CROSSRUN_ALLOW_USER_RLIB=1` to turn the isolation off.

## Inputs

**`manifest.tsv`** (see `manifest_template.tsv`): one row per IP sample.

- **Required columns:**
  - `sample_id`: the pipeline's `meta.id`, i.e. the samplesheet `sample` plus `_REP<replicate>` (e.g. `HSV_THP1_WT_24h_REP1`). If you leave off `_REP<n>`, setup accepts the ID when exactly one replicate matches, and warns.
  - `run_id`.
  - `results_dir`: that run's `--outdir`. Relative paths are resolved from the manifest's folder.
- **Metadata columns:** any names you like. List them in `config.yaml` under `categorical_covariates`, `numeric_covariates`, `series_var` (time-course ID, blank for stand-alone samples) and `time_var`.
- Spreadsheet placeholders (`N/A`, `NA`, `#N/A`, `null`, `-`) in any cell are treated as blank. `none` is kept as a real value, e.g. `virus = none` for mock samples.
- **Optional columns:**
  - `control_id`: used for read-level control comparisons and input-peak artifact flags.
  - `aligner`, `peak_type`: per-row overrides; leave blank (or `N/A`/`NA`) to auto-detect. By default (`auto`) both are detected per run from the `<aligner>/merged_library/macs3/<narrow_peak|broad_peak>` directories that exist. The pipeline's own defaults are `bwa` and `broad_peak`.
  - `peaks_file`, `xstreme_dir`, `bam_file`, `homer_file`: explicit paths, if the standard layout doesn't apply.
  - `with_allo`, `genome_build`: bookkeeping.

**`config.yaml`** sets the genome, annotations, motif databases and parameters. Every setting is commented. The main things to set:

- **Genome:** `genome_fasta` (indexed) and `chrom_sizes` (`<outdir>/genome/*.sizes` from any run).
- **Annotations:** `repeatmasker` and `censat`. Use the pipeline's `genome/annotation/*.bb`, or any BED, UCSC `rmsk.txt`, or RepeatMasker `.out`.
- **Null-model regions:** `include_bed` and/or `exclude_bed`. The pipeline's `*.include_regions.bed` is ideal, because it encodes the blacklist you actually filtered on.
- **`known_motif_dbs`:** for example, JASPAR 2024 CORE vertebrates plus HOCOMOCO v12, both in MEME format. The same files are applied to every sample. This deliberately replaces the per-run `SIM_MOTIF` columns, which are empty wherever `--xstreme_motif_db` was not set.

**Chromosome names.** All tables and figures use UCSC-style names (`chr1`…`chr22`, `chrX`, `chrY`, `chrM`), whatever convention the inputs use. Peaks, BigBeds, BEDs, the chrom-sizes file and the include/exclude BEDs are renamed on the way in. Where a tool must read the original FASTA or BAMs, the canonical names are translated back to that file's own naming. The map is built automatically:

1. alias files you list under `chrom_alias` (an NCBI `*_assembly_report.txt`, a UCSC `chromAlias.txt`, or any TSV with a `chr`-style column);
2. FASTA header descriptions such as `>NC_060925.1 ... chromosome 1`;
3. a built-in T2T-CHM13v2.0 table covering RefSeq `NC_0609xx.1` and GenBank `CP0682xx.2` accessions;
4. simple rules (`1` → `chr1`, `MT` → `chrM`).

Check `work/annot/chrom_naming_report.tsv`, `rename_log.tsv` (records renamed or unrecognised per input) and `chrom_overlap.tsv` (chromosomes shared between peaks and each annotation). Step 01 stops if an annotation shares no chromosomes with the peaks.

All samples in one analysis must share a genome build; setup stops if they don't. Keep Allo and non-Allo runs, or runs with different blacklists, in separate analyses if satellites are the focus. The filtering changes what is callable in centromeric DNA.

## Run

```bash
./run_all.sh config.yaml            # everything
./run_all.sh config.yaml 9 13       # only re-draw the R tables/figures
```

Steps are resumable. Finished per-sample outputs are reused; delete them from `<outdir>/work` to force a rerun.

| Step | Script | What it does |
|---|---|---|
| 0 | `bin/00_setup.py` | Validate config/manifest, locate every pipeline file, write `work/samples.tsv`, `work/layouts.tsv` (layouts found and used) and `work/missing_files.tsv` (what was looked for, and what exists nearby) |
| 1 | `steps/01_prepare_inputs.sh` | Chromosome alias table and canonical genome sizes; label BEDs for RepeatMasker (`class/family`) and CenSat (class); normalized peak BEDs; naming consistency checks |
| 2 | `bin/02_parse_xstreme.py` | Parse `xstreme.tsv`, `sea_out/sea.tsv`, `combined.meme` for every sample; pool all de novo motifs |
| 3 | `steps/03_sequences.sh` | Peak sequences (whole peak or summit ± flank), plus one genomic-background relocation per sample |
| 4 | `steps/04_motif_families.sh` | Tomtom all-vs-all, then clustering into cross-run motif families with one representative PWM each |
| 5 | `steps/05_motif_scan.sh` | FIMO scan of every sample (and its background) with every family; SEA of known motifs per sample; Tomtom vs known |
| 6 | `steps/06_repeat_enrichment.sh` | Observed vs shuffled peak counts per repeat/CenSat label |
| 7 | `steps/07_reads_in_repeats.sh` | Optional: IP vs control read fractions per label (featureCounts) |
| 8 | `steps/08_regions.sh` | Merged cross-sample regions, labels, nearest gene, input overlap, pairwise Jaccard, genome bins |
| 9–13 | `R/10_motifs.R` … `R/50_metadata.R` | Tables and figures for Q1–Q5 |
| 14 | `bin/build_report.py` | One self-contained HTML report: every figure with its meaning, axes and abbreviations, plus a data dictionary for every table column (no tables) |

`steps/optional_call_input_peaks.sh` calls MACS3 peaks on control BAMs. Point `regions.input_peaks_dir` at its output, and recurrent regions that are also "peaks" in inputs are flagged as likely artifacts.

## Outputs

**`<outdir>/report/crossrun_report.html`** is a single document containing every figure. Before each figure it explains what the figure means, what every axis, colour and symbol encodes, and the abbreviations used. It ends with a data dictionary covering every column of every table. Figure numbers match the figure file names. Open it in a browser; print to PDF for a paginated copy (each figure starts on a new page). Figure descriptions live in `report/figures.yaml` and column descriptions in `report/data_dictionary.tsv`. The build warns if a table column has no dictionary entry.

Run IDs (`batch_var`) are provenance only. They are never a heatmap key, a model term, a PERMANOVA stratum or a PCA colour.

All outputs go to `<outdir>/tables` and `<outdir>/figures`.

**Q1. Motif occurrence and repeat-likeness**

- `table01_master_motifs`: every de novo motif, with its family, best known match, and sample metadata.
- `table02_motif_prevalence`: family × sample. Contains % of peaks with ≥ 1 site, genomic-background %, fold enrichment, Fisher q, `detected`, sites per occupied peak, and whether XSTREME found it de novo in that sample. The `_matrix` version is wide.
- `table02_repeat_likeness`: one row per family. Contains:
  - mean and max sites per occupied peak
  - spacing mode and periodic fraction
  - PWM self-periodicity and trimer diversity
  - % of sites in RepeatMasker and CenSat labels, and the top label
  - optional best hit among repeat consensus sequences
  - a flag: `repeat (tandem)`, `repeat-embedded`, `low-complexity / periodic`, or `TF-like`
- Figures: `fig01` prevalence heatmap, `fig02` repeat-likeness scatter, `fig03` site-spacing histograms.

**Q2. Peaks over repeats, centromeres and satellites**

- `table03_repeat_enrichment`: observed, expected (mean of shuffles), SD, log2 O/E, z-score, empirical enrichment/depletion p and BH q, with the null settings recorded in each row.
- `table03_log2oe_*`: wide matrices.
- `table03b_censat_by_chromosome`: peaks per CenSat class per chromosome, summed over all samples.
- `table03c_repeat_enrichment_vs_universe`: this sample vs all merged regions across samples (a LOLA-style stricter null).
- `table03d_reads_in_repeats`: optional.
- Figures: `fig04` O/E heatmaps, `fig05` CenSat composition, `fig06` CenSat by chromosome, `fig07` read-level heatmaps.

**Q3. Similarity to known DNA-binding-protein motifs**

- `table04_known_motif_ranking`: per known motif, the number of samples where SEA finds it enriched (with median enrichment ratio), the number of samples with a Tomtom-matching de novo motif, the best E-value, and the matched families.
- `table04b_known_by_tf_class`: collapses redundant motifs by TF class. Supply `tf_class_table`, which `bin/fetch_jaspar_classes.py` can build for JASPAR IDs; otherwise motifs are grouped by TF name.
- `table04c_family_known_matches`: top 5 known matches per family.
- `table04d_sea_known_long`.
- Figures: `fig08` dot plot (a black ring marks samples where a de novo motif also matched), `fig09` class bar chart.

**Q4. Recurrent regions and overlap**

- `table05_recurrent_regions` (top N) and `table05_all_regions`. Per region: number of samples, distinct conditions, batches and levels of each condition variable, repeat/CenSat labels, nearest gene, input-peak overlap, and a `likely_artifact` flag.
- `table05b_jaccard_matrix`, `table05c_hotspot_bins`.
- Figures: `fig10` recurrence histogram, `fig11` Jaccard heatmap, `fig12` UpSet by `upset_group_var`, `fig13` genome-wide hotspot bins.

**Q5. Metadata effects**

- `table06_sample_features`: one row per sample. Contains peak load, FRiP, HOMER composition, repeat/CenSat log2 O/E, motif family prevalence, and known-motif SEA enrichment.
- `table06_metadata_effects`: each feature regressed on all covariates jointly (ordinary least squares, no run term). Peak-count features also adjust for log10 depth. Falls back to one covariate at a time if there are too few samples for a joint model. Every row has a readable `feature_label` (e.g. "Peak enrichment in Human satellite 3 (HSat3) (log2 O/E)") and `term_label` (e.g. "virus: VACV vs HSV-1"). Aliased terms are reported as not estimable. Optional `covariate_labels` in the config renames covariates in labels.
- `table06b_covariate_associations`: how strongly each pair of covariates is tied (Cramér's V, eta², or |Spearman ρ|), flagging pairs whose effects cannot be separated.
- `table06c_timecourse_trends`: per series, linear slope versus time (with its p- and q-value) and Spearman ρ. Needs ≥ 3 time points; no replicates needed. `timecourse.significance` (`q` or `p`) and `timecourse.alpha` in the config decide what counts as significant. With 3–4 time points, q < 0.05 is rarely reachable.
- `table06d_timecourse_series_differences`: tests whether the time slope differs between series (time × series interaction).
- `table07_motif_detection_vs_metadata`: Fisher tests (categorical, level vs rest) or Wilcoxon tests (numeric) for detection of each motif family and known motif.
- `table07b_permanova`: `adonis2` per feature block, unrestricted and with permutations blocked within run.
- Figures: `fig14` effects heatmap (top 30 features), `fig15_timecourse_<series>` one figure per series showing its significant trends (top 20), `fig16` PCA colored by each covariate.

## Methods notes

**Motif families.** Tomtom (Pearson, min overlap 5) is run all-vs-all on every de novo motif. The distance is d = min(1, E / `family_evalue`), and families come from average-linkage clustering cut at `family_cut_height`. Each family's representative is the member with the best XSTREME SEA p-value. Because every sample is rescanned with every representative, a family found de novo in only one run is still counted everywhere.

**Motif occurrence.** FIMO (p < `fimo_thresh`) reports every site, unlike MEME's ZOOPS `SITES` column. Palindromic sites counted on both strands are collapsed. "Detected" means enriched over one genomic relocation of the same windows (Fisher q and fold thresholds in the config).

**Expected counts.** Each sample's peaks are shuffled `n_shuffles` times with `bedtools shuffle`, preserving widths. The shuffles use `include_bed`/`exclude_bed`, and `-chrom` if `same_chrom`. With `null_model: gc_matched`, each peak takes the closest-GC relocation out of `gc_candidates` tries. Observed and shuffled peaks use the same overlap rule. The default `midpoint` rule counts a peak when its midpoint lies inside a feature; `fraction` counts it when at least `overlap_fraction` of the peak lies inside. The pipeline's own `class_summary.tsv` (≥ 80% of the feature covered) is not used here, because a peak can never cover 80% of a multi-kb satellite array. The empirical p-value floor is 1/(n_shuffles + 1), so use ≥ 1000 shuffles for final q-values across many labels.

**Interpretation.** Motif similarity to a known factor is not protein homology. For a non-sequence-specific DNA binder, known-motif hits indicate whose binding sites or which repeat structures the peaks sit on.

## Testing

```bash
bash tests/run_test.sh /tmp/crossrun_test
```

This builds a synthetic genome with planted satellites, a CENP-B-box array, and ETS sites, plus 10 fake pipeline runs. Satellite enrichment rises with time in the WT-construct series only. The whole bundle then runs on it, with stub `fimo`/`tomtom`/`sea` that do exact consensus matching. The stubs exercise the plumbing only; they are not substitutes for MEME Suite.
