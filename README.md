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
| 2 | `steps/02_motif_discovery.sh` | Motif discovery: rerun XSTREME per sample with control sequences from its control alignment (default), or reuse the pipeline's XSTREME; then parse `xstreme.tsv`, `sea_out/sea.tsv`, `combined.meme` and pool all de novo motifs |
| 3 | `steps/03_sequences.sh` | Peak sequences (whole peak or summit ± flank), plus one genomic-background relocation per sample |
| 4 | `steps/04_motif_families.sh` | Tomtom all-vs-all, then clustering into cross-run motif families with one representative PWM each |
| 5 | `steps/05_motif_scan.sh` | FIMO scan of every sample (and its background) with every family; SEA of known motifs per sample; Tomtom vs known |
| 6 | `steps/06_repeat_enrichment.sh` | Observed vs shuffled peak counts per repeat/CenSat label |
| 7 | `steps/07_reads_in_repeats.sh` | Optional: IP vs control read fractions per label (featureCounts) |
| 8 | `steps/08_regions.sh` | Merged cross-sample regions, labels, nearest gene, input overlap, pairwise Jaccard, genome bins |
| 9–13 | `R/10_motifs.R` … `R/50_metadata.R` | Tables and figures for Q1–Q5 |
| 14 | `steps/14_qc.sh` | QC: single-copy k-mer mappability (GenMap, k = insert size) turned into per-base mappability and summarised per CenSat class and per peak; insert-size distributions per CenSat class for every IP and control library |
| 15 | `R/60_qc.R` | QC tables and figures (report section "Quality control") |
| 16 | `bin/build_report.py` | One self-contained HTML report: every figure with its meaning, axes and abbreviations, plus a data dictionary for every table column (no tables) |

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
- `table06e_sample_clusters`: cluster of each sample from Ward hierarchical clustering of the scaled sample features, with the number of clusters chosen by mean silhouette width (2–8).
- Figures: `fig14` effects heatmap (top 30 features), `fig15_timecourse_<series>` one figure per series showing its significant trends (top 20), `fig16` PCA with a sample-cluster panel followed by one panel per covariate.

## Motif discovery with control-library sequences

With `motifs.discovery: control_bam` (the default), step 02 reruns XSTREME for every sample with `--n` control sequences built from that sample's control (input/IgG) alignment, instead of shuffled peak sequences:
1. Control-library fragments are sampled from the control BAM: properly paired read-1 records for paired-end data, read alignments for single-end data.
2. A window is centred on each fragment midpoint. Window widths are drawn from the sample's own peak windows (whole peaks, or summit ± flank if `seq_mode: summit`), so primary and control sequences have the same length distribution. That matters because XSTREME's p-values assume it.
3. Windows overlapping the sample's peaks are discarded, and `control_ratio` sets how many control windows there are per peak window.
4. The sequences are taken from the genome FASTA at those positions.

This background carries the control library's own biases (mappability, copy number, fragmentation, open chromatin), so motifs must be enriched beyond them. Samples without a control BAM fall back to shuffled controls; `xstreme_run_info` records which was used. The pipeline's `pipeline_maxw.txt` is reused for `--maxw`.

**Discovery engine** (`motifs.discovery_engine`):
- **`text`** (default) runs the tools XSTREME wraps directly, each in its text-only mode, so no HTML report is ever written and a MEME install whose HTML templates are broken still works:
  1. `fasta-get-markov -m 2` builds an order-2 background model from the control sequences, as XSTREME does.
  2. `streme --text` discovers motifs against the controls.
  3. `meme -text` (ZOOPS model, differential-enrichment objective, `-neg` controls, `meme_nmotifs` motifs) discovers motifs against the same controls.
  4. `sea --text` scores every discovered motif.

  The results are written in XSTREME's file layout, so all later steps are unchanged.
- **`xstreme`** runs XSTREME itself, which always writes an HTML report.

Step 02 first runs the chosen engine's tools on a tiny input and stops with a clear message if they fail. If MEME itself is broken, fix the install or set `meme_prefix` (e.g. `apptainer exec /path/to/meme.sif`); the prefix is then used for every MEME tool, and paths must be visible inside the container. Changing discovery changes the motifs, and steps 04–05 detect that and recompute families and scans. Set `motifs.discovery: pipeline` to reuse the pipeline's results instead.

## Quality control

- **Mappability.** GenMap computes, for every 150-mer (`qc.mappability_k`), whether it occurs once in the genome (both strands, no mismatches). Each base's mappability is the fraction of the k-mers covering it that are single-copy. These per-base values are summarised over every base of each CenSat class (`qc02_censat_mappability_stats`, QC figure 3) and over every peak, and peak mappability is related to peak score and CenSat class by linear models (`qc01_peak_mappability_models`, QC figures 1–2). Indexing the genome is a one-off (about an hour for a human genome); keep the index via `qc.genmap_index`, or skip GenMap entirely by supplying a k-mer bedGraph in `qc.mappability_bedgraph`. Processing needs about 3 GB of RAM for the largest chromosome. The per-base track is written to `work/qc/mappability_k150.bedgraph.gz`.
- **Insert sizes.** For every IP and control library, fragments whose midpoint lies in each CenSat class are counted by length, plus a random sample of fragments outside CenSat (`qc.insert_genome_fraction`). The output is one figure per IP sample comparing it with its control (QC figure 4), a summary of median shifts (QC figure 5), and `qc03_insert_size_stats`. Single-end libraries are reported and skipped.
- **Reuse.** All QC results reuse step 01's peaks and annotation and step 06's peak CenSat labels.

## Methods notes

**Motif families.** Tomtom (Pearson, min overlap 5) is run all-vs-all on every de novo motif. The distance is d = min(1, E / `family_evalue`), and families come from average-linkage clustering cut at `family_cut_height`. Each family's representative is the member with the best XSTREME SEA p-value. Because every sample is rescanned with every representative, a family found de novo in only one run is still counted everywhere.

**Motif occurrence and detection.** FIMO (p < `fimo_thresh`) reports every site, unlike MEME's ZOOPS `SITES` column. Palindromic sites counted on both strands are collapsed. A family counts as "detected" in a sample when it passes the q and fold thresholds (`detect_qvalue`, `detect_min_enrichment`) against the reference set by `motifs.detect_against`:
- **`genomic`** (default): peaks versus the same windows relocated at random in the genome. This is the stricter question: is the motif more common at the peaks than in DNA at large?
- **`shuffled`**: peaks versus shuffled copies of the same peaks, via SEA. This is the comparison XSTREME used to discover the motifs. Repeats and other sequence structure pass it easily, because shuffling destroys that structure.
- **`both`** or **`either`** combine the two.

Step 9 prints a one-line diagnosis to the log: median prevalence in peaks and background, median fold, and how many tests pass each criterion. Use it to tell whether nothing is enriched over the genome, or whether prevalence is saturated (most windows contain a site at all), which calls for `seq_mode: summit` or a smaller `fimo_thresh`.

Motif families enter the Q5 models if they are detected anywhere or found de novo in at least `q5_min_denovo_samples` samples, capped at the `q5_max_motif_features` most variable.

**Expected counts.** Each sample's peaks are shuffled `n_shuffles` times with `bedtools shuffle`, preserving widths. The shuffles use `include_bed`/`exclude_bed`, and `-chrom` if `same_chrom`. With `null_model: gc_matched`, each peak takes the closest-GC relocation out of `gc_candidates` tries. Observed and shuffled peaks use the same overlap rule. The default `midpoint` rule counts a peak when its midpoint lies inside a feature; `fraction` counts it when at least `overlap_fraction` of the peak lies inside. The pipeline's own `class_summary.tsv` (≥ 80% of the feature covered) is not used here, because a peak can never cover 80% of a multi-kb satellite array. The empirical p-value floor is 1/(n_shuffles + 1), so use ≥ 1000 shuffles for final q-values across many labels.

**Interpretation.** Motif similarity to a known factor is not protein homology. For a non-sequence-specific DNA binder, known-motif hits indicate whose binding sites or which repeat structures the peaks sit on.

## Testing

```bash
bash tests/run_test.sh /tmp/crossrun_test
```

This builds a synthetic genome with planted satellites, a CENP-B-box array, and ETS sites, plus 10 fake pipeline runs. Satellite enrichment rises with time in the WT-construct series only. The whole bundle then runs on it, with stub `fimo`/`tomtom`/`sea` that do exact consensus matching. The stubs exercise the plumbing only; they are not substitutes for MEME Suite.
