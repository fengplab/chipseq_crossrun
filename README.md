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
./run_all.sh config.yaml 16 18      # regression models and the report
./run_all.sh config.yaml 18 18      # only rebuild the report
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
| 16 | `steps/16_model_features.sh` | Per-window features for the regression models: GC- and width-matched background windows, motif sites and tandem arrays per window, largest-overlap CenSat/RepeatMasker class, IP and control read counts per window (featureCounts); runs step 7 if it was skipped |
| 17 | `R/70_models.R` | Regression models: motif family selection, confounding analysis, peak-level mixed logistic model per family, sample-level quasi-binomial / quasi-Poisson models against metadata, time-point contrasts |
| 18 | `bin/build_report.py` | One self-contained HTML report: every figure with its meaning, axes and abbreviations, plus a data dictionary for every table column (no tables) |

`steps/optional_call_input_peaks.sh` calls MACS3 peaks on control BAMs. Point `regions.input_peaks_dir` at its output, and recurrent regions that are also "peaks" in inputs are flagged as likely artifacts.

## Outputs

**`<outdir>/report/crossrun_report.html`** is a single document containing every figure. Before each figure it explains what the figure means, what every axis, colour and symbol encodes, and the abbreviations used. It ends with a data dictionary covering every column of every table. Figure numbers match the figure file names. Open it in a browser; print to PDF for a paginated copy (each figure starts on a new page). Figure descriptions live in `report/figures.yaml` and column descriptions in `report/data_dictionary.tsv`. The build warns if a table column has no dictionary entry.

Run IDs (`batch_var`) are never a heatmap key, a model term, part of the confounding analysis or a PCA colour. When `batch_var` is set, the regression models use it only to restrict permutations and for run-sensitivity refits; `batch_var: null` turns that off.

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
- `table06e_sample_clusters`: cluster of each sample from Ward hierarchical clustering of the scaled sample features, with the number of clusters chosen by mean silhouette width (2–8).
- Figure: `fig16` PCA with a sample-cluster panel followed by one panel per covariate.
- Covariate tests live in the regression models (steps 16–17, below). The earlier per-feature linear models, linear time trends, detection tests and PERMANOVA (`table06_metadata_effects`, `table06b`, `table06c`, `table06d`, `table07`, `table07b`, Figures 14–15) were removed; step 13 deletes their files from earlier runs.

## Regression models (steps 16–17)

Settings are in the `models:` block of the config. After upgrading from an earlier version of the bundle, rerun step 0 once (`./run_all.sh config.yaml 0 0`) so the new settings reach `work/config.env`, then run steps 16–18.

**What enters the models**

- Peaks: `models.peak_qvalue` keeps peaks with MACS q-value at or below the cutoff (column 9 of the peak files). At most `max_peaks_per_sample` of them (random subset) get background windows and enter the peak-level model; the sample-level models use all peaks that pass.
- Motif families: `models.family_stat` picks the family-wise statistic and `family_cutoff` the threshold (null = no cutoff); `max_families: null` removes the limit on how many families are modelled. The default, `sea_pvalue`, combines the family's per-sample SEA p-values into one (Cauchy combination). `evalue` uses the lowest member E-value times the number of members, `recurrence` the number of samples in which a member was discovered, and `pooled_q` a pooled peaks-versus-background test. `table08_family_selection` lists every statistic for every family, whichever one is used. `pooled_q` tests the same thing as the peak-level model, so selecting on it inflates "universal" calls.
- Scannability: step 05 computes, for each family, the best p-value FIMO could give a perfect match (`work/motifs/family_best_p.tsv`). Families whose best possible p-value is above `motifs.fimo_thresh` (typically 6–7 bp or degenerate motifs) can never get a FIMO site; `table08_family_selection` marks them `scannable = FALSE` and they are left out of the models.
- Background: each modelled peak gets `background_draws` windows of the same width at GC-matched positions.

**Peak-level model (is a motif repeat-driven or universal?)** One mixed logistic regression per family, over peaks and one background window per peak:

`motif present ~ peak + peak:(IP/input excess) + repeat class + peak:repeat class + log2 width + GC + peak:covariates + (1 + peak | sample)`

- *Repeat class* is the window's largest-overlap CenSat class, else RepeatMasker class, if it covers at least `class_min_fraction` of the window. Its coefficients are sequence effects: how much more often the motif occurs in that repeat than outside repeats, in peaks and background alike.
- *peak* is the enrichment in peaks outside repeats, for the average sample, at the sample's baseline IP/input ratio. The baseline is the `ratio_baseline_quantile` (default 10th percentile) of the ratio among the sample's own peaks, so the peak term absorbs the ratio a window needs to be called at all and the slope measures binding strength above it.
- *peak:repeat class* is extra enrichment in peaks inside a repeat class.
- Classification (`table09_motif_classification`, Figure 18), at `alpha` after BH correction across families: **repeat-specific binding** if the interaction terms are significant; otherwise **universal** if the peak term is significant and positive; otherwise **repeat-driven** if only the sequence effects are significant; otherwise **no association**. **Confined to repeats** means at least `confined_fraction` of the motif's windows lie in repeats and too few elsewhere to fit anything.
- A repeat class enters a family's model only if peaks and background each have at least `min_cell_count` windows with and without the motif. Classes failing that (the motif is always or never there) have no finite coefficient and are listed in `table09d_peak_model_excluded_classes`.
- A region that is a peak in several samples has the same sequence each time. Each peak is therefore weighted by 1 / (number of samples sharing its merged region), so a region counts once (`region_weighting`). `n_motif_regions` in table09 shows how many distinct regions carry the motif; when it is far below `n_motif_peaks`, a large pooled odds ratio comes from a few recurrent regions.
- Covariate effects on peak enrichment are in `table09c_peak_model_covariates`, with a second p-value from a refit in which the peak effect also varies by run.

**Sample-level models (what varies with metadata?)** One GLM per response against all covariates (`table11_covariate_tests`, `table11b_covariate_coefficients`, Figure 19):

| Response | Model | Offset |
|---|---|---|
| IP reads per repeat class | quasi-binomial | control library's fraction in the class |
| IP reads in peaks | quasi-binomial | control library's fraction in the same peaks |
| Share of peaks per repeat class | quasi-binomial | none |
| Tandem array score (sites in runs of ≥ 3 evenly spaced sites / all sites) | quasi-binomial | none |
| Motif density per kb in peaks | quasi-Poisson | kb of peak sequence × background density |
| Motif density per kb inside each repeat class | quasi-Poisson | kb of peak sequence in the class |

- Time is an ordered category (`time_bins`, shared across series), not a continuous variable.
- Each covariate gets an F test, BH-corrected within the type of response. Tests with p ≤ `perm_max_p` also get a permutation p-value, with labels shuffled within `batch_var` where it is set and the covariate varies within it, else within `block_fallback`, else freely (`perm_block` says which).
- When `batch_var` is set, `p_value_within_run` repeats the test with it as a blocking factor (NA, with a note, when the covariate is constant within batches). With `batch_var: null` the column is omitted.
- `table11c_stratified_tests` refits one covariate at a time inside each level of `stratify_by`.

**Confounding** (`table10*`, Figure 17). Cramér's V between every pair of covariates; which covariate levels are aliased and were dropped; and, for each covariate, the share of variance in the sample-level responses that it explains alone or shares with other covariates. Run is not part of this analysis.

**Time-point contrasts** (`table12_timepoint_contrasts`, Figure 20). Within each series, motif density (peaks relative to background) at each time point against the first, with bootstrap intervals over the windows of each library. A time point is marked as different only if q < `alpha` and the change is at least `min_log2_change` log2 units. Families with fewer than 5 motif sites in either library, or whose bootstrap does not vary, are not tested. With one library per time point these show that two libraries differ, not that the difference is reproducible.

**Heatmap rows.** Figures 4, 18 and 19 show, for every column, every row significant in that column, or that column's 5 rows with the smallest p-values when fewer are significant; the figure shows the union. Figures 18 and 19 are split into parts of at most 60 rows. Figure 1 shows the 5 families most enriched over background in each sample; Figure 6 every chromosome; Figure 7 (CenSat) every sample and every CenSat class; Figure 11 every sample.

If `repeats.reads_in_repeats` was off, step 16 counts reads per repeat class itself; rerun step 10 afterwards if you also want Figure 7 drawn from those counts.

**Run time.** Step 16 runs FIMO on the background windows and featureCounts on every IP and control BAM. Step 17 took about 3 minutes for 6 families at 30 samples × 5,000 peaks on 8 threads in testing; time grows with the number of families and windows.

## MEME Suite 5.5.4 and later

From version 5.5.4, FIMO (and SEA from 5.5.5) read FASTA headers of the form `name::chrom:start-end`, as written by `bedtools getfasta -name`, as genomic coordinates: they report the chromosome as the sequence name and genome positions as site positions. The bundle therefore writes headers that contain only the window name, passes `--no-pgc` to FIMO and SEA when the installed version supports it, and stops with an error if FIMO reports names that are not the scanned windows. Outputs of steps 03, 05 and 16 carry a format stamp, so a bundle update that changes sequence naming reruns them automatically. Results produced by earlier versions of the bundle with MEME 5.5.4+ (motif prevalence, spacing, sites in repeats, step 16–17 motif counts) should be regenerated from step 3.

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

- **Mappability.** GenMap computes, for every 150-mer (`qc.mappability_k`), whether it occurs once in the genome (both strands, no mismatches). Each base's mappability is the fraction of the k-mers covering it that are single-copy. These per-base values are summarised over every base of each CenSat class (`qc02_censat_mappability_stats`) and over every peak. Peak mappability is related to peak score and CenSat class by linear models (`qc01_peak_mappability_models`; per-class trends in `qc01c_mappability_score_by_censat`). QC figure 1 shows mappability against score per class, with a trend line only where significant. QC figure 2 shows violins of peak mappability per class, against each class's genome-wide median. Indexing the genome is a one-off (about an hour for a human genome); keep the index via `qc.genmap_index`, or skip GenMap entirely by supplying a k-mer bedGraph in `qc.mappability_bedgraph`. Processing needs about 3 GB of RAM for the largest chromosome. The per-base track is written to `work/qc/mappability_k150.bedgraph.gz`.
- **Insert sizes.** For every IP and control library, fragments whose midpoint lies in each CenSat class are counted by length, plus a random sample of fragments outside CenSat (`qc.insert_genome_fraction`). The output is one figure per IP sample comparing it with its control (QC figure 3), a summary of median shifts (QC figure 4), and `qc03_insert_size_stats`. Single-end libraries are reported and skipped.
- **Strand balance and read stacking** (QC figure 5, `qc05_peak_strand_stacking_by_censat`). For every peak, read 1 start positions give strand imbalance, start-site diversity relative to random placement, and the share of reads at the most-used start, compared between CenSat classes and peaks outside CenSat. Strand-biased, stacked peaks in repeats are likely alignment artifacts.
- **Position within arrays** (QC figure 6, `qc06_array_edge_distance_by_censat`). Distance from each summit to the nearest edge of its CenSat array, compared with random positions in the same arrays by Kolmogorov–Smirnov test. Pile-ups at array edges are a hallmark of mapping artifacts at unique/repeat boundaries.
- **Reuse.** All QC results reuse step 01's peaks and annotation and step 06's peak CenSat labels.

## Utilities

**Human-only JASPAR motifs** (`bin/filter_jaspar_species.py`). JASPAR's MEME downloads carry no species; this queries the JASPAR REST API for the matrices annotated to a species and keeps only those. A motif is also kept if another version of the same profile is annotated to the species. Motif blocks are copied verbatim, and a `.species_filter.tsv` records what was kept or dropped and why.
```bash
python3 bin/filter_jaspar_species.py --meme JASPAR2024_CORE_vertebrates_non-redundant.meme --tax-id 9606 \
    --out JASPAR2024_CORE_vertebrates_human.meme --write-ids jaspar_human_ids.txt
# offline (e.g. a node without internet): reuse the ID list
python3 bin/filter_jaspar_species.py --meme JASPAR...meme --ids-file jaspar_human_ids.txt --out JASPAR..._human.meme
```
Point `known_motif_dbs` at the filtered file and rerun from step 5. Step 5 notices that the databases changed, by path or content, and recomputes the known-motif SEA and Tomtom results; everything else is reused.

**IGV session for one sample** (`bin/make_igv_session.py`, standalone: Python standard library only). It writes `session.xml` with the genome FASTA, its GTF annotation, the CenSat and RepeatMasker bigBeds, the sample's peaks, and its reads (IGV adds the coverage track). First it checks that every file uses the same chromosome notation as the genome, reading the names from the FASTA index, BAM header, bigBed chromosome index, and the GTF and peak files, and reports each file's notation and how many names match. If notations differ it stops (exit 2) unless told how to resolve it:
- `--resolve alias` writes an IGV genome definition (`genome.json`) with an alias file, so IGV maps all names onto the FASTA's names. No files are copied.
- `--resolve chr-fasta` also writes a header-renamed copy of the FASTA, so IGV displays `chr1`, `chr2`, ….
```bash
python3 bin/make_igv_session.py --out igv/HSV_24h --crossrun-config config_human.yaml \
    --sample THP1_HSV1_8h_IP-cGAS_REP1 --gtf annotation.gtf --resolve chr-fasta
```
With `--crossrun-config`, the FASTA, bigBeds, BAM, peaks and the bundle's alias table come from the run; explicit options override them. A GTF that is already bgzipped and indexed (`.tbi` or `.csi` beside it) is used as-is, and its chromosome names are read from the index; the script warns if the index is older than the GTF. Otherwise the script writes a sorted, bgzipped, tabix-indexed copy so IGV loads it quickly (`--no-index-gtf` uses the GTF as given). The two large derived files are written once, beside the originals, and reused by every later session, so all sessions reference one central copy: `chr_<fasta name>` plus `.fai` next to the FASTA (with `--resolve chr-fasta`), and `sorted_<gtf name>.gz` plus `.tbi` next to the GTF. Only the file-name prefix changes, because IGV chooses the format from the extension; a GTF must therefore end in `.gtf`, `.gff` or `.gff3` (optionally `.gz`), and the script warns otherwise. A derived file is rebuilt when its original is newer. `--shared-dir DIR` puts both in another folder; if the original's folder is not writable and no `--shared-dir` is given, the copy goes into the session folder with a warning. Copies that earlier versions wrote into a session folder (`genome.chr.fa`, `<name>.gtf.sorted.gz`, `<name>.sorted.gtf.gz`) are moved to the central location on the next run for that folder; sessions in other folders that pointed at them need regenerating. `--rewrite-prefix OLD=NEW` (repeatable) rewrites paths in the session for a computer that mounts the files elsewhere, e.g. `--rewrite-prefix /home/bradlem9/isilon=/Volumes/fengplab`. Files are still read through their local paths for the notation check. Matching is by whole folder name, and trailing slashes don't matter. Any session path that matches no prefix is listed as a warning, since IGV on the other computer won't find it. Put `--out` somewhere the other computer can see, because the session, `genome.json` and the alias file are written there; the FASTA and GTF folders must be visible there too. Open the result in IGV Desktop (2.11 or later) with File > Open Session.

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

This builds a synthetic genome with planted satellites, a CENP-B-box array, and ETS sites, plus 10 fake pipeline runs. Satellite enrichment rises with time in the WT-construct series only. The whole bundle then runs on it, with stub `fimo`/`tomtom`/`sea` that do exact consensus matching. The stubs exercise the plumbing only; they are not substitutes for MEME Suite. The script then runs `R/70_models.R` on simulated windows (`tests/make_model_test_data.py`) with planted motif behaviours and covariate effects, and checks that the classifications and effects are recovered.
