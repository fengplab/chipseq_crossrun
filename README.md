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
./run_all.sh config.yaml 16 21      # bin table, binding models (steps 16-20) and the report
./run_all.sh config.yaml 21 21      # only rebuild the report
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
| 16 | `steps/16_bin_table.sh` | Bin table for the binding models at the tuned bin size: per-bin repeat class, satellite array and edge distance, repeat fraction, mappability, GC and N fraction; genome-wide FIMO sites per motif family; per-library fragment coverage (all reads, MAPQ ≥ `bins.unique_mapq`, two read-name halves, fragment length); drift check against the tuning record |
| 17 | `R/65_design.R` | Binding models, Stage 0: sample-level design, estimable and aliased covariate levels, covariate associations (the model fits follow in a later version) |
| 18 | `bin/build_report.py` | One self-contained HTML report: every figure with its meaning, axes and abbreviations, plus a data dictionary for every table column (no tables) |

### Profiling resource use

`CROSSRUN_PROFILE=1 ./run_all.sh config.yaml` samples every step's whole process tree every 5 s (`CROSSRUN_PROFILE_INTERVAL` to change) with `bin/profile_run.py` (Linux `/proc`, no extra packages): memory (PSS, so forked R workers are not double-counted), cores in use, and disk I/O. Per-step traces go to `<outdir>/logs/profile/stepNN.tsv`. At the end, `bin/profile_summary.py <outdir>` prints and writes `summary.tsv` with wall time, peak memory, peak and mean cores, CPU hours and efficiency (mean cores / `threads`) per step, and a suggested `--mem`, `--cpus-per-task` and `--time`. Profile a subset of steps the same way, for example `CROSSRUN_PROFILE=1 ./run_all.sh config.yaml 16 17`.

Parallelism: steps 2, 3, 5, 6, 9 (BAM counts), 14, 16 and 19 run one sample per thread (`xargs -P threads`), so threads beyond the number of samples (or libraries, for 9 and 14) do not help. Step 7 passes `threads` to featureCounts; step 16 runs FIMO (one family per thread) and the per-library signal (one BAM per thread). Steps 1, 4, 8, 10–13, 15 and 18 are single-threaded.

`steps/optional_call_input_peaks.sh` calls MACS3 peaks on control BAMs. Point `regions.input_peaks_dir` at its output, and recurrent regions that are also "peaks" in inputs are flagged as likely artifacts.

## Outputs

**`<outdir>/report/crossrun_report.html`** is a single document containing every figure. Before each figure it explains what the figure means, what every axis, colour and symbol encodes, and the abbreviations used. It ends with a data dictionary covering every column of every table. Figure numbers match the figure file names. Open it in a browser; print to PDF for a paginated copy (each figure starts on a new page). Figure descriptions live in `report/figures.yaml` and column descriptions in `report/data_dictionary.tsv`. The build warns if a table column has no dictionary entry.

Run IDs (`batch_var`) are never a heatmap key, a model term, part of the confounding analysis or a PCA colour. When `batch_var` is set, the binding models' Stage 0 reports its association with each covariate (`design_associations`); it is never a model term.

All outputs go to `<outdir>/tables` and `<outdir>/figures`.

**Viral genomes and chrM.** Reads aligned to the viral genomes in `contigs.viral` (display name: accession in the BAMs; lengths from `contigs.viral_fasta` or the BAM headers) and to host contigs in `contigs.host_exclude` (default chrM) are counted per library but never enter the peak-based analyses. `work/annot/genome.sizes` holds only the host nuclear chromosomes; `genome_all.sizes` and `contigs.tsv` list every contig with its group. Library sizes for IP/control ratios (steps 7, 16 and `plot_chrom_coverage.R`) are reads on the host nuclear chromosomes, so viral load does not shift them. Viral reads are counted at MAPQ ≥ `contigs.viral_mapq` (default 10); host contigs are not MAPQ-filtered. `tables/qc07_viral_mapq_summary` shows the MAPQ distribution before that filter, so you can choose the cutoff. Each BAM is read once (one streaming pass, cached in `work/contigs/`; any cutoff is applied later from the cached MAPQ table).

**Pipeline contig tables.** Pipeline versions that remove chrM and viral reads before writing the filtered BAMs report them in `<merged_library>/contig_reads/<id>.contig_reads.tsv` (counted before host filtering, duplicates removed), with the duplicate-marked viral alignments in `contig_reads/viral_bam/<id>*.bam`. Step 0 finds both for every IP and control (columns `contig_table`, `viral_bam`, `control_contig_table`, `control_viral_bam` in `work/samples.tsv`; override per sample in the manifest). When a library has the table, its chrM and viral reads are the table's `reads_dedup`, and the library splits into parts that add up to `total_reads_dedup`: host nuclear (`mapped_reads − duplicates −` the listed contigs, split over the chromosomes in proportion to the filtered BAM), chrM, each viral genome, and unaligned (`total_reads − mapped_reads`). Shares in table00c, Figure P5 and the pies are then of all reads, as in the pipeline's report. With `contigs.viral_mapq` > 0 the viral reads are recounted from the viral BAM at that cutoff, and qc07 comes from the viral BAM; set it to 0 to use the pipeline's counts unchanged (they are not MAPQ-filtered, and reads placed by Allo keep low MAPQ). Libraries without a table fall back to the filtered BAM, with a warning.

**Sample overview** (first section of the report; made at the start of step 9)

- `table00_peak_summary`: one row per sample with its metadata, number of peaks, base pairs and % of the genome covered, peak width (min, quartiles, median, mean, max), fold enrichment (same statistics), peaks with fold enrichment ≥ 2 and with q < 0.05 / 0.01, and the chromosome with the most peaks.
- `table00b_peaks_by_chromosome`: sample × chromosome, for every chromosome of the genome: peak count, % of the sample's peaks, peaks per Mb, and log2(share of peaks / share of genome).
- Figures: `pk01` peaks per sample (count and Mb covered), `pk02` peaks per chromosome (count, coloured by density relative to chromosome length), `pk03` peak width and `pk04` fold enrichment distributions per sample. They are numbered P1–P4 in the report.
- `table00c_contig_reads`: reads per sample, library and contig (host chromosomes, chrM, viral genomes), share of reads, relative coverage and log2 IP/control; `qc07_viral_mapq_summary`: viral-read MAPQ distribution before filtering. Figure `pk05` (P5): share of reads, relative coverage and IP enrichment for chrM and each viral genome, per sample and pooled.
- Fold enrichment, -log10 p and -log10 q are columns 7–9 of the narrowPeak/broadPeak files (for MACS, fold enrichment at the summit for narrow peaks, over the region for broad peaks). They are kept per peak in `work/peaks/<sample>.info.tsv`.

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
  To plot the spacing histogram of any family (not just the 20 in `fig03`): `Rscript R/plot_spacing.R <outdir> F001 [--period] [--max-bp 60] [--by-sample] [--samples a,b]`; `--help` lists all options.

**Extra figures** (step 22, written to `<outdir>/extra/`, not in the report; settings under `extra_plots`)

- `extra/coverage/<sample>/<contig>_<norm>.png`: coverage along every host chromosome, chrM and viral genome for each normalisation in `extra_plots.coverage_norms` (raw, cpm, genome, input). Viral genomes use the viral MAPQ cutoff and are drawn from the pipeline's viral BAM when there is one (library sizes still from the filtered BAM). Contigs other than the host nuclear chromosomes are skipped when the BAM used has no reads on them (e.g. chrM removed by the pipeline); the log lists them. Binned counts are cached in `work/coverage/`, and plots newer than their BAM are skipped on reruns.
- `extra/pies/`: per sample, reads per contig (IP and control side by side) with and without the viral genomes, each with and without an inner ring showing contig-length shares, plus peaks per chromosome with and without the ring; pooled over all samples, coverage as RPKM shares (no ring, with and without viruses) and peaks (with and without the ring). Every pie has a TSV of its numbers.
- `extra/insert_sizes/<sample>.png`: insert-size distributions by CenSat class (formerly Figure QC3; made in step 15).

**Coverage along a chromosome** (standalone; does not need a finished run)

`extra/make_supp_flag_endogenous.R --crossrun <outdir> --flag <id> --endo <id> --labels-file <tsv> [--motif-labels-file <tsv>]` makes the FLAG-cGAS versus endogenous supplementary figure: A pies of peaks per chromosome for the matched pair (table00b, via `extra/redraw_pies.R`), B pairwise Jaccard similarity of peak sets (table05b) with the pair outlined, C motif family prevalence for chosen samples (via `extra/plot_motif_prevalence.R`), each as .png/.pdf/.tsv at its printed size, and `supp_figure.png/.pdf` with the panels placed 1:1. `extra/compile_qc_table.py --crossrun <outdir> --labels-file <tsv>` writes a supplementary QC table (TSV and Excel with a data dictionary): one row per chosen IP sample and per input it uses, with sequencing (FastQC, Trim Galore), alignment (Bowtie2, samtools flagstat, Picard, deepTools fingerprint), composition (chrM, viral, unaligned) and peak columns, no thresholds or flags; `_sources.tsv` names the file behind every value and `_missing.tsv` what was not found. `extra/plot_motif_prevalence.R` draws motif family prevalence (`tables/table02_motif_prevalence.tsv`) for chosen samples only, as dots (size = % of peaks with a site, filled = detected) or a heatmap, at a fixed point size: `--labels-file` (the same sample_id/label file as the chromosome heatmap) picks the samples and their order, families are the report's top-5-per-sample rule applied to those samples (`--top`, `--max-families`) or `--families F001,F007`, and `--highlight` shades chosen columns. `bin/plot_family_logos.py --samples a,b` likewise ranks the logo families on chosen samples only. `extra/plot_chrom_heatmap.R` redraws the report's peaks-per-chromosome heatmap (`tables/table00b_peaks_by_chromosome.tsv`) for a figure panel, at the size it will be printed so its text comes out at a fixed point size (`--width 5.75 --font 9`, Arial by default via `--family`). It measures the laid-out tiles with the real font: `--count-font auto` (default) picks the largest count size up to `--count-max` (6 pt) that fits a tile, `--x-angle auto` turns the chromosome names 90 degrees when they don't fit upright, and fixed sizes that don't fit are reported; `--labels-file` (columns sample_id, label) sets the row names and also which samples are drawn and their order, top to bottom as in the file (`--samples` overrides the selection), row names can also be given with `--labels a=Label,...`; chromosomes in genome order (1–22, X, Y). By default a bold Total column (each sample's peaks), Total row (each chromosome's peaks over the samples drawn) and grand total are added on grey tiles outside the colour scale (`--no-totals` to leave them out; `--total-font`); chromosome totals too wide for their tile are turned 90 degrees in a taller Total row. Chromosomes in table00b and in the report heatmap are now also in genome order.

`R/plot_chrom_coverage.R` plots binned read coverage over a whole chromosome, a viral genome or a region, for one or more libraries, with optional annotation tracks: an ideogram from a cytoband file and CenSat arrays. Coverage comes from BAM files (via `samtools bedcov`), bigWig or bedGraph. BAM tracks can be shown as mean depth, per million reads, fold over the genome-wide mean, or log2 IP/control. `--crossrun <outdir> --samples a,b` takes the BAMs, controls, peaks, annotations and chromosome aliases from a run; add `--viral-bam` to plot a viral genome from the pipeline's viral BAMs. With many bins (a whole chromosome), bins narrower than a pixel are drawn faint; `--px-per-bin 2` widens the figure so each bin gets at least 2 pixels (PNG up to 32,767 pixels wide; `--dpi` sets the resolution); without `--bin` it picks the smallest bin (at least 10 bp) that fits and draws the figure at the full 32,767 pixels. The x-axis always labels both ends of the plotted region (to 0.01 of the unit, e.g. 248.39 Mb), with round-number ticks in between. `--label-pos left|center|right` moves each track's sample label along the top of its panel. `--label-top 10` numbers the 10 highest-scoring peaks in the plotted region in a band under each track's peak ticks (never on the coverage) (`--top-by score|signal|p|q`) and writes them to `<out>_top_peaks.tsv`; `Rscript R/plot_peak_sequences.R --top-peaks <that TSV> --crossrun <outdir>` then draws the genome sequence in a 120 bp window around each summit as a coloured grid (A yellow, C green, G blue, T pink-red; 3560 px wide by default) and writes the sequences to a TSV. To bring out repeats: `--repeat-units` boxes every copy of each row's tandem repeat unit and names it; `--highlight TTCGNNNNANNCGGG` outlines motif matches on either strand; `--shared-kmer 8` fades everything except k-mers shared by at least half the rows; `--align` aligns the rows with MAFFT (reverse-complemented rows marked *) and fades bases that differ from their column's consensus. With `--align`, motif matches and repeat copies whose bases are split by alignment gaps are not boxed (logged; TSV columns `highlight_gapped`, `repeat_copies_drawn`). They can be combined. `Rscript R/plot_chrom_coverage.R --help` lists every option, for example:

```bash
Rscript R/plot_chrom_coverage.R --chrom chr9 --crossrun crossrun_out_human_broad \
    --samples THP1_HSV1_8h_IP-cGAS,THP1_IAV_6h_IP-cGAS --norm genome --bands chm13v2.0_cytobands_allchrs.bed --cap 0.999
```

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
- Covariate effects on binding are estimated by the binding models (steps 16–17, below). The earlier per-feature linear models, linear time trends, detection tests and PERMANOVA (`table06_metadata_effects`, `table06b`, `table06c`, `table06d`, `table07`, `table07b`, Figures 14–15) were removed; step 13 deletes their files from earlier runs.

## Binding models (steps 16–20) and bin-size tuning

Binding is the mean per-base fragment coverage of a fixed-size genomic bin, as log2 IP over input (both per million fragments, plus ε = the 1st percentile of non-zero input). Bins tile every host nuclear chromosome from position 0 (BED, 0-based half-open). The old peak-window regression models (tables 08–12, Figures 17–20, `models:` block) were removed; a config that still has a `models:` block gets a warning and the block is ignored.

**Bin-size tuning (run once).** `tune/run_tuning.sh config.yaml` (on the cluster: `sbatch tune/run_tuning.sbatch config.yaml`, bigmem) computes coverage at `tuning.base` (10 bp) on `tuning.chromosomes` for every IP and control library and sums it exactly to each size in `tuning.sizes`. For each size it measures, separately for CenSat bins and other bins:
- signal retained, R(w): the covariance of binding between the two halves of the reads (split by a hash of the read name; noise is independent between halves, so the covariance estimates real signal variance), relative to the smallest size;
- reproducibility: the split-half correlation;
- optionally, with `tuning.reference_motifs`, each reference motif's within-stratum effect on binding relative to the smallest size.

The chosen size is the largest one for which it and every smaller size keep the median ratio ≥ 1 − `tuning.tolerance` in both groups (and for the motif effects). The record `<tuning.outdir>/bin_size.tsv` is reused by every later run: set `bins.tuning_record` to it (or set `bins.size` directly) and rerun step 0. Outputs: `bin_size.tsv`, `tuning_curves.tsv`, `tuning_summary.tsv`, `tuning_motif_effects.tsv`, `tuning_curves.png`.

**Drift check.** With `bins.drift_check: true`, step 16 recomputes R(chosen)/R(floor) on `tuning.drift_chromosome` for the current samples and writes `work/bins/drift_check.tsv`; a median below 1 − tolerance prints a warning to retune.

**Bin table (step 16).** `work/bins/w<size>/`: `manifest.tsv`, `classes.tsv`, `arrays.tsv`, per chromosome raw little-endian columns (`class.int16`, `array_id.int32`, `edge_dist.int32`, `repeat_frac.float32`, `mappability.float32`, `gc.float32`, `n_frac.float32`), `sites.tsv.gz` (family, chrom, bin, number of FIMO sites whose midpoint is in the bin; FIMO 1-based closed → 0-based midpoint), and one directory per library with `all`, `uniq` (MAPQ ≥ `bins.unique_mapq`; 255 = placed uniquely by Bowtie2 in `-k` mode), `h1`, `h2` and `fl` columns plus `library.json`. Fragments: proper pairs counted once from read 1 (`-f 66 -F 3852`), length ≤ `bins.insert_max`; single-end reads as reads (`-F 3844`). Annotation, FIMO and signal are cached and recomputed only when their inputs change. Step 16 needs the per-base mappability track from step 14 (`qc.write_mappability_track: true`) or `bins.mappability_bedgraph`.

**Stage 0 (step 17).** `R/65_design.R` builds the sample-level design from `categorical_covariates` and `time_var` binned at `binding_models.time_bins`, for samples with both an IP and a control BAM: `tables/design_samples.tsv`, `design_estimability.tsv`, `design_associations.tsv`, `work/models/estimable_terms.tsv`, `work/models/design_matrix.tsv`.

Step numbers changed when the models were added: the report is now step 21 (was 18) and the extra figures step 22 (was 19).

**Settings.** Everything is under `binding_models` in `config.yaml`; the bin size comes from `bins.size` or `bins.tuning_record` (step 0 writes it to `work/config.env`), so the models always run on the bin table of that size.
- `terms`: `auto` (every level step 17 marks estimable) or a list of covariates, e.g. `[virus, cgas_construct]`, which uses their estimable levels. A covariate with no estimable level (aliased, or one level only) stops steps 19–20 with step 17's reason.
- `families`: `detected` (detected in at least one sample in table02, at most `max_families`, ordered by samples detected then median % of peaks), `all`, or a list of family ids.
- `count_transform`: motif term per bin, `log2p1` = log2(1 + sites) (default), `raw`, `presence`.
- Bin filter: N fraction ≤ `max_n_frac`, finite GC and mappability, not in `work/annot/exclude.bed` (blacklist), inside `include.bed` if given (bin midpoint), and input per million fragments ≥ the `input_min_quantile` quantile of non-zero input. Counts at every filter are in `models_level1_info`.
- Strata: GC in `gc_bins` quantile bins × mappability cut at `mappability_breaks`. Clusters for standard errors: satellite arrays, and `nonarray_block` bp blocks elsewhere; a fit with fewer than `min_clusters` clusters gets heteroskedasticity-robust standard errors instead.
- `seed` (default 7) fixes the validation folds and the similarity subsample and permutations.

**Level 1 (step 18).** For each modelled sample and read version (`arms`: `all`, `uniq`), `R/70_level1.R` regresses binding on every family's motif term within strata (fixest `feols`):
- `per_class`: one fit per repeat class with at least `min_bins_per_class` bins;
- `pooled_adjusted`: all classes, strata = class × GC × mappability (bins compared within their class);
- `pooled_unadjusted`: all classes, strata = GC × mappability only. The difference to `pooled_adjusted` shows how much of a motif's apparent effect is its prevalence in highly bound repeats.
- Held-out arrays (`validation_arms`): in CenSat classes with at least `validation_folds` arrays, arrays are split into folds; `delta_mse` is the held-out error without the family minus with it.

With `level1_mode: slurm` (default), step 18 writes `work/models/level1_tasks_todo.tsv`, prints `sbatch --array=1-N models/level1.sbatch config.yaml <todo>` and stops (exit 2). Submit that from the bundle directory, wait for the array to finish, then rerun `./run_all.sh config.yaml 18`: finished tasks are skipped and the run carries on. A task reruns when its output is older than the bin table, the design matrix or the config. Size `models/level1.sbatch` (`--cpus-per-task`, `--mem`) from the benchmark below and set `level1_threads` to the CPUs. `level1_mode: inline` runs the tasks inside step 18 (small data, tests).

**Levels 2–3 (step 19).** `R/72_levels23.R` combines each family's per-sample effects with a random-effects meta-analysis (metafor `rma`, REML; fixed effect when k < 3), and regresses them on the design terms (meta-regression). Effects are then shrunk with ashr (normal mixture) within read version × model (× term), giving `shrunk`, `lfsr` and a rank. Tables `models_motif_effects`, `models_condition_effects`, `models_before_after`, `models_read_version_agreement`, `models_validation(_by_sample)`; figures mod01–mod05.

**Similarity (step 20).** `R/74_similarity.R` correlates bin-level binding between every pair of samples (bins passing the filters in all samples, CenSat and other bins separately, at most `similarity_max_bins`), scales r by each sample's split-half ceiling (Spearman–Brown), and tests the design terms on 1 − scaled r with PERMANOVA (vegan `adonis2`, marginal, `permutations`) and a one-way ICC per covariate. Figures mod06–mod07.

**Not implemented.** The fragment-length (`fl`) sensitivity version is computed by step 16 but not modelled.

**Tests.** `tests/run_model_test.sh` builds synthetic bins with planted effects (`tests/make_model_test_data.py`) and checks that the models recover them: a motif effect and its change with a condition, a motif acting only outside satellites, a null motif, and a motif that looks strong only because its sites sit in satellites (removed by `pooled_adjusted`). Effects are recovered at about 80% of the planted value because of ε in the binding definition (see above).

**Level-1 benchmark.** `sbatch tests/bench_level1.sbatch` times the within-stratum regression (fixest) at 50 families × 15 classes on synthetic bins and reports peak memory per bin.

## MEME Suite 5.5.4 and later

From version 5.5.4, FIMO (and SEA from 5.5.5) read FASTA headers of the form `name::chrom:start-end`, as written by `bedtools getfasta -name`, as genomic coordinates: they report the chromosome as the sequence name and genome positions as site positions. The bundle therefore writes headers that contain only the window name, passes `--no-pgc` to FIMO and SEA when the installed version supports it, and stops with an error if FIMO reports names that are not the scanned windows. Outputs of steps 03, 05 and 16 carry a format stamp, so a bundle update that changes sequence naming reruns them automatically. Results produced by earlier versions of the bundle with MEME 5.5.4+ (motif prevalence, spacing, sites in repeats, step 16–17 motif counts) should be regenerated from step 3.

## Motif discovery with control-library sequences

`motifs.min_width` (default 6) sets the shortest motif discovery may report (`--minw` for STREME, MEME and XSTREME). Motifs shorter than it from the pipeline's own XSTREME runs (`discovery: pipeline`) are dropped when the motifs are pooled. Changing it reruns discovery for every sample. Short motifs are worth excluding when they cannot reach `fimo_thresh`: a 6 bp motif's best possible FIMO p-value is usually above 1e-4, so it can never be scanned (`work/motifs/family_best_p.tsv` lists each family's best possible p-value).

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

Step 02 first runs the chosen engine's tools on a tiny input and stops with a clear message if they fail. If MEME itself is broken, fix the install or set `meme_prefix` (e.g. `singularity exec /path/to/meme.sif`); the prefix is then used for every MEME tool, and paths must be visible inside the container. Changing discovery changes the motifs, and steps 04–05 detect that and recompute families and scans. Set `motifs.discovery: pipeline` to reuse the pipeline's results instead.

## Quality control

- **Mappability.** GenMap computes, for every 150-mer (`qc.mappability_k`), whether it occurs once in the genome (both strands, no mismatches). Each base's mappability is the fraction of the k-mers covering it that are single-copy. These per-base values are summarised over every base of each CenSat class (`qc02_censat_mappability_stats`) and over every peak. Peak mappability is related to peak score and CenSat class by linear models (`qc01_peak_mappability_models`; per-class trends in `qc01c_mappability_score_by_censat`). QC figure 1 shows mappability against score per class, with a trend line only where significant. QC figure 2 shows violins of peak mappability per class, against each class's genome-wide median. Indexing the genome is a one-off (about an hour for a human genome); keep the index via `qc.genmap_index`, or skip GenMap entirely by supplying a k-mer bedGraph in `qc.mappability_bedgraph`. Processing needs about 3 GB of RAM for the largest chromosome. The per-base track is written to `work/qc/mappability_k150.bedgraph.gz`.
- **Insert sizes.** For every IP and control library, fragments whose midpoint lies in each CenSat class are counted by length, plus a random sample of fragments outside CenSat (`qc.insert_genome_fraction`). The output is one figure per IP sample comparing it with its control (QC figure 3), a summary of median shifts (QC figure 4), and `qc03_insert_size_stats`. Single-end libraries are reported and skipped.
- **Strand balance and read stacking** (QC figure 5, `qc05_peak_strand_stacking_by_censat`). For every peak, read 1 start positions give strand imbalance, start-site diversity relative to random placement, and the share of reads at the most-used start, compared between CenSat classes and peaks outside CenSat. Strand-biased, stacked peaks in repeats are likely alignment artifacts.
- **Position within arrays** (QC figure 6, `qc06_array_edge_distance_by_censat`). Distance from each summit to the nearest edge of its CenSat array, compared with random positions in the same arrays by Kolmogorov–Smirnov test. Pile-ups at array edges are a hallmark of mapping artifacts at unique/repeat boundaries.
- **Reuse.** All QC results reuse step 01's peaks and annotation and step 06's peak CenSat labels.

## Utilities

**Peak pies on their own.** `sbatch extra/peak_pies.sbatch config.yaml` (defq, 1 core, 4 GB) runs steps 0–1 and `Rscript R/plot_pies.R config.yaml --peaks-only`: pie charts of peaks per host chromosome for every sample and pooled, with and without the chromosome-length ring, plus `extra/pies/peaks_per_chromosome.tsv`. chrM and the viral genomes are not peak-calling contigs; step 1 leaves any peaks on them out and logs a warning. The same pies are drawn by step 22 in a full run.

**One-off pie comparisons.** `extra/compare_pies.R config.yaml comparisons.tsv` (or `sbatch extra/compare_pies.sbatch config.yaml comparisons.tsv`) draws pies side by side for hand-defined groups of samples, e.g. `extra/comparisons_cgas.tsv` (A549 vs THP1, THP1 Flag-cGAS vs endogenous, THP1 IAV vs HSV-1). Each row is `comparison`, `group`, `filter` (`column=value`, `;` = and, `|` = or). Output in `extra/pies/compare/`: peaks per host chromosome, IP coverage as RPKM shares (when step 9 has run) and a members table. Groups are pooled by summing; `--equal-weight` averages per-sample shares instead. Not part of `run_all.sh`.

**Motif family logos.** `python3 bin/plot_family_logos.py --crossrun <outdir> --top 5` draws the representative motif of the top families as sequence logos (bits, A green, C blue, G yellow, T red, numbered pink circles, as in the cGAS figure). Families are ranked from `table02_motif_prevalence` by the number of samples in which they are detected (`--rank-by detected`, the default; or `prevalence`, `denovo`, `sea`), optionally on `--samples` only; `--families F001,F007` draws chosen families. A TSV next to the figure lists the families and their statistics. Needs matplotlib and logomaker.

**Pies for a figure panel.** `Rscript extra/redraw_pies.R <pie>.tsv --ncol 4 --width 8 --height 1.75 --font 4.5` redraws any pie TSV written by step 22 or `extra/compare_pies.R` without a legend: slices of at least `--inside-min` % (default 6) are labelled inside, smaller ones outside with the chromosome and percentage, stacked so they never overlap.

**Pies pooled by metadata.** `extra_plots.pie_group_by: [virus, cell_type]` in the config (step 22), or `--group-by virus,cell_type` on `R/plot_pies.R` (second argument of `extra/peak_pies.sbatch`), adds one figure per manifest column with one pie per level, labelled with the number of samples: `peaks_by_<column>_{ring,noring}.png` (peaks summed over the level's samples) and, outside `--peaks-only`, `coverage_rpkm_by_<column>_{host,with_viruses}.png` (IP reads pooled over the level's samples, as RPKM shares; controls are left out because one input serves several samples). Samples with no value are pooled as `NA` and listed in the log; an unknown column stops with the list of manifest columns.

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

This builds a synthetic genome with planted satellites, a CENP-B-box array, and ETS sites, plus 10 fake pipeline runs. Satellite enrichment rises with time in the WT-construct series only. The whole bundle then runs on it, with stub `fimo`/`tomtom`/`sea` that do exact consensus matching. The stubs exercise the plumbing only; they are not substitutes for MEME Suite. Before the main run it runs the bin-size tuning at test scale (chromosomes chr1–chr2, sizes 50–1000 bp) and then checks the bin table (coverage halves add up, unique ≤ all, GC in 0–1, CenSat bins classed), the tuning record, the drift check and the Stage 0 design tables. `tests/test_retention.R` checks the tuning measure on a planted 400 bp signal.
