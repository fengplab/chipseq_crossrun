#!/usr/bin/env bash
# Smoke test on synthetic data with stub fimo/tomtom/sea (exact-consensus matching; NOT the real tools).
# Needs bedtools, python3 (pyyaml), R with the packages listed in environment.yml.
set -euo pipefail
T="$(cd "$(dirname "$0")" && pwd)"; ROOT="$T/.."
D="${1:-$T/test_data}"
python3 "$T/make_test_data.py" "$D"
sed -e "s#manifest: manifest.tsv#manifest: $D/manifest.tsv#" -e "s#outdir: crossrun_out#outdir: $D/out#" \
    -e "s#genome_build: T2T-CHM13v2.0#genome_build: test#" -e "s#/path/to/genome.fa.sizes#$D/genome.sizes#" \
    -e "s#/path/to/genome.fa#$D/genome.fa#" -e "s#include_bed: null#include_bed: $D/include.bed#" \
    -e "s#/path/to/repeatmasker.bb#$D/rmsk.bed#" -e "s#/path/to/censat.bb#$D/censat.bed#" \
    -e "s#genes_bed: null #genes_bed: $D/genes.bed#" -e "s#  - /path/to/JASPAR2024_CORE_vertebrates_non-redundant.meme#  - $D/known.meme#" \
    -e "/H12CORE_meme_format/d" -e "s#n_shuffles: 200#n_shuffles: 20#" "$ROOT/config.yaml" > "$D/config.yaml"
PATH="$T/fakebin:$PATH" "$ROOT/run_all.sh" "$D/config.yaml"
