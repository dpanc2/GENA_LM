#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${GENA_PYTHON:-python}"
REPO_ROOT="${GENA_REPO_ROOT:-$(cd "$SCRIPT_DIR/../../.." && pwd)}"
ASSET_ROOT="${GENA_ASSET_ROOT:-/home/jovyan/dpanc/benchmarking/GENA_LM}"
KOMAREK_ROOT="${KOMAREK_ROOT:-/home/jovyan/dpanc/komarek}"
GENOME="${GENA_GENOME:-$KOMAREK_ROOT/data/T2Thuman_AsteI2_plm_merged_ref.fa}"
CONDITION="${GENA_CONDITION:-HEK293}"
OUTPUT_DIR="${GENA_OUTPUT_DIR:-$KOMAREK_ROOT/atac_center50_${CONDITION}_predictions}"
LIMIT_ARGS=()

if [[ -n "${GENA_LIMIT:-}" ]]; then
    LIMIT_ARGS=(--limit "$GENA_LIMIT")
fi

SCRIPT="$REPO_ROOT/downstream_tasks/expression_prediction/atac_seq_predictions/predict_atac_center50.py"

run_one() {
    local bed="$1"
    local name="$2"
    if [[ ! -f "$bed" ]]; then
        echo "Skipping missing BED: $bed" >&2
        return
    fi
    "$PYTHON" "$SCRIPT" \
        --bed "$bed" \
        --output "$OUTPUT_DIR/${name}.center50.tsv" \
        --repo-root "$REPO_ROOT" \
        --asset-root "$ASSET_ROOT" \
        --genome "$GENOME" \
        --condition "$CONDITION" \
        "${LIMIT_ARGS[@]}"
}

mkdir -p "$OUTPUT_DIR"

# run_one "/home/jovyan/dpanc/komarek/data/motif_test.random_20000_intergenic_10000bp.bed" "intergenic"
# run_one "/home/jovyan/dpanc/komarek/data/Astel2_random_20000_intervals_100000bp.strand.sorted.bed" "mosquito"
# run_one "/home/dpanc/komarek/data/protein_coding_promoters_m5-5Kb.bed" "promoters"

run_one "$KOMAREK_ROOT/MPRA/log2ratio_results/komar_MPRA_10kb_plus_strand.bed" "MPRA_plus_strand"
run_one "$KOMAREK_ROOT/MPRA/log2ratio_results/komar_MPRA_10kb_minus_strand.bed" "MPRA_minus_strand"

