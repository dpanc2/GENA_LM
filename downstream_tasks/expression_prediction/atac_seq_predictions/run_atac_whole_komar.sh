#!/usr/bin/env bash
set -euo pipefail

GPU_ID="${1:-0}"
CONDITION="${2:-HCT116}"

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${GENA_PYTHON:-python}"
REPO_ROOT="${GENA_REPO_ROOT:-$(cd "$SCRIPT_DIR/../../.." && pwd)}"
ASSET_ROOT="${GENA_ASSET_ROOT:-/home/jovyan/dpanc/benchmarking/GENA_LM}"
KOMAREK_ROOT="${KOMAREK_ROOT:-/home/jovyan/dpanc/komarek}"
SCRIPT="$REPO_ROOT/downstream_tasks/expression_prediction/atac_seq_predictions/predict_atac_whole_komar.py"

BED="${GENA_BED:-$KOMAREK_ROOT/data/AsteI2_win10000_slide1000_whole_genome.bed}"
GENOME="${GENA_GENOME:-$KOMAREK_ROOT/data/AsteI2_V4.clean.fa}"
CHROM_SIZES="${GENA_CHROM_SIZES:-$KOMAREK_ROOT/data/AsteI2.clean.chrom.sizes}"

OUTPUT_DIR="${GENA_OUTPUT_DIR:-$KOMAREK_ROOT/atac_whole_${CONDITION}}"
BIGWIG="$OUTPUT_DIR/whole_AsteI2_${CONDITION}_ATAC.bw"

LIMIT_ARGS=()
if [[ -n "${GENA_LIMIT:-}" ]]; then
    LIMIT_ARGS=(--limit "$GENA_LIMIT")
fi

mkdir -p "$OUTPUT_DIR"

echo "GPU: $GPU_ID"
echo "Condition: $CONDITION"
echo "BED: $BED"
echo "bigWig: $BIGWIG"

CUDA_VISIBLE_DEVICES="$GPU_ID" "$PYTHON" "$SCRIPT" \
    --bed "$BED" \
    --bigwig "$BIGWIG" \
    --chrom-sizes "$CHROM_SIZES" \
    --repo-root "$REPO_ROOT" \
    --asset-root "$ASSET_ROOT" \
    --genome "$GENOME" \
    --condition "$CONDITION" \
    --device cuda:0 \
    --sequence-size 10000 \
    --center-bp 1000 \
    --chunk-size 20000 \
    "${LIMIT_ARGS[@]}"

echo "Finished: $BIGWIG"
