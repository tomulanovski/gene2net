#!/bin/bash
#SBATCH --job-name=diploidize
#SBATCH --array=1-21
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/diploidize_%A_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/diploidize_%A_%a.err
#SBATCH --time=2:00:00
#SBATCH --mem=4g
#SBATCH --cpus-per-task=1
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

# ============================================================================
# SUBMIT_DIPLOIDIZATION.SH
# ============================================================================
# Applies post-WGD diploidization (fractionation) to an existing config's gene
# trees + alignments, writing a new diploidized config.  Runs as a 21-task array
# (one network per task); each task loops over the replicates.  Must run AFTER
# AliSim (gene trees + alignments exist) and BEFORE the method-prep scripts.
#
# Usage (defaults: beta retention, mean 0.5, concentration 4, 5 replicates):
#   sbatch --export=CONFIG=conf_ils_low_10M submit_diploidization.sh
#
# Sweep mean retention into separate configs:
#   sbatch --export=CONFIG=conf_ils_low_10M,MEAN=0.3 submit_diploidization.sh
#   sbatch --export=CONFIG=conf_ils_low_10M,MEAN=0.7 submit_diploidization.sh
#
# Other distributions:
#   sbatch --export=CONFIG=conf_ils_low_10M,DIST=fixed,VALUE=0.4 submit_diploidization.sh
#   sbatch --export=CONFIG=conf_ils_low_10M,DIST=uniform,LOW=0.2,HIGH=0.8 submit_diploidization.sh
#
# Environment variables:
#   CONFIG          - input config name (required, e.g. conf_ils_low_10M)
#   OUT_CONFIG      - output config name (default: {CONFIG}_{tag}, tag from params)
#   NUM_REPLICATES  - replicates to process (default: 5)
#   SEED            - base random seed (default: 42)
#   DIST            - beta | uniform | fixed (default: beta)
#   MEAN, CONCENTRATION  - beta params (default MEAN=0.5, CONCENTRATION=4)
#   ALPHA, BETA          - beta params (override MEAN)
#   LOW, HIGH            - uniform params
#   VALUE                - fixed retention rate
# ============================================================================

set -eo pipefail
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

# ---------------------------------------------------------------------------
# Config / params
# ---------------------------------------------------------------------------
CONFIG="${CONFIG:?ERROR: CONFIG is required. Use --export=CONFIG=conf_name}"
NUM_REPLICATES="${NUM_REPLICATES:-5}"
SEED="${SEED:-42}"
DIST="${DIST:-beta}"

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"
CONDA_PATH="/groups/itay_mayrose/tomulanovski/miniconda3"
PYTHON_SCRIPT="/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/apply_diploidization.py"

networks=(
  "Bendiksby_2011" "Koenen_2020" "Brysting_2007" "Lawrence_2016"
  "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Liang_2019" "Popp_2005" "Wu_2015"
  "Liu_2023" "Ren_2024" "Marcussen_2011" "Marcussen_2012" "Sessa_2012b" "Zhao_2021"
  "Hori_2014" "Marcussen_2015" "Shahrestani_2015" "Morales-Briones_2021" "Soza_2014"
)

network_idx=$((SLURM_ARRAY_TASK_ID - 1))
network="${networks[$network_idx]}"

# Build distribution args for the Python CLI
DIST_ARGS=(--dist "$DIST")
case "$DIST" in
  beta)
    if [ -n "${ALPHA:-}" ] && [ -n "${BETA:-}" ]; then
      DIST_ARGS+=(--alpha "$ALPHA" --beta "$BETA")
    else
      DIST_ARGS+=(--mean "${MEAN:-0.5}" --concentration "${CONCENTRATION:-4}")
    fi
    ;;
  uniform)
    DIST_ARGS+=(--low "${LOW:?LOW required for DIST=uniform}" --high "${HIGH:?HIGH required for DIST=uniform}")
    ;;
  fixed)
    DIST_ARGS+=(--value "${VALUE:?VALUE required for DIST=fixed}")
    ;;
  *)
    echo "ERROR: unknown DIST=$DIST"; exit 1 ;;
esac

# Replicate list 1..NUM_REPLICATES
REPLICATES=()
for r in $(seq 1 "$NUM_REPLICATES"); do REPLICATES+=("$r"); done

OUT_CONFIG_ARGS=()
[ -n "${OUT_CONFIG:-}" ] && OUT_CONFIG_ARGS=(--out-config "$OUT_CONFIG")

echo "============================================================================"
echo "DIPLOIDIZATION - ${network}"
echo "Config: ${CONFIG}   Replicates: ${NUM_REPLICATES}   Seed: ${SEED}"
echo "Dist args: ${DIST_ARGS[*]}"
echo "Date: $(date)"
echo "============================================================================"

if [ ! -d "${BASE_DIR}/${network}/data/${CONFIG}" ]; then
    echo "ERROR: data dir not found: ${BASE_DIR}/${network}/data/${CONFIG}"
    exit 1
fi

# ---------------------------------------------------------------------------
# Activate conda + run
# ---------------------------------------------------------------------------
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate gene2net || { echo "ERROR: could not activate gene2net"; exit 1; }

python "$PYTHON_SCRIPT" "$CONFIG" \
    --base-dir "$BASE_DIR" \
    --networks "$network" \
    --replicates "${REPLICATES[@]}" \
    --seed "$SEED" \
    "${OUT_CONFIG_ARGS[@]}" \
    "${DIST_ARGS[@]}"

echo "Done: ${network}"
