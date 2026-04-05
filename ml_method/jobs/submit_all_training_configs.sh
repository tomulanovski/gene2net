#!/bin/bash
# ==============================================================================
# Submit SimPhy training data generation for ALL configurations.
#
# Runs SimPhy on the 2K MUL-trees for each ILS/dup-loss config.
# After SimPhy completes, run ASTRAL + packaging for each config.
#
# USAGE:
#   ./submit_all_training_configs.sh MUL_TREES_DIR [STEP]
#
# STEPS:
#   simphy   - Submit SimPhy jobs (run first)
#   astral   - Submit ASTRAL jobs (after SimPhy completes)
#   package  - Submit packaging jobs (after ASTRAL completes)
#   all      - Submit all three (with dependencies — NOT recommended,
#              better to verify each step)
#
# EXAMPLES:
#   ./submit_all_training_configs.sh /path/to/mul_trees simphy
#   ./submit_all_training_configs.sh /path/to/mul_trees astral
#   ./submit_all_training_configs.sh /path/to/mul_trees package
# ==============================================================================

MUL_TREES_DIR=${1:?Usage: submit_all_training_configs.sh MUL_TREES_DIR [STEP]}
STEP=${2:-simphy}

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/ml_method"
N_TREES=2000
START=0

# ============================================================================
# Configuration table: CONFIG_NAME  NE  DUP_RATE  LOSS_RATE
# ils_low (Ne=200K, no dup/loss) is already done
# ============================================================================

declare -a CONFIGS=(
    # ILS only (no dup/loss)
    "ils_medium        1000000  0     0"
    "ils_high          2000000  0     0"
    # Low ILS + dup/loss
    "dup_loss_low      200000   1e-9  1e-9"
    "dup_loss_medium   200000   1e-8  1e-8"
    "dup_loss_high     200000   1e-7  1e-7"
    # Medium ILS + dup/loss
    "dup_loss_low_ne1M      1000000  1e-9  1e-9"
    "dup_loss_medium_ne1M   1000000  1e-8  1e-8"
    "dup_loss_high_ne1M     1000000  1e-7  1e-7"
)

echo "============================================================"
echo "Training Data Generation — All Configs"
echo "============================================================"
echo "MUL-trees dir: ${MUL_TREES_DIR}"
echo "Trees: ${N_TREES} (indices ${START}-$((START + N_TREES - 1)))"
echo "Step: ${STEP}"
echo "Configs: ${#CONFIGS[@]} (+ ils_low already done)"
echo "============================================================"
echo ""

for config_line in "${CONFIGS[@]}"; do
    read -r CONFIG NE DUP LOSS <<< "$config_line"

    echo "--- ${CONFIG} (Ne=${NE}, dup=${DUP}, loss=${LOSS}) ---"

    case "$STEP" in
        simphy)
            "${BASE_DIR}/jobs/submit_simphy_training.sh" \
                "$N_TREES" "$MUL_TREES_DIR" "$START" \
                "$CONFIG" "$NE" 500 1 "$DUP" "$LOSS"
            ;;
        astral)
            "${BASE_DIR}/jobs/submit_astral_training.sh" \
                "$N_TREES" "$MUL_TREES_DIR" "$START" "$CONFIG"
            ;;
        package)
            "${BASE_DIR}/jobs/submit_package_training.sh" \
                "$N_TREES" "$MUL_TREES_DIR" "$START" "$CONFIG"
            ;;
        *)
            echo "Unknown step: ${STEP}. Use: simphy, astral, package"
            exit 1
            ;;
    esac

    echo ""
done

echo "============================================================"
echo "Submitted ${STEP} for all configs."
echo ""
echo "Next steps:"
echo "  1. Wait for ${STEP} jobs to complete (squeue -u \$USER)"
echo "  2. Run the next step:"
if [ "$STEP" = "simphy" ]; then
    echo "     ./submit_all_training_configs.sh ${MUL_TREES_DIR} astral"
elif [ "$STEP" = "astral" ]; then
    echo "     ./submit_all_training_configs.sh ${MUL_TREES_DIR} package"
elif [ "$STEP" = "package" ]; then
    echo "  Training data ready! Update training script to load all configs."
fi
echo "============================================================"
