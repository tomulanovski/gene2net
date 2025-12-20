#!/bin/bash

# ============================================================================
# DEBUG COMPARISON ERRORS
# ============================================================================
# Investigates why network comparisons fail with "division by zero"
#
# Usage:
#   ./debug_comparison_error.sh <network> <config> <method>
#
# Example:
#   ./debug_comparison_error.sh Diaz-Perez_2018 conf_ils_low_10M polyphest_p50
# ============================================================================

if [ $# -lt 3 ]; then
  echo "Usage: $0 <network> <config> <method>"
  echo "Example: $0 Diaz-Perez_2018 conf_ils_low_10M polyphest_p50"
  exit 1
fi

NETWORK="$1"
CONFIG="$2"
METHOD="$3"

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"
NETWORKS_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/networks"

# Determine result filename based on method
if [ "$METHOD" = "grampa" ]; then
  RESULT_FILE="grampa_result.tre"
elif [[ "$METHOD" == polyphest* ]]; then
  RESULT_FILE="polyphest_result.tre"
elif [ "$METHOD" = "padre" ]; then
  RESULT_FILE="padre_result.tre"
elif [ "$METHOD" = "mpsugar" ]; then
  RESULT_FILE="mpsugar_result.tre"
else
  echo "Unknown method: $METHOD"
  exit 1
fi

echo "========================================================================"
echo "DEBUGGING: $NETWORK / $CONFIG / $METHOD"
echo "========================================================================"
echo ""

# ============================================================================
# CHECK GROUND TRUTH
# ============================================================================

GT_FILE="${NETWORKS_DIR}/${NETWORK}.tre"

echo "--- GROUND TRUTH NETWORK ---"
echo "File: $GT_FILE"
echo ""

if [ ! -f "$GT_FILE" ]; then
  echo "ERROR: Ground truth file not found!"
  exit 1
fi

echo "Content:"
cat "$GT_FILE"
echo ""

# Count reticulations (# symbols)
ret_count=$(grep -o '#' "$GT_FILE" | wc -l)
echo "Number of reticulations (# symbols): $ret_count"

# Count species (taxa)
# Assuming Newick format - count commas + 1
species_count=$(($(grep -o ',' "$GT_FILE" | wc -l) + 1))
echo "Estimated number of species: $species_count"

echo ""

# ============================================================================
# CHECK INFERRED NETWORKS
# ============================================================================

echo "--- INFERRED NETWORKS ---"
echo ""

for rep in {1..5}; do
  INFERRED="${BASE_DIR}/${NETWORK}/results/${CONFIG}/${METHOD}/replicate_${rep}/${RESULT_FILE}"

  echo "Replicate $rep:"
  echo "  File: $INFERRED"

  if [ ! -f "$INFERRED" ]; then
    echo "  Status: FILE NOT FOUND"
    echo ""
    continue
  fi

  # Get file size
  size=$(stat -c%s "$INFERRED" 2>/dev/null || stat -f%z "$INFERRED" 2>/dev/null)
  echo "  Size: $size bytes"

  # Show content
  echo "  Content:"
  cat "$INFERRED" | sed 's/^/    /'
  echo ""

  # Count reticulations
  ret_count=$(grep -o '#' "$INFERRED" | wc -l)
  echo "  Reticulations (# symbols): $ret_count"

  # Count species
  species_count=$(($(grep -o ',' "$INFERRED" | wc -l) + 1))
  echo "  Estimated species: $species_count"

  # Check if it's an empty or trivial network
  if [ "$size" -lt 50 ]; then
    echo "  WARNING: Very small file - might be empty or incomplete!"
  fi

  if [ "$ret_count" -eq 0 ]; then
    echo "  WARNING: No reticulations - this is a plain tree, not a network!"
  fi

  echo ""
done

# ============================================================================
# TEST COMPARISON MANUALLY
# ============================================================================

echo "========================================================================"
echo "TESTING COMPARISON"
echo "========================================================================"
echo ""

echo "Trying to compare replicate 1 manually..."
echo ""

INFERRED="${BASE_DIR}/${NETWORK}/results/${CONFIG}/${METHOD}/replicate_1/${RESULT_FILE}"

if [ -f "$INFERRED" ]; then
  # Try to run the comparison
  cd /groups/itay_mayrose/tomulanovski/gene2net

  python -c "
from reticulate_tree import ReticulateTree
from simulations.scripts.compare_reticulations import pairwise_compare

# Load trees
print('Loading ground truth...')
try:
    gt = ReticulateTree('$GT_FILE')
    print(f'  Loaded: {gt.num_species} species')
except Exception as e:
    print(f'  ERROR: {e}')
    exit(1)

print('')
print('Loading inferred network...')
try:
    inferred = ReticulateTree('$INFERRED')
    print(f'  Loaded: {inferred.num_species} species')
except Exception as e:
    print(f'  ERROR: {e}')
    exit(1)

print('')
print('Attempting comparison...')
try:
    result = pairwise_compare(gt, inferred)
    print('  SUCCESS!')
    print('')
    print('Metrics:')
    for key, value in result.items():
        print(f'  {key}: {value}')
except ZeroDivisionError as e:
    print(f'  DIVISION BY ZERO ERROR!')
    print(f'  Details: {e}')
    print('')
    print('This usually means one of the networks has:')
    print('  - Zero reticulations')
    print('  - Zero polyploid species')
    print('  - Some other zero-count property')
    print('')
    print('Let me check the network properties...')
    print(f'  GT reticulations: {len(gt.reticulations) if hasattr(gt, \"reticulations\") else \"unknown\"}')
    print(f'  Inferred reticulations: {len(inferred.reticulations) if hasattr(inferred, \"reticulations\") else \"unknown\"}')
except Exception as e:
    print(f'  OTHER ERROR: {e}')
    import traceback
    traceback.print_exc()
"
else
  echo "Cannot test - inferred file not found: $INFERRED"
fi

echo ""
echo "========================================================================"
echo "RECOMMENDATIONS"
echo "========================================================================"
echo ""

echo "If the error is 'division by zero':"
echo "  1. Check if the ground truth has zero reticulations"
echo "  2. Check if the inferred network is empty or has zero reticulations"
echo "  3. Look at the comparison code to see which metric is dividing by zero"
echo ""
echo "Common causes:"
echo "  - Polyphest produced an empty or trivial network"
echo "  - Ground truth is actually a tree (no reticulations)"
echo "  - Species mismatch between ground truth and inferred"
echo ""
echo "To fix:"
echo "  - If Polyphest output is bad, re-run with different parameters"
echo "  - If ground truth is a tree, this comparison might not be valid"
echo "  - Check the compare_reticulations.py code for the division"
echo ""
