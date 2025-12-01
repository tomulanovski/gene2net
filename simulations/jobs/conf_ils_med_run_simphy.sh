#!/bin/bash
export SIMPHY_CONFIG="ils_medium_10M"
export SIMPHY_NE=1000000
export SIMPHY_DUP=0
export SIMPHY_LOSS=0
sbatch "/groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs/simphy_reusable.sh"