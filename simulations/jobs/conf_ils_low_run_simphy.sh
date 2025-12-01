#!/bin/bash
export SIMPHY_CONFIG="ils_low_10M"
export SIMPHY_NE=200000
export SIMPHY_DUP=0
export SIMPHY_LOSS=0
sbatch "/groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs/simphy_reusable.sh"