# CLAUDE.md

This file provides guidance to Claude Code when working with this repository.

## Overview

`gene2net` is a research pipeline for phylogenetic network inference from gene trees. It includes:
- Scripts for data preparation and format conversion
- SLURM job scripts for cluster computing
- Simulation pipelines using SimPhy
- Network reconstruction methods: GRAMPA, Polyphest, MPSUGAR, PADRE, ASTRAL

## Repository Structure

```
gene2net/
├── scripts/           # Data processing and format conversion
├── jobs/              # SLURM job scripts for analyses
└── simulations/
    ├── scripts/       # Simulation analysis tools
    ├── jobs/          # Simulation job scripts
    ├── networks/      # Network files (21 networks in Newick format)
    └── PIPELINE_STATUS_GUIDE.md  # Complete pipeline validation guide
```

## Key Tools

### Scripts Directory
- **`fix_substrings_grampa.py`** - Critical preprocessing for GRAMPA (fixes substring label issues)
- Format conversion tools for different phylogenetic methods
- Taxa mapping generators for MPAllopp/MPSUGAR

### Simulations Directory
- **`check_pipeline_status.py`** - Validates SimPhy simulations and all method inputs/outputs
- **`compare_nets.py`** - Network comparison pipeline (handles MUL-trees and networks)
- **Full documentation:** `simulations/PIPELINE_STATUS_GUIDE.md`

## Common Workflows

### Running Simulations

```bash
# 1. Configure and submit SimPhy
sbatch simulations/jobs/conf_ils_low_run_simphy.sh

# 2. Verify simulations completed (1000 trees × 5 replicates × 21 networks)
python simulations/scripts/check_pipeline_status.py conf_ils_low_10M --step simphy

# 3. Submit method pipelines
sbatch simulations/jobs/submit_grampa_pipeline.sh conf_ils_low_10M
sbatch simulations/jobs/submit_polyphest_pipeline.sh conf_ils_low_10M

# 4. Verify results (see PIPELINE_STATUS_GUIDE.md for details)
python simulations/scripts/check_pipeline_status.py conf_ils_low_10M --step run
```

### Running Analyses on Real Data

```bash
# 1. Prepare gene trees (CRITICAL for GRAMPA)
python scripts/fix_substrings_grampa.py input_trees.tre output_trees.tre

# 2. Create taxa mapping if needed
python scripts/create_taxa_map.py -i input_trees.nex -o taxamap.txt

# 3. Submit analysis
sbatch jobs/grampa_general.sh
```

## Critical Notes

### GRAMPA
**MUST run `fix_substrings_grampa.py` first!** GRAMPA crashes if one tip label is a substring of another (e.g., "Qlo" and "Qlo03"). The script adds suffixes to make labels unique.

### Polyphest
- Creates output directories with percentile: `polyphest_p{PERCENTILE}/`
- **When validating, always specify the percentile you used:**
  ```bash
  python scripts/check_pipeline_status.py conf_ils_low_10M --method polyphest --step run --percentile 50
  ```

### MPSUGAR
- Requires NEXUS format: `mpsugar_trees.nex`
- Requires JSON taxon map: `taxon_map.json`

## Pipeline Validation

The pipeline checker validates every stage of your workflow. **Complete guide:** `simulations/PIPELINE_STATUS_GUIDE.md`

**Quick commands:**
```bash
# Check SimPhy simulations
python simulations/scripts/check_pipeline_status.py conf_ils_low_10M --step simphy

# Check inputs are ready before running methods
python simulations/scripts/check_pipeline_status.py conf_ils_low_10M --step prep

# Check outputs after methods complete
python simulations/scripts/check_pipeline_status.py conf_ils_low_10M --step run --verbose
```

## Environment Setup

- **Cluster path:** `/groups/itay_mayrose/tomulanovski/gene2net/`
- **Conda environments:** `gene2net` (main), `polyphest` (Polyphest-specific)
- **Job scheduler:** SLURM

When adapting scripts, update absolute paths in job scripts to match your directories.

## Working with Claude Code

### Git Commands
- **IMPORTANT:** Do NOT execute git commands automatically
- Suggest commands but let the user run them manually
- This allows the user to review changes before committing

### Development Approach
- Keep solutions simple and focused
- Only modify what's necessary for the task
- Avoid adding unnecessary features or abstractions
- Trust that the user knows their research workflow
