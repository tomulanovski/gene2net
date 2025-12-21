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
- **`submit_all_methods.sh`** - Master orchestration script for running methods across multiple configs
- **`check_pipeline_status.py`** - Validates SimPhy simulations and all method inputs/outputs
- **`compare_nets.py`** - Network comparison pipeline (handles MUL-trees and networks)
- **`run_reticulation_stats.py`** - Generates network statistics including WGD event detection
- **`create_analysis_figures.py`** - Publication-quality figures and tables for method evaluation
- **`test_autopolyploid_detection.py`** - Test script for verifying autopolyploidization event counting
- **Full documentation:**
  - `simulations/md_files/COMPLETE_PIPELINE_GUIDE.md` - End-to-end workflow from SimPhy to summary
  - `simulations/md_files/METHODS_GUIDE.md` - Running phylogenetic network inference methods
  - `simulations/md_files/PIPELINE_STATUS_GUIDE.md` - Validation and status checking
  - `simulations/md_files/VISUALIZATION_GUIDE.md` - Publication figures and analysis (UPDATED 2025-12-20)

## Common Workflows

### Running Simulations

```bash
# 1. Configure and submit SimPhy
cd simulations/jobs
./submit_simphy.sh conf_ils_low_10M

# 2. Verify simulations completed (1000 trees × 5 replicates × 21 networks)
cd ../scripts
python check_pipeline_status.py conf_ils_low_10M --step simphy

# 3. Submit ALL method pipelines with master script (RECOMMENDED)
cd ../jobs
./submit_all_methods.sh conf_ils_low_10M

# OR submit methods individually
./submit_grampa_pipeline.sh conf_ils_low_10M
./submit_polyphest_pipeline.sh conf_ils_low_10M
./submit_padre_pipeline.sh conf_ils_low_10M
./submit_mpsugar_pipeline.sh conf_ils_low_10M

# 4. Verify results
cd ../scripts
python check_pipeline_status.py conf_ils_low_10M --step run --verbose

# 5. Post-process and analyze
python postprocess_results.py conf_ils_low_10M
python run_full_summary.py conf_ils_low_10M
```

**For complete workflow:** See `simulations/COMPLETE_PIPELINE_GUIDE.md`

### Publication Analysis and Figures

After completing simulations, generate comprehensive publication-quality figures and tables.

**Complete guide:** See `simulations/md_files/VISUALIZATION_GUIDE.md`

**Quick start:**
```bash
cd simulations/scripts

# Generate all analysis figures for one or more configurations
python create_analysis_figures.py --config conf_ils_low_10M
python create_analysis_figures.py --config conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M
```

**Output structure** (per configuration):
```
simulations/analysis/summary/{config}/
├── plots/
│   ├── 01_success_vs_reticulations.pdf/png  # Success rate vs H_Strict
│   ├── 02_success_vs_polyploids.pdf/png     # Success rate vs Num_Polyploids
│   ├── 03_success_vs_wgd.pdf/png            # Success rate vs Total_WGD
│   ├── 04_method_performance_overview.pdf/png
│   ├── 05_method_network_heatmap.pdf/png
│   ├── 06_reticulation_accuracy.pdf/png
│   └── 07_difficulty_correlations.pdf/png
└── tables/
    ├── summary_table.csv
    └── summary_simple.csv
```

**Key features:**
- Non-interactive matplotlib backend (no X11 required)
- Per-configuration organization for easy comparison across ILS levels
- Comprehensive evaluation: success rates, reticulation accuracy, polyploid identification
- Network characteristics from `simulations/networks/mul_tree_final_stats.csv`

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
- **Must run with each percentile separately** to compare p50, p70, p90:
  ```bash
  ./submit_all_methods.sh conf_ils_low_10M --methods polyphest --polyphest-percentile 50
  ./submit_all_methods.sh conf_ils_low_10M --methods polyphest --polyphest-percentile 70
  ./submit_all_methods.sh conf_ils_low_10M --methods polyphest --polyphest-percentile 90
  ```

### MPSUGAR
- Requires NEXUS format: `mpsugar_trees.nex`
- Requires JSON taxon map: `taxon_map.json`

### Network Statistics and WGD Events

**Ground truth statistics:** `simulations/networks/mul_tree_final_stats.csv`

**Important distinctions:**
- **Autopolyploidization EVENT**: A single WGD that duplicates an entire clade (creates identical sibling subtrees)
  - Example: Ding_2023 has ONE event (Rch clade duplication) → creates 7 autopolyploid species
  - Detected by finding identical siblings in MUL-tree before network folding
- **Allopolyploidization/Hybridization**: Reticulation events where different lineages merge (H_Strict count)
  - Detected as reticulations with 2+ parents after network folding
- **Total_WGD**: Sum of autopolyploidization events + allopolyploidization events (H_Strict)
- **Num_Polyploids**: Count of species appearing >1 time in the MUL-tree (NOT the number of events)

**Key CSV columns:**
- `H_Strict`: Number of reticulations (allopolyploidization/hybridization events)
- `Num_Polyploids`: Number of polyploid species (species with >1 copy)
- `Num_Autopolyploidization_Events`: Number of autopolyploidization events
- `Total_WGD`: Total whole genome duplication events

**Regenerate statistics:**
```bash
cd simulations/scripts
python run_reticulation_stats.py ../networks/ ../networks/mul_tree_final_stats.csv
```

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

### Debugging Failed Runs
- **`investigate_network.sh`** - Automated script to debug why specific networks failed
- **`LOG_SEARCH_GUIDE.md`** - Complete guide for searching and navigating log files
```bash
# Investigate why a network failed
simulations/scripts/investigate_network.sh Diaz-Perez_2018 conf_ils_low_10M
```

## Environment Setup

- **Cluster path:** `/groups/itay_mayrose/tomulanovski/gene2net/`
- **Conda environments:** `gene2net` (main), `polyphest` (Polyphest-specific)
- **Job scheduler:** SLURM

When adapting scripts, update absolute paths in job scripts to match your directories.

### Git Configuration
Analysis outputs should be excluded from version control:
```bash
# Add to .gitignore if not present
simulations/analysis/
**/__pycache__/
**/*.pyc
```

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

## Recent Updates (2025-12-20)

### Autopolyploidization Event Detection
**Problem:** Previous implementation counted each duplicated species as a separate autopolyploidization event, rather than counting the actual WGD events.

**Solution:** Modified `run_reticulation_stats.py` with `count_autopolyploid_events()` function that:
- Uses top-down (preorder) traversal to identify top-level duplicated clades
- Detects identical sibling subtrees (same parent) as ONE event
- Tracks descendants to avoid double-counting nested duplications
- Example: Ding_2023 has ONE autopolyploidization event → creates 7 autopolyploid species

**Key distinction:**
- **Autopolyploidization** (identical siblings) → simplified away in network folding → H_Strict = 0
- **Allopolyploidization** (identical non-siblings) → persist as reticulations → H_Strict > 0

### Publication Analysis Pipeline
Created `create_analysis_figures.py` for comprehensive method evaluation:

**Features:**
- Per-configuration output organization (`simulations/analysis/summary/{config}/plots/` and `tables/`)
- Non-interactive matplotlib backend (no X11 required for cluster use)
- Automatic path resolution relative to script location
- Uses `inferred_exists` boolean column from inventory CSV

**Core plots:**
1. Success rate vs. number of reticulations (H_Strict)
2. Success rate vs. number of polyploids
3. Success rate vs. total WGD events
4. Method performance overview
5. Method-network heatmap
6. Reticulation accuracy scatter plot
7. Difficulty correlation matrix

**Usage:**
```bash
cd simulations/scripts
python create_analysis_figures.py --config conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M
```

### Files Modified/Created
- `simulations/scripts/run_reticulation_stats.py` - Added autopolyploidization event detection
- `simulations/scripts/create_analysis_figures.py` - New comprehensive analysis script
- `simulations/scripts/test_autopolyploid_detection.py` - New test script for verification
- `simulations/networks/mul_tree_final_stats.csv` - Updated with new columns
- `simulations/md_files/VISUALIZATION_GUIDE.md` - Complete documentation for analysis workflow
