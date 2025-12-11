# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This repository (`gene2net`) is a research pipeline for phylogenetic network inference from gene trees. It contains scripts for preparing, running, and analyzing various phylogenetic network reconstruction methods (GRAMPA, Polyphest, MPAllopp, MPSUGAR, PADRE, Notung, ASTRAL) as well as simulation pipelines using SimPhy.

## Architecture

The repository is organized into three main directories:

### 1. `scripts/` - Data Processing and Format Conversion
Python scripts for data preparation and format conversion between different phylogenetic tools. Key categories:

- **Format conversion**: `grampa_to_polyphest_reformat.py`, `polyphest_to_grampa_reformat.py`, `fasta_to_nex.r`
- **Tree preprocessing**: `fix_substrings_grampa.py` (fixes tip labels where one is a substring of another, which crashes GRAMPA), `clean_tree_from_accessions.py`, `concatenate_trees.py`
- **Taxa mapping**: `create_taxa_map.py` (generates species-to-gene mappings for tools like MPAllopp), `create_taxa_table.py`
- **Copy number handling**: `add_copy_number.py`, `copies_smoothing.py` (for polyploidy analysis)
- **Network inference**: `MPSUGAR.py` (implements MP-SUGAR network inference using PhyNetPy)
- **Parameter inference**: `infer_gtr_parameters_distribution.py`, `infer_tree_heigh_distribution.py`, `infer_alignment_lengths.py`

### 2. `jobs/` - SLURM Job Scripts
Bash scripts for submitting analyses to a compute cluster (configured for SLURM). All scripts use conda environments and follow the pattern:
- Activate specific conda environment (e.g., `gene2net`, `polyphest`)
- Set input/output paths
- Run the tool with appropriate parameters
- Deactivate conda

Key job types:
- `grampa_general.sh`, `grampa_run.sh`, `grampa_iter.sh` - GRAMPA network inference
- `polyphest_general.sh`, `polyphest_run.sh` - Polyphest network inference
- `mpallop.sh`, `mpallop_general.sh` - MPAllopp polyploid network inference
- `mpsugar.sh`, `mpsugar_general.sh` - MPSUGAR network inference
- `padre_run.sh` - PADRE species tree inference
- `notung_for_dup.sh`, `notung_resolve.sh` - Notung reconciliation
- `Astral_4_run.sh`, `Astral_pro3_run.sh` - ASTRAL species tree inference

### 3. `simulations/` - Simulation Pipeline
Contains two subdirectories mirroring the main structure:

- **`simulations/scripts/`**: Simulation-specific analysis scripts
  - `compare_nets.py` - Comprehensive network comparison pipeline that handles both MUL-trees (standard Newick) and networks (extended Newick with #H markers)
  - `compare_reticulations.py` - Pairwise network comparison utilities
  - `add_nhx_tags.py`, `add_nhx_single.py` - Add NHX format tags for reconciliation
  - `aggregate_dup_loss.py`, `dup_loss_summary_to_csv.py` - Duplication/loss analysis
  - `extract_newick_from_nex.py`, `extract_sql.py` - Data extraction utilities

- **`simulations/jobs/`**: SLURM jobs for running simulations
  - `conf_ils_{high,med,low}_run_simphy.sh` - SimPhy simulations with different ILS levels
  - `calculate_dup_loss.sh`, `extract_dup_loss_reusable.sh` - Duplication/loss calculations
  - `calculate_rf_reusable.sh` - Robinson-Foulds distance calculations
  - `make_grampa_gene_trees.sh`, `make_polyphest_trees.sh`, `make_padre_trees.sh` - Generate input for different methods

## Common Development Workflows

### Running a phylogenetic analysis

1. **Prepare gene trees**: Ensure gene trees are in the correct format for your chosen method
   - For GRAMPA: Use `fix_substrings_grampa.py` to avoid substring label issues
   - For format conversion: Use `grampa_to_polyphest_reformat.py` or `polyphest_to_grampa_reformat.py`

2. **Create taxa mapping** (for methods requiring it like MPAllopp):
   ```bash
   python scripts/create_taxa_map.py -i input_trees.nex -o taxamap.txt
   ```

3. **Submit job**: Edit the appropriate job script in `jobs/` to set your input/output paths, then submit:
   ```bash
   sbatch jobs/grampa_general.sh
   ```

### Running simulations

1. Configure simulation parameters by editing the config script (e.g., `simulations/jobs/conf_ils_low_run_simphy.sh`)
2. Submit the simulation job which will call the reusable script
3. Use comparison scripts to analyze results:
   ```bash
   python simulations/scripts/compare_nets.py
   ```

### Network comparison pipeline

The `compare_nets.py` script provides a complete pipeline that:
- Automatically detects tree format (standard Newick vs extended Newick with #H markers)
- Handles both MUL-trees and networks
- Supports STRICT (exact) or POLYPHEST (fuzzy) matching for MUL-tree conversion
- Computes pairwise distances and generates comparison matrices
- Creates summary reports comparing methods to ground truth

## Tool-Specific Notes

### GRAMPA
- **Critical preprocessing**: Always use `fix_substrings_grampa.py` on gene trees before running GRAMPA
- Reason: GRAMPA crashes when building MUL-trees if one tip label is a substring of another (e.g., "Qlo" and "Qlo03")
- The script adds suffixes (default: "X") to problematic labels to make them unique

### Polyphest
- Requires consensus multiset file indicating gene copy numbers per species
- Uses percentile-based filtering (typical: 98th percentile)
- Can use near-isomorphic networks with threshold (typical: 0.2)

### MPSUGAR
- Uses the PhyNetPy library (`from PhyNetPy.MPSugar import MP_SUGAR`)
- Requires a taxon map dictionary mapping species to their gene copies
- Example in `scripts/MPSUGAR.py` shows format for Helianthus data

### MPAllopp
- Specifically designed for polyploid phylogenetic networks
- Requires taxa mapping in format: `<species1:gene1,gene2; species2:gene3,gene4>`
- Use `create_taxa_map.py` to generate this mapping automatically

## Conda Environments

The pipeline uses multiple conda environments:
- `gene2net` - Main environment for GRAMPA and general tools
- `polyphest` - Polyphest-specific environment
- Others as specified in individual job scripts

## File Formats

- **Gene trees**: Newick (.tre, .nwk) or NEXUS (.nex)
- **Networks**: Extended Newick with #H markers for reticulations
- **MUL-trees**: Standard Newick with duplicated tip labels
- **Taxa mappings**: Custom format for MPAllopp (generated by `create_taxa_map.py`)

## Path Conventions

The codebase uses absolute paths pointing to `/groups/itay_mayrose/tomulanovski/gene2net/`. When adapting scripts for new analyses:
1. Update paths in job scripts to point to your input/output directories
2. Ensure conda paths match your environment setup
3. Scripts in `scripts/` directory can be called with absolute paths or relative to the working directory

## Working with Claude Code

### Git Commands
- **IMPORTANT**: Do NOT execute git commands (add, commit, push, pull, etc.) automatically
- When git operations are needed, inform the user and let them execute the commands manually
- You can suggest the appropriate git commands, but the user will run them
