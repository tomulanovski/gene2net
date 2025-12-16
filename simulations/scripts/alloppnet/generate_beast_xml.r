#!/usr/bin/env Rscript
#
# generate_beast_xml.r - Generate BEAST XML for AlloppNET simulation analysis
#
# This script wraps the AlloppDT R scripts to generate BEAST XML files
# for AlloppNET phylogenetic network inference.
#
# Usage:
#   Rscript generate_beast_xml.r <input_dir> <output_dir>
#
# Arguments:
#   input_dir  - Directory containing NEXUS files and taxa_table.txt
#   output_dir - Directory where BEAST XML will be created
#
# Requirements:
#   - AlloppDT_5beastxml_toplevel.r (in same directory as this script)
#   - AlloppDT_6beastxml_bits.r (in same directory as this script)
#   - R package: ape
#
# Example:
#   Rscript generate_beast_xml.r \
#       /path/to/processed/config/alloppnet_input/replicate_1/ \
#       /path/to/results/config/alloppnet/replicate_1/

library(ape)

# Get script directory (where AlloppDT scripts should be located)
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
script_dir <- dirname(script_path)

# Source AlloppDT scripts from same directory
alloppdt_toplevel <- file.path(script_dir, "AlloppDT_5beastxml_toplevel.r")
alloppdt_bits <- file.path(script_dir, "AlloppDT_6beastxml_bits.r")

if (!file.exists(alloppdt_toplevel)) {
    stop("ERROR: AlloppDT_5beastxml_toplevel.r not found in ", script_dir, "\n",
         "Please copy AlloppDT scripts from cluster. See README.md for instructions.")
}

if (!file.exists(alloppdt_bits)) {
    stop("ERROR: AlloppDT_6beastxml_bits.r not found in ", script_dir, "\n",
         "Please copy AlloppDT scripts from cluster. See README.md for instructions.")
}

source(alloppdt_toplevel)
source(alloppdt_bits)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    cat("ERROR: Insufficient arguments\n\n")
    cat("Usage: Rscript generate_beast_xml.r <input_dir> <output_dir>\n\n")
    cat("Arguments:\n")
    cat("  input_dir  - Directory with NEXUS files and taxa_table.txt\n")
    cat("  output_dir - Directory where BEAST XML will be created\n\n")
    cat("Example:\n")
    cat("  Rscript generate_beast_xml.r \\\n")
    cat("      /path/to/processed/config/alloppnet_input/replicate_1/ \\\n")
    cat("      /path/to/results/config/alloppnet/replicate_1/\n")
    quit(status = 1)
}

input_dir <- args[1]
output_dir <- args[2]

cat("==============================================================================\n")
cat("AlloppNET BEAST XML Generation\n")
cat("==============================================================================\n")
cat("Input directory:  ", input_dir, "\n")
cat("Output directory: ", output_dir, "\n")
cat("==============================================================================\n\n")

# Verify input directory exists
if (!dir.exists(input_dir)) {
    stop("ERROR: Input directory not found: ", input_dir)
}

# Create output directory if needed
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory: ", output_dir, "\n")
}

# Find all .nex files in input directory
nex_files <- list.files(path = input_dir, pattern = "\\.nex$", full.names = FALSE)

if (length(nex_files) == 0) {
    stop("ERROR: No .nex files found in directory: ", input_dir)
}

cat("Found", length(nex_files), "NEXUS alignment files\n\n")

# Extract alignment names (remove .nex extension)
alignmentnames <- sub("\\.nex$", "", nex_files)

# BEAST configuration
# Note: fpath.taxatable should be JUST the filename, not full path
# AlloppDT will construct the full path using data.dpath + fpath.taxatable
fpath.taxatable <- "taxa_table.txt"

# Verify taxa table exists (construct full path for verification)
full_taxatable_path <- file.path(input_dir, fpath.taxatable)
if (!file.exists(full_taxatable_path)) {
    stop("ERROR: taxa_table.txt not found in input directory: ", input_dir)
}

cat("Taxa table: ", full_taxatable_path, "\n\n")

# BEAST XML output file
beastXMLfilename <- file.path(output_dir, "alloppnet.XML")

# BEAST runtime parameters
beast.chain.length <- "100000000"      # 100M iterations (~5 days)
beast.screen.logevery <- "10000"       # Screen log every 10K iterations
beast.params.logevery <- "1000"        # Parameter sampling every 1K
beast.gtrees.logevery <- "1000"        # Gene tree sampling every 1K
beast.multree.logevery <- "1000"       # Species tree sampling every 1K
beast.dbugtune.logevery <- "10000"     # Debug tuning every 10K

# BEAST output file names (all in output_dir)
sampledgtrees.fnamebase <- file.path(output_dir, "sampledgtrees")
sampledmultrees.fpath <- file.path(output_dir, "sampledmultrees.txt")
sampledparams.fpath <- file.path(output_dir, "sampledparams.txt")
DBUGTUNE.fpath <- file.path(output_dir, "DBUGTUNE.txt")

cat("BEAST configuration:\n")
cat("  Chain length:   ", beast.chain.length, " iterations\n")
cat("  Screen logging: every ", beast.screen.logevery, " iterations\n")
cat("  Sampling rate:  every ", beast.params.logevery, " iterations\n")
cat("\n")

# Create AlloppDT configuration list
beastAlloppDTxmlinfo <- list(
    data.dpath = input_dir,
    fpath.taxatable = fpath.taxatable,
    alignmentnames = alignmentnames,
    beastXMLfilename = beastXMLfilename,
    beast.chain.length = beast.chain.length,
    beast.screen.logevery = beast.screen.logevery,
    beast.params.logevery = beast.params.logevery,
    beast.gtrees.logevery = beast.gtrees.logevery,
    beast.multree.logevery = beast.multree.logevery,
    beast.dbugtune.logevery = beast.dbugtune.logevery,
    sampledgtrees.fpathbase = sampledgtrees.fnamebase,
    sampledmultrees.fpath = sampledmultrees.fpath,
    sampledparams.fpath = sampledparams.fpath,
    DBUGTUNE.fpath = DBUGTUNE.fpath
)

# Generate BEAST XML using AlloppDT
cat("Generating BEAST XML...\n")
tryCatch({
    make.beastxml.AlloppDT.forRealData(beastAlloppDTxmlinfo)
    cat("\n")
    cat("==============================================================================\n")
    cat("BEAST XML generated successfully!\n")
    cat("==============================================================================\n")
    cat("Output file: ", beastXMLfilename, "\n")
    cat("\n")
    cat("Expected BEAST output files:\n")
    cat("  - sampledmultrees.txt  (species networks)\n")
    cat("  - sampledgtrees*.txt   (gene trees)\n")
    cat("  - sampledparams.txt    (parameters)\n")
    cat("  - DBUGTUNE.txt         (debug info)\n")
    cat("==============================================================================\n")
}, error = function(e) {
    cat("\nERROR: Failed to generate BEAST XML\n")
    cat("Error message:", e$message, "\n")
    quit(status = 1)
})
