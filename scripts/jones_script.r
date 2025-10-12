#!/usr/bin/env Rscript
# General Jobes script to create the xml file for beast run
# Usage: Rscript "/groups/itay_mayrose/tomulanovski/gene2net/scripts/jones_script.r" path_to_directory
# takes all nex files in the directory and the taxa table and creates the xml
setwd("/groups/itay_mayrose/tomulanovski/gene2net/jones_examples/2013-05-15-manual (2)/")
library(ape)
source("AlloppDT_5beastxml_toplevel.r")
source("AlloppDT_6beastxml_bits.r")

#########################################################################
# Get directory containing .nex files from command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript script.R <nex_directory>")
}
nex_dir <- args[1]

# List all .nex files and remove extensions
nex_files <- list.files(path = nex_dir, pattern = "\\.nex$", full.names = FALSE)
if (length(nex_files) == 0) {
  stop(paste("No .nex files found in directory:", nex_dir))
}
alignmentnames <- sub("\\.nex$", "", nex_files)

#########################################################################
# directory for input and output
data.dpath <- nex_dir
fpath.taxatable <- "taxa_table.txt"

# BEAST information
# Extract parent folder name of nex_dir
parent_folder <- basename(dirname(nex_dir))
beastXMLfilename <- paste0(parent_folder, "_100mil.XML")

beast.chain.length <- "100000000"
beast.screen.logevery <- "10000"
beast.params.logevery <- "1000"
beast.gtrees.logevery <- "1000"
beast.multree.logevery <- "1000"
beast.dbugtune.logevery <- "10000"

# Print XML filename for confirmation
cat("BEAST XML filename will be:", beastXMLfilename, "\n")

# BEAST output file names
sampledgtrees.fnamebase <- "sampledgtrees"
sampledmultrees.fname <- "sampledmultrees.txt"
sampledparams.fname <- "sampledparams.txt"
DBUGTUNE.fname <- "DBUGTUNE.txt"

#########################################################################
beastAlloppDTxmlinfo <- list(
    data.dpath = data.dpath,
    fpath.taxatable = fpath.taxatable,
    alignmentnames = alignmentnames,
    beastXMLfilename = beastXMLfilename,
    beast.chain.length = beast.chain.length,
    beast.screen.logevery = beast.screen.logevery,
    beast.params.logevery = beast.params.logevery,
    beast.gtrees.logevery = beast.gtrees.logevery,
    beast.multree.logevery = beast.multree.logevery,
    beast.dbugtune.logevery = beast.dbugtune.logevery,
    sampledgtrees.fpathbase = paste(data.dpath, sampledgtrees.fnamebase, sep = ""),
    sampledmultrees.fpath = paste(data.dpath, sampledmultrees.fname, sep = ""),
    sampledparams.fpath = paste(data.dpath, sampledparams.fname, sep = ""), 
    DBUGTUNE.fpath = paste(data.dpath, DBUGTUNE.fname, sep = "")
)

make.beastxml.AlloppDT.forRealData(beastAlloppDTxmlinfo)
#########################################################################
