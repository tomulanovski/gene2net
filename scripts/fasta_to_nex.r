#!/usr/bin/env Rscript

# General FASTA to NEXUS converter
# Usage: Rscript "/groups/itay_mayrose/tomulanovski/gene2net/scripts/fasta_to_nex.r" input.fasta output.nex

library(ape)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments provided
if(length(args) != 2) {
  cat("\nUsage: Rscript fasta_to_nexus.r <input.fasta> <output.nex>\n")
  cat("Example: Rscript fasta_to_nexus.r sequences.fasta sequences.nex\n\n")
  quit(status = 1)
}

input_file <- args[1]
output_file <- args[2]

# Check if input file exists
if(!file.exists(input_file)) {
  cat("Error: Input file '", input_file, "' not found\n", sep="")
  quit(status = 1)
}

# Check if ape package is available
if(!require(ape, quietly = TRUE)) {
  cat("Error: R package 'ape' is required but not installed\n")
  cat("Install it with: install.packages('ape')\n")
  quit(status = 1)
}

cat("Converting FASTA to NEXUS...\n")
cat("Input file: ", input_file, "\n")
cat("Output file:", output_file, "\n")

# Read FASTA file
tryCatch({
  sequences <- read.FASTA(input_file)
}, error = function(e) {
  cat("Error reading FASTA file:", e$message, "\n")
  quit(status = 1)
})

# Write NEXUS file
tryCatch({
  write.nexus.data(sequences, output_file, format="DNA")
}, error = function(e) {
  cat("Error writing NEXUS file:", e$message, "\n")
  quit(status = 1)
})

# Print summary
cat("\nConversion successful!\n")
cat("Number of sequences:", length(sequences), "\n")
cat("Sequence names (first 5):\n")
print(head(names(sequences), 5))
cat("Output saved to:", output_file, "\n")