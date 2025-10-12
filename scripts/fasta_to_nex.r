#!/usr/bin/env Rscript

# General FASTA to NEXUS converter
# Usage: Rscript "/groups/itay_mayrose/tomulanovski/gene2net/scripts/fasta_to_nex.r" path_to_input_directory path_to_output_directory
# Converts all fasta files in the directory to NEXUS files


library(ape)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 2) {
  cat("\nUsage: Rscript fasta_to_nex.r <input_directory> <output_directory>\n")
  cat("Example: Rscript fasta_to_nex.r /path/to/fasta_files /path/to/nexus_output\n\n")
  quit(status = 1)
}

input_dir <- args[1]
output_dir <- args[2]

# Check if input directory exists
if(!dir.exists(input_dir)) {
  cat("Error: Input directory '", input_dir, "' not found\n", sep="")
  quit(status = 1)
}

# Create output directory if it doesn't exist
if(!dir.exists(output_dir)) {
  cat("Creating output directory:", output_dir, "\n")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if(!dir.exists(output_dir)) {
    cat("Error: Could not create output directory '", output_dir, "'\n", sep="")
    quit(status = 1)
  }
} else {
  cat("Using existing output directory:", output_dir, "\n")
}

# Get all .fasta files
fasta_files <- list.files(path = input_dir, pattern = "\\.fasta$", full.names = TRUE)

if(length(fasta_files) == 0) {
  cat("No .fasta files found in the directory.\n")
  quit(status = 1)
}

cat("Found", length(fasta_files), "FASTA files to convert.\n\n")

# Loop through each FASTA file
successful_conversions <- 0
failed_conversions <- 0

for(fasta_file in fasta_files) {
  # Get the base filename without path and extension
  base_filename <- tools::file_path_sans_ext(basename(fasta_file))
  
  # Create output file path in the output directory
  output_file <- file.path(output_dir, paste0(base_filename, ".nex"))
  
  cat("Converting:", basename(fasta_file), "->", basename(output_file), "\n")
  
  tryCatch({
    # Read FASTA file
    sequences <- read.FASTA(fasta_file)
    
    # Write NEXUS file
    write.nexus.data(sequences, output_file, format="DNA")
    
    cat("  Success! Number of sequences:", length(sequences), "\n\n")
    successful_conversions <- successful_conversions + 1
  }, error = function(e) {
    cat("  Error processing file:", e$message, "\n\n")
    failed_conversions <- failed_conversions + 1
  })
}

cat("Conversion completed!\n")
cat("Successful conversions:", successful_conversions, "\n")
cat("Failed conversions:", failed_conversions, "\n")
cat("Output directory:", output_dir, "\n")