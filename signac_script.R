#!/usr/bin/env Rscript

# Run Signac's variant identification and export only the tables required by
# the Python-vs-Signac comparison. Usage is documented in README.md.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: signac_script.R <input_dir> <low_coverage_threshold> <output_dir>")
}

input_dir <- args[[1]]
low_coverage_threshold <- as.integer(args[[2]])
output_dir <- args[[3]]
if (is.na(low_coverage_threshold)) {
  stop("low_coverage_threshold must be an integer")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({
  library(Matrix)
  library(Signac)
})

# Load each allele's forward/reverse counts and the mitochondrial reference.
alleles <- c("A", "C", "T", "G")
counts_by_allele <- setNames(lapply(alleles, function(base) {
  path <- file.path(input_dir, paste0("mgatk.", base, ".txt"))
  read.csv(
    path,
    header = FALSE,
    col.names = c("pos", "cell", "forw", "rev")
  )
}), alleles)

reference <- read.csv(
  file.path(input_dir, "chrM_refAllele.txt"),
  header = FALSE,
  col.names = c("pos", "ref")
)
reference$pos <- as.integer(reference$pos)
reference$ref <- toupper(reference$ref)

cells <- unique(unlist(
  lapply(counts_by_allele, function(data) data$cell),
  use.names = FALSE
))
max_position <- max(reference$pos)

# Build one position-by-cell sparse matrix for a base and strand.
make_strand_matrix <- function(base, strand) {
  data <- counts_by_allele[[base]]
  count_column <- if (strand == "fwd") "forw" else "rev"
  cell_index <- match(data$cell, cells)
  if (anyNA(cell_index)) {
    stop("A count-file barcode is missing from the cell list")
  }

  Matrix::sparseMatrix(
    i = as.integer(data$pos),
    j = cell_index,
    x = as.numeric(data[[count_column]]),
    dims = c(max_position, length(cells)),
    dimnames = list(
      paste(base, seq_len(max_position), strand, sep = "-"),
      cells
    )
  )
}

# Keep the row ordering expected by Signac::ReadMGATK and AlleleFreq.
counts <- do.call(rbind, list(
  make_strand_matrix("A", "fwd"),
  make_strand_matrix("C", "fwd"),
  make_strand_matrix("T", "fwd"),
  make_strand_matrix("G", "fwd"),
  make_strand_matrix("A", "rev"),
  make_strand_matrix("C", "rev"),
  make_strand_matrix("T", "rev"),
  make_strand_matrix("G", "rev")
))

# IdentifyVariants evaluates substitutions and provides Signac's filter metrics.
candidate_stats <- Signac::IdentifyVariants(
  object = counts,
  refallele = reference,
  stabilize_variance = TRUE,
  low_coverage_threshold = low_coverage_threshold,
  verbose = FALSE
)

# Apply the Signac vignette's high-confidence thresholds.
high_confidence <- subset(
  candidate_stats,
  n_cells_conf_detected >= 5 &
    strand_correlation >= 0.65 &
    vmr > 0.01
)

# AlleleFreq requires at least two variants in this matrix method. Temporarily
# include an extra identified variant for a one-call result, then discard it.
frequency_variants <- unique(as.character(high_confidence$variant))
frequency_path <- file.path(output_dir, "signac_allele_freq.csv")
if (length(frequency_variants) == 1) {
  padding_variant <- setdiff(as.character(candidate_stats$variant), frequency_variants)[1]
  if (!is.na(padding_variant)) {
    frequency_matrix <- Signac::AlleleFreq(
      object = counts,
      variants = c(frequency_variants, padding_variant)
    )
    frequency_matrix <- frequency_matrix[frequency_variants, , drop = FALSE]
  } else {
    frequency_matrix <- Signac::AlleleFreq(object = counts, variants = frequency_variants)
  }
} else if (length(frequency_variants) > 1) {
  frequency_matrix <- Signac::AlleleFreq(object = counts, variants = frequency_variants)
  frequency_matrix <- frequency_matrix[frequency_variants, , drop = FALSE]
} else {
  frequency_matrix <- Matrix::Matrix(
    0,
    nrow = 0,
    ncol = length(cells),
    sparse = TRUE,
    dimnames = list(character(), cells)
  )
}

# Save cells x variants to align with Python's AnnData VAF matrix.
frequency_frame <- data.frame(
  cell = cells,
  as.matrix(t(frequency_matrix)),
  check.names = FALSE
)
write.csv(frequency_frame, frequency_path, row.names = FALSE)
write.csv(
  candidate_stats,
  file.path(output_dir, "signac_candidate_stats.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(variant = as.character(high_confidence$variant)),
  file.path(output_dir, "signac_high_confidence.csv"),
  row.names = FALSE
)

cat("Signac candidates evaluated:", nrow(candidate_stats), "\n")
cat("Signac high-confidence variants:", nrow(high_confidence), "\n")
cat("Cells in Signac allele-frequency matrix:", length(cells), "\n")
