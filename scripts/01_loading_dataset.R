# ============================================================
# HU LAB REPRODUCTION PROJECT
# load_dataset.R — LOAD, VALIDATE, AND EXTRACT DATASETS
#
# Integrates:
#   01_load_dataset.R
#   02_extract_HER2_samples.R
#
# Normal reference:
#   GSE113197 — Ind5, Ind6, Ind7
#
# HER2+/ER+ samples:
#   CID3586 — Naive
#   CID4066 — Treated
#
# CID3586 and CID4066 are extracted directly from the full
# GSE176078 sparse expression matrix.
#
# No dense conversion of the full HER2 matrix is performed.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
})

# ------------------------------------------------------------
# 0. PROJECT PATHS
# ------------------------------------------------------------

project_dir <- getwd()

raw_dir <- file.path(project_dir, "data", "raw")
gse176078_dir <- file.path(raw_dir, "GSE176078")

processed_her2_dir <- file.path(
  project_dir,
  "data",
  "processed",
  "HER2"
)

validation_dir <- file.path(
  project_dir,
  "results",
  "01_dataset_validation"
)

extraction_dir <- file.path(
  project_dir,
  "results",
  "02_HER2_extraction"
)

dir.create(
  processed_her2_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  validation_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  extraction_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

message("==============================================")
message("load_dataset.R")
message("LOAD, VALIDATE, AND EXTRACT DATASETS")
message("==============================================")
message("Project directory: ", project_dir)
message("")

# ------------------------------------------------------------
# 1. HELPER FUNCTION — CHECK FILE
# ------------------------------------------------------------

check_file <- function(path) {

  if (!file.exists(path)) {
    stop(
      "Required file not found:\n",
      path
    )
  }

  message("Found: ", path)
}

# ============================================================
# PART A — NORMAL BREAST REFERENCE DATA
# ============================================================

# ------------------------------------------------------------
# 2. REQUIRED NORMAL FILES
# ------------------------------------------------------------

normal_files <- c(
  file.path(
    raw_dir,
    "GSM3099847_Ind5_Expression_Matrix.txt"
  ),
  file.path(
    raw_dir,
    "GSM3099848_Ind6_Expression_Matrix.txt"
  ),
  file.path(
    raw_dir,
    "GSM3099849_Ind7_Expression_Matrix.txt"
  ),
  file.path(
    raw_dir,
    "Normal_celltype.txt"
  )
)

message("Checking normal breast files...")

for (f in normal_files) {
  check_file(f)
}

# ------------------------------------------------------------
# 3. LOAD NORMAL BREAST METADATA
# ------------------------------------------------------------

message("")
message("==============================================")
message("LOADING NORMAL BREAST METADATA")
message("==============================================")

normal_annotation <- read.table(
  file.path(
    raw_dir,
    "Normal_celltype.txt"
  ),
  header = TRUE,
  stringsAsFactors = FALSE
)

message(
  "Normal annotation rows: ",
  nrow(normal_annotation)
)

message(
  "Normal annotation columns: ",
  paste(
    colnames(normal_annotation),
    collapse = ", "
  )
)

print(head(normal_annotation))

# ------------------------------------------------------------
# 4. LOAD NORMAL BREAST EXPRESSION MATRICES
# ------------------------------------------------------------

load_normal_matrix <- function(
  file_path,
  sample_name
) {

  message("")
  message("----------------------------------------------")
  message("Loading ", sample_name)
  message("----------------------------------------------")

  message(
    "File: ",
    basename(file_path)
  )

  # Expression matrices are tab-delimited text
  # with genes in rows and cells in columns.

  mat <- fread(
    file_path,
    sep = "\t",
    header = TRUE,
    data.table = FALSE,
    check.names = FALSE
  )

  message(
    "Raw dimensions: ",
    nrow(mat),
    " genes × ",
    ncol(mat) - 1,
    " cells"
  )

  gene_names <- mat[[1]]

  expression_data <- as.matrix(
    mat[, -1, drop = FALSE]
  )

  rownames(expression_data) <- gene_names

  # Convert to numeric matrix
  expression_data <- matrix(
    as.numeric(expression_data),
    nrow = nrow(expression_data),
    ncol = ncol(expression_data),
    dimnames = dimnames(expression_data)
  )

  # Check duplicate genes
  duplicate_genes <- sum(
    duplicated(rownames(expression_data))
  )

  message(
    "Duplicated gene names: ",
    duplicate_genes
  )

  if (duplicate_genes > 0) {

    rownames(expression_data) <- make.unique(
      rownames(expression_data)
    )

    message(
      "Duplicate gene names were made unique."
    )
  }

  # Check duplicate cell barcodes
  duplicate_cells <- sum(
    duplicated(colnames(expression_data))
  )

  message(
    "Duplicated cell barcodes: ",
    duplicate_cells
  )

  if (duplicate_cells > 0) {
    stop(
      "Duplicate cell barcodes detected in ",
      sample_name
    )
  }

  message(
    "Final dimensions: ",
    nrow(expression_data),
    " genes × ",
    ncol(expression_data),
    " cells"
  )

  return(expression_data)
}

normal_matrix_files <- c(
  Ind5 = file.path(
    raw_dir,
    "GSM3099847_Ind5_Expression_Matrix.txt"
  ),
  Ind6 = file.path(
    raw_dir,
    "GSM3099848_Ind6_Expression_Matrix.txt"
  ),
  Ind7 = file.path(
    raw_dir,
    "GSM3099849_Ind7_Expression_Matrix.txt"
  )
)

normal_matrices <- list()

for (sample_name in names(normal_matrix_files)) {

  normal_matrices[[sample_name]] <-
    load_normal_matrix(
      normal_matrix_files[[sample_name]],
      sample_name
    )
}

# ------------------------------------------------------------
# 5. VALIDATE NORMAL CELL BARCODES
# ------------------------------------------------------------

message("")
message("==============================================")
message("NORMAL BREAST BARCODE VALIDATION")
message("==============================================")

normal_barcode_summary <- data.frame()

for (sample_name in names(normal_matrices)) {

  expression_barcodes <-
    colnames(normal_matrices[[sample_name]])

  annotation_barcodes <- rownames(
    normal_annotation[
      normal_annotation$sample == sample_name,
      ,
      drop = FALSE
    ]
  )

  matched <-
    sum(
      expression_barcodes %in% annotation_barcodes
    )

  message(
    sample_name,
    ": ",
    length(expression_barcodes),
    " expression cells"
  )

  message(
    sample_name,
    ": ",
    length(annotation_barcodes),
    " annotation cells"
  )

  message(
    sample_name,
    ": ",
    matched,
    " matched cells"
  )

  normal_barcode_summary <- rbind(
    normal_barcode_summary,
    data.frame(
      sample = sample_name,
      expression_cells = length(
        expression_barcodes
      ),
      annotation_cells = length(
        annotation_barcodes
      ),
      matched_cells = matched
    )
  )
}

print(normal_barcode_summary)

write.csv(
  normal_barcode_summary,
  file.path(
    validation_dir,
    "normal_barcode_validation.csv"
  ),
  row.names = FALSE
)

# ============================================================
# PART B — HER2+/ER+ CID3586 AND CID4066
# ============================================================

# ------------------------------------------------------------
# 6. REQUIRED GSE176078 FILES
# ------------------------------------------------------------

message("")
message("==============================================")
message("CHECKING HER2+ SOURCE FILES")
message("==============================================")

matrix_file <- file.path(
  gse176078_dir,
  "count_matrix_sparse.mtx"
)

gene_file <- file.path(
  gse176078_dir,
  "count_matrix_genes.tsv"
)

barcode_file <- file.path(
  gse176078_dir,
  "count_matrix_barcodes.tsv"
)

metadata_file <- file.path(
  gse176078_dir,
  "metadata.csv"
)

her2_source_files <- c(
  matrix_file,
  gene_file,
  barcode_file,
  metadata_file
)

for (f in her2_source_files) {
  check_file(f)
}

# ------------------------------------------------------------
# 7. LOAD HER2+ GENES AND BARCODES
# ------------------------------------------------------------

message("")
message("==============================================")
message("LOADING HER2+ MATRIX DIMENSIONS")
message("==============================================")

genes <- readLines(
  gene_file,
  warn = FALSE
)

barcodes <- readLines(
  barcode_file,
  warn = FALSE
)

message(
  "Genes: ",
  length(genes)
)

message(
  "Barcodes: ",
  length(barcodes)
)

# ------------------------------------------------------------
# 8. BASIC FEATURE/BARCODE QC
# ------------------------------------------------------------

message("")
message("==============================================")
message("FEATURE / BARCODE QC")
message("==============================================")

if (anyDuplicated(genes) > 0) {
  stop(
    "Duplicated gene identifiers detected."
  )
}

if (anyDuplicated(barcodes) > 0) {
  stop(
    "Duplicated cell barcodes detected."
  )
}

message("Duplicated genes: 0")
message("Duplicated barcodes: 0")

# ------------------------------------------------------------
# 9. LOAD HER2+ METADATA
# ------------------------------------------------------------

message("")
message("==============================================")
message("LOADING HER2+ METADATA")
message("==============================================")

metadata <- read.csv(
  metadata_file,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

message(
  "Metadata rows: ",
  nrow(metadata)
)

message(
  "Metadata columns: ",
  paste(
    colnames(metadata),
    collapse = ", "
  )
)

# ------------------------------------------------------------
# 10. VERIFY MATRIX / BARCODE / METADATA DIMENSIONS
# ------------------------------------------------------------

message("")
message("==============================================")
message("STRUCTURAL VALIDATION")
message("==============================================")

# Read MatrixMarket header without loading the matrix.

con <- file(
  matrix_file,
  open = "r"
)

header_line <- NULL

repeat {

  line <- readLines(
    con,
    n = 1
  )

  if (length(line) == 0) {
    break
  }

  if (!startsWith(line, "%")) {
    header_line <- line
    break
  }
}

close(con)

if (is.null(header_line)) {
  stop(
    "Could not find MatrixMarket dimension line."
  )
}

matrix_dimensions <- as.integer(
  strsplit(
    trimws(header_line),
    "\\s+"
  )[[1]]
)

if (length(matrix_dimensions) != 3) {
  stop(
    "Could not parse MatrixMarket dimensions."
  )
}

n_genes_matrix <- matrix_dimensions[1]
n_cells_matrix <- matrix_dimensions[2]
n_nonzero <- matrix_dimensions[3]

message(
  "Matrix genes: ",
  n_genes_matrix
)

message(
  "Matrix cells: ",
  n_cells_matrix
)

message(
  "Matrix non-zero entries: ",
  n_nonzero
)

if (
  n_genes_matrix != length(genes)
) {

  stop(
    "Gene count mismatch: matrix has ",
    n_genes_matrix,
    " rows but genes.tsv has ",
    length(genes),
    " genes."
  )
}

if (
  n_cells_matrix != length(barcodes)
) {

  stop(
    "Cell count mismatch: matrix has ",
    n_cells_matrix,
    " columns but barcode file has ",
    length(barcodes),
    " barcodes."
  )
}

if (
  nrow(metadata) != length(barcodes)
) {

  stop(
    "Metadata/barcode count mismatch: metadata has ",
    nrow(metadata),
    " rows but barcode file has ",
    length(barcodes),
    " cells."
  )
}

message("Matrix dimensions: PASSED")
message("Gene mapping: PASSED")
message("Barcode mapping: PASSED")
message("Metadata dimensions: PASSED")

# ------------------------------------------------------------
# 11. VERIFY BARCODE ORDER AGAINST METADATA
# ------------------------------------------------------------

message("")
message("==============================================")
message("BARCODE / METADATA ORDER VALIDATION")
message("==============================================")

metadata_barcodes <- rownames(metadata)

same_order <- identical(
  barcodes,
  metadata_barcodes
)

if (!same_order) {

  message(
    "Barcode order is not identical to metadata row order."
  )

  matched_positions <- match(
    barcodes,
    metadata_barcodes
  )

  if (anyNA(matched_positions)) {
    stop(
      "Some matrix barcodes are missing from metadata."
    )
  }

  message(
    "All barcodes are present in metadata."
  )

  message(
    "Metadata will be reordered to matrix barcode order."
  )

  metadata <- metadata[
    barcodes,
    ,
    drop = FALSE
  ]

} else {

  message(
    "Barcode order matches metadata order: PASSED"
  )
}

# ------------------------------------------------------------
# 12. IDENTIFY HER2+ CELLS
# ------------------------------------------------------------

message("")
message("==============================================")
message("IDENTIFYING HER2+ CELLS")
message("==============================================")

target_samples <- c(
  "CID3586",
  "CID4066"
)

if (
  !"orig.ident" %in%
  colnames(metadata)
) {

  stop(
    "Metadata does not contain 'orig.ident'."
  )
}

selected_cells <- barcodes[
  metadata$orig.ident %in%
    target_samples
]

message(
  "Target samples: ",
  paste(
    target_samples,
    collapse = ", "
  )
)

message(
  "Total HER2+ cells selected: ",
  length(selected_cells)
)

# ------------------------------------------------------------
# 13. SAMPLE-SPECIFIC CELL COUNTS
# ------------------------------------------------------------

message("")
message("==============================================")
message("HER2+ SAMPLE CELL COUNTS")
message("==============================================")

sample_counts <- as.data.frame(
  table(
    metadata[
      selected_cells,
      "orig.ident"
    ]
  )
)

colnames(sample_counts) <- c(
  "sample",
  "cells"
)

print(sample_counts)

expected_counts <- data.frame(
  sample = c(
    "CID3586",
    "CID4066"
  ),
  expected_cells = c(
    6178,
    5309
  )
)

sample_counts <- merge(
  expected_counts,
  sample_counts,
  by = "sample",
  all = TRUE
)

if (
  any(
    sample_counts$expected_cells !=
      sample_counts$cells
  )
) {

  stop(
    "HER2 sample cell counts do not match expected values."
  )
}

message("")
message(
  "CID3586: 6,178 cells — PASSED"
)

message(
  "CID4066: 5,309 cells — PASSED"
)

if (
  length(selected_cells) != 11487
) {

  stop(
    "Expected 11,487 HER2+ cells but found ",
    length(selected_cells)
  )
}

message(
  "Total HER2+ cells: 11,487 — PASSED"
)

# ------------------------------------------------------------
# 14. HER2+ CELLTYPE VALIDATION
# ------------------------------------------------------------

message("")
message("==============================================")
message("HER2+ CELLTYPE DISTRIBUTION")
message("==============================================")

if (
  "celltype_major" %in%
  colnames(metadata)
) {

  message("CID3586:")

  print(
    sort(
      table(
        metadata[
          selected_cells,
          "celltype_major"
        ][
          metadata[
            selected_cells,
            "orig.ident"
          ] == "CID3586"
        ]
      ),
      decreasing = TRUE
    )
  )

  message("")
  message("CID4066:")

  print(
    sort(
      table(
        metadata[
          selected_cells,
          "celltype_major"
        ][
          metadata[
            selected_cells,
            "orig.ident"
          ] == "CID4066"
        ]
      ),
      decreasing = TRUE
    )
  )
}

# ------------------------------------------------------------
# 15. LOAD FULL MATRIX AS SPARSE MATRIX
# ------------------------------------------------------------

message("")
message("==============================================")
message("LOADING SPARSE HER2+ SOURCE MATRIX")
message("==============================================")

message(
  "Reading MatrixMarket file..."
)

message(
  "This may take some time because the source matrix",
  " contains 177,994,136 non-zero entries."
)

full_matrix <- readMM(
  matrix_file
)

message(
  "Sparse matrix loaded."
)

message(
  "Dimensions: ",
  nrow(full_matrix),
  " genes × ",
  ncol(full_matrix),
  " cells"
)

message(
  "Class: ",
  class(full_matrix)[1]
)

# ------------------------------------------------------------
# 16. FINAL MATRIX STRUCTURE VALIDATION
# ------------------------------------------------------------

if (
  nrow(full_matrix) !=
    length(genes)
) {

  stop(
    "Loaded matrix row count does not match genes."
  )
}

if (
  ncol(full_matrix) !=
    length(barcodes)
) {

  stop(
    "Loaded matrix column count does not match barcodes."
  )
}

rownames(full_matrix) <- make.unique(
  genes
)

colnames(full_matrix) <- barcodes

message(
  "Matrix dimensions and identifiers: PASSED"
)

# ------------------------------------------------------------
# 17. EXTRACT HER2+ CELLS
# ------------------------------------------------------------

message("")
message("==============================================")
message("EXTRACTING HER2+ EXPRESSION")
message("==============================================")

her2_matrix <- full_matrix[
  ,
  selected_cells,
  drop = FALSE
]

message(
  "Extracted matrix: ",
  nrow(her2_matrix),
  " genes × ",
  ncol(her2_matrix),
  " cells"
)

# ------------------------------------------------------------
# 18. SPLIT CID3586 / CID4066
# ------------------------------------------------------------

message("")
message("==============================================")
message("SPLITTING HER2+ SAMPLES")
message("==============================================")

cid3586_cells <- selected_cells[
  metadata[
    selected_cells,
    "orig.ident"
  ] == "CID3586"
]

cid4066_cells <- selected_cells[
  metadata[
    selected_cells,
    "orig.ident"
  ] == "CID4066"
]

cid3586_matrix <- her2_matrix[
  ,
  cid3586_cells,
  drop = FALSE
]

cid4066_matrix <- her2_matrix[
  ,
  cid4066_cells,
  drop = FALSE
]

message(
  "CID3586 matrix: ",
  nrow(cid3586_matrix),
  " genes × ",
  ncol(cid3586_matrix),
  " cells"
)

message(
  "CID4066 matrix: ",
  nrow(cid4066_matrix),
  " genes × ",
  ncol(cid4066_matrix),
  " cells"
)

# ------------------------------------------------------------
# 19. CREATE ANALYSIS-READY HER2+ METADATA
# ------------------------------------------------------------

her2_metadata <- metadata[
  selected_cells,
  ,
  drop = FALSE
]

her2_metadata$sample <- her2_metadata$orig.ident

her2_metadata$condition <- ifelse(
  her2_metadata$orig.ident == "CID3586",
  "Naive",
  ifelse(
    her2_metadata$orig.ident == "CID4066",
    "Treated",
    NA
  )
)

her2_metadata$disease_status <- "HER2_ER"

# ------------------------------------------------------------
# 20. SAVE SPARSE MATRICES AND METADATA
# ------------------------------------------------------------

message("")
message("==============================================")
message("SAVING HER2+ DATA")
message("==============================================")

saveRDS(
  cid3586_matrix,
  file.path(
    processed_her2_dir,
    "CID3586_expression_sparse.rds"
  ),
  compress = FALSE
)

saveRDS(
  cid4066_matrix,
  file.path(
    processed_her2_dir,
    "CID4066_expression_sparse.rds"
  ),
  compress = FALSE
)

saveRDS(
  her2_matrix,
  file.path(
    processed_her2_dir,
    "HER2_combined_expression_sparse.rds"
  ),
  compress = FALSE
)

saveRDS(
  her2_metadata,
  file.path(
    processed_her2_dir,
    "HER2_combined_metadata.rds"
  )
)

# ------------------------------------------------------------
# 21. SAVE VALIDATION TABLES
# ------------------------------------------------------------

write.csv(
  sample_counts,
  file.path(
    extraction_dir,
    "HER2_sample_cell_counts.csv"
  ),
  row.names = FALSE
)

validation_summary <- data.frame(

  matrix_genes = nrow(full_matrix),

  matrix_cells = ncol(full_matrix),

  matrix_nonzero = n_nonzero,

  CID3586_cells = ncol(cid3586_matrix),

  CID4066_cells = ncol(cid4066_matrix),

  HER2_total_cells = ncol(her2_matrix),

  genes_unique = !anyDuplicated(genes),

  barcodes_unique = !anyDuplicated(barcodes),

  matrix_gene_match = (
    nrow(full_matrix) ==
      length(genes)
  ),

  matrix_barcode_match = (
    ncol(full_matrix) ==
      length(barcodes)
  ),

  metadata_barcode_match = (
    nrow(metadata) ==
      length(barcodes)
  )
)

write.csv(
  validation_summary,
  file.path(
    extraction_dir,
    "HER2_extraction_QC_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  her2_metadata,
  file.path(
    extraction_dir,
    "HER2_combined_metadata.csv"
  ),
  row.names = TRUE
)

# ============================================================
# PART C — FINAL SUMMARY
# ============================================================

message("")
message("==============================================")
message("load_dataset.R COMPLETE")
message("==============================================")

message(
  "Normal datasets loaded: ",
  length(normal_matrices)
)

message(
  "Normal samples: ",
  paste(
    names(normal_matrices),
    collapse = ", "
  )
)

message(
  "HER2+ datasets: CID3586, CID4066"
)

message(
  "CID3586: ",
  ncol(cid3586_matrix),
  " cells — Naive"
)

message(
  "CID4066: ",
  ncol(cid4066_matrix),
  " cells — Treated"
)

message(
  "Total HER2+ cells: ",
  ncol(her2_matrix)
)

message("")
message(
  "Processed HER2+ matrices saved to:"
)

message(
  processed_her2_dir
)

message("")
message(
  "Validation results saved to:"
)

message(
  validation_dir
)

message(
  extraction_dir
)

message("")
message(
  "The full HER2 source matrix remains sparse."
)

# ------------------------------------------------------------
# 22. MEMORY CLEANUP
# ------------------------------------------------------------

rm(
  full_matrix,
  her2_matrix
)

gc()

message("")
message(
  "load_dataset.R finished successfully."
)