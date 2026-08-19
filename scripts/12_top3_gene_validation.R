
# ============================================================
# 12 — TOP-3 GENE-LEVEL VALIDATION
#
# Summarize directional DEG evidence supporting the final
# top-three computational drug candidates.
#
# Candidates:
#   1. Neratinib
#   2. Vorinostat
#   3. Bortezomib
#
# This script documents the existing validation results.
# No new DEG analysis is performed.
# ============================================================

library(dplyr)

input_dir <-
  "/mnt/lincs/analysis/13_top3_gene_validation"

output_dir <-
  "/mnt/lincs/analysis/13_top3_gene_validation"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# TOP-3 CANDIDATES
# ------------------------------------------------------------

top3_drugs <- c(
  "neratinib",
  "vorinostat",
  "bortezomib"
)

# ------------------------------------------------------------
# LOAD GENE OVERLAP TABLE
# ------------------------------------------------------------

overlap_file <- file.path(
  input_dir,
  "top3_gene_overlap.csv"
)

overlap <- read.csv(
  overlap_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "\nGene-overlap records:",
  nrow(overlap),
  "\n"
)

cat(
  "\nColumns:\n"
)

print(
  colnames(overlap)
)

# ------------------------------------------------------------
# LOAD EXISTING SUMMARY
# ------------------------------------------------------------

summary_file <- file.path(
  input_dir,
  "top3_gene_validation_summary.txt"
)

if (
  file.exists(summary_file)
) {

  cat(
    "\n============================================================\n"
  )

  cat(
    "EXISTING GENE-LEVEL VALIDATION SUMMARY\n"
  )

  cat(
    "============================================================\n\n"
  )

  cat(
    paste(
      readLines(summary_file),
      collapse = "\n"
    )
  )

  cat("\n")

}

# ------------------------------------------------------------
# LOAD DIRECTIONAL DEG LISTS
# ------------------------------------------------------------

gene_files <- list(
  neratinib =
    "neratinib_directional_DEGs.txt",

  vorinostat =
    "vorinostat_directional_DEGs.txt",

  bortezomib =
    "bortezomib_directional_DEGs.txt"
)

gene_counts <- data.frame(
  Drug = character(),
  Number.of.genes = integer(),
  stringsAsFactors = FALSE
)

for (
  drug in names(gene_files)
) {

  file_path <- file.path(
    input_dir,
    gene_files[[drug]]
  )

  if (
    file.exists(file_path)
  ) {

    genes <- readLines(
      file_path
    )

    genes <- genes[
      nzchar(
        trimws(genes)
      )
    ]

    gene_counts <- rbind(
      gene_counts,
      data.frame(
        Drug = drug,
        Number.of.genes = length(
          unique(genes)
        ),
        stringsAsFactors = FALSE
      )
    )

  }

}

# ------------------------------------------------------------
# SAVE GENE COUNT SUMMARY
# ------------------------------------------------------------

write.csv(
  gene_counts,
  file.path(
    output_dir,
    "top3_gene_counts.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# PRINT GENE COUNTS
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "TOP-3 GENE-LEVEL VALIDATION\n"
)

cat(
  "============================================================\n"
)

print(
  gene_counts
)

# ------------------------------------------------------------
# PRINT OVERLAP TABLE
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "TOP-3 GENE OVERLAP\n"
)

cat(
  "============================================================\n"
)

print(
  overlap
)

# ------------------------------------------------------------
# COMPLETE
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "TOP-3 GENE VALIDATION COMPLETE\n"
)

cat(
  "============================================================\n"
)

