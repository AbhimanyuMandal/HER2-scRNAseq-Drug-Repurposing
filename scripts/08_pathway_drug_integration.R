
# ============================================================
# 08 — PATHWAY–DRUG INTEGRATION
#
# Integrate:
#   1. Treatment-state pathway enrichment
#   2. Cell-type-specific drug evidence
#   3. Treatment-contrast therapeutic drug scores
#
# Comparison:
#   CID4066 vs CID3586
# ============================================================

library(dplyr)

# ------------------------------------------------------------
# DIRECTORIES
# ------------------------------------------------------------

pathway_dir <- "/mnt/lincs/analysis/07_pathway_analysis"

drug_dir <- "/mnt/lincs/analysis/06_drug_repurposing"

output_dir <- "/mnt/lincs/analysis/08_pathway_drug_integration"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 1. LOAD TREATMENT PATHWAY RESULTS
# ============================================================

pathway_files <- list.files(
  pathway_dir,
  pattern = "^CID4066_vs_CID3586_.*_pathways\\.csv$",
  full.names = TRUE
)

if (length(pathway_files) == 0) {
  stop("No CID4066 vs CID3586 pathway files found.")
}

pathway_results <- lapply(
  pathway_files,
  function(f) {

    x <- read.csv(
      f,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    x$Source.file <- basename(f)

    x
  }
)

pathway_data <- bind_rows(
  pathway_results
)

cat(
  "\nTreatment pathway records:",
  nrow(pathway_data),
  "\n"
)

# ============================================================
# 2. EXTRACT CELL TYPE AND DIRECTION
# ============================================================

pathway_data$Cell.type <- sub(
  "^CID4066_vs_CID3586_(.*)_(UP|DOWN)_pathways\\.csv$",
  "\\1",
  pathway_data$Source.file
)

pathway_data$Direction <- sub(
  "^CID4066_vs_CID3586_.*_(UP|DOWN)_pathways\\.csv$",
  "\\1",
  pathway_data$Source.file
)

# ============================================================
# 3. STANDARDIZE PATHWAY NAME
# ============================================================

if ("name" %in% colnames(pathway_data)) {

  pathway_data$Pathway <-
    pathway_data$name

} else if ("term_name" %in% colnames(pathway_data)) {

  pathway_data$Pathway <-
    pathway_data$term_name

} else if ("Pathway" %in% colnames(pathway_data)) {

  pathway_data$Pathway <-
    pathway_data$Pathway

} else {

  pathway_data$Pathway <-
    pathway_data$term_id
}

# Make sure p-value is numeric

if ("p_value" %in% colnames(pathway_data)) {

  pathway_data$p_value <-
    as.numeric(pathway_data$p_value)

}

# ============================================================
# 4. LOAD CELL-TYPE DRUG RESULTS
# ============================================================

celltype_file <- file.path(
  drug_dir,
  "CID4066_vs_CID3586_celltype_drugs.csv"
)

celltype_drugs <- read.csv(
  celltype_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "Cell-type drug records:",
  nrow(celltype_drugs),
  "\n"
)

cat(
  "\nCell-type drug columns:\n"
)

print(
  colnames(celltype_drugs)
)

if (!"Drug" %in% colnames(celltype_drugs)) {

  stop(
    "Drug column not found in cell-type drug results."
  )

}

if (!"Cell.type" %in% colnames(celltype_drugs)) {

  stop(
    "Cell.type column not found in cell-type drug results."
  )

}

# ============================================================
# 5. LOAD AUTHORITATIVE THERAPEUTIC SCORES
# ============================================================

score_file <- file.path(
  drug_dir,
  "CID4066_vs_CID3586_final_drugs.rds"
)

drug_scores <- readRDS(
  score_file
)

drug_scores <- as.data.frame(
  drug_scores
)

drug_scores$Drug <- rownames(
  drug_scores
)

drug_scores <- drug_scores %>%

  select(
    Drug,
    Drug.therapeutic.score,
    P.value,
    FDR
  )

cat(
  "\nTreatment-contrast drug scores:\n"
)

print(
  drug_scores
)

# ============================================================
# 6. INTEGRATE PATHWAYS WITH CELL-TYPE DRUG EVIDENCE
# ============================================================

integration <- pathway_data %>%

  inner_join(
    celltype_drugs,
    by = "Cell.type"
  )

# ============================================================
# 7. ADD AUTHORITATIVE DRUG SCORES
# ============================================================

integration <- integration %>%

  select(
    -any_of(
      c(
        "Drug.therapeutic.score",
        "P.value",
        "FDR",
        "Drug.FDR"
      )
    )
  ) %>%

  left_join(
    drug_scores,
    by = "Drug"
  )

# ============================================================
# 8. KEEP RELEVANT COLUMNS
# ============================================================

preferred_columns <- c(
  "Pathway",
  "term_id",
  "source",
  "p_value",
  "intersection_size",
  "query_size",
  "Cell.type",
  "Direction",
  "Drug",
  "Cell.type.coverage",
  "Drug.coverage",
  "Drug.therapeutic.score",
  "P.value",
  "FDR"
)

keep_columns <- intersect(
  preferred_columns,
  colnames(integration)
)

integration <- integration[
  ,
  keep_columns,
  drop = FALSE
]

# Rename drug-level statistical columns for clarity

if ("P.value" %in% colnames(integration)) {

  colnames(integration)[
    colnames(integration) == "P.value"
  ] <- "Drug.P.value"

}

if ("FDR" %in% colnames(integration)) {

  colnames(integration)[
    colnames(integration) == "FDR"
  ] <- "Drug.FDR"

}

# ============================================================
# 9. SAVE COMPLETE PATHWAY–DRUG TABLE
# ============================================================

write.csv(
  integration,
  file.path(
    output_dir,
    "pathway_drug_integration_treatment.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 10. DRUG SUMMARY
# ============================================================

drug_summary <- integration %>%

  group_by(
    Drug
  ) %>%

  summarise(

    Cell.types =
      paste(
        unique(Cell.type),
        collapse = "; "
      ),

    Number.of.celltypes =
      n_distinct(Cell.type),

    Number.of.pathways =
      n_distinct(Pathway),

    Best.pathway.p =
      min(
        p_value,
        na.rm = TRUE
      ),

    Drug.therapeutic.score =
      first(
        Drug.therapeutic.score
      ),

    Drug.P.value =
      first(
        Drug.P.value
      ),

    Drug.FDR =
      first(
        Drug.FDR
      ),

    .groups = "drop"

  ) %>%

  arrange(
    desc(
      Drug.therapeutic.score
    )
  )

write.csv(
  drug_summary,
  file.path(
    output_dir,
    "treatment_drug_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 11. DRUG–CELL-TYPE ASSOCIATION
# ============================================================

celltype_summary <- integration %>%

  group_by(
    Drug,
    Cell.type
  ) %>%

  summarise(

    Number.of.pathways =
      n_distinct(Pathway),

    Best.pathway.p =
      min(
        p_value,
        na.rm = TRUE
      ),

    Drug.therapeutic.score =
      first(
        Drug.therapeutic.score
      ),

    .groups = "drop"

  )

write.csv(
  celltype_summary,
  file.path(
    output_dir,
    "treatment_drug_celltype_association.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 12. PATHWAY SUMMARY
# ============================================================

pathway_summary <- integration %>%

  group_by(
    Pathway,
    Direction,
    Cell.type
  ) %>%

  summarise(

    Number.of.drugs =
      n_distinct(Drug),

    Best.p.value =
      min(
        p_value,
        na.rm = TRUE
      ),

    .groups = "drop"

  )

write.csv(
  pathway_summary,
  file.path(
    output_dir,
    "treatment_pathway_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 13. SUMMARY
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "PATHWAY–DRUG INTEGRATION COMPLETE\n"
)

cat(
  "============================================\n"
)

cat(
  "\nIntegrated records:",
  nrow(integration),
  "\n"
)

cat(
  "Drugs represented:",
  n_distinct(integration$Drug),
  "\n"
)

cat(
  "Cell types represented:",
  n_distinct(integration$Cell.type),
  "\n"
)

cat(
  "\nDrug summary:\n"
)

print(
  drug_summary
)

cat(
  "\nSaved:\n",
  file.path(
    output_dir,
    "pathway_drug_integration_treatment.csv"
  ),
  "\n",
  file.path(
    output_dir,
    "treatment_drug_summary.csv"
  ),
  "\n",
  file.path(
    output_dir,
    "treatment_drug_celltype_association.csv"
  ),
  "\n",
  file.path(
    output_dir,
    "treatment_pathway_summary.csv"
  ),
  "\n"
)

