
# ============================================================
# 06 — COMPARATIVE DRUG ANALYSIS
#
# Compare drug candidates identified from:
#
#   1. CID4066 vs Normal
#   2. CID3586 vs Normal
#   3. CID4066 vs CID3586
#
# Purpose:
#   Identify shared disease-state candidates and drugs that
#   are specifically associated with the treatment-state
#   transition.
# ============================================================

library(dplyr)

# ------------------------------------------------------------
# DIRECTORIES
# ------------------------------------------------------------

drug_dir <- "/mnt/lincs/analysis/06_drug_repurposing"

output_dir <- "/mnt/lincs/analysis/06_drug_repurposing"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# INPUT FILES
# ------------------------------------------------------------

files <- list(

  CID4066_vs_Normal =
    "CID4066_vs_Normal_final_drugs.rds",

  CID3586_vs_Normal =
    "CID3586_vs_Normal_final_drugs.rds",

  CID4066_vs_CID3586 =
    "CID4066_vs_CID3586_final_drugs.rds"
)

# ------------------------------------------------------------
# LOAD RESULTS
# ------------------------------------------------------------

results <- lapply(
  files,
  function(f) {
    readRDS(
      file.path(
        drug_dir,
        f
      )
    )
  }
)

# ------------------------------------------------------------
# CONVERT TO COMPARABLE TABLES
# ------------------------------------------------------------

make_table <- function(
    x,
    comparison
) {

  x <- as.data.frame(x)

  x$Drug <- rownames(x)

  x$Comparison <- comparison

  x
}

tables <- Map(
  make_table,
  results,
  names(results)
)

all_drugs <- bind_rows(
  tables
)

# ------------------------------------------------------------
# FINAL DRUG COUNTS
# ------------------------------------------------------------

cat("\n============================================\n")
cat("FINAL DRUG COUNTS\n")
cat("============================================\n")

print(
  table(
    all_drugs$Comparison
  )
)

# ------------------------------------------------------------
# DRUG PRESENCE MATRIX
# ------------------------------------------------------------

presence <- all_drugs %>%

  select(
    Drug,
    Comparison
  ) %>%

  distinct() %>%

  mutate(
    Present = 1
  ) %>%

  tidyr::pivot_wider(
    names_from = Comparison,
    values_from = Present,
    values_fill = 0
  )

cat("\n============================================\n")
cat("DRUG PRESENCE MATRIX\n")
cat("============================================\n")

print(presence)

write.csv(
  presence,
  file.path(
    output_dir,
    "comparative_drug_presence.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# SHARED ACROSS ALL THREE
# ------------------------------------------------------------

shared_all <- presence %>%

  filter(
    CID4066_vs_Normal == 1,
    CID3586_vs_Normal == 1,
    CID4066_vs_CID3586 == 1
  )

cat("\n============================================\n")
cat("SHARED ACROSS ALL THREE\n")
cat("============================================\n")

print(shared_all)

# ------------------------------------------------------------
# TREATMENT-CONTRAST ONLY
#
# Candidate appears in CID4066 vs CID3586 but not in either
# disease-vs-normal comparison.
# ------------------------------------------------------------

treatment_only <- presence %>%

  filter(
    CID4066_vs_CID3586 == 1,
    CID4066_vs_Normal == 0,
    CID3586_vs_Normal == 0
  )

cat("\n============================================\n")
cat("TREATMENT-CONTRAST ONLY\n")
cat("============================================\n")

print(treatment_only)

# ------------------------------------------------------------
# SHARED BY BOTH DISEASE-vs-NORMAL ANALYSES
# ------------------------------------------------------------

shared_disease <- presence %>%

  filter(
    CID4066_vs_Normal == 1,
    CID3586_vs_Normal == 1
  )

cat("\n============================================\n")
cat("SHARED BY BOTH DISEASE-vs-NORMAL ANALYSES\n")
cat("============================================\n")

print(shared_disease)

# ------------------------------------------------------------
# SCORE COMPARISON
# ------------------------------------------------------------

score_tables <- lapply(
  tables,
  function(x) {

    x %>%

      select(
        Drug,
        Drug.therapeutic.score,
        P.value,
        FDR
      )

  }
)

score_comparison <- Reduce(
  function(x, y) {

    full_join(
      x,
      y,
      by = "Drug",
      suffix = c(
        "",
        ""
      )
    )

  },
  score_tables
)

# The Reduce approach above can produce duplicate column names.
# Construct the score table explicitly for clarity.

score_comparison <- all_drugs %>%

  select(
    Drug,
    Comparison,
    Drug.therapeutic.score,
    P.value,
    FDR
  ) %>%

  tidyr::pivot_wider(
    names_from = Comparison,
    values_from = c(
      Drug.therapeutic.score,
      P.value,
      FDR
    )
  )

cat("\n============================================\n")
cat("SCORE COMPARISON\n")
cat("============================================\n")

print(
  score_comparison
)

write.csv(
  score_comparison,
  file.path(
    output_dir,
    "comparative_drug_scores.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# SUMMARY
# ------------------------------------------------------------

cat("\n============================================\n")
cat("COMPARATIVE ANALYSIS COMPLETE\n")
cat("============================================\n")

cat(
  "\nShared across all three:",
  nrow(shared_all),
  "\n"
)

cat(
  "Treatment-contrast-specific:",
  nrow(treatment_only),
  "\n"
)

cat(
  "Shared by both disease-vs-normal analyses:",
  nrow(shared_disease),
  "\n"
)

cat(
  "\nSaved:\n",
  file.path(
    output_dir,
    "comparative_drug_presence.csv"
  ),
  "\n",
  file.path(
    output_dir,
    "comparative_drug_scores.csv"
  ),
  "\n"
)

