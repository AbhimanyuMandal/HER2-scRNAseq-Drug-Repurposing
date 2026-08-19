
# ============================================================
# 11 — TOP-3 PATHWAY VALIDATION
#
# Summarize pathway-level evidence supporting the final
# top-three computational drug candidates.
#
# Candidates:
#   1. Neratinib
#   2. Vorinostat
#   3. Bortezomib
#
# This script does not perform new enrichment analysis.
# It organizes the pathway evidence generated previously.
# ============================================================

library(dplyr)

input_dir <-
  "/mnt/lincs/analysis/12_top3_validation"

output_dir <-
  "/mnt/lincs/analysis/12_top3_validation"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# LOAD PATHWAY THEME SUMMARY
# ------------------------------------------------------------

theme_file <- file.path(
  input_dir,
  "top3_pathway_theme_summary.csv"
)

themes <- read.csv(
  theme_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "\nTop-3 pathway theme records:",
  nrow(themes),
  "\n"
)

cat(
  "\nColumns:\n"
)

print(
  colnames(themes)
)

# ------------------------------------------------------------
# LOAD STRONGEST PATHWAYS
# ------------------------------------------------------------

strong_file <- file.path(
  input_dir,
  "top3_strongest_pathways.csv"
)

strongest <- read.csv(
  strong_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "\nStrongest pathway records:",
  nrow(strongest),
  "\n"
)

# ------------------------------------------------------------
# STANDARDIZE DRUG NAMES
# ------------------------------------------------------------

top3_drugs <- c(
  "neratinib",
  "vorinostat",
  "bortezomib"
)

themes <- themes %>%
  filter(
    Drug %in% top3_drugs
  )

strongest <- strongest %>%
  filter(
    Drug %in% top3_drugs
  )

# ------------------------------------------------------------
# SAVE CLEANED TABLES
# ------------------------------------------------------------

write.csv(
  themes,
  file.path(
    output_dir,
    "top3_pathway_theme_summary_clean.csv"
  ),
  row.names = FALSE
)

write.csv(
  strongest,
  file.path(
    output_dir,
    "top3_strongest_pathways_clean.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# PRINT SUMMARY
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "TOP-3 PATHWAY VALIDATION\n"
)

cat(
  "============================================================\n"
)

for (drug in top3_drugs) {

  cat(
    "\n--------------------------------------------\n"
  )

  cat(
    toupper(drug),
    "\n"
  )

  cat(
    "--------------------------------------------\n"
  )

  x <- themes %>%
    filter(
      Drug == drug
    )

  if (nrow(x) > 0) {

    print(x)

  } else {

    cat(
      "No pathway-theme records found.\n"
    )

  }
}

# ------------------------------------------------------------
# HUMAN-READABLE SUMMARY
# ------------------------------------------------------------

summary_file <- file.path(
  output_dir,
  "top3_biological_summary.txt"
)

con <- file(
  summary_file,
  open = "wt"
)

writeLines(
  c(
    "TOP-3 BIOLOGICAL PATHWAY VALIDATION",
    "===================================",
    "",
    "Final computational candidates:",
    "1. Neratinib",
    "2. Vorinostat",
    "3. Bortezomib",
    "",
    "Pathway-level evidence was summarized from the",
    "treatment-state pathway enrichment results.",
    ""
  ),
  con
)

for (drug in top3_drugs) {

  x <- themes %>%
    filter(
      Drug == drug
    )

  writeLines(
    paste0(
      drug,
      ":"
    ),
    con
  )

  if (nrow(x) > 0) {

    for (theme in unique(x$Theme)) {

      writeLines(
        paste0(
          "  - ",
          theme
        ),
        con
      )

    }

  }

  writeLines(
    "",
    con
  )
}

close(con)

cat(
  "\n============================================================\n"
)

cat(
  "TOP-3 PATHWAY VALIDATION COMPLETE\n"
)

cat(
  "============================================================\n"
)

