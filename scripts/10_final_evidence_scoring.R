
# ============================================================
# 10 — FINAL EVIDENCE SCORING
#
# Integrate:
#   • therapeutic drug score
#   • statistical evidence
#   • candidate classification
#   • cell-type evidence
#   • non-redundant pathway evidence
#
# The non-redundant pathway table is used as the authoritative
# source for pathway-theme evidence.
# ============================================================

library(dplyr)

# ------------------------------------------------------------
# DIRECTORIES
# ------------------------------------------------------------

drug_dir <- "/mnt/lincs/analysis/06_drug_repurposing"

validation_dir <-
  "/mnt/lincs/analysis/10_comparative_validation"

ranking_dir <-
  "/mnt/lincs/analysis/11_final_evidence_scoring"

output_dir <-
  "/mnt/lincs/analysis/11_final_evidence_scoring"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 1. LOAD CANDIDATE CLASSIFICATION
# ============================================================

candidates <- read.csv(
  file.path(
    validation_dir,
    "comparative_candidate_classes.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Only candidates appearing in the treatment contrast

candidates <- candidates %>%
  filter(
    CID4066_vs_CID3586 == 1
  )

cat(
  "\nTreatment/transition candidates:",
  nrow(candidates),
  "\n"
)

# ============================================================
# 2. LOAD AUTHORITATIVE DRUG SCORES
# ============================================================

drug_data <- readRDS(
  file.path(
    drug_dir,
    "CID4066_vs_CID3586_final_drugs.rds"
  )
)

drug_data <- as.data.frame(
  drug_data
)

drug_data$Drug <-
  rownames(drug_data)

drug_data <- drug_data %>%
  select(
    Drug,
    Drug.therapeutic.score,
    P.value,
    FDR
  )

# ============================================================
# 3. LOAD ORIGINAL NON-REDUNDANT PATHWAY THEMES
# ============================================================

theme_file <- file.path(
  ranking_dir,
  "nonredundant_pathway_themes.csv"
)

themes <- read.csv(
  theme_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "Non-redundant pathway records:",
  nrow(themes),
  "\n"
)

# ============================================================
# 4. COLLAPSE REDUNDANT THEMES ACROSS CELL TYPES
#
# A biological theme is counted once per drug even if the same
# theme occurs in multiple cell types.
#
# This is important for Neratinib:
#
#   L1.1:
#     Metabolic
#     Other
#     Translation
#
#   L1.2:
#     Translation
#
# The second Translation signal is redundant and therefore
# does not create a fourth theme.
# ============================================================

drug_theme_summary <- themes %>%

  group_by(
    Drug
  ) %>%

  summarise(

    Number.of.celltypes =
      n_distinct(Cell.type),

    Number.of.themes =
      n_distinct(Theme),

    Strong.pathway.themes =
      n_distinct(
        Theme[
          Best.pathway.p < 0.05
        ]
      ),

    Themes =
      paste(
        sort(
          unique(Theme)
        ),
        collapse = "; "
      ),

    Cell.types =
      paste(
        sort(
          unique(Cell.type)
        ),
        collapse = "; "
      ),

    .groups = "drop"
  )

# ------------------------------------------------------------
# Save the summary actually used for scoring
# ------------------------------------------------------------

write.csv(
  drug_theme_summary,
  file.path(
    output_dir,
    "drug_pathway_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 5. MERGE CANDIDATES + DRUG SCORES + PATHWAY EVIDENCE
# ============================================================

evidence <- candidates %>%

  select(
    Drug,
    Candidate.class
  ) %>%

  left_join(
    drug_data,
    by = "Drug"
  ) %>%

  left_join(
    drug_theme_summary,
    by = "Drug"
  )

# ============================================================
# 6. CHECK FOR MISSING PATHWAY EVIDENCE
# ============================================================

evidence$Number.of.celltypes[
  is.na(evidence$Number.of.celltypes)
] <- 0

evidence$Number.of.themes[
  is.na(evidence$Number.of.themes)
] <- 0

evidence$Strong.pathway.themes[
  is.na(evidence$Strong.pathway.themes)
] <- 0

# ============================================================
# 7. NORMALIZE THERAPEUTIC SCORE
# ============================================================

score_range <- range(
  evidence$Drug.therapeutic.score,
  na.rm = TRUE
)

if (
  diff(score_range) > 0
) {

  evidence$Therapeutic.component <-
    (
      evidence$Drug.therapeutic.score -
      score_range[1]
    ) /
    diff(score_range)

} else {

  evidence$Therapeutic.component <- 1

}

# ============================================================
# 8. NORMALIZE STATISTICAL EVIDENCE
# ============================================================

evidence$Statistical.component <-
  -log10(
    pmax(
      evidence$FDR,
      1e-300
    )
  )

stat_range <- range(
  evidence$Statistical.component,
  na.rm = TRUE
)

if (
  diff(stat_range) > 0
) {

  evidence$Statistical.component <-
    (
      evidence$Statistical.component -
      stat_range[1]
    ) /
    diff(stat_range)

} else {

  evidence$Statistical.component <- 1

}

# ============================================================
# 9. CELL-TYPE EVIDENCE
# ============================================================

evidence$Celltype.component <-
  pmin(
    evidence$Number.of.celltypes / 2,
    1
  )

# ============================================================
# 10. PATHWAY EVIDENCE
# ============================================================

evidence$Pathway.component <-
  pmin(
    evidence$Strong.pathway.themes / 4,
    1
  )

# ============================================================
# 11. FINAL EVIDENCE SCORE
#
# Weighted integration:
#
#   Therapeutic score    40%
#   Statistical evidence 25%
#   Cell-type evidence   15%
#   Pathway evidence     20%
# ============================================================

evidence$Final.evidence.score <-

  0.40 *
  evidence$Therapeutic.component +

  0.25 *
  evidence$Statistical.component +

  0.15 *
  evidence$Celltype.component +

  0.20 *
  evidence$Pathway.component

# ============================================================
# 12. FINAL RANKING
# ============================================================

evidence <- evidence %>%

  arrange(
    desc(Final.evidence.score)
  ) %>%

  mutate(
    Final.rank = row_number()
  )

# ============================================================
# 13. SAVE FINAL EVIDENCE TABLE
# ============================================================

write.csv(
  evidence,
  file.path(
    output_dir,
    "final_drug_evidence_scores.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 14. TOP 3
# ============================================================

top3 <- evidence %>%

  slice_head(
    n = 3
  )

write.csv(
  top3,
  file.path(
    output_dir,
    "final_top3_drugs.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 15. PRINT FINAL RANKING
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "FINAL EVIDENCE-BASED DRUG RANKING\n"
)

cat(
  "============================================================\n"
)

print(
  evidence %>%

    select(
      Final.rank,
      Drug,
      Candidate.class,
      Drug.therapeutic.score,
      FDR,
      Number.of.celltypes,
      Number.of.themes,
      Strong.pathway.themes,
      Themes,
      Final.evidence.score
    )
)

cat(
  "\n============================================================\n"
)

cat(
  "TOP 3 CANDIDATES\n"
)

cat(
  "============================================================\n"
)

print(
  top3 %>%

    select(
      Final.rank,
      Drug,
      Candidate.class,
      Drug.therapeutic.score,
      FDR,
      Number.of.celltypes,
      Number.of.themes,
      Themes,
      Cell.types,
      Final.evidence.score
    )
)

# ============================================================
# 16. WRITE HUMAN-READABLE SUMMARY
# ============================================================

summary_file <- file.path(
  output_dir,
  "final_candidate_summary.txt"
)

con <- file(
  summary_file,
  open = "wt"
)

writeLines(
  c(
    "FINAL COMPUTATIONAL DRUG-REPURPOSING CANDIDATES",
    "==============================================",
    ""
  ),
  con
)

for (i in seq_len(nrow(top3))) {

  d <- top3[i, ]

  writeLines(
    c(
      paste0(
        "Rank ",
        d$Final.rank,
        ": ",
        d$Drug
      ),

      paste0(
        "Candidate class: ",
        d$Candidate.class
      ),

      paste0(
        "Therapeutic score: ",
        sprintf(
          "%.4f",
          d$Drug.therapeutic.score
        )
      ),

      paste0(
        "Drug FDR: ",
        format(
          d$FDR,
          scientific = TRUE
        )
      ),

      paste0(
        "Cell types: ",
        d$Cell.types
      ),

      paste0(
        "Pathway themes: ",
        d$Themes
      ),

      paste0(
        "Final evidence score: ",
        sprintf(
          "%.4f",
          d$Final.evidence.score
        )
      ),

      ""
    ),
    con
  )
}

close(con)

# ============================================================
# COMPLETE
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "FINAL EVIDENCE SCORING COMPLETE\n"
)

cat(
  "============================================================\n"
)

cat(
  "\nOutput directory:\n",
  output_dir,
  "\n"
)

