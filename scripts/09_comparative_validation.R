
# ============================================================
# 09 — COMPARATIVE CANDIDATE VALIDATION
#
# Classify drug candidates according to their occurrence across:
#
#   CID4066 vs Normal
#   CID3586 vs Normal
#   CID4066 vs CID3586
#
# Candidate classes:
#
#   1. Shared disease-state candidates
#   2. Treatment-specific candidates
#   3. Transition-associated candidates
#
# This step does not perform new drug scoring.
# It reorganizes and interprets the existing results.
# ============================================================

library(dplyr)

# ------------------------------------------------------------
# DIRECTORIES
# ------------------------------------------------------------

drug_dir <- "/mnt/lincs/analysis/06_drug_repurposing"

pathway_dir <- "/mnt/lincs/analysis/08_pathway_drug_integration"

output_dir <- "/mnt/lincs/analysis/10_comparative_validation"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 1. LOAD COMPARATIVE DRUG PRESENCE
# ============================================================

presence_file <- file.path(
  drug_dir,
  "comparative_drug_presence.csv"
)

presence <- read.csv(
  presence_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "\nCandidate presence table:",
  nrow(presence),
  "drugs\n"
)

# ============================================================
# 2. LOAD DRUG SCORES
# ============================================================

score_file <- file.path(
  drug_dir,
  "comparative_drug_scores.csv"
)

scores <- read.csv(
  score_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ============================================================
# 3. SHARED DISEASE-STATE CANDIDATES
#
# Present in both:
#
#   CID4066 vs Normal
#   CID3586 vs Normal
#
# but absent from the treatment contrast.
# ============================================================

shared_disease <- presence %>%

  filter(
    CID4066_vs_Normal == 1,
    CID3586_vs_Normal == 1,
    CID4066_vs_CID3586 == 0
  )

write.csv(
  shared_disease,
  file.path(
    output_dir,
    "shared_disease_state_candidates.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 4. TREATMENT-SPECIFIC CANDIDATES
#
# Present only in CID4066 vs CID3586.
# ============================================================

treatment_specific <- presence %>%

  filter(
    CID4066_vs_CID3586 == 1,
    CID4066_vs_Normal == 0,
    CID3586_vs_Normal == 0
  )

write.csv(
  treatment_specific,
  file.path(
    output_dir,
    "treatment_specific_candidates.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 5. TRANSITION-ASSOCIATED CANDIDATES
#
# Drugs present in the treatment contrast that are also
# associated with at least one disease-vs-normal comparison.
#
# These candidates may reflect a treatment-associated
# cellular transition rather than a purely treatment-specific
# signal.
# ============================================================

transition_associated <- presence %>%

  filter(
    CID4066_vs_CID3586 == 1,
    (
      CID4066_vs_Normal == 1 |
      CID3586_vs_Normal == 1
    )
  )

write.csv(
  transition_associated,
  file.path(
    output_dir,
    "transition_associated_candidates.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 6. COMBINED CANDIDATE CLASSIFICATION
# ============================================================

candidate_classes <- presence %>%

  mutate(

    Candidate.class = case_when(

      CID4066_vs_CID3586 == 1 &
      CID4066_vs_Normal == 0 &
      CID3586_vs_Normal == 0
        ~ "Treatment-contrast-specific",

      CID4066_vs_CID3586 == 1 &
      (
        CID4066_vs_Normal == 1 |
        CID3586_vs_Normal == 1
      )
        ~ "Transition-associated",

      CID4066_vs_Normal == 1 &
      CID3586_vs_Normal == 1 &
      CID4066_vs_CID3586 == 0
        ~ "Shared disease-state",

      TRUE
        ~ "Other"

    )

  )

# ------------------------------------------------------------
# Add treatment therapeutic scores
# ------------------------------------------------------------

treatment_scores <- readRDS(
  file.path(
    drug_dir,
    "CID4066_vs_CID3586_final_drugs.rds"
  )
)

treatment_scores <- as.data.frame(
  treatment_scores
)

treatment_scores$Drug <-
  rownames(treatment_scores)

treatment_scores <- treatment_scores %>%

  select(
    Drug,
    Drug.therapeutic.score,
    P.value,
    FDR
  )

candidate_classes <- candidate_classes %>%

  left_join(
    treatment_scores,
    by = "Drug"
  )

write.csv(
  candidate_classes,
  file.path(
    output_dir,
    "comparative_candidate_classes.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 7. TREATMENT-SPECIFIC PATHWAY EVIDENCE
# ============================================================

pathway_file <- file.path(
  pathway_dir,
  "pathway_drug_integration_treatment.csv"
)

if (file.exists(pathway_file)) {

  pathway_drug <- read.csv(
    pathway_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Keep only treatment-specific candidates

  treatment_drug_names <-
    treatment_specific$Drug

  treatment_pathway <- pathway_drug %>%

    filter(
      Drug %in% treatment_drug_names
    )

  write.csv(
    treatment_pathway,
    file.path(
      output_dir,
      "treatment_specific_pathway_evidence.csv"
    ),
    row.names = FALSE
  )

} else {

  cat(
    "\nPathway integration file not found; skipping pathway subset.\n"
  )

}

# ============================================================
# 8. PRINT RESULTS
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "COMPARATIVE CANDIDATE VALIDATION\n"
)

cat(
  "============================================\n"
)

cat(
  "\nShared disease-state candidates:",
  nrow(shared_disease),
  "\n"
)

print(
  shared_disease
)

cat(
  "\nTreatment-specific candidates:",
  nrow(treatment_specific),
  "\n"
)

print(
  treatment_specific
)

cat(
  "\nTransition-associated candidates:",
  nrow(transition_associated),
  "\n"
)

print(
  transition_associated
)

cat(
  "\nComplete candidate classification:\n"
)

print(
  candidate_classes
)

# ============================================================
# COMPLETE
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "COMPARATIVE VALIDATION COMPLETE\n"
)

cat(
  "============================================\n"
)

cat(
  "\nOutput directory:\n",
  output_dir,
  "\n"
)

