
# ============================================================
# 13 — GENERATE PORTFOLIO FIGURES
#
# Purpose:
#   Generate the final figures documenting the reconstructed
#   HER2 scRNA-seq drug-repurposing analysis.
#
# Figures:
#   1. Analysis workflow
#   2. Drug landscape
#   3. Final drug ranking
#   4. Top-3 biological evidence
#
# No new biological analysis is performed.
# ============================================================

library(ggplot2)
library(dplyr)

# ============================================================
# DIRECTORIES
# ============================================================

ranking_dir <-
  "/mnt/lincs/analysis/11_final_evidence_scoring"

validation_dir <-
  "/mnt/lincs/analysis/12_top3_validation"

out_dir <-
  "/mnt/lincs/analysis/14_figures"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# FIGURE 1 — WORKFLOW
# ============================================================

workflow <- data.frame(
  x = 1,
  y = 8:1,
  label = c(
    "scRNA-seq breast cancer data",
    "Cell-type-specific analysis",
    "HER2/ER treatment-state comparisons",
    "Cell-type-specific DEG analysis",
    "LINCS drug repurposing",
    "Comparative drug prioritization",
    "Pathway + gene-level validation",
    "Final candidate drugs"
  )
)

p1 <- ggplot(
  workflow,
  aes(
    x = x,
    y = y
  )
) +

  geom_label(
    aes(
      label = label
    ),
    size = 4,
    linewidth = 0.4
  ) +

  scale_x_continuous(
    limits = c(
      0.5,
      1.5
    )
  ) +

  scale_y_continuous(
    breaks = NULL
  ) +

  labs(
    title =
      "HER2 scRNA-seq Drug Repurposing Workflow"
  ) +

  theme_void() +

  theme(
    plot.title =
      element_text(
        hjust = 0.5,
        face = "bold",
        size = 16
      )
  )

ggsave(
  file.path(
    out_dir,
    "Figure1_workflow.png"
  ),
  p1,
  width = 8,
  height = 10,
  dpi = 300
)

ggsave(
  file.path(
    out_dir,
    "Figure1_workflow.pdf"
  ),
  p1,
  width = 8,
  height = 10
)

# ============================================================
# FIGURE 2 — DRUG LANDSCAPE
# ============================================================

candidate_file <- file.path(
  ranking_dir,
  "final_drug_evidence_scores.csv"
)

ranking <- read.csv(
  candidate_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ranking <- ranking %>%
  arrange(
    Final.rank
  )

p2 <- ggplot(
  ranking,
  aes(
    x = reorder(
      Drug,
      Drug.therapeutic.score
    ),
    y = Drug.therapeutic.score
  )
) +

  geom_col(
    width = 0.7
  ) +

  coord_flip() +

  labs(
    title =
      "Therapeutic Drug Landscape",

    x =
      "Drug",

    y =
      "Therapeutic score"
  ) +

  theme_minimal(
    base_size = 12
  ) +

  theme(
    plot.title =
      element_text(
        face = "bold"
      )
  )

ggsave(
  file.path(
    out_dir,
    "Figure2_drug_landscape.png"
  ),
  p2,
  width = 9,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(
    out_dir,
    "Figure2_drug_landscape.pdf"
  ),
  p2,
  width = 9,
  height = 6
)

# ============================================================
# FIGURE 3 — FINAL DRUG RANKING
# ============================================================

top_ranking <- ranking %>%
  arrange(
    Final.rank
  ) %>%
  head(7)

p3 <- ggplot(
  top_ranking,
  aes(
    x = reorder(
      Drug,
      Final.evidence.score
    ),
    y = Final.evidence.score
  )
) +

  geom_col(
    width = 0.7
  ) +

  coord_flip() +

  labs(
    title =
      "Final Evidence-Based Drug Ranking",

    x =
      "Drug",

    y =
      "Final evidence score"
  ) +

  theme_minimal(
    base_size = 12
  ) +

  theme(
    plot.title =
      element_text(
        face = "bold"
      )
  )

ggsave(
  file.path(
    out_dir,
    "Figure3_final_drug_ranking.png"
  ),
  p3,
  width = 9,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(
    out_dir,
    "Figure3_final_drug_ranking.pdf"
  ),
  p3,
  width = 9,
  height = 6
)

# ============================================================
# FIGURE 4 — TOP-3 BIOLOGICAL EVIDENCE
# ============================================================

theme_file <- file.path(
  validation_dir,
  "top3_pathway_theme_summary.csv"
)

themes <- read.csv(
  theme_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

themes$Drug <- factor(
  themes$Drug,
  levels = c(
    "neratinib",
    "vorinostat",
    "bortezomib"
  )
)

# Count pathway themes per drug

theme_counts <- themes %>%
  group_by(
    Drug,
    Theme
  ) %>%
  summarise(
    Pathway.count =
      sum(
        Pathway.count,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

p4 <- ggplot(
  theme_counts,
  aes(
    x = Drug,
    y = Pathway.count,
    fill = Theme
  )
) +

  geom_col(
    position = "stack"
  ) +

  labs(
    title =
      "Top-3 Candidates: Biological Pathway Evidence",

    x =
      "Drug",

    y =
      "Number of enriched pathways",

    fill =
      "Pathway theme"
  ) +

  theme_minimal(
    base_size = 12
  ) +

  theme(
    plot.title =
      element_text(
        face = "bold"
      ),
    axis.text.x =
      element_text(
        angle = 30,
        hjust = 1
      )
  )

ggsave(
  file.path(
    out_dir,
    "Figure4_top3_biological_evidence.png"
  ),
  p4,
  width = 10,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(
    out_dir,
    "Figure4_top3_biological_evidence.pdf"
  ),
  p4,
  width = 10,
  height = 6
)

# ============================================================
# COMPLETE
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "FIGURE GENERATION COMPLETE\n"
)

cat(
  "============================================================\n"
)

cat(
  "\nFigures saved to:\n",
  out_dir,
  "\n\n"
)

print(
  list.files(
    out_dir,
    pattern = "\\.(png|pdf)$",
    full.names = FALSE
  )
)

