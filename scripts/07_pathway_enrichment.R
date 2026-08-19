
# ============================================================
# 07 — PATHWAY ENRICHMENT
#
# Perform cell-type-specific pathway enrichment for the
# significant genes identified by the LIMMA analysis.
#
# Comparisons:
#   CID4066 vs Normal
#   CID3586 vs Normal
#   CID4066 vs CID3586
#
# Significant genes:
#   adjusted P-value < 0.05
#
# Enrichment:
#   g:Profiler
# ============================================================

library(gprofiler2)

# ------------------------------------------------------------
# DIRECTORIES
# ------------------------------------------------------------

deg_dir <- "/mnt/lincs/analysis/03_DEG"

output_dir <- "/mnt/lincs/analysis/07_pathway_analysis"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# INPUT FILES
# ------------------------------------------------------------

comparisons <- list(

  CID4066_vs_Normal =
    "HER2_ER_CID4066_vs_Normal_limma.rds",

  CID3586_vs_Normal =
    "HER2_ER_CID3586_vs_Normal_limma.rds",

  CID4066_vs_CID3586 =
    "HER2_ER_CID4066_vs_CID3586_limma.rds"
)

# ------------------------------------------------------------
# HELPER FUNCTION
# ------------------------------------------------------------

run_gprofiler <- function(
    genes,
    comparison,
    cell_type,
    direction
) {

  if (length(genes) < 5) {

    cat(
      "Skipping:",
      comparison,
      cell_type,
      direction,
      "- too few genes\n"
    )

    return(NULL)
  }

  cat(
    "Running g:Profiler for",
    direction,
    "genes...\n"
  )

  result <- gost(
    query = genes,
    organism = "hsapiens",
    sources = c(
      "GO:BP",
      "GO:MF",
      "GO:CC",
      "REAC",
      "KEGG"
    ),
    user_threshold = 0.05,
    correction_method = "g_SCS",
    significant = TRUE
  )

  if (is.null(result$result)) {

    cat(
      "No significant pathways found.\n"
    )

    return(NULL)
  }

  out <- result$result

  # ----------------------------------------------------------
  # Convert list-columns to character
  #
  # g:Profiler can return list columns. Base R's write.csv()
  # cannot directly serialize these.
  # ----------------------------------------------------------

  out[] <- lapply(
    out,
    function(x) {

      if (is.list(x)) {

        vapply(
          x,
          function(y) {

            paste(
              y,
              collapse = ";"
            )

          },
          character(1)
        )

      } else {

        x
      }
    }
  )

  filename <- paste0(
    comparison,
    "_",
    cell_type,
    "_",
    direction,
    "_pathways.csv"
  )

  # Make filename filesystem-safe
  filename <- gsub(
    "[^A-Za-z0-9_.-]",
    "_",
    filename
  )

  write.csv(
    out,
    file.path(
      output_dir,
      filename
    ),
    row.names = FALSE
  )

  cat(
    "Saved:",
    filename,
    "\n"
  )

  out
}

# ------------------------------------------------------------
# RUN ALL COMPARISONS
# ------------------------------------------------------------

for (comparison in names(comparisons)) {

  cat("\n============================================\n")
  cat("COMPARISON:", comparison, "\n")
  cat("============================================\n")

  x <- readRDS(
    file.path(
      deg_dir,
      comparisons[[comparison]]
    )
  )

  for (cell_type in names(x)) {

    cat(
      "\n--------------------------------------------\n"
    )

    cat(
      "CELL TYPE:",
      cell_type,
      "\n"
    )

    cat(
      "--------------------------------------------\n"
    )

    df <- x[[cell_type]]

    # --------------------------------------------------------
    # Significant genes
    # --------------------------------------------------------

    sig <- df[
      !is.na(df$adj.P.Val) &
      df$adj.P.Val < 0.05,
      ,
      drop = FALSE
    ]

    up <- rownames(
      sig[
        sig$score > 0,
        ,
        drop = FALSE
      ]
    )

    down <- rownames(
      sig[
        sig$score < 0,
        ,
        drop = FALSE
      ]
    )

    cat(
      "Significant genes:",
      nrow(sig),
      "\n"
    )

    cat(
      "UP:",
      length(up),
      "\n"
    )

    cat(
      "DOWN:",
      length(down),
      "\n"
    )

    # --------------------------------------------------------
    # UP-regulated genes
    # --------------------------------------------------------

    run_gprofiler(
      genes = up,
      comparison = comparison,
      cell_type = cell_type,
      direction = "UP"
    )

    # --------------------------------------------------------
    # DOWN-regulated genes
    # --------------------------------------------------------

    run_gprofiler(
      genes = down,
      comparison = comparison,
      cell_type = cell_type,
      direction = "DOWN"
    )
  }
}

# ------------------------------------------------------------
# COMPLETE
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "PATHWAY ENRICHMENT COMPLETE\n"
)

cat(
  "============================================\n"
)

cat(
  "\nOutput directory:\n",
  output_dir,
  "\n"
)

