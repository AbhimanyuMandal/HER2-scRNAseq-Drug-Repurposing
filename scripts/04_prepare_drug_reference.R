# ============================================================
# 04 — PREPARE LINCS DRUG REFERENCE
# ============================================================
# Uses the Asgard PrepareReference() helper to construct the
# breast-tissue drug reference used by the original workflow.
# Large LINCS files are external inputs and are not committed to GitHub.
# ============================================================

suppressPackageStartupMessages(library(Asgard))

LINCS_DIR <- "data/lincs"
OUTPUT_DIR <- "results/drug_reference"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

required_files <- file.path(LINCS_DIR, c(
  "GSE92742_Broad_LINCS_cell_info.txt",
  "GSE92742_Broad_LINCS_gene_info.txt",
  "GSE70138_Broad_LINCS_sig_info.txt",
  "GSE92742_Broad_LINCS_sig_info.txt",
  "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx",
  "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
))

missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  stop(
    "Missing LINCS reference files:\n",
    paste(missing, collapse = "\n")
  )
}

PrepareReference(
  cell.info = required_files[1],
  gene.info = required_files[2],
  GSE70138.sig.info = required_files[3],
  GSE92742.sig.info = required_files[4],
  GSE70138.gctx = required_files[5],
  GSE92742.gctx = required_files[6],
  Output.Dir = paste0(OUTPUT_DIR, "/")
)

message("LINCS drug reference preparation complete.")
