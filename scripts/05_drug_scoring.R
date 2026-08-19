# ============================================================
# 05 — DRUG REPURPOSING AND DRUG SCORING
# ============================================================
# Retains the original HuLab workflow:
#   1. Generate drug-repurposing candidates from cell-type DE profiles
#   2. Calculate drug therapeutic scores
#   3. Select high-scoring candidates
#   4. Identify cell-type-specific drugs
#   5. Evaluate two-drug combinations
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Asgard)
})

INPUT_OBJECT <- "results/alignment/combined_aligned_object.rds"
INPUT_GENES <- "results/differential_expression/HER2_ER_genelist_limma.rds"
DRUG_REFERENCE_DIR <- "results/drug_reference"
OUTPUT_DIR <- "results/drug_repurposing"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

SC.merge <- readRDS(INPUT_OBJECT)
Gene.list <- readRDS(INPUT_GENES)

Case <- "CID4066"
Tissue <- "breast"

# The two Broad LINCS resources are external inputs.
GSE92742.gctx.path <- file.path(
  "data/lincs",
  "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
)
GSE70138.gctx.path <- file.path(
  "data/lincs",
  "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx"
)

# -----------------------------
# Load breast drug reference
# -----------------------------
my_gene_info <- read.table(
  file.path(DRUG_REFERENCE_DIR, "breast_gene_info.txt"),
  sep = "\t",
  header = TRUE,
  quote = ""
)

my_drug_info <- read.table(
  file.path(DRUG_REFERENCE_DIR, "breast_drug_info.txt"),
  sep = "\t",
  header = TRUE,
  quote = ""
)

cmap.ref.profiles <- GetDrugRef(
  drug.response.path = file.path(DRUG_REFERENCE_DIR, "breast_rankMatrix.txt"),
  probe.to.genes = my_gene_info,
  drug.info = my_drug_info
)

# -----------------------------
# Drug repurposing
# -----------------------------
Drug.ident.res <- GetDrug(
  gene.data = Gene.list,
  drug.ref.profiles = cmap.ref.profiles,
  repurposing.unit = "drug",
  connectivity = "negative",
  drug.type = "FDA"
)

saveRDS(
  Drug.ident.res,
  file.path(OUTPUT_DIR, "HER2_ER_drugs_FDA.rds")
)

# -----------------------------
# Drug scoring
# -----------------------------
Drug.score <- DrugScore(
  SC.integrated = SC.merge,
  Gene.data = Gene.list,
  Cell.type = NULL,
  Drug.data = Drug.ident.res,
  FDA.drug.only = TRUE,
  Case = Case,
  Tissue = Tissue,
  GSE92742.gctx = GSE92742.gctx.path,
  GSE70138.gctx = GSE70138.gctx.path
)

saveRDS(
  Drug.score,
  file.path(OUTPUT_DIR, "HER2_ER_drugscore_FDA.rds")
)

# -----------------------------
# High-scoring candidates
# -----------------------------
score_cutoff <- quantile(
  Drug.score$Drug.therapeutic.score,
  0.99,
  na.rm = TRUE
)

Final.drugs <- subset(
  Drug.score,
  Drug.therapeutic.score > score_cutoff & FDR < 0.05
)

saveRDS(
  Final.drugs,
  file.path(OUTPUT_DIR, "HER2_ER_final_drugs.rds")
)

# -----------------------------
# Cell-type-specific drugs
# -----------------------------
celltype_drugs <- TopDrug(
  SC.integrated = SC.merge,
  Drug.data = Drug.ident.res,
  Drug.FDR = 0.19,
  FDA.drug.only = TRUE,
  Case = Case
)

saveRDS(
  celltype_drugs,
  file.path(OUTPUT_DIR, "HER2_ER_celltype_drugs.rds")
)

# -----------------------------
# Drug combinations
# -----------------------------
Drug.combinations <- DrugCombination(
  SC.integrated = SC.merge,
  Gene.data = Gene.list,
  Drug.data = Drug.ident.res,
  Drug.FDR = 0.2,
  FDA.drug.only = TRUE,
  Combined.drugs = 2,
  Case = Case,
  Tissue = Tissue,
  GSE92742.gctx = GSE92742.gctx.path,
  GSE70138.gctx = GSE70138.gctx.path
)

Final.combinations <- TopCombination(
  Drug.combination = Drug.combinations,
  Combination.FDR = 0.05,
  Min.combination.score = 1
)

saveRDS(
  Final.combinations,
  file.path(OUTPUT_DIR, "HER2_ER_final_drug_combinations.rds")
)

message("Drug repurposing and scoring workflow complete.")
