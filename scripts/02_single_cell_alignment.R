# ============================================================
# HU LAB REPRODUCTION PROJECT
# SCRIPT 02 — SINGLE-CELL ALIGNMENT
#
# Reproduces the original HuLab SC-alignment workflow while
# including BOTH HER2+/ER+ samples:
#
#   CID3586 = Naive
#   CID4066 = Treated
#
# Reference:
#   GSE113197 normal breast epithelial samples
#   Ind5, Ind6, Ind7
#
# Original HuLab workflow:
#   NormalizeData
#   FindVariableFeatures (5,000)
#   FindIntegrationAnchors (dims 1:30)
#   IntegrateData (dims 1:30)
#   ScaleData
#   PCA (30 PCs)
#   UMAP (dims 1:30)
#   FindTransferAnchors
#   TransferData
#   MapQuery
#   merge reference + query
#   ScaleData
#   PCA
#   FindNeighbors
#   FindClusters
#   UMAP
#
# Extension:
#   Both CID3586 and CID4066 are used as the query dataset.
#   Their original sample identities and Naive/Treated status
#   are retained in metadata.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
})

# ------------------------------------------------------------
# 0. PROJECT PATHS
# ------------------------------------------------------------

project_dir <- getwd()

raw_dir <- file.path(
  project_dir,
  "data",
  "raw"
)

her2_processed_dir <- file.path(
  project_dir,
  "data",
  "processed",
  "HER2"
)

result_dir <- file.path(
  project_dir,
  "results",
  "02_single_cell_alignment"
)

object_dir <- file.path(
  project_dir,
  "data",
  "processed",
  "Seurat"
)

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  object_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

message("==============================================")
message("SCRIPT 02 — SINGLE-CELL ALIGNMENT")
message("==============================================")
message("Project directory: ", project_dir)
message("")

# ------------------------------------------------------------
# 1. CHECK REQUIRED INPUTS
# ------------------------------------------------------------

required_files <- c(
  file.path(
    raw_dir,
    "GSM3099847_Ind5_Expression_Matrix.txt"
  ),
  file.path(
    raw_dir,
    "GSM3099848_Ind6_Expression_Matrix.txt"
  ),
  file.path(
    raw_dir,
    "GSM3099849_Ind7_Expression_Matrix.txt"
  ),
  file.path(
    raw_dir,
    "Normal_celltype.txt"
  ),
  file.path(
    her2_processed_dir,
    "CID3586_expression_sparse.rds"
  ),
  file.path(
    her2_processed_dir,
    "CID4066_expression_sparse.rds"
  ),
  file.path(
    her2_processed_dir,
    "HER2_combined_metadata.rds"
  )
)

for (f in required_files) {

  if (!file.exists(f)) {
    stop(
      "Required input not found:\n",
      f
    )
  }

  message("Found: ", f)
}

# ============================================================
# PART A — CONSTRUCT THE NORMAL REFERENCE
# ============================================================

# ------------------------------------------------------------
# 2. LOAD NORMAL ANNOTATION
# ------------------------------------------------------------

message("")
message("==============================================")
message("LOADING NORMAL REFERENCE")
message("==============================================")

normal_annotation <- read.table(
  file.path(
    raw_dir,
    "Normal_celltype.txt"
  ),
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Original HuLab workflow uses only these four epithelial
# populations for the reference.
reference_celltypes <- c(
  "Luminal_L2_epithelial_cells",
  "Luminal_L1.1_epithelial_cells",
  "Luminal_L1.2_epithelial_cells",
  "Basal_epithelial_cells"
)

# ------------------------------------------------------------
# 3. HELPER FUNCTION — CREATE NORMAL SEURAT OBJECT
# ------------------------------------------------------------

create_normal_object <- function(
  sample_name,
  matrix_file,
  object_name
) {

  message("")
  message("----------------------------------------------")
  message(
    "Creating reference object: ",
    sample_name
  )
  message("----------------------------------------------")

  mat <- fread(
    matrix_file,
    sep = "\t",
    header = TRUE,
    data.table = FALSE,
    check.names = FALSE
  )

  gene_names <- mat[[1]]

  counts <- as.matrix(
    mat[, -1, drop = FALSE]
  )

  rownames(counts) <- make.unique(
    gene_names
  )

  counts <- matrix(
    as.numeric(counts),
    nrow = nrow(counts),
    ncol = ncol(counts),
    dimnames = dimnames(counts)
  )

  annotation <- subset(
    normal_annotation,
    sample == sample_name &
      celltype %in% reference_celltypes
  )

  common_cells <- intersect(
    colnames(counts),
    rownames(annotation)
  )

  message(
    "Expression cells: ",
    ncol(counts)
  )

  message(
    "Annotated epithelial cells: ",
    nrow(annotation)
  )

  message(
    "Matched epithelial cells: ",
    length(common_cells)
  )

  if (length(common_cells) == 0) {
    stop(
      "No common cells found for ",
      sample_name
    )
  }

  counts <- counts[
    ,
    common_cells,
    drop = FALSE
  ]

  annotation <- annotation[
    common_cells,
    ,
    drop = FALSE
  ]

  metadata <- data.frame(
    celltype = annotation$celltype,
    sample = sample_name,
    type = "Normal",
    row.names = common_cells,
    stringsAsFactors = FALSE
  )

  obj <- CreateSeuratObject(
    counts = counts,
    project = "Epithelial",
    min.cells = 3,
    min.features = 200,
    meta.data = metadata
  )

  message(
    object_name,
    " created: ",
    ncol(obj),
    " cells × ",
    nrow(obj),
    " genes"
  )

  return(obj)
}

# ------------------------------------------------------------
# 4. CREATE IND5 / IND6 / IND7 REFERENCE OBJECTS
# ------------------------------------------------------------

Epithelial2 <- create_normal_object(
  sample_name = "Ind5",
  matrix_file = file.path(
    raw_dir,
    "GSM3099847_Ind5_Expression_Matrix.txt"
  ),
  object_name = "Epithelial2"
)

Epithelial3 <- create_normal_object(
  sample_name = "Ind6",
  matrix_file = file.path(
    raw_dir,
    "GSM3099848_Ind6_Expression_Matrix.txt"
  ),
  object_name = "Epithelial3"
)

Epithelial4 <- create_normal_object(
  sample_name = "Ind7",
  matrix_file = file.path(
    raw_dir,
    "GSM3099849_Ind7_Expression_Matrix.txt"
  ),
  object_name = "Epithelial4"
)

# ============================================================
# PART B — CREATE CID3586 + CID4066 QUERY OBJECT
# ============================================================

# ------------------------------------------------------------
# 5. LOAD HER2+ SPARSE MATRICES AND METADATA
# ------------------------------------------------------------

message("")
message("==============================================")
message("LOADING HER2+ QUERY DATA")
message("==============================================")

cid3586_counts <- readRDS(
  file.path(
    her2_processed_dir,
    "CID3586_expression_sparse.rds"
  )
)

cid4066_counts <- readRDS(
  file.path(
    her2_processed_dir,
    "CID4066_expression_sparse.rds"
  )
)

her2_metadata <- readRDS(
  file.path(
    her2_processed_dir,
    "HER2_combined_metadata.rds"
  )
)

message(
  "CID3586: ",
  ncol(cid3586_counts),
  " cells"
)

message(
  "CID4066: ",
  ncol(cid4066_counts),
  " cells"
)

# ------------------------------------------------------------
# 6. ADD ANALYSIS METADATA
# ------------------------------------------------------------

# Ensure metadata matches each matrix.

metadata_3586 <- her2_metadata[
  colnames(cid3586_counts),
  ,
  drop = FALSE
]

metadata_4066 <- her2_metadata[
  colnames(cid4066_counts),
  ,
  drop = FALSE
]

metadata_3586$sample <- "CID3586"
metadata_3586$condition <- "Naive"
metadata_3586$disease_status <- "HER2_ER"

metadata_4066$sample <- "CID4066"
metadata_4066$condition <- "Treated"
metadata_4066$disease_status <- "HER2_ER"

# ------------------------------------------------------------
# 7. CREATE SEURAT OBJECTS
# ------------------------------------------------------------

message("")
message("Creating CID3586 query object...")

CID3586 <- CreateSeuratObject(
  counts = cid3586_counts,
  project = "HER2_ER",
  min.cells = 3,
  min.features = 200,
  meta.data = metadata_3586
)

message(
  "CID3586 Seurat object: ",
  ncol(CID3586),
  " cells"
)

message("")
message("Creating CID4066 query object...")

CID4066 <- CreateSeuratObject(
  counts = cid4066_counts,
  project = "HER2_ER",
  min.cells = 3,
  min.features = 200,
  meta.data = metadata_4066
)

message(
  "CID4066 Seurat object: ",
  ncol(CID4066),
  " cells"
)

# ------------------------------------------------------------
# 8. COMBINE THE TWO QUERY SAMPLES
# ------------------------------------------------------------

message("")
message("==============================================")
message("COMBINING HER2+ QUERY SAMPLES")
message("==============================================")

SC.query <- merge(
  CID3586,
  CID4066
)

message(
  "Combined query: ",
  ncol(SC.query),
  " cells × ",
  nrow(SC.query),
  " genes"
)

message("")
message("Query sample distribution:")
print(
  table(
    SC.query$sample
  )
)

message("")
message("Query condition distribution:")
print(
  table(
    SC.query$condition
  )
)

# ============================================================
# PART C — ORIGINAL HULAB SC-ALIGNMENT
# ============================================================

# ------------------------------------------------------------
# 9. CREATE SC.list
# ------------------------------------------------------------

SC.list <- list(
  HER2_ER = SC.query,
  Epithelial2 = Epithelial2,
  Epithelial3 = Epithelial3,
  Epithelial4 = Epithelial4
)

CellCycle <- TRUE

anchor.features <- 5000

# ------------------------------------------------------------
# 10. NORMALIZE + FIND VARIABLE FEATURES
# ------------------------------------------------------------

message("")
message("==============================================")
message("NORMALIZATION AND VARIABLE FEATURES")
message("==============================================")

for (i in 1:length(SC.list)) {

  message(
    "Processing dataset ",
    i,
    " / ",
    length(SC.list)
  )

  SC.list[[i]] <- NormalizeData(
    SC.list[[i]],
    verbose = FALSE
  )

  SC.list[[i]] <- FindVariableFeatures(
    SC.list[[i]],
    selection.method = "vst",
    nfeatures = anchor.features,
    verbose = FALSE
  )
}

# Keep the named objects synchronized with SC.list.

SC.query <- SC.list[["HER2_ER"]]
Epithelial2 <- SC.list[["Epithelial2"]]
Epithelial3 <- SC.list[["Epithelial3"]]
Epithelial4 <- SC.list[["Epithelial4"]]

# ------------------------------------------------------------
# 11. INTEGRATE NORMAL REFERENCE
# ------------------------------------------------------------

message("")
message("==============================================")
message("REFERENCE INTEGRATION")
message("==============================================")

SC.reference <- list(
  Epithelial2 = Epithelial2,
  Epithelial3 = Epithelial3,
  Epithelial4 = Epithelial4
)

SC.anchors <- FindIntegrationAnchors(
  object.list = SC.reference,
  dims = 1:30
)

SC.integrated <- IntegrateData(
  anchorset = SC.anchors,
  dims = 1:30
)

message(
  "Integrated reference cells: ",
  ncol(SC.integrated)
)

# ------------------------------------------------------------
# 12. SCALE / PCA / UMAP REFERENCE
# ------------------------------------------------------------

message("")
message("==============================================")
message("REFERENCE PCA / UMAP")
message("==============================================")

DefaultAssay(
  SC.integrated
) <- "integrated"

SC.integrated <- ScaleData(
  SC.integrated,
  verbose = FALSE
)

SC.integrated <- RunPCA(
  SC.integrated,
  npcs = 30,
  verbose = FALSE
)

SC.integrated <- RunUMAP(
  SC.integrated,
  reduction = "pca",
  dims = 1:30,
  verbose = FALSE
)

# ------------------------------------------------------------
# 13. REFERENCE VISUALIZATION
# ------------------------------------------------------------

p_reference <- DimPlot(
  SC.integrated,
  reduction = "umap",
  split.by = "sample",
  group.by = "celltype"
) +
  ggtitle(
    "Reference annotations"
  )

ggsave(
  filename = file.path(
    result_dir,
    "reference_annotations.pdf"
  ),
  plot = p_reference,
  width = 12,
  height = 8
)

# ------------------------------------------------------------
# 14. FIND TRANSFER ANCHORS FOR BOTH HER2+ SAMPLES
# ------------------------------------------------------------

message("")
message("==============================================")
message("QUERY LABEL TRANSFER")
message("==============================================")

SC.query <- JoinLayers(
  SC.query,
  assay = "RNA"
)


SC.anchors <- FindTransferAnchors(
  reference = SC.integrated,
  query = SC.query,
  dims = 1:30,
  reference.reduction = "pca"
)

# ------------------------------------------------------------
# 15. TRANSFER REFERENCE CELL-TYPE LABELS
# ------------------------------------------------------------

predictions <- TransferData(
  anchorset = SC.anchors,
  refdata = SC.integrated$celltype,
  dims = 1:30
)

SC.query <- AddMetaData(
  SC.query,
  metadata = predictions
)

message("")
message("Transferred label distribution:")
print(
  table(
    SC.query$predicted.id
  )
)

# ------------------------------------------------------------
# 16. MAP QUERY ONTO REFERENCE UMAP
# ------------------------------------------------------------

SC.integrated <- RunUMAP(
  SC.integrated,
  dims = 1:30,
  reduction = "pca",
  return.model = TRUE
)

SC.query <- MapQuery(
  anchorset = SC.anchors,
  reference = SC.integrated,
  query = SC.query,
  refdata = list(
    celltype = "celltype"
  ),
  reference.reduction = "pca",
  reduction.model = "umap"
)

# ------------------------------------------------------------
# 17. QUERY VISUALIZATION
# ------------------------------------------------------------

if ("celltype" %in% colnames(SC.query@meta.data)) {

  # Already present

} else if ("predicted.celltype" %in% colnames(SC.query@meta.data)) {

  SC.query$celltype <- SC.query$predicted.celltype

} else if ("predicted.id" %in% colnames(SC.query@meta.data)) {

  SC.query$celltype <- SC.query$predicted.id

} else {

  stop(
    "No transferred cell-type label found in SC.query"
  )
}

message("")
message("Transferred label distribution:")
print(
  table(SC.query$celltype)
)

p_query <- DimPlot(
  SC.query,
  reduction = "ref.umap",
  split.by = "sample",
  group.by = "celltype"
)

# ------------------------------------------------------------
# 18. CLEAN QUERY METADATA
# ------------------------------------------------------------

# Original HuLab script checks prediction agreement and then
# removes temporary transfer-score columns.

if (
  all(
    c(
      "predicted.id",
      "predicted.celltype"
    ) %in%
      colnames(SC.query@meta.data)
  )
) {

  SC.query$prediction.match <-
    SC.query$predicted.id ==
    SC.query$predicted.celltype

  message("")
  message("Prediction match:")
  print(
    table(
      SC.query$prediction.match
    )
  )
}

temporary_columns <- c(
  "prediction.match",
  "predicted.celltype.score",
  "prediction.score.max",
  "predicted.id",
  "prediction.score.Basal_epithelial_cells",
  "prediction.score.Luminal_L1.1_epithelial_cells",
  "prediction.score.Luminal_L1.2_epithelial_cells",
  "prediction.score.Luminal_L2_epithelial_cells"
)

temporary_columns <- intersect(
  temporary_columns,
  colnames(SC.query@meta.data)
)

if (length(temporary_columns) > 0) {

  SC.query@meta.data[
    temporary_columns
  ] <- NULL
}

if (
  "predicted.celltype" %in%
  colnames(SC.query@meta.data)
) {

  colnames(
    SC.query@meta.data
  )[
    colnames(
      SC.query@meta.data
    ) == "predicted.celltype"
  ] <- "celltype"
}

# ============================================================
# PART D — MERGE REFERENCE + BOTH QUERY SAMPLES
# ============================================================

# ------------------------------------------------------------
# 19. MERGE
# ------------------------------------------------------------

message("")
message("==============================================")
message("MERGING REFERENCE AND QUERY")
message("==============================================")

SC.merge <- merge(
  SC.integrated,
  SC.query
)

message(
  "Merged object: ",
  ncol(SC.merge),
  " cells × ",
  nrow(SC.merge),
  " genes"
)

# ------------------------------------------------------------
# 20. DOWNSTREAM SCALING / PCA / CLUSTERING / UMAP
# ------------------------------------------------------------

DefaultAssay(
  SC.merge
) <- "integrated"

VariableFeatures(
  SC.merge
) <- unique(
  c(
    VariableFeatures(SC.integrated),
    VariableFeatures(SC.query)
  )
)

SC.merge <- ScaleData(
  SC.merge,
  verbose = FALSE
)

SC.merge <- RunPCA(
  SC.merge,
  verbose = FALSE
)

SC.merge <- FindNeighbors(
  SC.merge,
  dims = 1:30
)

SC.merge <- FindClusters(
  SC.merge,
  verbose = FALSE
)

SC.merge <- RunUMAP(
  SC.merge,
  dims = 1:30
)

# ------------------------------------------------------------
# 21. STANDARDIZE NORMAL SAMPLE NAMES
# ------------------------------------------------------------

sample <- SC.merge@meta.data$sample

sample[
  which(
    sample == "Ind5"
  )
] <- "Normal1"

sample[
  which(
    sample == "Ind6"
  )
] <- "Normal2"

sample[
  which(
    sample == "Ind7"
  )
] <- "Normal3"

SC.merge@meta.data$sample <- sample

# ============================================================
# PART E — SAVE RESULTS
# ============================================================

# ------------------------------------------------------------
# 22. SAVE SEURAT OBJECTS
# ------------------------------------------------------------

message("")
message("==============================================")
message("SAVING SEURAT OBJECTS")
message("==============================================")

saveRDS(
  SC.integrated,
  file.path(
    object_dir,
    "SC_integrated_reference.rds"
  ),
  compress = FALSE
)

saveRDS(
  SC.query,
  file.path(
    object_dir,
    "SC_query_HER2_CID3586_CID4066.rds"
  ),
  compress = FALSE
)

saveRDS(
  SC.merge,
  file.path(
    object_dir,
    "SC_merge_HuLab_HER2_both.rds"
  ),
  compress = FALSE
)

saveRDS(
  SC.anchors,
  file.path(
    object_dir,
    "SC_transfer_anchors_HER2_both.rds"
  ),
  compress = FALSE
)

# ------------------------------------------------------------
# 23. SAVE METADATA AND BASIC QC
# ------------------------------------------------------------

write.csv(
  SC.merge@meta.data,
  file.path(
    result_dir,
    "SC_merge_metadata.csv"
  ),
  row.names = TRUE
)

write.csv(
  as.data.frame(
    table(
      SC.merge$sample
    )
  ),
  file.path(
    result_dir,
    "sample_cell_counts.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(
    table(
      SC.merge$sample,
      SC.merge$celltype
    )
  ),
  file.path(
    result_dir,
    "sample_celltype_counts.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(
    table(
      SC.query$sample,
      SC.query$celltype
    )
  ),
  file.path(
    result_dir,
    "HER2_query_celltype_counts.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 24. FINAL MERGED VISUALIZATION
# ------------------------------------------------------------

p_merged <- DimPlot(
  SC.merge,
  reduction = "umap",
  split.by = "sample",
  group.by = "celltype"
) +
  ggtitle(
    "HuLab reference + HER2 query alignment"
  )

ggsave(
  filename = file.path(
    result_dir,
    "SC_merge_celltype_umap.pdf"
  ),
  plot = p_merged,
  width = 16,
  height = 10
)

p_merged_condition <- DimPlot(
  SC.merge,
  reduction = "umap",
  group.by = "condition"
) +
  ggtitle(
    "HER2 condition / reference distribution"
  )

ggsave(
  filename = file.path(
    result_dir,
    "SC_merge_condition_umap.pdf"
  ),
  plot = p_merged_condition,
  width = 10,
  height = 8
)

# ============================================================
# 25. FINAL SUMMARY
# ============================================================

message("")
message("==============================================")
message("SCRIPT 02 COMPLETE")
message("==============================================")

message(
  "Reference samples: Ind5, Ind6, Ind7"
)

message(
  "Reference cells: ",
  ncol(SC.integrated)
)

message(
  "Query samples: CID3586 (Naive), CID4066 (Treated)"
)

message(
  "Query cells: ",
  ncol(SC.query)
)

message(
  "Merged cells: ",
  ncol(SC.merge)
)

message("")
message("Reference cell types:")
print(
  sort(
    table(
      SC.integrated$celltype
    ),
    decreasing = TRUE
  )
)

message("")
message("Query cell types after transfer:")
print(
  sort(
    table(
      SC.query$celltype
    ),
    decreasing = TRUE
  )
)

message("")
message("Merged sample counts:")
print(
  table(
    SC.merge$sample
  )
)

message("")
message(
  "Seurat objects saved to:"
)

message(
  object_dir
)

message("")
message(
  "Alignment results saved to:"
)

message(
  result_dir
)

message("")
message(
  "Next analysis should use SC.merge_HuLab_HER2_both.rds"
)

message("")
message(
  "SC-alignment completed successfully."
)