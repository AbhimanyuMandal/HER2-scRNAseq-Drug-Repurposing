# ============================================================
# HU LAB REPRODUCTION PROJECT
# SCRIPT 03 — DEG ANALYSIS (LIMMA)
#
# Faithful reconstruction of the original HuLab DEG workflow,
# extended to include both HER2+/ER+ samples:
#
#   CID3586 = Naive
#   CID4066 = Treated
#
# Primary method in original HuLab workflow:
#   LIMMA
#
# Comparisons:
#   1. CID4066 vs Normal1 + Normal2 + Normal3
#      -> original HuLab comparison
#
#   2. CID3586 vs Normal1 + Normal2 + Normal3
#      -> new extension
#
#   3. CID4066 vs CID3586
#      -> new treated-vs-naive comparison
#
# DEG is performed separately for each transferred epithelial
# cell type.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(limma)
})

# ------------------------------------------------------------
# 0. PROJECT PATHS
# ------------------------------------------------------------

project_dir <- getwd()

seurat_file <- file.path(
  project_dir,
  "data",
  "processed",
  "Seurat",
  "SC_merge_HuLab_HER2_both.rds"
)

result_dir <- file.path(
  project_dir,
  "results",
  "03_DEG"
)

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. LOAD MERGED OBJECT
# ------------------------------------------------------------

message("==============================================")
message("SCRIPT 03 — DEG ANALYSIS (LIMMA)")
message("==============================================")
message("Project directory: ", project_dir)
message("")

if (!file.exists(seurat_file)) {
  stop(
    "SC.merge object not found:\n",
    seurat_file
  )
}

SC.merge <- readRDS(seurat_file)

message(
  "Loaded SC.merge: ",
  ncol(SC.merge),
  " cells × ",
  nrow(SC.merge),
  " genes"
)

# ------------------------------------------------------------
# 2. BASIC METADATA VALIDATION
# ------------------------------------------------------------

required_metadata <- c(
  "sample",
  "celltype"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(SC.merge@meta.data)
)

if (length(missing_metadata) > 0) {
  stop(
    "Required metadata missing: ",
    paste(
      missing_metadata,
      collapse = ", "
    )
  )
}

message("")
message("Sample counts:")
print(
  table(
    SC.merge$sample
  )
)

message("")
message("Cell-type counts:")
print(
  table(
    SC.merge$celltype
  )
)

# ------------------------------------------------------------
# 3. USE RNA ASSAY
# ------------------------------------------------------------

if (!"RNA" %in% names(SC.merge@assays)) {
  stop("RNA assay not found in SC.merge.")
}

DefaultAssay(SC.merge) <- "RNA"

SC.merge <- JoinLayers(
  SC.merge,
  assay = "RNA"
)

# ------------------------------------------------------------
# 4. HELPER FUNCTION TO ACCESS RNA EXPRESSION
#
# Original HuLab script used:
#   c_cells@assays$RNA@data
#
# Seurat v5 stores assay data in layers. Use LayerData()
# where available while retaining the original normalized-data
# logic.
# ------------------------------------------------------------

get_rna_data <- function(seurat_object) {

  assay <- seurat_object[["RNA"]]

  if ("data" %in% Layers(assay)) {

    expr <- LayerData(
      seurat_object,
      assay = "RNA",
      layer = "data"
    )

  } else {

    stop(
      "RNA assay does not contain a 'data' layer. ",
      "Available layers: ",
      paste(
        Layers(assay),
        collapse = ", "
       )
    )
  }

  return(
    as.matrix(expr)
  )
}


# ------------------------------------------------------------
# 5. LIMMA FUNCTION
#
# This reproduces the original HuLab logic:
#
#   - subset one cell type
#   - select Case and Control cells
#   - use normalized RNA data
#   - remove genes expressed in <3 cells
#   - model Case vs Control
#   - empirical Bayes
#   - BY multiple-testing correction
#   - p.value threshold = 1e-7
#   - score = LIMMA t-statistic
#
# The only deliberate extension is that case/control sample
# names are supplied as function arguments.
# ------------------------------------------------------------

run_limma_deg <- function(
  seurat_object,
  case_samples,
  control_samples,
  comparison_name,
  output_prefix
) {

  message("")
  message("==============================================")
  message(
    "COMPARISON: ",
    comparison_name
  )
  message("==============================================")

  message(
    "Case: ",
    paste(
      case_samples,
      collapse = ", "
    )
  )

  message(
    "Control: ",
    paste(
      control_samples,
      collapse = ", "
    )
  )

  set.seed(123456)

  min.cells <- 3

  gene_list <- list()

  summary_rows <- list()

  celltypes <- unique(
    seurat_object@meta.data$celltype
  )

  for (celltype_name in celltypes) {

    message("")
    message(
      "Cell type: ",
      celltype_name
    )

    # --------------------------------------------------------
    # Subset cell type
    # --------------------------------------------------------

    cell_ids <- rownames(
      seurat_object@meta.data[
        seurat_object@meta.data$celltype ==
          celltype_name,
        ,
        drop = FALSE
      ]
    )

    c_cells <- subset(
      seurat_object,
      cells = cell_ids
    )

    samples <- c_cells@meta.data

    control_cells <- rownames(
      subset(
        samples,
        sample %in% control_samples
      )
    )

    case_cells <- rownames(
      subset(
        samples,
        sample %in% case_samples
      )
    )

    message(
      "Case cells: ",
      length(case_cells)
    )

    message(
      "Control cells: ",
      length(control_cells)
    )

    # Original HuLab threshold:
    # both groups must contain > min.cells cells.
    if (
      length(control_cells) <= min.cells ||
      length(case_cells) <= min.cells
    ) {

      message(
        "SKIPPED — insufficient cells."
      )

      summary_rows[[celltype_name]] <- data.frame(
        comparison = comparison_name,
        celltype = celltype_name,
        case_cells = length(case_cells),
        control_cells = length(control_cells),
        genes_tested = 0,
        genes_significant = 0,
        status = "SKIPPED",
        stringsAsFactors = FALSE
      )

      next
    }

    # --------------------------------------------------------
    # Extract normalized RNA expression
    # --------------------------------------------------------

    expr <- get_rna_data(
      c_cells
    )

    # Make sure expression columns match Seurat cells.
    common_case <- intersect(
      case_cells,
      colnames(expr)
    )

    common_control <- intersect(
      control_cells,
      colnames(expr)
    )

    if (
      length(common_case) <= min.cells ||
      length(common_control) <= min.cells
    ) {

      message(
        "SKIPPED — insufficient cells after expression matching."
      )

      summary_rows[[celltype_name]] <- data.frame(
        comparison = comparison_name,
        celltype = celltype_name,
        case_cells = length(common_case),
        control_cells = length(common_control),
        genes_tested = 0,
        genes_significant = 0,
        status = "SKIPPED_AFTER_EXPRESSION_MATCH",
        stringsAsFactors = FALSE
      )

      next
    }

    selected_cells <- c(
      common_case,
      common_control
    )

    expr <- expr[
      ,
      selected_cells,
      drop = FALSE
    ]

    # --------------------------------------------------------
    # Remove genes expressed in <3 cells
    # --------------------------------------------------------

    bad <- which(
      rowSums(
        expr > 0
      ) < 3
    )

    if (length(bad) > 0) {

      expr <- expr[
        -bad,
        ,
        drop = FALSE
      ]
    }

    genes_tested <- nrow(expr)

    # --------------------------------------------------------
    # Build sample information
    # --------------------------------------------------------

    new_sample <- data.frame(
      Samples = selected_cells,
      type = c(
        rep(
          "Case",
          length(common_case)
        ),
        rep(
          "Control",
          length(common_control)
        )
      ),
      stringsAsFactors = FALSE
    )

    rownames(new_sample) <- paste(
      new_sample$Samples,
      seq_len(
        nrow(new_sample)
      ),
      sep = "_"
    )

    # --------------------------------------------------------
    # LIMMA
    # --------------------------------------------------------

    mm <- model.matrix(
      ~ 0 + type,
      data = new_sample
    )

    colnames(mm) <- c(
      "Case",
      "Control"
    )

    fit <- lmFit(
      expr,
      mm
    )

    contr <- makeContrasts(
      Case - Control,
      levels = colnames(
        coef(fit)
      )
    )

    tmp <- contrasts.fit(
      fit,
      contrasts = contr
    )

    tmp <- eBayes(
      tmp
    )

    # Original HuLab threshold:
    # p.value = 1e-7
    C_data <- topTable(
      tmp,
      adjust.method = "BY",
      sort.by = "P",
      n = nrow(tmp),
      p.value = 1e-7
    )

    # --------------------------------------------------------
    # Drug-repurposing format
    # --------------------------------------------------------

    C_data_for_drug <- data.frame(
      row.names = rownames(C_data),
      score = C_data$t,
      adj.P.Val = C_data$adj.P.Val,
      P.Value = C_data$P.Value
    )

    gene_list[[celltype_name]] <-
      C_data_for_drug

    summary_rows[[celltype_name]] <- data.frame(
      comparison = comparison_name,
      celltype = celltype_name,
      case_cells = length(common_case),
      control_cells = length(common_control),
      genes_tested = genes_tested,
      genes_significant = nrow(C_data_for_drug),
      status = "COMPLETED",
      stringsAsFactors = FALSE
    )

    message(
      "Genes tested: ",
      genes_tested
    )

    message(
      "Genes passing p.value <= 1e-7: ",
      nrow(C_data_for_drug)
    )
  }

  names(gene_list) <- names(
    gene_list
  )

  summary_df <- do.call(
    rbind,
    summary_rows
  )

  # ----------------------------------------------------------
  # SAVE
  # ----------------------------------------------------------

  gene_file <- file.path(
    result_dir,
    paste0(
      output_prefix,
      "_limma.rds"
    )
  )

  summary_file <- file.path(
    result_dir,
    paste0(
      output_prefix,
      "_summary.csv"
    )
  )

  saveRDS(
    gene_list,
    file = gene_file
  )

  write.csv(
    summary_df,
    file = summary_file,
    row.names = FALSE
  )

  message("")
  message(
    "Saved: ",
    gene_file
  )

  message(
    "Saved: ",
    summary_file
  )

  return(
    gene_list
  )
}

# ============================================================
# 6. ANALYSIS 1 — ORIGINAL HULAB COMPARISON
#
# CID4066 vs Normal1 + Normal2 + Normal3
# ============================================================

Gene.list.CID4066 <- run_limma_deg(
  seurat_object = SC.merge,
  case_samples = "CID4066",
  control_samples = c(
    "Normal1",
    "Normal2",
    "Normal3"
  ),
  comparison_name =
    "CID4066_vs_Normal",
  output_prefix =
    "HER2_ER_CID4066_vs_Normal"
)

# ============================================================
# 7. ANALYSIS 2 — CID3586 NAIVE VS NORMAL
# ============================================================

Gene.list.CID3586 <- run_limma_deg(
  seurat_object = SC.merge,
  case_samples = "CID3586",
  control_samples = c(
    "Normal1",
    "Normal2",
    "Normal3"
  ),
  comparison_name =
    "CID3586_vs_Normal",
  output_prefix =
    "HER2_ER_CID3586_vs_Normal"
)

# ============================================================
# 8. ANALYSIS 3 — TREATED VS NAIVE
#
# CID4066 vs CID3586
# ============================================================

Gene.list.Treated_vs_Naive <- run_limma_deg(
  seurat_object = SC.merge,
  case_samples = "CID4066",
  control_samples = "CID3586",
  comparison_name =
    "CID4066_vs_CID3586",
  output_prefix =
    "HER2_ER_CID4066_vs_CID3586"
)

# ============================================================
# 9. FINAL SUMMARY
# ============================================================

message("")
message("==============================================")
message("SCRIPT 03 COMPLETE")
message("==============================================")

message("")
message("Three LIMMA comparisons completed:")
message(
  "1. CID4066 vs Normal1 + Normal2 + Normal3"
)
message(
  "2. CID3586 vs Normal1 + Normal2 + Normal3"
)
message(
  "3. CID4066 vs CID3586"
)

message("")
message("Results directory:")
message(
  result_dir
)

message("")
message(
  "The CID4066-vs-Normal result reproduces the ",
  "original HuLab DEG comparison."
)

message(
  "CID3586-vs-Normal and CID4066-vs-CID3586 ",
  "are the two extensions for the dual-sample analysis."
)

message("")
message(
  "Next stage: DrugScore / Asgard analysis."
)