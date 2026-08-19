# Analysis Overview

## Project Context

This repository contains a computational single-cell RNA-seq drug-repurposing
analysis originally conducted during an internship in 2023.

The analysis was subsequently reconstructed and organized to document the
original workflow, scripts, intermediate analyses, validation steps, and
final results.

## Analysis Workflow

The analysis consisted of the following major stages:

1. Dataset loading and preprocessing
2. Single-cell reference/query alignment
3. Cell-type annotation and comparison
4. Differential expression analysis
5. Drug signature scoring
6. Comparative drug analysis across HER2 conditions
7. Pathway enrichment analysis
8. Drug–pathway evidence integration
9. Comparative candidate validation
10. Final evidence-based drug ranking
11. Top-candidate pathway validation
12. Top-candidate gene-level validation
13. Figure generation

## Script Organization

| Script | Description |
|---|---|
| `01_loading_dataset.R` | Loads and prepares the datasets and metadata |
| `02_single_cell_alignment.R` | Performs reference/query alignment and cell-type annotation |
| `03_DEG_limma.R` | Performs differential expression analysis |
| `04_prepare_drug_reference.R` | Prepares the drug reference data |
| `05_drug_scoring.R` | Calculates drug therapeutic scores |
| `06_drug_comparison.R` | Compares drug signatures across conditions |
| `07_pathway_enrichment.R` | Performs pathway enrichment analysis |
| `08_pathway_drug_integration.R` | Integrates pathway and drug evidence |
| `09_comparative_validation.R` | Classifies and validates drug candidates |
| `10_final_evidence_scoring.R` | Calculates the final evidence-based ranking |
| `11_top3_pathway_validation.R` | Validates pathway evidence for the top candidates |
| `12_top3_gene_validation.R` | Performs gene-level validation of top candidates |
| `13_generate_figures.R` | Generates the final project figures |

## Conceptual Workflow

Dataset
→ Single-cell alignment
→ Cell-type annotation
→ Differential expression
→ Drug signature scoring
→ Pathway enrichment
→ Drug–pathway integration
→ Comparative validation
→ Evidence scoring
→ Top candidate identification
→ Pathway and gene validation

## Final Candidate Prioritization

The final analysis prioritized candidates using multiple evidence layers,
including therapeutic drug scores, statistical significance, cell-type
specificity, pathway evidence, and gene-level validation.

The final analysis identified:

1. Neratinib
2. Vorinostat
3. Bortezomib

as the top three computationally prioritized candidates.

These rankings represent computational prioritization and should not be
interpreted as clinical efficacy or therapeutic recommendation.
