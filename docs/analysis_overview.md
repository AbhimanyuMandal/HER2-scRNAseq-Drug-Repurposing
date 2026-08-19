# Analysis Workflow

The analysis was performed as a multi-stage computational drug
repurposing workflow.

## Stage 1 — Dataset preparation

Normal breast epithelial and HER2-positive breast cancer scRNA-seq
datasets were processed and integrated.

## Stage 2 — Differential expression

Differentially expressed genes were identified using limma-based
comparisons between:

- HER2-positive breast cancer and normal breast
- CID4066 and normal
- CID3586 and normal
- CID4066 and CID3586

## Stage 3 — Drug signature scoring

Disease-associated transcriptional signatures were compared with
LINCS perturbational signatures to identify candidate drugs.

## Stage 4 — Pathway analysis

Candidate drug-associated transcriptional changes were evaluated
using pathway enrichment analysis.

## Stage 5 — Comparative validation

Candidates were classified according to their relationship with
the disease and treatment/transition contrasts.

## Stage 6 — Evidence-based ranking

Candidate drugs were ranked using integrated evidence from:

- therapeutic/drug signature score
- statistical significance
- cell-type specificity
- pathway-level evidence
- non-redundant biological pathway themes

## Stage 7 — Top candidate validation

The top three candidates were subjected to additional pathway and
gene-level validation.

## Stage 8 — Visualization

Final results were summarized using four publication-style figures.
