<p align="center">
  <img src="assets/drug_repurposing_banner.png" alt="HER2 Luminal Breast Cancer Drug Repurposing Banner" width="100%">
</p>

<h1 align="center">
HER2-Positive Breast Cancer scRNA-seq Drug Repurposing
</h1>

<p align="center">

Single-cell transcriptomic analysis and LINCS-based computational drug repurposing in HER2-positive breast cancer

</p>

<p align="center">
  
<img src="https://img.shields.io/badge/R-4.x-blue?logo=r" />
<img src="https://img.shields.io/badge/Seurat-scRNA--seq-blueviolet" />
<img src="https://img.shields.io/badge/LINCS-Connectivity%20Analysis-orange" />
<img src="https://img.shields.io/badge/Data-NCBI%20GEO-blue" />
<img src="https://img.shields.io/badge/Drug%20Repurposing-Computational-red" />
<img src="https://img.shields.io/badge/License-MIT-brightgreen" />
</p>
---

## Overview

This project investigated **computational drug repurposing opportunities in HER2-positive breast cancer** using single-cell RNA sequencing (scRNA-seq) and large-scale perturbational gene-expression signatures from the **LINCS L1000** resource.

The analysis integrated normal breast epithelial cells with HER2-positive breast cancer samples and compared transcriptional states across disease and treatment-related contrasts. Disease-associated transcriptional signatures were then compared with LINCS drug perturbation profiles to identify candidate compounds whose transcriptional effects could potentially reverse or modulate disease-associated states.

The workflow subsequently integrated **drug-level statistical evidence, cell-type specificity, pathway enrichment, gene-level evidence, and comparative validation** to prioritize candidate drugs.

### Key computational candidates

| Rank | Drug | Candidate class | Final evidence score |
|---:|---|---|---:|
| 1 | **Neratinib** | Transition-associated | **0.9212** |
| 2 | **Vorinostat** | Treatment-contrast-specific | **0.7058** |
| 3 | **Bortezomib** | Treatment-contrast-specific | **0.4425** |

These candidates represent **computational drug-repurposing hypotheses** and require experimental and clinical validation.

---

## Project Context

This analysis was originally conducted during a **2023 research internship** involving computational analysis of breast cancer single-cell transcriptomic data.

The original analysis was not fully documented at the time. The present repository consolidates the analysis scripts, intermediate/final result tables, validation outputs, and figures into a structured and reproducible project.

> **Note:** The repository documents and organizes the original research workflow; it does not imply that the underlying biological experiments were performed as part of this repository.

---

# Research Question

Can single-cell transcriptional states in HER2-positive breast cancer be used together with LINCS perturbational signatures to identify and prioritize candidate drugs for computational repurposing?

The analysis was designed to address this question through several levels of evidence:

1. Identify disease-associated transcriptional changes.
2. Compare disease signatures against drug perturbation signatures.
3. Identify candidate compounds with potentially relevant transcriptional effects.
4. Compare candidates across HER2-positive breast cancer states.
5. Evaluate pathway-level biological evidence.
6. Validate the strongest candidates using gene- and pathway-level evidence.
7. Produce an integrated evidence-based drug ranking.

---

# Datasets

The project integrates publicly available transcriptomic datasets.

## Normal breast tissue

**NCBI GEO:** GSE113197  
**Relevant subseries:** GSE113196

Normal adult breast epithelial samples used:

| Sample | GEO accession |
|---|---|
| Individual 5 | GSM3099847 |
| Individual 6 | GSM3099848 |
| Individual 7 | GSM3099849 |

## HER2-positive breast cancer

**NCBI GEO:** GSE176078

Samples used:

| Sample | GEO accession |
|---|---|
| CID3586 | GSM5354513 |
| CID4066 | GSM5354521 |

The analysis uses the CID3586 and CID4066 samples to investigate transcriptional differences between HER2-positive breast cancer states.

## LINCS perturbational signatures

Two LINCS L1000 datasets were used:

- **GSE70138**
- **GSE92742**

These datasets provide large-scale perturbational gene-expression signatures used for computational drug-response and connectivity analysis.

### Data availability

The original datasets are publicly available through NCBI GEO and the LINCS resource.

The raw datasets are **not redistributed in this repository** because of their size. Users should download the source datasets directly from their respective repositories.

See [`docs/data_sources.md`](docs/data_sources.md) for the complete dataset description.

---

# Analysis Workflow

```text
Public scRNA-seq datasets
        │
        ▼
Dataset loading & preprocessing
        │
        ▼
Single-cell alignment / integration
        │
        ▼
Cell-type-specific transcriptional analysis
        │
        ▼
Differential expression analysis
        │
        ├───────────────┐
        ▼               ▼
Disease signatures   HER2 state comparison
        │               │
        └───────┬───────┘
                ▼
       LINCS drug signatures
                │
                ▼
        Drug signature scoring
                │
                ▼
      Candidate drug landscape
                │
                ▼
     Comparative validation
                │
                ▼
       Pathway enrichment
                │
                ▼
      Evidence integration
                │
                ▼
      Top candidate validation
                │
                ▼
      Final drug prioritization
