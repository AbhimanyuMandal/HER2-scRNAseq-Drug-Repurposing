# Data Sources

This project integrates publicly available transcriptomic datasets
for normal breast tissue, HER2-positive breast cancer, and drug
perturbation signatures.

## 1. Normal breast single-cell RNA-seq

**GEO:** GSE113197  
**Relevant subseries:** GSE113196

Normal adult breast epithelial samples used:

| Sample | GEO accession |
|---|---|
| Individual 5 | GSM3099847 |
| Individual 6 | GSM3099848 |
| Individual 7 | GSM3099849 |

The samples are part of the adult human breast epithelial single-cell
RNA-seq dataset generated using the 10x Genomics platform.

## 2. HER2-positive breast cancer single-cell RNA-seq

**GEO:** GSE176078

Samples used:

| Sample | GEO accession |
|---|---|
| CID3586 | GSM5354513 |
| CID4066 | GSM5354521 |

These samples are from the single-cell and spatially resolved atlas
of human breast cancers.

## 3. LINCS L1000 perturbation data

Two GEO datasets from the Broad Institute LINCS/Connectivity Map
project were used:

- **GSE70138** — LINCS Phase II
- **GSE92742** — LINCS Phase I

These datasets provide L1000 perturbational gene-expression profiles
used for drug-response/signature analysis.

## Data availability

All primary datasets are publicly available through GEO.

The raw/large-scale datasets are not redistributed in this repository.
Users should download the datasets directly from the corresponding
repositories and follow the preprocessing steps implemented in the
analysis scripts.
