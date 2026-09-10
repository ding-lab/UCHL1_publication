# Functional and spatial multimodal dissection reveals a UCHL1-driven malignant program in clear cell renal cell carcinoma

## Overview

This repository contains the custom code, analytical pipelines, and data processing scripts used in **Li et al. (under review)** to investigate the UCHL1-driven malignant program in ccRCC. We integrate single-nucleus RNA/ATAC-seq, high-resolution spatial transcriptomics (Xenium, Xenium 5K), CODEX imaging, and deep ubiquitinomic profiling to characterize ccRCC heterogeneity. The scripts provided here reproduce the bioinformatic analyses, including consensus clustering, spatial neighborhood delineation, and multi-omic pathway integrations.

## Repo structure

UCHL1_publication/
├── README.md
├── LICENSE
├── envs/
│   ├── UCHL1_python_env.yml               # Conda environment for Python
│   └── UCHL1_R_env.yml                    # Conda environment for R
├── Data/
│   └── README.md                          # Instructions for data access
└── Scripts/
    ├── Preprocessing/                     # Shared processing steps (e.g., Harmony, BANKSY)
    ├── Fig1_study_design_clinical_cohort/ # Figure 1
    ├── Fig2_functional_assays/            # Figure 2
    ├── Fig3_NFkB_signaling/               # Figure 3
    ├── Fig4_metabolic_rewiring/           # Figure 4
    ├── Fig5_spatial_TME/                  # Figure 5
    └── Fig6_therapeutic_models/           # Figure 6

## System Requirements

### Hardware requirements
The preprocessing and spatial integration scripts (e.g., Harmony integration, BANKSY spatial clustering) require a high-performance computing (HPC) cluster or a standard workstation with at least:
* 64+ GB RAM
* 8+ CPU cores
Standard downstream statistical analyses and plotting can be performed on a standard laptop with 16 GB RAM.

### Software requirements
* **OS**: Linux (Ubuntu 20.04/22.04 or CentOS 7) or macOS.
* **Languages**: R (v4.3.0+) and Python (v3.11.6+).

## Installation Guide
We recommend using [Anaconda](https://www.anaconda.com/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) to manage the environments.

1. Clone this repository:
```bash
git clone https://github.com/ding-lab/UCHL1_publication.git
cd UCHL1_publication
```

2. Create and activate the Python environment (install time: ~5 minutes):
```bash
conda env create -f envs/UCHL1_python_env.yml
conda activate UCHL1_python
```

3. Create and activate the R environment (install time: ~10-15 minutes):
```bash
conda env create -f envs/UCHL1_R_env.yml
conda activate UCHL1_R
```

4. Packages not available via Conda, such as CellChat v1.6.1 and scider v1.7.3, must be installed via devtools::install_github() in R. The complete list can be found below:
```bash
R Packages:
AnnotationDbi v1.68.0
argparse v2.2.5
AUCell v1.28.0
BiocParallel v1.40.0
biomaRt v2.62.1
broom v1.0.8
BSgenome.Hsapiens.UCSC.hg38 v1.4.5
car v3.1-3
caret v6.0-94
CellChat v1.6.1
ChIPseeker v1.42.1
chromVAR v1.12.0
circlize v0.4.16
clinfun v1.1.5
clusterProfiler v4.14.6
ComplexHeatmap v2.22.0
ConsensusClusterPlus v1.74.0
cowplot v1.1.3
cutpointr v1.2.1
data.table v1.17.4
DelayedMatrixStats v1.28.1
DESeq2 v1.46.0
dplyr v1.1.4
edgeR v4.4.2
emmeans v1.11.2-8
EnhancedVolcano v1.24.0
ensembldb v2.30.0
EnsDb.Hsapiens.v86 v2.99.0
EnsDb.Hsapiens.v100 v0.0.1
factoextra v1.0.7
FactoMineR v2.12
forcats v1.0.0
future v1.49.0
genefilter v1.88.0
GenomeInfoDb v1.42.0
GenomicRanges v1.58.0
ggalluvial v0.12.5
ggplot2 v3.5.2
ggpubr v0.6.0
ggrastr v1.0.2
ggrepel v0.9.6
gplots v3.2.0
ggVennDiagram v1.5.7
grid v4.4.3
gridExtra v2.3
GSEABase v1.68.0
GSVA v2.0.7
harmony v1.2.3
inferCNV v0.99.7
JASPAR2020 v0.99.10
limma v3.56.2
maftools v2.14.0
magrittr v2.0.3
MAST v1.32.0
Matrix v1.7-3
matrixStats v1.5.0
mclust v6.1.1
monocle3 v1.3.1
motifmatchr v1.28.0
msigdbr v25.1.1
optparse v1.7.5
org.Hs.eg.db v3.20.0
patchwork v1.3.0
pheatmap v1.0.13
purrr v1.0.4
qs v0.27.3
RColorBrewer v1.1-3
readr v2.1.5
reticulate v1.42.0
rstatix v0.7.2
scales v1.4.0
scCustomize v3.0.1
scider v1.7.3
Seurat v5.3.1
Signac v1.3.0
SpatialExperiment v1.12.0
SummarizedExperiment v1.36.0
survival v3.7.0
survminer v0.5.0
TCGAbiolinks v2.34.1
TFBSTools v1.44.0
tibble v3.2.1
tidyr v1.3.1
tidyverse v2.0.0
viridis v0.6.5

Python Packages:
Python v3.11.6
anndata v0.11.4
banksy_py v1.3.3
harmony-pytorch v0.0.10
matplotlib v3.10.3
numpy v2.3.0
pandas v2.2.3
scanpy v1.9.5
scikit-image v0.25.2
scipy v1.16.2
Scrublet v0.2.3
scikit-learn v0.24.2
spatialdata-io v0.7.1
statsmodels v0.15.0
```

## Bioinformatics Tools
The scripts and additional information for reproducing the analyses and figures are available in the GitHub public repository (https://github.com/ding-lab/UCHL1_publication). We used bioinformatics tools in this work, as detailed in the Methods section. All tools written and/or published by the authors are freely available at our public GitHub repository (https://github.com/ding-lab/). Other tools used here are available from their associated authors.

## Data Access
All data are available in the paper and its Supplementary Information. Source data are provided with this paper. WES, RNA-seq, snRNA-seq, snATAC-seq, Xenium, Xenium 5K, and CODEX imaging data are part of the HTAN dbGaP study accession phs002371.v1.p1 and can be accessed through the HTAN Data Coordinating Center (DCC) Portal (https://data.humantumoratlas.org/) under the HTAN WUSTL Atlas. CPTAC: CPTAC snRNA-seq data can be accessed through NIH GDC Data Portal using the CPTAC Case ID (https://portal.gdc.cancer.gov/cases/). GRCh38 references used for scRNA-seq and snRNA-seq (refdata-gex-GRCh38-2020-A) are freely available from the 10x Genomics website (https://support.10xgenomics.com). The reference GRCh38 genome (GRCh38.d1.vd1.fa.tar.gz) used for WES read alignment is available from GDC (https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files). The MSigDB hallmark gene sets are available from the GSEA website (https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp). Due to GitHub file size limits, all processed data objects (Seurat .rds, Scanpy .h5ad, and DIA-MS matrices) required to run these scripts will be hosted on Zenodo.

## Instructions for Use (Reproducibility)
The `Scripts` directory is organized by figure to facilitate direct reproducibility of the manuscript's panels.

* **`Scripts/Preprocessing/`**: Scripts for initial QC, Scrublet doublet removal, and Harmony batch correction for snRNA-seq and snATAC-seq, and preprocessing of spatial assays.
* **`Scripts/Fig1_study_design_clinical_cohort/`**: Gaussian Mixture Modeling (GMM) for UCHL1-high classification and survival analyses.
* **`Scripts/Fig2_functional_assays/`**: Processing scripts for *in vitro* functional assays, including cell line CCK-8 proliferation curves, CRISPR-KO variant confirmation, and multi-omic (proteomic, phosphoproteomic, and ubiquitomic) consensus clustering and differential expression.
* **`Scripts/Fig3_NFkB_signaling/`**: Differential expression, ssGSEA scoring, and integration of ubiquitination data for IKKβ/NF-κB axis.
* **`Scripts/Fig4_metabolic_rewiring/`**: Pseudotime trajectory inference (Monocle3) and HIF1A ubiquitin correlation analysis.
* **`Scripts/Fig5_spatial_TME/`**: BANKSY spatial niche generation, grid-density contouring, and CellChat intercellular communication across Xenium and Xenium 5K datasets.
* **`Scripts/Fig6_therapeutic_models/`**: Analysis pipelines for evaluating *in vivo* therapeutic efficacy across cell-line-derived xenograft (CDX) and patient-derived xenograft (PDX) models, including combination drug treatments (e.g., IMP-1710 with sunitinib or cabozantinib).

Each folder contains a standalone script (e.g., `Fig5.R`) that specifies the detailed analyses.

---

**Manuscript**: Functional and spatial multimodal dissection reveals a UCHL1-driven malignant program in clear cell renal cell carcinoma

**Status**: Under Review

**Contact**: Yize Li at [yize.li@wustl.edu]