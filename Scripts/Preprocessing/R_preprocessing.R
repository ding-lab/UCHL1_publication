################################################################################
### Master R Pipeline: ccRCC Multi-Omics, Spatial & Proteogenomic Analysis
### Associated Manuscript: Functional and spatial multimodal dissection reveals 
### a UCHL1-driven malignant program in clear cell renal cell carcinoma
################################################################################

################################################################################
### 1. LOAD LIBRARIES & CONFIGURE ENVIRONMENT
################################################################################
# Core & Data Manipulation
library(tidyverse)         # v2.0.0 (dplyr, tidyr, readr, purrr, tibble, forcats, etc.)
library(data.table)        # v1.17.4
library(Matrix)            # v1.7-3

# Single-Cell & Spatial Core
library(Seurat)            # v5.3.1
library(Signac)            # v1.3.0
library(SpatialExperiment) # v1.12.0
library(harmony)           # v1.2.3

# High-Performance I/O & Parallelization
library(qs)                # v0.27.3 (Saving/loading of intermediate results)
library(future)            # v1.49.0
library(BiocParallel)      # v1.40.0
plan("multisession", workers = 30) 

# Specialized Preprocessing & Orthogonal Validation
library(inferCNV)          # v0.99.7
library(chromVAR)          # v1.12.0
library(motifmatchr)       # v1.28.0
library(JASPAR2020)        # v0.99.10
library(BSgenome.Hsapiens.UCSC.hg38) # v1.4.5

# Downstream Tools
library(monocle3)          # v1.3.1
library(CellChat)          # v1.6.1
library(scider)            # v1.7.3
library(clusterProfiler)   # v4.14.6
library(GSVA)              # v2.0.7
library(limma)             # v3.56.2
library(survival)          # v3.7.0
library(mclust)            # v6.1.1
library(ConsensusClusterPlus) # v1.74.0

# Visualization
library(ComplexHeatmap)    # v2.22.0
library(ggplot2)           # v3.5.2
library(ggpubr)            # v0.6.0
library(patchwork)         # v1.3.0
library(ggalluvial)        # v0.12.5
library(ggrastr)           # v1.0.2


################################################################################
### 2. AHEAD-OF-TIME R PREPROCESSING (Orthogonal Annotations)
################################################################################

# A. inferCNV (Validate malignant vs. normal compartments)
# infercnv_obj <- CreateInfercnvObject(raw_counts_matrix, annotations_file, ref_group_names=c("PT_cells", "Endothelial"))
# infercnv_obj <- run(infercnv_obj, cutoff=0.1, out_dir="intermediate_data/inferCNV_out")

# B. snATAC-seq Motif Matrix Generation (Map JASPAR motifs to accessibility peaks)
# motif_matrix <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
# atac_obj <- AddMotifs(atac_obj, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = motif_matrix)
# atac_obj <- RunChromVAR(atac_obj, genome = BSgenome.Hsapiens.UCSC.hg38)
# qsave(atac_obj, "intermediate_data/atac_chromvar_prep.qs") 


################################################################################
### 3. CORE INTEGRATION (Seurat & Harmony)
################################################################################

# Load filtered objects 
# seurat_merged <- qread("intermediate_data/filtered_rna_merged.qs")

# Normalize & Batch Correct
# seurat_merged <- NormalizeData(seurat_merged, scale.factor = 10000) %>%
#   FindVariableFeatures(nfeatures = 3000) %>%
#   ScaleData(vars.to.regress = c("percent.mito", "nCount_RNA", "S.Score", "G2M.Score")) %>%
#   RunPCA(npcs = 50) %>%
#   RunHarmony(group.by.vars = "sample_id", dims.use = 1:50) %>%
#   FindNeighbors(reduction = "harmony", dims = 1:50) %>%
#   FindClusters(resolution = 1.0)

# Save post-integration intermediate
# qsave(seurat_merged, "intermediate_data/seurat_harmony_integrated.qs")


################################################################################
### 4. DOWNSTREAM ANALYSIS BLOCKS
################################################################################

# 4A. Trajectory Inference (Monocle3)
# cds <- as.cell_data_set(seurat_merged) 

# 4B. Tumor Classification via Gaussian Mixture Modeling (mclust)
# gmm_model <- Mclust(uchl1_fractions, G = 1:3) 

# 4C. Cell Communication (CellChat)
# cellchat <- createCellChat(object = seurat_merged) 

# 4D. Pathway Enrichment (clusterProfiler / GSVA)
# enrich_res <- enricher(gene = sig_genes, TERM2GENE = msig_db) 

# 4E. Survival & Clinical Association Analysis (survival / survminer)
# fit <- survfit(Surv(OS_time, OS_status) ~ group, data = surv_data)


################################################################################
### 5. VISUALIZATION EXPORT (ggplot2, ComplexHeatmap, ggpubr)
################################################################################

# Execute dedicated visualization scripts per figure panel.
# Note: Use ggrastr::rasterize() for heavily populated spatial / UMAP plots,
# ggalluvial for pathway flow mapping, and ComplexHeatmap for bulk proteomics 
# consensus clustering matrices.