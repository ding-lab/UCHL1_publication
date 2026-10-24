################################################################################
### Figure 5: UCHL1-high tumor cells are spatially embedded within an inflamed 
### yet immunosuppressive TME.
################################################################################
# a, TME spatial delineation using Xenium (UMAP of integrated 19 ccRCC samples).
# b, Tumor classification by UCHL1 expression in Xenium (Gaussian mixture modeling).
# c, TME intertumoral heterogeneity in U-high vs. U-low sample groups (Boxplots).
# d, Co-enrichment of U-high state and immune infiltration in Xenium (Pie charts).
# e, Distinct spatial patterns of U-high vs. U-low regions in Xenium (Tumor density).
# f, Spatially-resolved cell type distribution across niches in Xenium (BANKSY).
# g, TME intratumoral heterogeneity in U-high vs. U-low regions (Transcript mappings).
# h, Spatial protein abundance of TME markers in U-high regions (CODEX).
# i, Significant cell-cell interactions between U-high and TME (Xenium 5K CellChat).
# j, Mechanistic crosstalk in Xenium 5K (Ligand-receptor pairs).
# k, Schematic demonstrating spatial embedding of UCHL1-high tumor cells.
################################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(survival)
library(gridExtra)
library(scales)
library(rhdf5)
library(readr)
library(tidyverse)
library(qs)
library(stringr)
library(tidyr)
library(patchwork)
library(mclust)
library(NMF)
library(CellChat)

################################################################################
### Fig. 5a: TME spatial delineation using Xenium (UMAP)
################################################################################
# Refer to run_all_477_ccRCC_x19.sh for MakeSlimSeurat -> Sketch+Harmony and the 
# scripts create_slim_object_layer_v6.R and SlimSeurat__run_sketch_harmony.R, 
# which generated the output file of 4-5_integrated_harmony.qs.

obj <- qread("./ccRCC_UCHL1_manuscript/Data/4-5_integrated_harmony.qs") # 4,301,606 cells
Idents(obj) <- obj@meta.data$seurat_clusters_full

ALL_MARKERS_TO_PLOT = c("UCHL1", "MKI67", "CA9", # TUMOR
                        "LRP2",
                        "CD19", "MS4A1", "CD79A", "MZB1",  # B_CELL_PLASMA
                        "CLEC10A", "CD1C", # DC
                        "CD68", "CD163", "MS4A4A", # TAM_MACROPHAGE
                        "CD3D","CD3E","TRAC","IL7R", # CD4+T
                        "CD8A",# EXHAUSTED_CD8_T
                        "GZMK","GZMB","GNLY","FGFBP2","KLRD1", # NK_CELLS
                        "KIT", "MS4A2", # MAST_CELLS
                        "THBS2", "PDGFRA", "C7", # CAF_MYOFIBROBLAST
                        "PDGFRB", "THY1", # VSM/P
                        "VWF","CD34") # ENDOTHELIAL

# Extended Data validation plot
DotPlot(obj, features = ALL_MARKERS_TO_PLOT, cols = c("blue", "red")) + RotatedAxis()

# Clusters 12, 19, 23 excluded due to doublet pattern or low quality
xenium.sub <- subset(obj, idents = c("12", "19", "23"), invert = TRUE) # 4,065,399 cells

xenium.sub@meta.data$celltype_braod <- plyr::mapvalues(as.character(Idents(xenium.sub)),
                                                       from=c("1","2","3","4","5","6","7","8","9","10","11","13","14","15","16","17","18","20","21","22","24","25","26","27","28"),
                                                       to=c("TUMOR_low","VEC","TUMOR_low","VSM/P","TUMOR_low","IMM_CD4_T","IMM_TAM","IMM_CD8_T","TUMOR_low","myCAF","IMM_TAM","TUMOR_low","TUMOR_high","IMM_B","IMM_TAM","IMM_NK","VEC","TUMOR_high","VEC","IMM_TAM","IMM_Mast","IMM_DC","TUMOR_low","IMM_CD8_T","PT"))

Idents(xenium.sub) <- xenium.sub@meta.data$celltype_braod
xenium.sub@active.ident <- factor(xenium.sub@active.ident, 
                                  levels=c("VEC", "VSM/P", "myCAF", "IMM_Mast", 
                                           "IMM_NK", "IMM_CD8_T", "IMM_CD4_T", 
                                           "IMM_TAM", "IMM_DC", "IMM_B", "PT", 
                                           "TUMOR_low", "TUMOR_high"))

celltype_color <- c("TUMOR_high"="#8e183d","TUMOR_low"="#f6c970","PT"="#1CFFCE",
                    "IMM_B"="#FB9A99","IMM_DC"="#CFA3CB","IMM_TAM"="#FF7F00",
                    "IMM_CD4_T"="#c8facc","IMM_CD8_T"="#E7298A","IMM_NK"="#E12027",
                    "IMM_Mast"="#999999","myCAF"="#66A61E","VSM/P"="#2834b5","VEC"="#377EB8")

# Generate UMAPs
DimPlot(object = xenium.sub, reduction = "full.umap.sketch", label=FALSE, label.size=6, cols = celltype_color)
DimPlot(object = xenium.sub, reduction = "umap.rerun.full", label=FALSE, label.size=6, cols = celltype_color)


################################################################################
### Fig. 5b: Tumor classification by UCHL1 expression (GMM)
################################################################################
tumor <- subset(xenium.sub, idents = c("TUMOR_low", "TUMOR_high"))
meta <- tumor@meta.data

df <- meta %>%
  group_by(orig.ident, celltype_concised) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(fraction = n / sum(n))

ratio_df <- df %>%
  group_by(orig.ident, celltype_concised) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(fraction = n / sum(n)) %>%
  ungroup() %>%
  select(orig.ident, celltype_concised, fraction) %>%
  pivot_wider(names_from = celltype_concised, values_from = fraction) %>%
  mutate(ratio = TUMOR_high / (TUMOR_high + TUMOR_low))

ratio_clean <- ratio_df %>% filter(!is.na(ratio) & !is.nan(ratio) & is.finite(ratio))

# Gaussian Mixture Modeling (GMM) classification
mod_all <- Mclust(ratio_clean$ratio)
mod <- Mclust(ratio_clean$ratio, G = 1:3)
means <- mod$parameters$mean
cutoff <- mean(means) 

# Left panel showing distribution and data-driven cutoff
ggplot(ratio_clean, aes(x = ratio)) +
  geom_density(fill = "#8E063B", alpha = 1, color = "black") +
  geom_vline(aes(xintercept = cutoff), linetype = "dashed", color = "blue", linewidth = 1) +
  labs(x = "Fraction of UCHL1 High Cluster", y = "Density", title = paste0("Mixture Model Cutoff at ", round(cutoff, 2))) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

ratio_clean <- ratio_clean %>%
  mutate(UCHL1_group = ifelse(ratio >= cutoff, "UCHL1_high_sample", "UCHL1_low_sample"))

# Sample mapping based on calculated classifications
xenium.sub@meta.data$category <- plyr::mapvalues(as.character(xenium.sub@meta.data$orig.ident),
                                                 from=c("WUR238-A3Us2_1","C3L-01287-11Us1_1","C3N-01200-11Us1_1","C3L-00970-13U1","HT282N1-S1H2A1Fp1U1","HT293N1-S1H2A1FpU1",
                                                        "HT338N2-S1H4Fp1U1","HT609N1-S1H2Fp1U1","HT918N1-S1H1Fp1U1","HT282N1-S1H3A4U2","HT338N2-H1FP1U1","HT522N1-Th1K4Fp1U1",
                                                        "HT609N1-S1H1Fp1U1","C3L-00982-11U2","C3L-01676-13U2","C3L-01929-11U2","C3L-02551-11U3","HT941N1-S1H1Fp1U1",
                                                        "HT973N1-S1H1A2U1"),
                                                 to=c("Low","High","High","High","Low","Low",
                                                      "High","Low","Low","Low","High","Low",
                                                      "Low","High","High","High","High","High",
                                                      "Low"))


################################################################################
### Fig. 5c: TME intertumoral heterogeneity in U-high vs. U-low (Boxplots)
################################################################################
meta_data <- xenium.sub@meta.data
Idents(xenium.sub) <- xenium.sub@meta.data$celltype_concised 
meta_data$CellType <- Idents(xenium.sub)

meta_data_filtered <- meta_data %>%
  filter(!CellType %in% c("TUMOR_high", "TUMOR_low", "PT"))

prop_table <- meta_data_filtered %>%
  group_by(orig.ident, category, CellType) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(TotalCells = sum(Count), Proportion = Count / TotalCells) %>%
  ungroup()

# Filter for the specifically requested cell types in legend 5c
prop_table_sub <- prop_table[which(prop_table$CellType %in% c("IMM_B", "IMM_TAM", "IMM_CD4_T", "IMM_CD8_T", "VSM/P", "VEC")),]

ggplot(prop_table_sub, aes(x = category, y = Proportion, fill = category)) +
  geom_boxplot(outlier.shape = NA, alpha = 1) +
  geom_jitter(width = 0.2, size = 1.5, aes(color = category)) +
  labs(x = "Group", y = "Proportion") +
  scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
  theme_minimal() +
  stat_compare_means(size = 3) + # Wilcoxon rank-sum test
  theme(legend.position = "none", strip.text = element_text(face = "bold")) +
  scale_fill_manual(values = c("Low" = "#f6c970", "High" = "#8e183d")) +
  scale_color_manual(values = c("High" = "#8e183d", "Low" = "#f6c970")) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2) +
  facet_wrap(~ CellType, scales = "free_y")


################################################################################
### Fig. 5d: Co-enrichment of U-high state and immune infiltration (Pie Charts)
################################################################################
# Group cell types into IMM (Immune) and STR (Stromal)
Idents(xenium.sub) <- xenium.sub@meta.data$celltype_concised
xenium.sub@meta.data$celltype_high <- plyr::mapvalues(as.character(Idents(xenium.sub)),
                                                  from=c("IMM_TAM","IMM_B","myCAF","TUMOR_high","TUMOR_low","IMM_CD8_T","IMM_Mast","IMM_CD4_T","VEC","VSM/P","IMM_NK","IMM_DC","PT"),
                                                  to=c("IMM","IMM","STR","TUMOR_high","TUMOR_low","IMM","IMM","IMM","STR","STR","IMM","IMM","PT"))

qsave(xenium.sub, "./ccRCC_UCHL1_manuscript/Data/477_ccRCC_x19_integrated_harmony_annotated.qs")
# ImageDimPlots and Pie charts generated iteratively per sample grouping.


################################################################################
### Fig. 5e: Distinct spatial patterns of U-high vs. U-low regions
################################################################################
# Tumor density spatial mapping (upper panel) and periphery vs core categorization 
# (bottom panel) calculated via regional distance functions.


################################################################################
### Fig. 5f: Spatially-resolved cell type distribution across niches (BANKSY)
################################################################################
file_list <- list.files(pattern = "./ccRCC_UCHL1_manuscript/Data/*_leiden_obs.csv", full.names = TRUE)
df_list <- lapply(file_list, read.csv)
combined_df <- do.call(rbind, df_list)

annotation <- combined_df
annotation$barcode2 <- sub("^([A-Za-z]+-[0-9]+)_(.*)$", "\\2_\\1", annotation$barcode)

idx <- match(rownames(xenium.sub@meta.data), annotation$barcode2)
xenium.sub@meta.data$Banksy_leiden <- annotation$group[idx]
xenium.sub@meta.data$Banksy_secuer <- annotation$secuer_cluster[idx]

valid_cells <- rownames(xenium.sub@meta.data)[!is.na(xenium.sub@meta.data$Banksy_secuer)]
xenium.sub_banksy <- subset(xenium.sub, cells = valid_cells) 

meta_data <- xenium.sub_banksy@meta.data
meta_data$CellType <- xenium.sub_banksy@meta.data$celltype_concised
meta_data$Group <- as.factor(meta_data$Banksy_secuer)

tumor_prop_per_neighborhood <- meta_data %>%
  dplyr::group_by(Banksy_secuer) %>%
  dplyr::summarize(tumor_fraction = mean(CellType %in% c("TUMOR_high", "TUMOR_low"), na.rm = TRUE))

meta_data$tumor_neighborhood_fraction <- tumor_prop_per_neighborhood$tumor_fraction[match(meta_data$Banksy_secuer, tumor_prop_per_neighborhood$Banksy_secuer)]
meta_data$TUMOR_neighbor <- ifelse(meta_data$tumor_neighborhood_fraction > 0.2, "Yes", "No")

# Identify N7 vs Other Tumor Ns
meta_data_filtered <- meta_data %>%
  filter(!CellType %in% c("PT", "TUMOR_high", "TUMOR_low"), TUMOR_neighbor == "Yes") %>%
  mutate(Plot_Group = ifelse(Banksy_secuer == 7, "N7", "Other")) %>%
  mutate(Plot_Group = factor(Plot_Group, levels = c("N7", "Other")))

# Paired Wilcoxon Signed-Rank Test plotting code executed for stacked bar charts.


################################################################################
### Fig. 5g: TME intratumoral heterogeneity in U-high vs. U-low regions
################################################################################
# Xenium spatial transcript colors as saved from Xenium Explorer
gene_color_fig.5g <- c("CD163"="#d1c829","TGFB1"="#f16566","IL7R"="#8dc640",
                       "LAG3"="#d3569f","PDGFRB"="#f89935","PECAM1"="#5982c2",
                       "EDNRB"="#69c395","EDN1"="#47bce5","UCHL1"="#8e183d",
                       "NFKB1"="#acacac")

# Run correlation check for spatial colocalization
neighborhood_avg <- meta_data %>%
  dplyr::group_by(Banksy_secuer) %>%
  dplyr::summarize(
    UCHL1 = mean(UCHL1, na.rm = TRUE),
    NFKB1 = mean(NFKB1, na.rm = TRUE),
    CD68 = mean(CD68, na.rm = TRUE)
  )

# Xenium Explorer view imaging being incorporated into the figure panel


################################################################################
### Fig. 5h: Spatial protein abundance of TME markers (CODEX)
################################################################################
# CODEX spatial images processed via QuPath-0.7.0-x64. Color code:
gene_color_fig.5h <- c("DAPI"="#3953a4","CA9"="#6fccdd","UCHL1"="#f7ec13",
                       "CD68"="#f57f20","CD4"="#e02628","CD8"="#e52c8b",
                       "CD20"="#f6999a","CD31"="#fcfcfc","SMA"="#6abd45")


################################################################################
### Fig. 5i & 5j: Cell-cell interactions and Mechanistic crosstalk (CellChat)
################################################################################
# Subsampling dataset for CellChat execution
total_cells <- colnames(xenium.sub)
set.seed(123)
half_cells <- sample(total_cells, size = floor(length(total_cells) / 2)) 
seurat_half <- subset(xenium.sub, cells = half_cells)
DefaultAssay(seurat_half) <- "Xenium"
seurat_half <- JoinLayers(seurat_half)

cellchat <- createCellChat(object = seurat_half, group.by = "celltype_braod")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 15)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Network Heatmap (Fig. 5i Left Panel)
pdf("./ccRCC_UCHL1_manuscript/Data/signaling_role_heatmap_outgoing.pdf", width = 20, height = 20)
cellchat_healthy <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_heatmap(cellchat_healthy, pattern = "outgoing")
dev.off()

# Mechanistic crosstalk specific LR pairs (Fig. 5j)
pathways.show <- c("FN1", "SPP1", "PECAM1", "LAMININ", "CD45", "EDN", "KIT",
                   "CDH", "CCL", "CD86", "SELL", "CXCL", "CLEC", "EPHB")

pdf("./ccRCC_UCHL1_manuscript/Data/PAR_top_LR_pair.pdf", width = 10, height = 10)
for (variable in pathways.show) {
  pairLR <- extractEnrichedLR(cellchat, signaling = variable, geneLR.return = FALSE)
  netVisual_individual(
    cellchat,
    signaling    = variable,
    pairLR.use   = pairLR[1, ], 
    layout       = "circle"
  )
}
dev.off()


################################################################################
### Fig. 5k: Schematic
################################################################################
# Schematic created in BioRender demonstrating spatial embedding of UCHL1-high 
# tumor cells within an inflamed yet immunosuppressive TME.