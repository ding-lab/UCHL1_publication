################################################################################
### Figure 1. Study design for functional and spatial multimodal dissection, 
### revealing a UCHL1-driven malignant state in ccRCC that is significantly 
### associated with worse prognosis.
################################################################################

# a, Schematic, created with BioRender, summarizing the study design, from 
#    multi-omics and PTM-based delineation to spatial dynamics and TME interaction 
#    of UCHL1, and from in vitro to in vivo to evaluate the therapeutic strategy, 
#    across datasets of human samples, cell lines, PDX, and CDX models. 
# b, Cell type annotation and tumor-intrinsic UCHL1 classification in the human 
#    snRNA-seq atlas. Left panel: the uniform manifold approximation and projection 
#    (UMAP) displays cell clusters that were subsequently annotated by cell types, 
#    consisting of 476,527 nuclei from 71 samples; right panel: the UMAP shows the 
#    tumor subclusters colored by refined tumor-intrinsic UCHL1 classification. 
# c, snRNA-seq Gaussian mixture modeling stratifies tumors into distinct UCHL1 
#    subgroups using the mixture model cutoff of 0.31, as noted in the left panel. 
#    13 samples were classified as UCHL1-high tumors (U-high), and 54 were 
#    classified as UCHL1-low tumors (U-low). Among them, 7 samples were labeled 
#    as UCHL1-ITH samples (i.e., #1-7), with considerable U-high and U-low content 
#    in the pie charts. 
# d, Refined UCHL1 tumor stratification in snRNA-seq and bulk multi-omic concordance. 
#    Boxplots show the significantly higher UCHL1 quantifications in the U-high group. 
# e, Elevated UCHL1 levels predict poorer overall survival in ccRCC patients. The 
#    Kaplan-Meier plots indicate the association between overall survival and UCHL1 
#    gene expression in this study cohort, the CPTAC ccRCC cohort, and the TCGA KIRC 
#    cohort, respectively, and between overall survival and UCHL1 protein abundance 
#    in the CPTAC ccRCC cohort. Log-rank test p = 0.0022, 0.016, 0.049, and 0.00027, 
#    respectively. 

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(survival)
library(survminer)
library(gridExtra)
library(scales)
library(tidyverse)
library(qs)
library(mclust)
library(patchwork)

################################################################################
### Figure 1a: Schematic
################################################################################
# Schematic created in BioRender. Reference final manuscript for composite image.


################################################################################
### Figure 1b: Cell type annotation and tumor-intrinsic UCHL1 classification
################################################################################
# Load annotated snRNA-seq atlas
xenium.sub <- qread("./ccRCC_UCHL1_manuscript/Data/ccRCC_NAT_snMultiome_snRNA_x71_annotated.qs")
Idents(xenium.sub) <- xenium.sub@meta.data$Cell_type_detailed

# Define cell type color palette
sn_celltype_color <- c("TUMOR_ccRCC"="#A5772B", "EPI_PC"="#1CFFCE", "EPI_PEC"="#FF7F00",
                       "EPI_POD"="#15616d", "EPI_PT"="#2834b5", "EPI_TAL"="#E6AB02",
                       "IMM_B"="#FB9A99", "IMM_DC"="#CFA3CB", "IMM_MAC"="#E12027",
                       "IMM_MAST"="#999999", "IMM_PL"="#377EB8", "IMM_T"="#E7298A",
                       "FIB"="#66A61E", "EC"="#F7E1A0")

xenium.sub@active.ident <- factor(xenium.sub@active.ident, 
                                  levels=c("TUMOR_ccRCC", "EPI_PC", "EPI_PEC",
                                           "EPI_POD", "EPI_PT", "EPI_TAL",
                                           "IMM_B", "IMM_DC", "IMM_MAC",
                                           "IMM_MAST", "IMM_PL", "IMM_T",
                                           "FIB", "EC"))

# Fig 1b Left Panel: UMAP showing all cell types
DimPlot(object = xenium.sub, reduction = "sct.umap", label = FALSE, label.size = 6, cols = sn_celltype_color) + 
  labs(title = "Cell Type Annotation (snRNA-seq Atlas)")

# Subset and process tumor populations for UCHL1 classification
ccRCC_tumor <- subset(xenium.sub, idents = c("TUMOR_ccRCC"))
ccRCC_tumor <- NormalizeData(ccRCC_tumor)
ccRCC_tumor <- FindVariableFeatures(ccRCC_tumor)
ccRCC_tumor <- ScaleData(ccRCC_tumor)
ccRCC_tumor <- RunPCA(ccRCC_tumor)
ccRCC_tumor <- FindNeighbors(ccRCC_tumor, dims = 1:30, reduction = "pca")
ccRCC_tumor <- FindClusters(ccRCC_tumor, resolution = 0.5, cluster.name = "recluster_PC")
ccRCC_tumor <- RunUMAP(ccRCC_tumor, dims = 1:30, reduction = "pca", reduction.name = "umap.integrated.recluster_PC")

# Map UCHL1 clusters
Idents(ccRCC_tumor) <- ccRCC_tumor@meta.data$UCHL1_cluster
uchl1_colors <- c("UCHL1_low_cluster" = "#F6C971", "UCHL1_high_cluster" = "#8E063B")

# Fig 1b Right Panel: UMAP showing tumor subclusters
DimPlot(object = ccRCC_tumor, reduction = "umap.integrated.recluster_PC", label = FALSE, label.size = 6, cols = uchl1_colors) + 
  labs(title = "Tumor-intrinsic UCHL1 Classification")


################################################################################
### Figure 1c: snRNA-seq Gaussian mixture modeling & Pie Charts
################################################################################
# Calculate cluster fractions for mixture modeling
meta <- ccRCC_tumor@meta.data
df <- meta %>%
  group_by(Sample_ID, UCHL1_cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Sample_ID) %>%
  mutate(fraction = n / sum(n)) %>%
  ungroup()

ratio_df <- df %>%
  select(Sample_ID, UCHL1_cluster, fraction) %>%
  pivot_wider(names_from = UCHL1_cluster, values_from = fraction, values_fill = 0) %>%
  mutate(ratio = UCHL1_high_cluster / (UCHL1_high_cluster + UCHL1_low_cluster))

ratio_clean <- ratio_df %>% filter(!is.na(ratio) & !is.nan(ratio) & is.finite(ratio))

# Gaussian Mixture Model (GMM)
mod <- Mclust(ratio_clean$ratio, G = 1:3)
cutoff <- mean(mod$parameters$mean)

# Fig 1c Left Panel: GMM Density Plot
ggplot(ratio_clean, aes(x = ratio)) +
  geom_density(fill = "#8E063B", alpha = 1, color = "black") +
  geom_vline(aes(xintercept = cutoff), linetype = "dashed", color = "blue", linewidth = 1) +
  labs(x = "Fraction of UCHL1 High Cluster", y = "Density", 
       title = paste0("Mixture Model Cutoff at ", round(cutoff, 2))) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

# Fig 1c Right Panel: Pie Charts for ITH and Subgroups
sample_stats <- meta %>%
  group_by(Sample_ID) %>%
  summarise(high_fraction = mean(UCHL1_cluster == "UCHL1_high_cluster")) %>%
  arrange(desc(high_fraction))

plot_list <- list()
for (i in seq_len(nrow(sample_stats))) {
  s <- sample_stats$Sample_ID[i]
  frac <- sample_stats$high_fraction[i]
  
  df_pie <- meta %>%
    filter(Sample_ID == s) %>%
    count(UCHL1_cluster) %>%
    mutate(Percent = n / sum(n) * 100)
  
  p <- ggplot(df_pie, aes(x = "", y = Percent, fill = UCHL1_cluster)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = uchl1_colors, drop = FALSE) +
    theme_void() +
    theme(legend.position = "none", plot.title = element_text(size = 9, hjust = 0.5)) +
    labs(title = paste0(s, "\nHigh=", round(frac * 100, 1), "%"))
  
  plot_list[[s]] <- p
}

pdf("./ccRCC_UCHL1_manuscript/Data/Figure_1c_piecharts.pdf", width = 30, height = 21)
ggarrange(plotlist = plot_list, ncol = 10, nrow = 7)
dev.off()


################################################################################
### Figure 1d: Refined UCHL1 tumor stratification boxplots
################################################################################
# Load integrated pseudobulk and multi-omic data
expr_matrix <- read_tsv('./ccRCC_UCHL1_manuscript/Data/UCHL1_pseudobulk_ssGSEA_t_annotated_v10.txt') %>% as.data.frame()
my_comparisons <- list(c("UCHL1_high_sample", "UCHL1_low_sample")) 

p_1d <- ggplot(expr_matrix, aes(x = UCHL1_group, y = UCHL1_avg_exp_scaled, fill = UCHL1_group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_fill_manual(values = c("UCHL1_low_sample" = "#f6c970", "UCHL1_high_sample" = "#8e183d")) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format") + 
  labs(title = "snRNA-seq Multi-omic Concordance", x = "UCHL1 Group", y = "UCHL1 Avg Exp (Scaled)") +
  theme_minimal() + theme(legend.position = "none")

print(p_1d)


################################################################################
### Figure 1e: Elevated UCHL1 levels predict poorer overall survival (Kaplan-Meier)
################################################################################
# Load survival datasets
surv_data_study <- read_tsv('./ccRCC_UCHL1_manuscript/Data/Study_Survival_Data.txt') %>% as.data.frame()
surv_data_cptac <- read_tsv('./ccRCC_UCHL1_manuscript/Data/CPTAC_Survival_Data.txt') %>% as.data.frame()
surv_data_tcga  <- read_tsv('./ccRCC_UCHL1_manuscript/Data/TCGA_Survival_Data.txt') %>% as.data.frame()

# Standardized KM plot function
plot_km <- function(data, time_col, status_col, group_col, title) {
  form <- as.formula(paste("Surv(", time_col, ",", status_col, ") ~", group_col))
  fit <- surv_fit(form, data = data)
  
  ggsurvplot(fit, data = data, pval = TRUE, pval.method = TRUE,
             palette = c("#8E063B", "#F6C971"), 
             title = title, xlab = "Time (Days)", ylab = "Overall Survival",
             legend.title = "UCHL1 Level",
             legend.labs = c("High", "Low"),
             risk.table = TRUE, theme = theme_classic())
}

# 1. Study Cohort (Gene Expression)
km_study <- plot_km(surv_data_study, "OS_time", "OS_status", "UCHL1_group", "Study Cohort (snRNA)")

# 2. CPTAC ccRCC Cohort (Gene Expression)
km_cptac_rna <- plot_km(surv_data_cptac, "OS_time", "OS_status", "UCHL1_rna_group", "CPTAC Cohort (RNA)")

# 3. TCGA KIRC Cohort (Gene Expression)
km_tcga <- plot_km(surv_data_tcga, "OS_time", "OS_status", "UCHL1_group", "TCGA KIRC Cohort (RNA)")

# 4. CPTAC ccRCC Cohort (Protein Abundance)
km_cptac_pro <- plot_km(surv_data_cptac, "OS_time", "OS_status", "UCHL1_protein_group", "CPTAC Cohort (Protein)")

# Arrange and export KM plots
pdf("./ccRCC_UCHL1_manuscript/Data/Figure_1e_Survival.pdf", width = 20, height = 5)
arrange_ggsurvplots(list(km_study, km_cptac_rna, km_tcga, km_cptac_pro), 
                    print = TRUE, ncol = 4, nrow = 1)
dev.off()