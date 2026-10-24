################################################################################
### Figure 3. UCHL1-mediated NF-κB activation modulates tumor-promoting 
### signaling pathways, including inflammatory and EMT programs in ccRCC.
################################################################################
# a, Tumor-intrinsic DEGs and altered hallmarks in human snRNA-seq.
# b, Upregulated NF-κB and EMT in UCHL1-high tumor nuclei (Dot plot).
# c, Enriched hallmarks in UCHL1-high tumors based on ssGSEA.
# d, Correlation of NF-κB, EMT, and IMM hallmarks in UCHL1-high tumors.
# e, Spatial expression of UCHL1 in human Xenium data (Regions of interest).
# f, Expression profiling by tumor states in Xenium (ModuleScore dot plot).
# g, Spatial protein abundance of EMT and inflammation markers using CODEX.
# h, Elevated UCHL1-associated gene signature predicts poorer overall survival 
#    in ccRCC patients (Quartile and Consensus Clustering approaches).
# i, NF-κB members correlated with UCHL1 protein abundance in cell line data.
# j, Schematic showing UCHL1-mediated NF-κB activation created with BioRender.
################################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(survival)
library(survminer)
library(gridExtra)
library(scales)
library(rhdf5)
library(readr)
library(tidyverse)
library(qs)
library(stringr)
library(tidyr)
library(patchwork)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)
library(ConsensusClusterPlus)

################################################################################
### Fig. 3a: Tumor-intrinsic DEGs and altered hallmarks (Volcano & Pathways)
################################################################################
# Using logistic regression (LR) model with latent variables
Idents(ccRCC_tumor) <- ccRCC_tumor@meta.data$UCHL1_exp
panc.markers2 <- FindMarkers(ccRCC_tumor,
                             ident.1 = "UCHL1_high",
                             ident.2 = "UCHL1_low",
                             test.use = "LR",  
                             logfc.threshold = 0.5,
                             min.pct = 0.1,
                             min.diff.pct = 0.1,
                             latent.vars = c("Sex", "Tumor_grade"))

# Filter significant DEGs
panc.markers3 <- panc.markers2[which(abs(panc.markers2$avg_log2FC) > 0.5 & panc.markers2$p_val_adj < 0.05),]
panc.markers3$gene <- rownames(panc.markers3)

# 1. Left Panel: Volcano Plot
panc.markers3$Significance <- "Not Sig"
panc.markers3$Significance[panc.markers3$avg_log2FC > 1 & panc.markers3$p_val_adj < 0.05] <- "Upregulated"
panc.markers3$Significance[panc.markers3$avg_log2FC < -1 & panc.markers3$p_val_adj < 0.05] <- "Downregulated"

genes_to_label <- c("UCHL1", "NFKB1", "RELA", "CCL2", "SERPINE1") 

volcano_plot <- ggplot(panc.markers3, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "#8e183d", "Downregulated" = "#377EB8", "Not Sig" = "grey")) +
  geom_text_repel(data = subset(panc.markers3, gene %in% genes_to_label),
                  aes(label = gene), size = 4, color = "black", box.padding = 0.5) +
  theme_classic() +
  labs(title = "DEGs: UCHL1-High vs Low", x = "Log2 Fold Change", y = "-Log10(Adj P-value)")

# 2. Right Panel: Pathway Enrichment
gene_list_up <- unique(rownames(panc.markers3[which(panc.markers3$avg_log2FC > 1),]))
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
term2gene <- hallmark[, c("gs_name", "gene_symbol")]

enrich_res <- enricher(gene_list_up, TERM2GENE = term2gene)
df_enrich <- as.data.frame(enrich_res)

pathway_plot <- ggplot(df_enrich[1:10, ], aes(x = reorder(Description, Count), y = Count, fill = p.adjust)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue") +
  coord_flip() +
  theme_classic() +
  labs(title = "Enriched Hallmarks in UCHL1-High", x = "Pathway", y = "Gene Count")

print(volcano_plot + pathway_plot)


################################################################################
### Fig. 3b: Upregulated NF-κB and EMT in UCHL1-high tumor nuclei (Dot Plot)
################################################################################
cc_sub <- readRDS("./ccRCC_UCHL1_manuscript/Data/ccRCC_snMultiome_snRNA_x72_Apr2025_merged_RNA_only_annotated_UCHL1.rds")

sn <- subset(cc_sub, idents = c("TUMOR_ccRCC_UCHL1_high_cluster", "TUMOR_ccRCC_UCHL1_low_cluster", "EPI_PT"))
sn@active.ident <- factor(sn@active.ident,
                          levels=c("EPI_PT", "TUMOR_ccRCC_UCHL1_low_cluster", "TUMOR_ccRCC_UCHL1_high_cluster"))

hallmark_panel <- c("UCHL1", "NFKB1", "RELA", "TNF", "IL6", "CCL2", 
                    "VIM", "SNAI1", "TWIST1", "SERPINE1", "CDH2", "FN1")

pdf("./ccRCC_UCHL1_manuscript/Data/Figure3b_hallmark_dotplot.pdf", width = 8, height = 6)
DotPlot(sn, features = hallmark_panel, group.by = "ident", cols = c("blue", "red")) + 
  RotatedAxis() + 
  coord_flip() +
  labs(title = "NF-kB, EMT, and IMM Hallmarks")
dev.off()


################################################################################
### Fig. 3c: Enriched Hallmarks ssGSEA (Violin Plots)
################################################################################
df_gsea <- read_tsv('./ccRCC_UCHL1_manuscript/Data/UCHL1_pseudobulk_ssGSEA_t_annotated_v10.txt') %>% as.data.frame()

plot_ssgsea <- function(df, feature, title) {
  ggviolin(df, "UCHL1_group_v2", feature, fill = "UCHL1_group_v2", 
           alpha = 1, width = 0.8, outlier.shape = NA, size = 0.1, add = "mean_sd") + 
    stat_compare_means() + theme_classic() + theme(legend.position = "none") + 
    scale_fill_manual(values = c("UCHL1_high_sample"="#8e183d", "UCHL1_low_sample"="#f6c970")) +
    labs(title = title, x = "", y = "ssGSEA Score")
}

A <- plot_ssgsea(df_gsea, "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "NF-kB Signaling")
B <- plot_ssgsea(df_gsea, "HALLMARK_INFLAMMATORY_RESPONSE", "Inflammatory Response")
C <- plot_ssgsea(df_gsea, "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "EMT")

ggarrange(A, B, C, ncol = 3, nrow = 1)


################################################################################
### Fig. 3d: Correlation of NF-κB, EMT, and IMM Hallmarks
################################################################################
df_corr <- read.delim("./ccRCC_UCHL1_manuscript/Data/UCHL1_pseudobulk_ssGSEA_t_annotated_summary_plot_v4_score.txt", header = TRUE, stringsAsFactors = FALSE)

plot_corr <- function(df, x_var, y_var) {
  ggplot(df, aes_string(x = x_var, y = y_var, color = "Group")) + 
    geom_point(size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, linetype = "dashed") + 
    stat_cor() +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
    scale_color_manual(values = c("UCHL1_high_sample"="#8e183d", "UCHL1_low_sample"="#f6c970")) +
    theme_bw(base_size = 13) +
    theme(strip.text = element_text(size = 14, face = "bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

corr_1 <- plot_corr(df_corr, "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE")
corr_2 <- plot_corr(df_corr, "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")

ggarrange(corr_1, corr_2, ncol = 2, nrow = 1, common.legend = TRUE)


################################################################################
### Fig. 3e & 3f: Xenium Spatial Expression and Profiling
################################################################################
xenium_obj <- readRDS("./ccRCC_UCHL1_manuscript/Data/Xenium_Spatial_Annotated.rds")

# 3e: Spatial DimPlot and FeaturePlots
pdf("./ccRCC_UCHL1_manuscript/Data/Figure3e_Xenium_Spatial.pdf", width = 15, height = 5)
ImageDimPlot(xenium_obj, group.by = "Cell_type", size = 0.5) 
ImageFeaturePlot(xenium_obj, features = c("UCHL1", "NFKB1", "RELA", "CCL2", "SERPINE1"), size = 0.5, ncol = 5)
dev.off()

# 3f: Dot plot for moduleScore in Xenium
pdf("./ccRCC_UCHL1_manuscript/Data/Figure3f_Xenium_DotPlot.pdf", width = 6, height = 4)
DotPlot(xenium_obj, features = hallmark_panel, group.by = "UCHL1_tumor_state") + 
  coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Xenium Explorer view imaging being incorporated into the figure panel


################################################################################
### Fig. 3g: CODEX Spatial Protein Abundance
################################################################################
# Rendered externally using QuPath / CODEX MAIA. See raw images in Supplementary Data.


################################################################################
### Fig. 3h: Elevated UCHL1-associated gene signature predicts poorer OS
################################################################################
# Load Data
raw_data <- read_tsv("/Users/Yize/Desktop/input_final.txt", col_names = TRUE)
clinical_df <- raw_data %>%
  rename_with(~ "OS_time", contains("Overall survival")) %>%
  rename_with(~ "OS_status", contains("Survival status")) %>%
  mutate(
    OS_years = OS_time / 365.25,
    Grade_Binary = case_when(
      Tumor_grade %in% c("G1", "G2") ~ "Low Grade",
      Tumor_grade %in% c("G3", "G4") ~ "High Grade",
      TRUE ~ NA_character_
    )
  )

df_export <- clinical_df
valid_surv_idx <- !is.na(clinical_df$OS_time) & !is.na(clinical_df$OS_status)

# Target 5-gene signature
target_genes <- c("UCHL1", "NFKB1", "VIM", "FN1", "ICAM1")
plot_title <- "UCHL1 | NFKB1 | VIM | FN1 | ICAM1"

expr_subset <- clinical_df %>% 
  select(all_of(target_genes)) %>% 
  mutate_all(~ as.numeric(as.character(.)))

valid_idx <- complete.cases(expr_subset) & valid_surv_idx
df_clean <- clinical_df[valid_idx, , drop = FALSE]
expr_clean <- expr_subset[valid_idx, , drop = FALSE]

# Strategy 1: Quartile-Based Analysis
z_scored_q <- scale(expr_clean)
df_export$Integrated_Signature_Score <- NA
df_export$Integrated_Signature_Score[valid_idx] <- rowMeans(z_scored_q, na.rm = TRUE)

q25 <- quantile(df_export$Integrated_Signature_Score, 0.25, na.rm = TRUE)
q75 <- quantile(df_export$Integrated_Signature_Score, 0.75, na.rm = TRUE)

df_export <- df_export %>%
  mutate(
    Quartile_Classification = case_when(
      is.na(Integrated_Signature_Score) ~ NA_character_,
      Integrated_Signature_Score >= q75 ~ "High",
      Integrated_Signature_Score <= q25 ~ "Low",
      TRUE ~ "Intermediate"
    )
  )

df_quartile <- df_export %>% filter(Quartile_Classification %in% c("High", "Low") & valid_surv_idx)
df_quartile$Score_Group <- factor(df_quartile$Quartile_Classification, levels = c("High", "Low"))

if(nrow(df_quartile) > 0) {
  fit_km_q <- survfit(Surv(OS_years, OS_status == 1) ~ Score_Group, data = df_quartile)
  km_plot_q <- ggsurvplot(fit_km_q, data = df_quartile, pval = TRUE, conf.int = TRUE,
                          palette = c("#8e183d", "#f6c970"), risk.table = FALSE)
  km_plot_q$table <- km_plot_q$table + theme(axis.line = element_blank())
  print(km_plot_q)
}

# Strategy 2: Consensus-Cluster-Based Analysis
z_scored_c <- scale(expr_clean)
if(any(is.na(z_scored_c))) {
  z_scored_c[is.na(z_scored_c)] <- 0
}

expr_matrix <- t(z_scored_c)
rownames(expr_matrix) <- target_genes
colnames(expr_matrix) <- paste0("Sample_", 1:ncol(expr_matrix))

max_k_search <- 6
ccp_res <- ConsensusClusterPlus(expr_matrix, maxK = max_k_search, reps = 500, pItem = 0.8, pFeature = 1,
                                clusterAlg = "hc", innerLinkage = "ward.D2", distance = "pearson", 
                                title = tempdir(), plot = NULL, verbose = FALSE)

pac_scores <- numeric(max_k_search)
for(k in 2:max_k_search) {
  cm <- ccp_res[[k]]$consensusMatrix
  pac_scores[k] <- sum(cm > 0.1 & cm < 0.9) / length(cm)
}

opt_k <- which.min(pac_scores[-1]) + 1 
df_export$Consensus_Cluster <- NA
df_export$Consensus_Cluster[valid_idx] <- ccp_res[[opt_k]]$consensusClass

cluster_means <- df_export %>%
  filter(!is.na(Consensus_Cluster)) %>%
  group_by(Consensus_Cluster) %>%
  summarise(Mean_UCHL1 = mean(UCHL1, na.rm = TRUE), .groups = 'drop')

high_cluster <- cluster_means$Consensus_Cluster[which.max(cluster_means$Mean_UCHL1)]
low_cluster <- cluster_means$Consensus_Cluster[which.min(cluster_means$Mean_UCHL1)]

df_export <- df_export %>%
  mutate(
    Consensus_Comparison_Group = case_when(
      is.na(Consensus_Cluster) ~ NA_character_,
      Consensus_Cluster == high_cluster ~ "High",
      Consensus_Cluster == low_cluster ~ "Low",
      TRUE ~ "Intermediate/Excluded"
    )
  )

df_cluster <- df_export %>% filter(Consensus_Comparison_Group %in% c("High", "Low") & valid_surv_idx)
df_cluster$Score_Group <- factor(df_cluster$Consensus_Comparison_Group, levels = c("High", "Low"))

if(length(unique(df_cluster$Score_Group)) == 2) {
  fit_km_c <- survfit(Surv(OS_years, OS_status == 1) ~ Score_Group, data = df_cluster)
  km_plot_c <- ggsurvplot(fit_km_c, data = df_cluster, pval = TRUE, conf.int = TRUE,
                          palette = c("#8e183d", "#f6c970"), risk.table = FALSE)
  km_plot_c$table <- km_plot_c$table + theme(axis.line = element_blank())
  print(km_plot_c)
}


################################################################################
### Fig. 3i: Cell Line Correlation (Proteomics/Phosphoproteomics)
################################################################################
cell_line_df <- read_tsv("./ccRCC_UCHL1_manuscript/Data/CellLine_Proteomics.txt") %>% as.data.frame()

p_3i <- ggplot(cell_line_df, aes(x = UCHL1_Protein, y = TAB2_Protein)) +
  geom_point(alpha = 0.7, color = "#8e183d") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor() + theme_classic() + labs(title = "Cell Line: UCHL1 vs TAB2")
print(p_3i)


################################################################################
### Fig. 3j: Schematic showing UCHL1-mediated NF-κB activation
################################################################################
# Schematic created in BioRender. Reference final manuscript for composite image.