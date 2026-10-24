################################################################################
### Figure 4. UCHL1 metabolically rewires a hypoxia-glycolysis tumor state 
### that emerges during human ccRCC progression.
################################################################################
# a, Distinct ssGSEA hallmarks by UCHL1 group based on human snRNA-seq (Radar plot).
# b, Differentially expressed metabolites relevant to glycolysis and FAO.
# c, Gene sets for hypoxia, glycolysis, and FAO among U-high, U-low, and PT (Dot plot).
# d, Significant pathways based on DEPs of C2/C9 compared to NC.
# e, UCHL1 stabilizes HIFs based on the cell line ubiquitinome and proteome (Volcano).
# f, Trajectory inference indicates a UCHL1-high state emerging during progression.
# g, UCHL1 expression in liver metastases vs paired primary ccRCCs (Boxplots).
# h, U-high and U-low tumor content distribution in paired primary and metastatic ccRCCs.
# i, Spatial dynamics of UCHL1 and metabolic landscape in Xenium 5K data (SP476).
# j, Schematic demonstrating how UCHL1 metabolically rewires a hypoxia-glycolysis state.
################################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(fmsb)
library(readr)
library(tidyverse)
library(ggrepel)
library(clusterProfiler)
library(monocle3)
library(patchwork)

################################################################################
### Fig. 4a: Radar Plot of Distinct ssGSEA Hallmarks
################################################################################
avg_scores <- read_tsv('./ccRCC_UCHL1_manuscript/Data/UCHL1_pseudobulk_ssGSEA_t_annotated_summary_plot_v4.txt') %>% as.data.frame()

selected_pathways <- c(
  "HALLMARK_GLYCOLYSIS", "HALLMARK_FATTY_ACID_METABOLISM", 
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_HYPOXIA",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_INFLAMMATORY_RESPONSE"
)

pathway_cols <- colnames(avg_scores)[(colnames(avg_scores) %in% selected_pathways)]

group_means <- avg_scores %>%
  group_by(Group) %>%
  summarise(across(all_of(pathway_cols), mean))

max_vals <- apply(group_means[,-1], 2, max)
min_vals <- apply(group_means[,-1], 2, min)

radar_data <- rbind(max = max_vals, min = min_vals, group_means[,-1])
rownames(radar_data) <- c("max", "min", group_means$Group)

line_colors <- c("#f6c970", "#8e183d")
fill_colors <- adjustcolor(line_colors, alpha.f = 0.3)

pdf("./ccRCC_UCHL1_manuscript/Data/Figure4a_RadarPlot.pdf", width = 6, height = 6)
radarchart(radar_data, pcol = line_colors, pfcol = fill_colors, plwd = 2,
           cglcol = "grey", cglty = 1, axislabcol = "grey", vlcex = 0.8)
legend("topright", legend = group_means$Group, col = line_colors, lty = 1, lwd = 2, bty = "n")
dev.off()


################################################################################
### Fig. 4b: Differentially Expressed Metabolites
################################################################################
df_metab <- read.delim(file = "./ccRCC_UCHL1_manuscript/Data/Metabolite_UCHL1_v3_sub4.txt", header = TRUE, stringsAsFactors = FALSE)

plot_metab <- function(df, feature, title) {
  ggplot(df, aes_string(x = "Group", y = feature, fill = "Group")) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
    stat_compare_means() + theme_classic() + theme(legend.position = "none") + 
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
    scale_fill_manual(values = c("UCHL1_high"="#8e183d", "UCHL1_low"="#f6c970", "NAT"="#82bfd9")) +
    labs(title = title, x = "", y = "Abundance")
}

p1 <- plot_metab(df_metab, "Succinate", "Succinate (Hypoxia)")
p2 <- plot_metab(df_metab, "Lactic.acid", "Lactic Acid (Glycolysis)")
p3 <- plot_metab(df_metab, "Glucose", "Glucose (Glycolysis)")
p4 <- plot_metab(df_metab, "Acetyl.L.carnitine", "Acetyl-L-carnitine (FAO)")

pdf("./ccRCC_UCHL1_manuscript/Data/Figure4b_Metabolites.pdf", width = 10, height = 8)
ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
dev.off()


################################################################################
### Fig. 4c: Dot Plot for Hypoxia, Glycolysis, and FAO (snRNA-seq)
################################################################################
cc_ith <- subset(cc_sub, idents = c("C3N-01200-03", "C3L-00813-01", "C3L-00010-01"))
Idents(cc_ith) <- factor(cc_ith@meta.data$Cell_type_detailed_UCHL1, 
                         levels = c("EPI_PT", "TUMOR_ccRCC_UCHL1_low_cluster", "TUMOR_ccRCC_UCHL1_high_cluster"))

metabolic_genes <- c(
  "HIF1A", "EPAS1", "VEGFA", "CA9",       # Hypoxia
  "SLC2A1", "HK2", "PKM", "LDHA",         # Glycolysis
  "CPT1A", "ACADM", "PPARA", "FABP4"      # FAO
)

pdf("./ccRCC_UCHL1_manuscript/Data/Figure4c_Metabolism_DotPlot.pdf", width = 8, height = 6)
DotPlot(cc_ith, features = metabolic_genes, group.by = "ident", cols = c("blue", "red")) + 
  RotatedAxis() + coord_flip()
dev.off()


################################################################################
### Fig. 4d: Significant Pathways (Cell Line Proteomics UCHL1-KO)
################################################################################
enrich_res_DEP <- enricher(DEP_gene_list, TERM2GENE = term2gene)

pdf("./ccRCC_UCHL1_manuscript/Data/Figure4d_Pathway_DotPlot.pdf", width = 7, height = 6)
dotplot(enrich_res_DEP, showCategory = 10, font.size = 10) + ggtitle("DEPs: UCHL1-KO vs NC")
dev.off()


################################################################################
### Fig. 4e: Ubiquitinome Volcano Plots (HIF1A & EPAS1)
################################################################################
c2_deu <- read_tsv('./ccRCC_UCHL1_manuscript/Data/DEU_Caki.1_C2_vs_NC_log2.tsv') %>% as.data.frame()
c9_deu <- read_tsv('./ccRCC_UCHL1_manuscript/Data/DEU_Caki.1_C9_vs_NC_log2.tsv') %>% as.data.frame()

plot_volcano_ub <- function(df, title) {
  df$logp <- -log10(df$FDR_Limma_trend_test)
  df$group <- ifelse(df$logFC > 0.5 & df$FDR_Limma_trend_test < 0.2, "Up",
                     ifelse(df$logFC < -0.5 & df$FDR_Limma_trend_test < 0.2, "Down", "NotSig"))
  
  ggplot(df, aes(x = logFC, y = logp, color = group)) +
    geom_point(alpha = 0.7, size = 1.8) +
    scale_color_manual(values = c("Up" = "#f6c970", "Down" = "#8e183d", "NotSig" = "gray")) +
    geom_text_repel(data = subset(df, Protein %in% c("HIF1A_ub1", "HIF1A_ub2", "HIF1A_ub3", "EPAS1_ub1")), 
                    aes(label = Protein), color = "black", size = 3.5, box.padding = 0.8) +
    theme_minimal() +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10(FDR)") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.2), linetype = "dashed", color = "black")
}

pdf("./ccRCC_UCHL1_manuscript/Data/Figure4e_Ubiquitinome.pdf", width = 10, height = 5)
ggarrange(plot_volcano_ub(c2_deu, "C2 vs NC"), plot_volcano_ub(c9_deu, "C9 vs NC"), ncol = 2)
dev.off()


################################################################################
### Fig. 4f: Trajectory Inference (Monocle3)
################################################################################
cds_sub <- readRDS("./ccRCC_UCHL1_manuscript/Data/cds_learngraph1_light.rds")

# Left Panel: Cell Type Distribution UMAP
p_traj_cell <- plot_cells(cds_sub, color_cells_by = "Cell_type_detailed_UCHL1_v2", 
                          label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE) + 
               scale_color_manual(values = celltype_cell_type_YL_detailed_UCHL1)

# Right Panel: Pseudotime UMAP
p_traj_pseudo <- plot_cells(cds_sub, color_cells_by = "pseudotime", 
                            label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

pdf("./ccRCC_UCHL1_manuscript/Data/Figure4f_Trajectory.pdf", width = 12, height = 5)
p_traj_cell | p_traj_pseudo
dev.off()


################################################################################
### Fig. 4g & 4h: Primary vs. Liver Metastases
################################################################################
df_mets <- read_tsv('./ccRCC_UCHL1_manuscript/Data/Primary_Metastasis_Paired.txt') %>% as.data.frame()
df_mets_counts <- read_tsv('./ccRCC_UCHL1_manuscript/Data/Primary_Metastasis_Counts.txt') %>% as.data.frame()

# 4g Boxplots (Paired Wilcoxon)
p_met_uchl1 <- ggpaired(df_mets, x = "Site", y = "UCHL1_Expression", fill = "Site", 
                        line.color = "gray", line.size = 0.4) +
               stat_compare_means(paired = TRUE) + theme_classic()

# 4h Chi-square Stacked Bar Plot
p_met_dist <- ggplot(df_mets_counts, aes(x = Patient, y = Count, fill = UCHL1_group)) +
              geom_bar(position = "fill", stat = "identity") +
              scale_fill_manual(values = c("UCHL1_high"="#8e183d", "UCHL1_low"="#f6c970")) +
              labs(y = "Proportion of Tumor Cells", x = "Paired Samples") + theme_minimal()

pdf("./ccRCC_UCHL1_manuscript/Data/Figure4gh_Metastasis.pdf", width = 10, height = 5)
p_met_uchl1 | p_met_dist
dev.off()


################################################################################
### Fig. 4i: Xenium 5K Spatial Dynamics (SP476)
################################################################################
xenium_pri <- readRDS("./ccRCC_UCHL1_manuscript/Data/SP476_Primary_Xenium.rds")
xenium_met <- readRDS("./ccRCC_UCHL1_manuscript/Data/SP476_LiverMet_Xenium.rds")

pdf("./ccRCC_UCHL1_manuscript/Data/Figure4i_Xenium_SP476.pdf", width = 15, height = 8)
p_xenium_pri <- ImageFeaturePlot(xenium_pri, features = c("UCHL1", "HK2", "HIF1A", "NFKB1"), size = 0.5)
p_xenium_met <- ImageFeaturePlot(xenium_met, features = c("UCHL1", "HK2", "HIF1A", "NFKB1"), size = 0.5)
print(p_xenium_pri / p_xenium_met)
dev.off()

# Xenium Explorer view imaging being incorporated into the figure panel


################################################################################
### Fig. 4j: Schematic (BioRender)
################################################################################
# Schematic created in BioRender demonstrating how UCHL1 metabolically rewires 
# a hypoxia-glycolysis tumor state. Reference final manuscript for composite image.