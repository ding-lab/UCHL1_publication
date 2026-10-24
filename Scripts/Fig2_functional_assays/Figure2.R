################################################################################
### Figure 2. UCHL1 CRISPR knockout, downregulation, and inhibition markedly 
### suppress ccRCC cell proliferation, as evidenced by in vitro assays.
################################################################################
# a, UCHL1 expression among RCC cell lines (log2 TMM-adjusted CPM).
# b, UCHL1 CRISPR-mediated knockout (KO) efficiency in Caki-1 cells (Lollipop plot).
# c, Western blot showing GAPDH and UCHL1 protein levels across Caki-1 clones.
# d, UCHL1-KO reduces proliferation (CCK-8 absorbance, two-way ANOVA).
# e, UCHL1 shRNA knockdown (KD) efficiency in 786-O (Proliferation and qRT-PCR).
# f, Identification of IMP-1710 as a highly potent UCHL1 inhibitor (CCK-8 viability).
# g, Cell viability of Caki-1 and 786-O treated cells with IMP-1710 over 96 hours.
# h, Schematic of cell-line-derived proteomics and PTM characterization.
# i, Cell line sample correlations at the protein, phosphoprotein, and ubiquitination levels.
# j, Protein abundances of UCHL1 and VHL across ccRCC cell lines.
# k, High correlations of C2 and C9 at the log2-scale protein level.
# l, Venn diagrams showing differentially expressed proteins (DEPs) by UCHL1 status.
################################################################################

library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(pheatmap)
library(eulerr)
library(tidyr)
library(limma)

################################################################################
### Data Loading & Preprocessing (Proteomics)
################################################################################
expr_matrix <- read.delim("./ccRCC_UCHL1_manuscript/Data/UCHL1_global_protein_Jul2025.txt", stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
cancer_s <- read.delim("./ccRCC_UCHL1_manuscript/Data/Protein_meta_v3_plot_IMP_updated.txt", stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)

pro <- as.matrix(expr_matrix[, c(1:14, 16:23)])
pro <- log2(pro + 1)
cancer_s <- cancer_s[c(1:10, 12:23), ]
cancer_s$sample <- rownames(cancer_s)


################################################################################
### Fig. 2a: UCHL1 expression among RCC cell lines (RNA-seq)
################################################################################
rna_cpm <- read.csv("./ccRCC_UCHL1_manuscript/Data/RNA_CPM_matrix.csv")
df_2a <- data.frame(CellLine = colnames(rna_cpm), UCHL1_log2CPM = log2(as.numeric(rna_cpm["UCHL1", ]) + 1))
ggplot(df_2a, aes(x = reorder(CellLine, -UCHL1_log2CPM), y = UCHL1_log2CPM)) +
  geom_bar(stat="identity", fill="steelblue") + theme_classic() +
  labs(x = "RCC Cell Lines", y = "UCHL1 Expression (log2 TMM-CPM)")


################################################################################
### Fig. 2b: UCHL1 CRISPR KO efficiency in Caki-1 cells
################################################################################
laml.maf = 'tumor_only_rerun.mutect2.maf' 
vc_nonsyn = c('Frame_Shift_Del','Missense_Mutation','Nonsense_Mutation','Frame_Shift_Ins','In_Frame_Ins','Splice_Site','In_Frame_Del','Nonstop_Mutation')
laml = read.maf(maf = laml.maf, removeDuplicatedVariants = FALSE, vc_nonSyn = vc_nonsyn)

#Shows sample summry
getSampleSummary(laml)
write.csv(getSampleSummary(laml), "./mutation_summary_sample.csv")

#Shows gene summary
getGeneSummary(laml)
write.csv(getGeneSummary(laml), "./mutation_summary_gene.csv")

oncoplot(
  maf = laml, colors = vc_cols, genes = driver_list, sortByAnnotation = TRUE, drawRowBar = TRUE, drawColBar = TRUE, showTumorSampleBarcodes = TRUE, sortByMutation = FALSE, 
  draw_titv = FALSE, removeNonMutated = FALSE, genesToIgnore = c("TTN", "MUC16", "USH2A"))

lollipopPlot(
  maf = laml,
  gene = 'UCHL1',
  AACol = 'HGVSp',
  showMutationRate = TRUE
  #labelPos = 882
)


################################################################################
### Fig. 2c: Western blot of GAPDH and UCHL1
################################################################################
# Processed externally (ImageJ/FIJI).


################################################################################
### Fig. 2d - 2g: Proliferation Assays, qRT-PCR, and Inhibitor Viability
################################################################################
# CCK-8 absorbance tracking, two-way ANOVA statistics, and viability metrics are processed via dedicated Python workflow scripts.


################################################################################
### Fig. 2h: Schematic of cell-line-derived proteomics
################################################################################
# Schematic created in BioRender.


################################################################################
### Fig. 2i: Cell line sample correlations (Protein Level)
################################################################################
data_matrix_clean <- pro[complete.cases(pro), ]
data_matrix_clean <- data_matrix_clean[apply(data_matrix_clean, 1, sd) > 0, ]
cor_matrix <- cor(data_matrix_clean, method = "pearson")

pdf("./ccRCC_UCHL1_manuscript/Data/Figure2i_Correlation_Heatmap.pdf", width = 8, height = 8)
pheatmap(cor_matrix, 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Cell Line Protein Correlations")
dev.off()


################################################################################
### Fig. 2j: Protein abundances of UCHL1 and VHL across cell lines
################################################################################
df_2j <- data.frame(
  Sample = colnames(pro),
  UCHL1 = pro["UCHL1", ],
  VHL = pro["VHL", ]
) %>% pivot_longer(cols = c("UCHL1", "VHL"), names_to = "Gene", values_to = "Log2_Abundance")

pdf("./ccRCC_UCHL1_manuscript/Data/Figure2j_UCHL1_VHL_Bar.pdf", width = 10, height = 5)
ggplot(df_2j, aes(x = Sample, y = Log2_Abundance, fill = Gene)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~Gene, scales = "free_y", ncol = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Protein Abundance (log2)", x = "Cell Lines")
dev.off()


################################################################################
### Fig. 2k: High correlations of CRISPR KO Clones C2 and C9
################################################################################
df_2k <- as.data.frame(pro)

pdf("./ccRCC_UCHL1_manuscript/Data/Figure2k_C2_C9_Correlation.pdf", width = 5, height = 5)
ggscatter(df_2k, x = "Caki.1_C2", y = "Caki.1_C9", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Caki-1 C2 (log2 protein)", ylab = "Caki-1 C9 (log2 protein)")
dev.off()


################################################################################
### Fig. 2l: Venn diagrams showing Differentially Expressed Proteins (DEPs)
################################################################################
# Load pre-calculated DEP lists generated from Cell_DE_final_2.R
ko1_down <- read.table("./ccRCC_UCHL1_manuscript/Data/ko_down_intersect_clean_FC0.7_FDR0.05.txt", stringsAsFactors = FALSE)$V1
ko2_down <- read.table("./ccRCC_UCHL1_manuscript/Data/ko2_down_intersect_clean_FC0.7_FDR0.05.txt", stringsAsFactors = FALSE)$V1

gene_sets <- list(
  `KO1 (C2) Down` = ko1_down,
  `KO2 (C9) Down` = ko2_down 
)
venn_fit <- euler(gene_sets)

pdf("./ccRCC_UCHL1_manuscript/Data/Figure2l_Venn_DEP.pdf", width = 5, height = 5)
plot(venn_fit, quantities = TRUE, fills = c("#FF9999", "#99CCFF"), alpha = 0.7,
     main = "Differentially Expressed Proteins by UCHL1 Status")
dev.off()