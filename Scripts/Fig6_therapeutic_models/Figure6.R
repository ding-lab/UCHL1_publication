################################################################################
### Figure 6. Therapeutic targeting of UCHL1 with inhibitors consistently 
### restrained tumor growth in ccRCC cell lines, CDX, and PDX models.
################################################################################
# a, IMP-1710 dismantles UCHL1-mediated proliferative and inflammatory networks (Sankey).
# b, CCK-8 proliferation assay of 786-O and Caki-1 cells treated with IMP-1710.
# c, UCHL1 knockout suppresses Caki-1 tumor growth in vivo.
# d, Pharmacologic targeting (IMP-1710 + sunitinib) suppresses Caki-1 CDX growth.
# e, Ligand and RTK expressions in human snRNA-seq (Dot plot).
# f, UCHL1-inhibitor treatment is effective in PDX models (Schematic and Plot).
################################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(ggalluvial)
library(patchwork)

################################################################################
### Fig. 6a: IMP-1710 Downregulated DEPs (Sankey Plot)
################################################################################
enrich_res <- readRDS("./ccRCC_UCHL1_manuscript/Data/enrich_res_IMP1710_down.rds")

# Extract top terms for visualization
top_terms <- enrich_res@result[1:10, ]
df_sankey <- do.call(rbind, lapply(1:nrow(top_terms), function(i) {
  data.frame(
    Pathway = rep(top_terms$Description[i], length(unlist(strsplit(top_terms$geneID[i], "/")))),
    Gene = unlist(strsplit(top_terms$geneID[i], "/"))
  )
}))

pdf("./ccRCC_UCHL1_manuscript/Data/Figure6a_Sankey_IMP1710_Down.pdf", width = 12, height = 25)
ggplot(df_sankey, aes(axis1 = Gene, axis2 = Pathway)) +
  geom_alluvium(aes(fill = Pathway), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.5) +
  theme_minimal() +
  ggtitle("Sankey: IMP-1710 Downregulated Hallmark Pathways and Genes") +
  scale_x_discrete(limits = c("Genes", "Pathways"), expand = c(.05, .05))
dev.off()


################################################################################
### Fig. 6b: CCK-8 proliferation assay (In Vitro IMP-1710 treatment)
################################################################################
# CCK-8 absorbance, blank subtraction, and two-way repeated-measures ANOVA 
# processed via dedicated Python workflow scripts.


################################################################################
### Fig. 6c: UCHL1 knockout suppresses tumor growth in vivo
################################################################################
# Subcutaneous tumor growth tracking, caliper measurement normalization, and 
# pairwise two-way ANOVA statistics processed via dedicated Python workflow scripts.


################################################################################
### Fig. 6d: Pharmacologic targeting suppresses CDX growth in vivo
################################################################################
# CDX relative tumor growth analysis (IMP-1710, sunitinib, combination), 
# AUC, and linear mixed-effects models processed via dedicated Python workflow scripts.


################################################################################
### Fig. 6e: Ligand and RTK expressions in human snRNA-seq (Dot Plot)
################################################################################
cc_sub <- readRDS("./ccRCC_UCHL1_manuscript/Data/ccRCC_snMultiome_snRNA_x72_Apr2025_merged_RNA_only_annotated_UCHL1.rds")

# Subset target populations: Tumor, myCAF (FIB), and VEC (EC)
ccc_sub <- subset(cc_sub, idents = c("TUMOR_ccRCC_UCHL1_high_cluster", 
                                     "TUMOR_ccRCC_UCHL1_low_cluster", 
                                     "FIB", "EC"))

# Refine Ident labels for clarity
ccc_sub@meta.data$Plot_Group <- plyr::mapvalues(
  as.character(ccc_sub@active.ident),
  from=c("TUMOR_ccRCC_UCHL1_high_cluster", "TUMOR_ccRCC_UCHL1_low_cluster", "FIB", "EC"),
  to=c("Tumor (U-high)", "Tumor (U-low)", "myCAF", "VEC")
)
Idents(ccc_sub) <- factor(ccc_sub@meta.data$Plot_Group, 
                          levels = c("VEC", "myCAF", "Tumor (U-low)", "Tumor (U-high)"))

# Define Ligand and RTK panel (PDGF and VEGF families)
rtk_ligand_panel <- c("PDGFA", "PDGFB", "PDGFC", "PDGFRA", "PDGFRB", 
                      "VEGFA", "VEGFB", "VEGFC", "FLT1", "KDR", "FLT4", "MET", "AXL")

pdf("./ccRCC_UCHL1_manuscript/Data/Figure6e_RTK_Ligand_DotPlot.pdf", width = 8, height = 5)
DotPlot(ccc_sub, features = rtk_ligand_panel, group.by = "Plot_Group", cols = c("blue", "red")) + 
  RotatedAxis() + 
  labs(title = "Ligand and RTK expressions (Tumor, myCAF, VEC)", x = "Genes", y = "Cell Types")
dev.off()


################################################################################
### Fig. 6f: UCHL1-inhibitor treatment is effective in the PDX models
################################################################################
# Left Panel: Schematic of experimental workflow created in BioRender.
# Right Panel: PDX relative tumor growth analysis (IMP-1710, cabozantinib, combination), 
# AUC, and linear mixed-effects models processed via dedicated Python workflow scripts.