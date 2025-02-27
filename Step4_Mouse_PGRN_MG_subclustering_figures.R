
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)

setwd("/Users/lifan/Desktop/data_analysis/mouse_pgrn/data_analysis/integration_with_Axl/")

PGRN <- readRDS(file = "Mouse_FTD_integrated_PCA_0.1_no141516.rds")

# Relevel object@ident
PGRN$Condition <- factor(x = PGRN$Condition, levels = c("Ctrl", "PGRN_KO","MerTK_KO","PGRN_MetTK_KO","Axl_KO","PGRN_Axl_KO"))
PGRN$Sample_Name <- factor(x = PGRN$Sample_Name, levels = c("Ctrl_1","Ctrl_2","Ctrl_3","Ctrl_4",
                                                            "PGRN_KO_1","PGRN_KO_2","PGRN_KO_3","PGRN_KO_4",
                                                            "MerTK_KO_1","MerTK_KO_2","MerTK_KO_3","MerTK_KO_4",
                                                            "PGRN_MetTK_KO_1","PGRN_MetTK_KO_2","PGRN_MetTK_KO_3","PGRN_MetTK_KO_4",
                                                            "Axl_KO_1","Axl_KO_2","Axl_KO_3","Axl_KO_4",
                                                            "PGRN_Axl_KO_1","PGRN_Axl_KO_2","PGRN_Axl_KO_3","PGRN_Axl_KO_4"))

# MG subclutering
MG <- subset(PGRN, idents = "4")
DefaultAssay(MG) <- 'integrated'
MG <- ScaleData(MG, verbose = FALSE)
MG <- RunPCA(MG, features = VariableFeatures(object = MG), verbose = FALSE)
ElbowPlot(MG)
MG <- FindNeighbors(MG, dims = 1:8)
MG <- FindClusters(MG, resolution = 0.1)
MG <- RunUMAP(MG, dims = 1: 8)
# rename cluster
Idents(MG) <- "seurat_clusters"
n <- dim(table(MG@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
MG@active.ident <- plyr::mapvalues(x = MG@active.ident, from = current.cluster.ids, to = new.cluster.ids)
MG@active.ident <- factor(MG@active.ident, levels=1:n)
saveRDS(MG, file = 'PGRN_MG_reclusted_res0.1_local.rds')

MG <- readRDS(file = "PGRN_MG_reclusted_res0.1_local.rds")
#pdf("PGRN_MG_umap.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
#dev.off()
#pdf("PGRN_MG_umap_Condition.pdf", width=9, height=3)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
#dev.off()
#pdf("PGRN_MG_umap_Sample.pdf", width=8, height=7.5)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
#dev.off()
DefaultAssay(MG) <- 'RNA'
PGRN_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(PGRN_MG_markers, "PGRN_MG_markers.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "MG_subcluster_cell_counts.csv")

VlnPlot(MG, features = c("P2ry12","Cx3cr1","Csf1r","Apoe"), pt.size = 0, ncol = 4)

Idents(MG) <- "Condition"
VlnPlot(MG, features = c("P2ry12","Cx3cr1","Csf1r","Apoe"), pt.size = 0, ncol = 4)
VlnPlot(MG, features = c("Grn","Axl","Mertk"), pt.size = 0.1, ncol = 3)

Idents(PGRN) <- "seurat_clusters"
VlnPlot(PGRN, features = c("P2ry12","Cx3cr1","Csf1r","Apoe"), pt.size = 0, ncol = 4)
VlnPlot(PGRN, features = c("Grn","Axl","Mertk"), pt.size = 0, ncol = 3)

MG <- RenameIdents(MG, `1` = "Homeostatic", `2`="DAM", `3`="Oligodendrocytes", `4`="mt- high", `5`="Macrophage", `6`="T cell")

pdf("PGRN_MG_umap_annotation.pdf", width=5.3, height=3)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()

pdf("PGRN_MG_umap_Condition_annotation.pdf", width=11, height=2.5)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = F, ncol = 6)
dev.off()


DefaultAssay(MG) <- "RNA"
pdf("PGRN_DotPlot_blue_red.pdf", width=8, height=3.5)
DotPlot(MG, features = c("P2ry12","Cx3cr1","Selplg",
                               "Arhgap24","Cacna1a","Apoe",
                               "Plp1","Mbp","St18",
                               "mt-Co3","mt-Atp6","mt-Co2",
                               "F13a1","Mrc1","Dab2",
                               "Skap1","Themis","Grap2"
                               
)) + RotatedAxis() + scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred")
dev.off()

df <- FindMarkers(MG, ident.1 = "DAM", ident.2 = "Homeostatic", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "MG_DAM_vs_Homeo_DEGs.csv")

# June 22, 2023 remove non-mg cluster
MG <- readRDS(file = "PGRN_MG_reclusted_res0.1_local.rds")
MG <- subset(MG, idents=c("1","2"))

setwd("/Users/lifan/Desktop/data_analysis/mouse_pgrn/data_analysis/integration_with_Axl/MG_local")
# Figure 4D
pdf("PGRN_MG_umap.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
# Figure 4E
pdf("PGRN_MG_umap_Condition.pdf", width=9, height=2.3)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 6)
dev.off()

DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)

saveRDS(MG, file = 'PGRN_MG_reclusted_res0.1_local_clean.rds')


pdf("PGRN_MG_umap_Sample.pdf", width=8, height=7.5)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()

# Figure 4F
prop.table(table(Idents(MG)))
table <- as.matrix(table(Idents(MG), MG$Condition))
sum=colSums(table)
table<-rbind(table[1,]/sum, table[2,]/sum)
barplot(as.matrix(table))

write.csv(table, "PGRN_MG_cell_ratio.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "MG_subcluster_cell_counts.csv")

# Figure 4H
DefaultAssay(MG) <- "RNA"
pdf("PGRN_MG_featurePlot.pdf", width=6, height=5)
FeaturePlot(MG, features = c("P2ry12","Cx3cr1","Apoe","Gpnmb"), ncol = 2)
dev.off()

# Figure 4J
pdf("PGRN_MG_DAM_DotPlot_by_Condition.pdf", width=9, height=3)
DotPlot(MG, features = c("Apoe","Trem2","Neat1","Vps13c", "Cd68","Csf2ra","Lyz2","Dtnbp1","Aplp2","Hif1a","Dgkz","Myo1e","Ctsz","Ctsa","Rap2b","Cd83","Galnt3","Csf1")) + scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred") + RotatedAxis()
dev.off()

# Figure 4G volcano plot
setwd("/Users/lifan/Desktop/data_analysis/mouse_pgrn/data_analysis/integration_with_Axl/")

marker <- read.csv(file = "MG_DAM_vs_Homeo_DEGs.csv", header=T,row.names =1)

marker$colours <- c("NC")
marker$colours[marker$avg_log2FC >= 0.25 & marker$p_val_adj <= 0.05] <- c("UP")
marker$colours[marker$avg_log2FC <= -0.25 & marker$p_val_adj <= 0.05] <- c("DN")

# Selected genes to highlight
genes_select_mature <- c("Apobec1","Csmd3","Gpnmb","Igf1","Nav3","Myo5a","Lgals3","Ctsb","Arhgap24","Ifi207","Apoe")
genes_to_plot_mature <- marker[row.names(marker) %in% genes_select_mature, ]
genes_to_plot_mature$Cluster <- "Mature"

genes_select_immature <- c("Siglech","Siglech","Selplg","Selplg","Selplg","P2ry12","Csf1r","Sall1")
genes_to_plot_immature <- marker[row.names(marker)  %in% genes_select_immature, ]
genes_to_plot_immature$Cluster <- c("Immature")

genes_to_plot <- rbind(genes_to_plot_mature, genes_to_plot_immature)

# Set color palette
my_color <- c("#2B8CBE", "#D7301F", "skyblue","seashell3", "plum1")
my_color_1 <- c("#D7301F","Grey", "#2B8CBE")

# Plot volcano plot

ggplot() + 
  geom_point(data=marker, aes(x=avg_log2FC, y=-log10(p_val_adj), colour=colours),
             shape=19, alpha=1, size=1) +
  scale_color_manual(values = my_color_1,
                     name="DEGs",
                     breaks=rev(names(table(marker$colours))),
                     labels=rev(names(table(marker$colours)))) +
  geom_point(data=genes_to_plot,
             aes(x=avg_log2FC, y=-log10(p_val_adj)),
             shape=19, alpha=1, size=3) +
  geom_text_repel(data=genes_to_plot,
                  aes(x=avg_log2FC, y=-log10(p_val_adj), label = row.names(genes_to_plot)), 
                  color="black", fontface = 'bold',size = 5, box.padding = 0.5,
                  point.padding = 0.5, segment.size=0.25, segment.colour="black") +
  ylab("-Log10[FDR]") + xlab("Log2FC") +
  ggtitle("Cluster 2")+
  theme_bw()+
  theme(panel.grid.major.x  = element_blank(),
        panel.grid.major.y  = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.text.x = element_text(colour = "black", size=15),
        axis.text.y = element_text(colour = "black", size=15),
        axis.title.x = element_text(colour = "black", size=15),
        axis.title.y = element_text(colour = "black", size=15),
        plot.title = element_text(size = 15, face = "bold")) +
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks=seq(-4, 4, 1), limits=c(-4, 4))+
  NoLegend()

ggsave("MG_DAM_vs_Homeo_volcanoPlot_new_axl_1.pdf", plot = last_plot(), device = "pdf",
       scale = 0.6, width = 9, height = 9, units = c("in"),
       dpi = 600, limitsize = FALSE)








