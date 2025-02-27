
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)

setwd("/Users/lifan/Desktop/data_analysis/Human_PGRN/integration_2023")
MG <- readRDS(file = "Human_FTD_MG_subset.rds")
DefaultAssay(MG) <- 'integrated'
MG <- ScaleData(MG, verbose = FALSE)
MG <- RunPCA(MG, features = VariableFeatures(object = MG), verbose = FALSE)
ElbowPlot(MG)
MG <- FindNeighbors(MG, dims = 1:13)
MG <- FindClusters(MG, resolution = 0.15)
MG <- RunUMAP(MG, dims = 1: 13)
# rMGame cluster
Idents(MG) <- "seurat_clusters"
n <- dim(table(MG@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
MG@active.ident <- plyr::mapvalues(x = MG@active.ident, from = current.cluster.ids, to = new.cluster.ids)
MG@active.ident <- factor(MG@active.ident, levels=1:n)

DimPlot(MG, reduction = 'umap', label = T)

DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 2)

# discard tiny clusters 6 and 7, 141 and 55 cells respectively. 
MG <- subset(MG, idents= c("1","2","3","4","5"))

# Figure 2A
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 2)

write.csv(table(MG$seurat_clusters, MG$orig.ident), "Human_PGRN_MG_subcluster_cell_counts_res0.15.csv")

DefaultAssay(MG) <- "RNA"
VlnPlot(MG, features = c("P2RY12","C3","SKAP1","MRC1","MT-CO1"))
VlnPlot(MG, features = c("P2RY12","CD74","MBP","ST18"), pt.size = 0)
FeaturePlot(MG, features = c("TYROBP","CD74","MBP","ST18"), label = T)
FeaturePlot(MG, features = c("P2RY12","C3","SKAP1","MRC1","MT-CO1"), label = T)

markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0, only.pos = F)
markers <- markers[markers$p_val_adj < 0.05,]
write.csv(markers, "MG_markers_res0.15_all.csv")


markers <- FindMarkers(MG, logfc.threshold = 0.1, test.use = "MAST", min.pct = 0, only.pos = F, ident.1 = "2")
write.csv(markers, "MG2_vs_all_others_DEGs_res0.15.csv")
markers <- FindMarkers(MG, logfc.threshold = 0.1, test.use = "MAST", min.pct = 0, only.pos = F, ident.1 = "3")
write.csv(markers, "MG3_vs_all_others_DEGs_res0.15.csv")
markers <- FindMarkers(MG, logfc.threshold = 0.1, test.use = "MAST", min.pct = 0, only.pos = F, ident.1 = "5")
write.csv(markers, "MG5_vs_all_others_DEGs_res0.15.csv")
markers <- FindMarkers(MG, logfc.threshold = 0.1, test.use = "MAST", min.pct = 0, only.pos = F, ident.1 = "1")
write.csv(markers, "MG1_vs_all_others_DEGs_res0.15.csv")
markers <- FindMarkers(MG, logfc.threshold = 0.1, test.use = "MAST", min.pct = 0, only.pos = F, ident.1 = "4")
write.csv(markers, "MG4_vs_all_others_DEGs_res0.15.csv")


# Figure 2C
setwd("/Users/lifan/Desktop/data_analysis/Human_PGRN/integration_2023")
marker <- read.csv(file = "MG2_vs_all_others_DEGs_res0.15.csv", header=T,row.names =1)
marker$colours <- c("NC")
marker$colours[marker$avg_log2FC >= 0.1 & marker$p_val_adj <= 0.05] <- c("UP")
marker$colours[marker$avg_log2FC <= -0.1 & marker$p_val_adj <= 0.05] <- c("DN")

# Selected genes to highlight
genes_select_mature <- c("KCNIP4","NRG3","DPP10","LRRTM4","LRP1B","NRXN1","FAM155A","CSMD1","ATRNL1","RBFOX1","LSAMP")
genes_to_plot_mature <- marker[row.names(marker) %in% genes_select_mature, ]
genes_to_plot_mature$Cluster <- "Mature"

genes_select_immature <- c("CD163","TMEM163","DPYD","FMN1","PTPRG","CPM","C9orf84","VSIG4","MERTK","APOE")
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
  ggtitle("MG2_vs_all_others_DEGs")+
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
  scale_x_continuous(breaks=seq(-2, 2, 0.5), limits=c(-2, 2))+
  NoLegend()

ggsave("MG2_vs_all_others_DEGs_res0.15_VolcanoPlot.pdf", plot = last_plot(), device = "pdf",
       scale = 0.5, width = 10, height = 10, units = c("in"),
       dpi = 600, limitsize = FALSE)

# Figure 2D
marker <- read.csv(file = "MG3_vs_all_others_DEGs_res0.15.csv", header=T,row.names =1)
marker$colours <- c("NC")
marker$colours[marker$avg_log2FC >= 0.1 & marker$p_val_adj <= 0.05] <- c("UP")
marker$colours[marker$avg_log2FC <= -0.1 & marker$p_val_adj <= 0.05] <- c("DN")

# Selected genes to highlight
genes_select_mature <- c("MTRNR2L12","DHFR","SERPINE1","FLT1","GFAP","SLC27A6","PRPF38B","CLU","RGS1","SRRM2","PADI2")
genes_to_plot_mature <- marker[row.names(marker) %in% genes_select_mature, ]
genes_to_plot_mature$Cluster <- "Mature"

genes_select_immature <- c("DLEU7","GRID2","PDK4","CCDC26","SYNDIG1","TANC1","FCGBP","TBC1D4","FAM177B","CD74")
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
  ggtitle("MG3_vs_all_others_DEGs")+
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
  scale_x_continuous(breaks=seq(-2, 2, 0.5), limits=c(-2, 2))+
  NoLegend()

ggsave("MG3_vs_all_others_DEGs_res0.15_VolcanoPlot.pdf", plot = last_plot(), device = "pdf",
       scale = 0.5, width = 10, height = 10, units = c("in"),
       dpi = 600, limitsize = FALSE)

# Figure 2E
marker <- read.csv(file = "MG5_vs_all_others_DEGs_res0.15.csv", header=T,row.names =1)
marker$colours <- c("NC")
marker$colours[marker$avg_log2FC >= 0.1 & marker$p_val_adj <= 0.05] <- c("UP")
marker$colours[marker$avg_log2FC <= -0.1 & marker$p_val_adj <= 0.05] <- c("DN")

# Selected genes to highlight
genes_select_mature <- c("IL1RAPL1","PCDH9","ST18","CTNNA3","RNF220","PPP2R2B","SLC24A2","FRMD5","MAGI2","CADM2","NCAM2")
genes_to_plot_mature <- marker[row.names(marker) %in% genes_select_mature, ]
genes_to_plot_mature$Cluster <- "Mature"

genes_select_immature <- c("CCDC26","P2RY12","CPED1","SORL1","A2M","CX3CR1","SYNDIG1","IL6ST","KCNQ3","ATP8B4")
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
  ggtitle("MG5_vs_all_others_DEGs")+
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
  scale_x_continuous(breaks=seq(-2, 4, 1), limits=c(-2, 4))+
  NoLegend()

ggsave("MG5_vs_all_others_DEGs_res0.15_VolcanoPlot.pdf", plot = last_plot(), device = "pdf",
       scale = 2, width = 10, height = 10, units = c("in"),
       dpi = 600, limitsize = FALSE)





