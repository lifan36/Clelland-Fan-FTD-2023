#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)
setwd("/athena/ganlab/scratch/lif4001/Human_PGRN/data_analysis/integration_2023")

Human_FTD_integrated <- readRDS("Human_FTD_integrated_PCA_0.1.rds")

# remove cluster 15 (doublets).
Human_FTD_integrated <- subset(Human_FTD_integrated, idents = c("15"), invert = T)

# Relevel object@ident
Human_FTD_integrated$Condition <- factor(x = Human_FTD_integrated$Condition, levels = c("No-pathology","FTD"))
Human_FTD_integrated$orig.ident <- factor(x = Human_FTD_integrated$orig.ident, levels = c("Non_1","Non_2","Non_3","Non_4","Non_5","Non_6","Non_7","Non_8",
                                                                                          "FTD_1","FTD_3","FTD_4","FTD_6","FTD_7","FTD_8","FTD_9" ))



Human_FTD_integrated <- RenameIdents(Human_FTD_integrated,
                                     `0` = "oligodendrocytes", `1`="astrocytes", `2`="excitatory neurons", `3`="microglia",
                                     `4`="OPCs", `5`="inhibitory neurons", `6`="inhibitory neurons", `7`="excitatory neurons",
                                     `8`="excitatory neurons", `9`="excitatory neurons", `10`="endothelial cells", `11`="excitatory neurons",
                                     `12`="excitatory neurons", `13`="inhibitory neurons", `14`="inhibitory neurons"
)

Human_FTD_integrated$celltype.orig.ident <- paste(Idents(Human_FTD_integrated), Human_FTD_integrated$orig.ident, sep = "_")
Human_FTD_integrated$celltype <- Idents(Human_FTD_integrated)

saveRDS(Human_FTD_integrated, file = "Human_FTD_integrated_Annotation.rds")

# Figure 1B
Idents(Human_FTD_integrated) <- "celltype"
DefaultAssay(Human_FTD_integrated) <- 'RNA'
pdf("Human_FTD_integrated_umap_annotation_noLabel.pdf", width=6, height=4)
DimPlot(Human_FTD_integrated, reduction = 'umap', label = F)
dev.off()


# Figure 1D
pdf("annotation_DotPlot.pdf", width=10.5, height=3.2)
DotPlot(data, features = c("PLP1", "MBP", "MOBP","AQP4","GFAP","SLC17A7", "CAMK2A", "NRGN","CD74","CSF1R","C3","PDGFRA","VCAN","GAD1", "GAD2",
                           "EBF1","IGFBP7","FLT1")) + RotatedAxis()
dev.off()


Cluster_EN <- subset(Human_FTD_integrated, idents = "excitatory neurons")
Cluster_IN <- subset(Human_FTD_integrated, idents = "inhibitory neurons")
Cluster_MG <- subset(Human_FTD_integrated, idents = "microglia")
Cluster_AST <- subset(Human_FTD_integrated, idents = "astrocytes")
Cluster_OL <- subset(Human_FTD_integrated, idents = "oligodendrocytes")
Cluster_OPC <- subset(Human_FTD_integrated, idents = "OPCs")
Cluster_EC <- subset(Human_FTD_integrated, idents = "endothelial cells")

saveRDS(Cluster_EN, file = "Human_FTD_EN_subset.rds")
saveRDS(Cluster_IN, file = "Human_FTD_IN_subset.rds")
saveRDS(Cluster_MG, file = "Human_FTD_MG_subset.rds")
saveRDS(Cluster_AST, file = "Human_FTD_AST_subset.rds")
saveRDS(Cluster_OL, file = "Human_FTD_OL_subset.rds")
saveRDS(Cluster_OPC, file = "Human_FTD_OPC_subset.rds")
saveRDS(Cluster_EC, file = "Human_FTD_EC_subset.rds")

# Figure 1C
data <- Human_FTD_integrated
# calculate ratio of each genotype in each cell type cluster
a<-as.data.frame(table(data$Condition,data$celltype))
colnames(a)<-c("clusters","cell.type","cell.no")
agg<-aggregate(cell.no~clusters,a,sum)
a$cluster.total <- agg$cell.no[match(a$clusters,agg$clusters)]
a$ratio<-a$cell.no/a$cluster.total

ggplot(a,aes(x=clusters, y=ratio, fill=cell.type))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Genotype")+
  ylab("Cell type ratio per genotype") + RotatedAxis()

ggsave("genotype_celltype_distribution.pdf",plot=last_plot(),path="/athena/ganlab/scratch/lif4001/Human_PGRN/data_analysis/integration_2023",
       width=4,height=4,units="in")


data <- Human_FTD_integrated
# calculate ratio of each sample in each cell type cluster
a<-as.data.frame(table(data$orig.ident,data$celltype))
colnames(a)<-c("clusters","cell.type","cell.no")
agg<-aggregate(cell.no~clusters,a,sum)
a$cluster.total <- agg$cell.no[match(a$clusters,agg$clusters)]
a$ratio<-a$cell.no/a$cluster.total

ggplot(a,aes(x=clusters, y=ratio, fill=cell.type))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Sample")+
  ylab("Cell type ratio per sample") + RotatedAxis()

ggsave("sample_celltype_distribution.pdf",plot=last_plot(),path="/athena/ganlab/scratch/lif4001/Human_PGRN/data_analysis/integration_2023",
       width=6,height=4,units="in")


Idents(Cluster_EN) <- "Condition"
Idents(Cluster_IN) <- "Condition"
Idents(Cluster_MG) <- "Condition"
Idents(Cluster_AST) <- "Condition"
Idents(Cluster_OL) <- "Condition"
Idents(Cluster_OPC) <- "Condition"
Idents(Cluster_EC) <- "Condition"


EN_FTD_vs_Non_DEGs <- FindMarkers(Cluster_EN, ident.1 = "FTD", ident.2 = "No-pathology", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
write.csv(EN_FTD_vs_Non_DEGs, "EN_FTD_vs_Non_DEGs.csv")
IN_FTD_vs_Non_DEGs <- FindMarkers(Cluster_IN, ident.1 = "FTD", ident.2 = "No-pathology", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
write.csv(IN_FTD_vs_Non_DEGs, "IN_FTD_vs_Non_DEGs.csv")
MG_FTD_vs_Non_DEGs <- FindMarkers(Cluster_MG, ident.1 = "FTD", ident.2 = "No-pathology", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
write.csv(MG_FTD_vs_Non_DEGs, "MG_FTD_vs_Non_DEGs.csv")
AST_FTD_vs_Non_DEGs <- FindMarkers(Cluster_AST, ident.1 = "FTD", ident.2 = "No-pathology", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
write.csv(AST_FTD_vs_Non_DEGs, "AST_FTD_vs_Non_DEGs.csv")
OL_FTD_vs_Non_DEGs <- FindMarkers(Cluster_OL, ident.1 = "FTD", ident.2 = "No-pathology", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
write.csv(OL_FTD_vs_Non_DEGs, "OL_FTD_vs_Non_DEGs.csv")
OPC_FTD_vs_Non_DEGs <- FindMarkers(Cluster_OPC, ident.1 = "FTD", ident.2 = "No-pathology", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
write.csv(OPC_FTD_vs_Non_DEGs, "OPC_FTD_vs_Non_DEGs.csv")
EC_FTD_vs_Non_DEGs <- FindMarkers(Cluster_EC, ident.1 = "FTD", ident.2 = "No-pathology", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
write.csv(EC_FTD_vs_Non_DEGs, "EC_FTD_vs_Non_DEGs.csv")


# Figure 1H
Idents(Human_FTD_integrated) <- "celltype"
pdf("Human_FTD_integrated_VlnPlot_Tardbp_celltype_split_by_Condition.pdf", width=8, height=6)
VlnPlot(Human_FTD_integrated, features = c("MERTK","AXL"), pt.size = 0.1, ncol = 2, split.by = "Condition")
VlnPlot(Human_FTD_integrated, features = c("MERTK","AXL"), pt.size = 0, ncol = 2, split.by = "Condition")
dev.off()


# Figure 1F:
setwd("/Users/lifan/Desktop/data_analysis/Human_PGRN/integration_2023")

marker <- read.csv(file = "MG_FTD_vs_Non_DEGs.csv", header=T,row.names =1)

# Identify DEGs for cluster 5 markers
marker$colours <- c("NC")
marker$colours[marker$avg_log2FC >= 0.1 & marker$p_val_adj <= 0.05] <- c("UP")
marker$colours[marker$avg_log2FC <= -0.1 & marker$p_val_adj <= 0.05] <- c("DN")

# Selected genes to highlight
genes_select_mature <- c("ARPC1B","ARPC1B","CRK","DOCK1","FCGR1A","FCGR2A","FCGR3A/FCGR3B","FGR","FYB1","HCK","LCP2","NCK2","PIK3R1","PLD3","PRKCE","PRKCH","PRKD3","PTK2B","PXN","TLN1","VASP","VAV3","CD163","DPYD","TMEM163","FMN1","FKBP5","F13A1","CCDC26","GRID2","SLC2A3","C5orf17","FOS","FRMD4A","SYNDIG1","SORL1","LRMDA","MERTK","AXL")
genes_to_plot_mature <- marker[row.names(marker) %in% genes_select_mature, ]
genes_to_plot_mature$Cluster <- "Mature"

genes_select_immature <- c("ADCY7","AKT3","APOE","CACNB2","CACNB4","CAMK2G","CDH10","CDH12","CDH13","CDH18","CDH23","CDH9","CLASP2","CNTNAP2","CREBBP","DAB1","EPHA6","EPHB2","GRIA2","GRIA4","GRIN2A","GRM5","GRM7","GUCY1A1","ITPR1","KALRN","LYN","MARCKS","NLGN1","NRXN1","NRXN3","NTRK2","PAK1","PIK3CD","RAPGEF1","SOS1","SYT1","SYT17","TIAM1")
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
  ggtitle("MG_GRN+/-_vs_Ctrl")+
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

ggsave("MG_FTD_vs_Non_volcanoPlot.pdf", plot = last_plot(), device = "pdf",
       scale = 1.8, width = 19, height = 19, units = c("in"),
       dpi = 600, limitsize = FALSE)




