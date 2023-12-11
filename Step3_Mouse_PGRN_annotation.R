#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)

setwd("/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/integration_with_Axl")
# remove cluster 14 (no significant marker genes), cluster 15 and 16 (doublets)
PGRN <- subset(PGRN, idents=c("14","15","16"), invert = T)
saveRDS(PGRN, file = 'Mouse_FTD_integrated_PCA_0.1_no141516.rds')

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

a<-as.data.frame(table(PGRN$Condition,PGRN$seurat_clusters))
colnames(a)<-c("clusters","seurat_clusters","cell.no")
agg<-aggregate(cell.no~clusters,a,sum)
a$cluster.total <- agg$cell.no[match(a$clusters,agg$clusters)]
a$ratio<-a$cell.no/a$cluster.total


ggplot(a,aes(x=clusters, y=ratio, fill=seurat_clusters))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Genotype")+
  ylab("Cell type ratio per genotype") + RotatedAxis()

ggsave("genotype_seurat_clusters_distribution_new.pdf",plot=last_plot(),
       width=4,height=5,units="in")

View(table(PGRN$seurat_clusters))


pdf("annotation_1.pdf", width=11, height=4.5)
DotPlot(PGRN, features = c("Plp1", "Mbp", "Mobp","Gpc5", "Plpp3","Pla2g7", "Snap25", "Cdh18","Sgcz","Kcnip4","Ntng1","Kcnc2",
                             "Cx3cr1", "P2ry12", "Csf1r","Vcan", "Pdgfra", "mt-Co3","mt-Atp6","Ttr", "Htr2c","Kcnq5","Dlgap2","Meis2",
                             "Nwd2","Kctd8","Scube1","Hs3st4","Foxp2","Lypd6b","Cped1","Lama1","Atp13a5","Ebf1","Kcnmb2","Kcnip1")) + RotatedAxis()
dev.off()

pdf("umap_by_condition.pdf", width=11.5, height=2.6)
DimPlot(PGRN, reduction = "umap", split.by = "Condition", ncol = 6, label = T)
dev.off()


Idents(PGRN) <- "seurat_clusters"
DefaultAssay(PGRN) <- "RNA"
VlnPlot(PGRN, features = c("Mertk"), split.by = "Condition")
VlnPlot(PGRN, features = c("Mertk"), split.by = "Condition", pt.size = 0)


