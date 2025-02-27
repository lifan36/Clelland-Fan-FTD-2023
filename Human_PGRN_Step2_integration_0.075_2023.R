
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)

#load in data from Cell Ranger or other counts data ====

#for loading Cell Ranger counts:
setwd("/athena/ganlab/scratch/lif4001/Human_PGRN/data_analysis/DF_2ndRound")
Human_control_1 <- readRDS(file = "Human_control_1_singlets_PCA.rds")
Human_control_2 <- readRDS(file = "Human_control_2_singlets_PCA.rds")
Human_control_3 <- readRDS(file = "Human_control_3_singlets_PCA.rds")
Human_control_4 <- readRDS(file = "Human_control_4_singlets_PCA.rds")
Human_control_5 <- readRDS(file = "Human_control_5_singlets_PCA.rds")
Human_control_6 <- readRDS(file = "Human_control_6_singlets_PCA.rds")
Human_control_7 <- readRDS(file = "Human_control_7_singlets_PCA.rds")
Human_control_8 <- readRDS(file = "Human_control_8_singlets_PCA.rds")
Human_PGRN_1 <- readRDS(file = "Human_PGRN_1_singlets_PCA.rds")
Human_PGRN_3 <- readRDS(file = "Human_PGRN_3_singlets_PCA.rds")
Human_PGRN_4 <- readRDS(file = "Human_PGRN_4_singlets_PCA.rds")
Human_PGRN_6 <- readRDS(file = "Human_PGRN_6_singlets_PCA.rds")
Human_PGRN_7 <- readRDS(file = "Human_PGRN_7_singlets_PCA.rds")
Human_PGRN_8 <- readRDS(file = "Human_PGRN_8_singlets_PCA.rds")
Human_PGRN_9 <- readRDS(file = "Human_PGRN_9_singlets_PCA.rds")

setwd("/athena/ganlab/scratch/lif4001/Human_PGRN/data_analysis/integration_2023")

Ctrl <- c(Human_control_1, Human_control_2, Human_control_3, Human_control_4, Human_control_5,Human_control_6, Human_control_7, Human_control_8)
anchors_Ctrl <- FindIntegrationAnchors(object.list = Ctrl, dims = 1:30)
Ctrl_integrated <- IntegrateData(anchorset = anchors_Ctrl, dims = 1:30)
rm(Human_control_1, Human_control_2, Human_control_3, Human_control_4, Human_control_5,Human_control_6, Human_control_7, Human_control_8, Ctrl)

FTD <- c(Human_PGRN_1, Human_PGRN_3,Human_PGRN_4,Human_PGRN_6,Human_PGRN_7,Human_PGRN_8,Human_PGRN_9)
anchors_FTD <- FindIntegrationAnchors(object.list = FTD, dims = 1:30)
FTD_integrated <- IntegrateData(anchorset = anchors_FTD, dims = 1:30)
rm(Human_PGRN_1, Human_PGRN_3,Human_PGRN_4,Human_PGRN_6,Human_PGRN_7,Human_PGRN_8,Human_PGRN_9, FTD)

Human_FTD <- c(Ctrl_integrated, FTD_integrated)
anchors_Human_FTD <- FindIntegrationAnchors(object.list = Human_FTD, dims = 1:30)
Human_FTD_integrated <- IntegrateData(anchorset = anchors_Human_FTD, dims = 1:30)
rm(Ctrl_integrated, FTD_integrated, Human_FTD)

#saveRDS(Human_FTD_integrated, file = "Human_FTD_integrated.rds")

#Human_FTD_integrated <- readRDS("Human_FTD_integrated.rds")
DefaultAssay(Human_FTD_integrated) <- 'integrated'

#Human_FTD_integrated <- NormalizeData(Human_FTD_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
#Human_FTD_integrated <- FindVariableFeatures(Human_FTD_integrated, selection.method = "vst", nfeatures = 3000)

Human_FTD_integrated <- ScaleData(Human_FTD_integrated, verbose = FALSE)
Human_FTD_integrated <- RunPCA(Human_FTD_integrated, features = VariableFeatures(object = Human_FTD_integrated), verbose = FALSE)

Human_FTD_integrated <- FindNeighbors(Human_FTD_integrated, dims = 1:20)
Human_FTD_integrated <- FindClusters(Human_FTD_integrated, resolution = 0.075)
Human_FTD_integrated <- RunUMAP(Human_FTD_integrated, dims = 1: 20)

str(Human_FTD_integrated)

DefaultAssay(Human_FTD_integrated) <- 'RNA'
Human_FTD_integrated <- NormalizeData(Human_FTD_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
Human_FTD_integrated <- ScaleData(Human_FTD_integrated, features = rownames(Human_FTD_integrated))

pdf("Human_FTD_integrated_umap.pdf", width=6, height=4)
DimPlot(Human_FTD_integrated, reduction = 'umap', label = T)
dev.off()
pdf("Human_FTD_integrated_umap_split_individual.pdf", width=13, height=3)
DimPlot(Human_FTD_integrated, reduction = "umap", split.by = "orig.ident", label = T, ncol = 8)
dev.off()
pdf("Human_FTD_integrated_umap_split_Condition.pdf", width=8.5, height=4)
DimPlot(Human_FTD_integrated, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()

saveRDS(Human_FTD_integrated, file = 'Human_FTD_integrated_PCA_0.1.rds')

#Human_FTD_integrated <- readRDS("Human_FTD_integrated_PCA_0.1.rds")

DefaultAssay(Human_FTD_integrated) <- 'RNA'
pdf("Human_FTD_integrated_umap_test.pdf", width=4, height=3)
DimPlot(Human_FTD_integrated, reduction = 'umap', label = T)
dev.off()


#Add marker genes

pdf("Human_FTD_integrated_annotation_combine.pdf", width=12, height=6)
sig_all<-c("SYT1","SNAP25","GRIN1","SLC17A7", "CAMK2A", "NRGN","GAD1", "GAD2","PLP1", "MBP", "MOBP","AQP4","GFAP", 
           "CD74","CSF1R","C3","PDGFRA","VCAN","EBF1","IGFBP7","FLT1","CLDN5")
markers.to.plot <- as.matrix(sig_all)
DotPlot(object = Human_FTD_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()
dev.off()

DefaultAssay(Human_FTD_integrated) <- 'RNA'

FTD_markers <- FindAllMarkers(Human_FTD_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1, test.use = "MAST")
write.csv(FTD_markers, "FTD_markers.csv")

write.csv(table(Human_FTD_integrated$seurat_clusters, Human_FTD_integrated$orig.ident), "cell_counts_cluster_sample.csv")
