
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)

#load in data from Cell Ranger or other counts data ====

#for loading Cell Ranger counts:
setwd("/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/DF_2ndRound")
Ctrl_1 <- readRDS(file = "Ctrl_1_singlets_PCA.rds")
Ctrl_2 <- readRDS(file = "Ctrl_2_singlets_PCA.rds")
Ctrl_3 <- readRDS(file = "Ctrl_3_singlets_PCA.rds")
Ctrl_4 <- readRDS(file = "Ctrl_4_singlets_PCA.rds")

PGRN_KO_1 <- readRDS(file = "PGRN_KO_1_singlets_PCA.rds")
PGRN_KO_2 <- readRDS(file = "PGRN_KO_2_singlets_PCA.rds")
PGRN_KO_3 <- readRDS(file = "PGRN_KO_3_singlets_PCA.rds")
PGRN_KO_4 <- readRDS(file = "PGRN_KO_4_singlets_PCA.rds")

MerTK_KO_1 <- readRDS(file = "MerTK_KO_1_singlets_PCA.rds")
MerTK_KO_2 <- readRDS(file = "MerTK_KO_2_singlets_PCA.rds")
MerTK_KO_3 <- readRDS(file = "MerTK_KO_3_singlets_PCA.rds")
MerTK_KO_4 <- readRDS(file = "MerTK_KO_4_singlets_PCA.rds")

Double_KO_1 <- readRDS(file = "Double_KO_1_singlets_PCA.rds")
Double_KO_2 <- readRDS(file = "Double_KO_2_singlets_PCA.rds")
Double_KO_3 <- readRDS(file = "Double_KO_3_singlets_PCA.rds")
Double_KO_4 <- readRDS(file = "Double_KO_4_singlets_PCA.rds")

Double_KO_1@active.ident

Idents(Double_KO_1) <- "Sample_Name"
Double_KO_1 <- RenameIdents(Double_KO_1, `Double_KO_1`="PGRN_MetTK_KO_1")
Double_KO_1$Sample_Name <- Idents(Double_KO_1)
Idents(Double_KO_1) <- "Condition"
Double_KO_1 <- RenameIdents(Double_KO_1, `Double_KO`="PGRN_MetTK_KO")
Double_KO_1$Condition <- Idents(Double_KO_1)

Idents(Double_KO_2) <- "Sample_Name"
Double_KO_2 <- RenameIdents(Double_KO_2, `Double_KO_2`="PGRN_MetTK_KO_2")
Double_KO_2$Sample_Name <- Idents(Double_KO_2)
Idents(Double_KO_2) <- "Condition"
Double_KO_2 <- RenameIdents(Double_KO_2, `Double_KO`="PGRN_MetTK_KO")
Double_KO_2$Condition <- Idents(Double_KO_2)

Idents(Double_KO_3) <- "Sample_Name"
Double_KO_3 <- RenameIdents(Double_KO_3, `Double_KO_3`="PGRN_MetTK_KO_3")
Double_KO_3$Sample_Name <- Idents(Double_KO_3)
Idents(Double_KO_3) <- "Condition"
Double_KO_3 <- RenameIdents(Double_KO_3, `Double_KO`="PGRN_MetTK_KO")
Double_KO_3$Condition <- Idents(Double_KO_3)

Idents(Double_KO_4) <- "Sample_Name"
Double_KO_4 <- RenameIdents(Double_KO_4, `Double_KO_4`="PGRN_MetTK_KO_4")
Double_KO_4$Sample_Name <- Idents(Double_KO_4)
Idents(Double_KO_4) <- "Condition"
Double_KO_4 <- RenameIdents(Double_KO_4, `Double_KO`="PGRN_MetTK_KO")
Double_KO_4$Condition <- Idents(Double_KO_4)

Axl_KO_1 <- readRDS(file = "Axl_KO_1_singlets_PCA.rds")
Axl_KO_2 <- readRDS(file = "Axl_KO_2_singlets_PCA.rds")
Axl_KO_3 <- readRDS(file = "Axl_KO_3_singlets_PCA.rds")
Axl_KO_4 <- readRDS(file = "Axl_KO_4_singlets_PCA.rds")

PGRN_Axl_KO_1 <- readRDS(file = "PGRN_Axl_KO_1_singlets_PCA.rds")
PGRN_Axl_KO_2 <- readRDS(file = "PGRN_Axl_KO_2_singlets_PCA.rds")
PGRN_Axl_KO_3 <- readRDS(file = "PGRN_Axl_KO_3_singlets_PCA.rds")
PGRN_Axl_KO_4 <- readRDS(file = "PGRN_Axl_KO_4_singlets_PCA.rds")

setwd("/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/integration_with_Axl")
Ctrl <- c(Ctrl_1, Ctrl_2,Ctrl_3,Ctrl_4)
anchors_Ctrl <- FindIntegrationAnchors(object.list = Ctrl, dims = 1:30)
Ctrl_integrated <- IntegrateData(anchorset = anchors_Ctrl, dims = 1:30)
rm(Ctrl_1, Ctrl_2,Ctrl_3,Ctrl_4, Ctrl)

PGRN_KO <- c(PGRN_KO_1, PGRN_KO_2,PGRN_KO_3,PGRN_KO_4)
anchors_PGRN_KO <- FindIntegrationAnchors(object.list = PGRN_KO, dims = 1:30)
PGRN_KO_integrated <- IntegrateData(anchorset = anchors_PGRN_KO, dims = 1:30)
rm(PGRN_KO_1, PGRN_KO_2,PGRN_KO_3,PGRN_KO_4, PGRN_KO)

MerTK_KO <- c(MerTK_KO_1, MerTK_KO_2,MerTK_KO_3,MerTK_KO_4)
anchors_MerTK_KO <- FindIntegrationAnchors(object.list = MerTK_KO, dims = 1:30)
MerTK_KO_integrated <- IntegrateData(anchorset = anchors_MerTK_KO, dims = 1:30)
rm(MerTK_KO_1, MerTK_KO_2,MerTK_KO_3,MerTK_KO_4, MerTK_KO)

Double_KO <- c(Double_KO_1, Double_KO_2,Double_KO_3,Double_KO_4)
anchors_Double_KO <- FindIntegrationAnchors(object.list = Double_KO, dims = 1:30)
Double_KO_integrated <- IntegrateData(anchorset = anchors_Double_KO, dims = 1:30)
rm(Double_KO_1, Double_KO_2,Double_KO_3,Double_KO_4, Double_KO)

Axl_KO <- c(Axl_KO_1, Axl_KO_2,Axl_KO_3,Axl_KO_4)
anchors_Axl_KO <- FindIntegrationAnchors(object.list = Axl_KO, dims = 1:30)
Axl_KO_integrated <- IntegrateData(anchorset = anchors_Axl_KO, dims = 1:30)
rm(Axl_KO_1, Axl_KO_2,Axl_KO_3,Axl_KO_4, Axl_KO)

PGRN_Axl_KO <- c(PGRN_Axl_KO_1, PGRN_Axl_KO_2,PGRN_Axl_KO_3,PGRN_Axl_KO_4)
anchors_PGRN_Axl_KO <- FindIntegrationAnchors(object.list = PGRN_Axl_KO, dims = 1:30)
PGRN_Axl_KO_integrated <- IntegrateData(anchorset = anchors_PGRN_Axl_KO, dims = 1:30)
rm(PGRN_Axl_KO_1, PGRN_Axl_KO_2,PGRN_Axl_KO_3,PGRN_Axl_KO_4, PGRN_Axl_KO)

Mouse_FTD <- c(Ctrl_integrated, PGRN_KO_integrated,MerTK_KO_integrated, Double_KO_integrated,Axl_KO_integrated, PGRN_Axl_KO_integrated)
anchors_Mouse_FTD <- FindIntegrationAnchors(object.list = Mouse_FTD, dims = 1:30)
Mouse_FTD_integrated <- IntegrateData(anchorset = anchors_Mouse_FTD, dims = 1:30)
rm(Ctrl_integrated, PGRN_KO_integrated,MerTK_KO_integrated, Double_KO_integrated,Axl_KO_integrated, PGRN_Axl_KO_integrated, Mouse_FTD)


pdf("Mouse_FTD_QC.pdf", width=12, height=4)
Idents(Mouse_FTD_integrated) <- "orig.ident"
VlnPlot(object = Mouse_FTD_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
Idents(Mouse_FTD_integrated) <- "Condition"
VlnPlot(object = Mouse_FTD_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

saveRDS(Mouse_FTD_integrated, file = "Mouse_FTD_integrated.rds")

DefaultAssay(Mouse_FTD_integrated) <- 'integrated'

# Mouse_FTD_integrated <- NormalizeData(Mouse_FTD_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
# Mouse_FTD_integrated <- FindVariableFeatures(Mouse_FTD_integrated, selection.method = "vst", nfeatures = 3000)

Mouse_FTD_integrated <- ScaleData(Mouse_FTD_integrated, verbose = FALSE)
Mouse_FTD_integrated <- RunPCA(Mouse_FTD_integrated, features = VariableFeatures(object = Mouse_FTD_integrated), verbose = FALSE)

Mouse_FTD_integrated <- FindNeighbors(Mouse_FTD_integrated, dims = 1:15)
Mouse_FTD_integrated <- FindClusters(Mouse_FTD_integrated, resolution = 0.1)
Mouse_FTD_integrated <- RunUMAP(Mouse_FTD_integrated, dims = 1: 15)

DefaultAssay(Mouse_FTD_integrated) <- 'RNA'
Mouse_FTD_integrated <- NormalizeData(Mouse_FTD_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
Mouse_FTD_integrated <- ScaleData(Mouse_FTD_integrated, features = rownames(Mouse_FTD_integrated))

# Relevel object@ident
#Mouse_FTD_integrated$Condition <- factor(x = Mouse_FTD_integrated$Condition, levels = c("Ctrl", "PGRN_KO","MerTK_KO","Axl_KO","PGRN_MerTK_KO","PGRN_Axl_KO"))
#Mouse_FTD_integrated$Sample_Name <- factor(x = Mouse_FTD_integrated$Sample_Name, levels = c("Ctrl_1","Ctrl_2","Ctrl_3","Ctrl_4",
#                                                                                        "PGRN_KO_1","PGRN_KO_2","PGRN_KO_3","PGRN_KO_4",
#                                                                                            "MerTK_KO_1","MerTK_KO_2","MerTK_KO_3","MerTK_KO_4",
#                                                                                            "Axl_KO_1","Axl_KO_2","Axl_KO_3","Axl_KO_4",
#                                                                                            "PGRN_MerTK_KO_1","PGRN_MerTK_KO_2","PGRN_MerTK_KO_3","PGRN_MerTK_KO_4",
#                                                                                            "PGRN_Axl_KO_1","PGRN_Axl_KO_2","PGRN_Axl_KO_3","PGRN_Axl_KO_4"))


pdf("Mouse_FTD_integrated_umap.pdf", width=5, height=4)
DimPlot(Mouse_FTD_integrated, reduction = 'umap', label = T)
dev.off()
pdf("Mouse_FTD_integrated_umap_split_individual.pdf", width=16, height=6)
DimPlot(Mouse_FTD_integrated, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 8)
dev.off()
pdf("Mouse_FTD_integrated_umap_split_Condition.pdf", width=6.5, height=9)
DimPlot(Mouse_FTD_integrated, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()

saveRDS(Mouse_FTD_integrated, file = 'Mouse_FTD_integrated_PCA_0.1.rds')

DefaultAssay(Mouse_FTD_integrated) <- 'RNA'
pdf("Mouse_FTD_umap_test.pdf", width=4, height=3)
DimPlot(Mouse_FTD_integrated, reduction = 'umap', label = T)
dev.off()

DefaultAssay(Mouse_FTD_integrated) <- 'RNA'

Mouse_FTD_markers <- FindAllMarkers(Mouse_FTD_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1, test.use = "MAST")
write.csv(Mouse_FTD_markers, "Mouse_FTD_markers.csv")

Mouse_FTD_markers <- read.csv(file = "Mouse_FTD_markers.csv", header=T,row.names =1)
top5 <- Mouse_FTD_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
top5$gene <- as.character(top5$gene)
pdf("Mouse_FTD_HeatMapTop5_0.1_new.pdf", width=24, height=16)
DoHeatmap(Mouse_FTD_integrated, features = top5$gene) + NoLegend()
dev.off()

DefaultAssay(Mouse_FTD_integrated) <- 'RNA'
pdf("Mouse_FTD_umap_test.pdf", width=8, height=6)
DimPlot(Mouse_FTD_integrated, reduction = 'umap', label = T)
dev.off()
#Add marker genes

sig_EN<-c("Snap25","Slc17a7", "Nrgn","Gad1", "Gad2","Plp1", "Mbp", "Mobp", "Clu", "Aldoc", "Pla2g7","Cx3cr1", "P2ry12", "Csf1r",
          "Pdgfra", "Vcan","Vtn", "Igfbp7")
markers.to.plot <- as.matrix(sig_EN)
pdf("Mouse_FTD_annotation_combine.pdf", width=10, height=5)
DotPlot(object = Mouse_FTD_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()
dev.off()

