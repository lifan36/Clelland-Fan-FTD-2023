#set working directory ====
setwd("/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/DF_2ndRound")
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(SoupX)

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/cellranger/LG180_1601/outs')
sc = autoEstCont(sc)
Ctrl_1.counts = adjustCounts(sc)
Ctrl_1 <- CreateSeuratObject(counts = Ctrl_1.counts, project = "LG180_1601", min.cells = 3, min.features = 200)
Ctrl_1[["Genotype"]] = c('C57BL/6J')
Ctrl_1[["Sample_Name"]] = c('Ctrl_1')
Ctrl_1[["Condition"]] = c('Ctrl')
Ctrl_1[["Sex"]] = c('M')
rm(Ctrl_1.counts)
#vizualize QC metrics and filtering====
Ctrl_1[["percent.mt"]] <- PercentageFeatureSet(object = Ctrl_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- Ctrl_1

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Ctrl_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Ctrl_1_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Ctrl_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "Ctrl_1_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.076*10594) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_805", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Ctrl_1_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Ctrl_1_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_805" #visualizing the singlet vs doublet cells
pdf("Ctrl_1_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Ctrl_1_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Ctrl_1_singlets.rds")
singlets<-readRDS("Ctrl_1_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Ctrl_1_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Ctrl_1_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Ctrl_1_singlets_PCA.rds")
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/cellranger/LG180_1603/outs')
sc = autoEstCont(sc)
Ctrl_2.counts = adjustCounts(sc)
Ctrl_2 <- CreateSeuratObject(counts = Ctrl_2.counts, project = "LG180_1603", min.cells = 3, min.features = 200)
Ctrl_2[["Genotype"]] = c('C57BL/6J')
Ctrl_2[["Sample_Name"]] = c('Ctrl_2')
Ctrl_2[["Condition"]] = c('Ctrl')
Ctrl_2[["Sex"]] = c('M')
rm(Ctrl_2.counts)
#vizualize QC metrics and filtering====
Ctrl_2[["percent.mt"]] <- PercentageFeatureSet(object = Ctrl_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- Ctrl_2

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Ctrl_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Ctrl_2_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Ctrl_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "Ctrl_2_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.069*9823) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.02, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.02, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.02_678", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Ctrl_2_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Ctrl_2_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.02_678" #visualizing the singlet vs doublet cells
pdf("Ctrl_2_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Ctrl_2_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Ctrl_2_singlets.rds")
singlets<-readRDS("Ctrl_2_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Ctrl_2_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Ctrl_2_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Ctrl_2_singlets_PCA.rds")
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/cellranger/LG180_1653/outs')
sc = autoEstCont(sc)
Ctrl_3.counts = adjustCounts(sc)
Ctrl_3 <- CreateSeuratObject(counts = Ctrl_3.counts, project = "LG180_1653", min.cells = 3, min.features = 200)
Ctrl_3[["Genotype"]] = c('C57BL/6J')
Ctrl_3[["Sample_Name"]] = c('Ctrl_3')
Ctrl_3[["Condition"]] = c('Ctrl')
Ctrl_3[["Sex"]] = c('F')
rm(Ctrl_3.counts)
#vizualize QC metrics and filtering====
Ctrl_3[["percent.mt"]] <- PercentageFeatureSet(object = Ctrl_3, pattern = "^mt-") #recognize mitochondrial transcripts
all <- Ctrl_3

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Ctrl_3_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Ctrl_3_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Ctrl_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "Ctrl_3_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7589) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_410", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Ctrl_3_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Ctrl_3_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_410" #visualizing the singlet vs doublet cells
pdf("Ctrl_3_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Ctrl_3_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Ctrl_3_singlets.rds")
singlets<-readRDS("Ctrl_3_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Ctrl_3_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Ctrl_3_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Ctrl_3_singlets_PCA.rds")
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/cellranger/LG180_1655/outs')
sc = autoEstCont(sc)
Ctrl_4.counts = adjustCounts(sc)
Ctrl_4 <- CreateSeuratObject(counts = Ctrl_4.counts, project = "LG180_1655", min.cells = 3, min.features = 200)
Ctrl_4[["Genotype"]] = c('C57BL/6J')
Ctrl_4[["Sample_Name"]] = c('Ctrl_4')
Ctrl_4[["Condition"]] = c('Ctrl')
Ctrl_4[["Sex"]] = c('F')
rm(Ctrl_4.counts)
#vizualize QC metrics and filtering====
Ctrl_4[["percent.mt"]] <- PercentageFeatureSet(object = Ctrl_4, pattern = "^mt-") #recognize mitochondrial transcripts
all <- Ctrl_4

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Ctrl_4_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Ctrl_4_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Ctrl_4_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "Ctrl_4_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8164) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_498", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Ctrl_4_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Ctrl_4_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_498" #visualizing the singlet vs doublet cells
pdf("Ctrl_4_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Ctrl_4_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Ctrl_4_singlets.rds")
singlets<-readRDS("Ctrl_4_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Ctrl_4_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Ctrl_4_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Ctrl_4_singlets_PCA.rds")
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/cellranger/LG332_174/outs')
sc = autoEstCont(sc)
PGRN_KO_1.counts = adjustCounts(sc)
PGRN_KO_1 <- CreateSeuratObject(counts = PGRN_KO_1.counts, project = "LG332_174", min.cells = 3, min.features = 200)
PGRN_KO_1[["Genotype"]] = c('PGRN KO: -/-')
PGRN_KO_1[["Sample_Name"]] = c('PGRN_KO_1')
PGRN_KO_1[["Condition"]] = c('PGRN_KO')
PGRN_KO_1[["Sex"]] = c('M')
rm(PGRN_KO_1.counts)
#vizualize QC metrics and filtering====
PGRN_KO_1[["percent.mt"]] <- PercentageFeatureSet(object = PGRN_KO_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- PGRN_KO_1

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("PGRN_KO_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("PGRN_KO_1_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("PGRN_KO_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "PGRN_KO_1_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8411) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_513", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("PGRN_KO_1_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("PGRN_KO_1_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_513" #visualizing the singlet vs doublet cells
pdf("PGRN_KO_1_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"PGRN_KO_1_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"PGRN_KO_1_singlets.rds")
singlets<-readRDS("PGRN_KO_1_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("PGRN_KO_1_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("PGRN_KO_1_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"PGRN_KO_1_singlets_PCA.rds")
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/cellranger/LG332_218/outs')
sc = autoEstCont(sc)
PGRN_KO_2.counts = adjustCounts(sc)
PGRN_KO_2 <- CreateSeuratObject(counts = PGRN_KO_2.counts, project = "LG332_218", min.cells = 3, min.features = 200)
PGRN_KO_2[["Genotype"]] = c('PGRN KO: -/-')
PGRN_KO_2[["Sample_Name"]] = c('PGRN_KO_2')
PGRN_KO_2[["Condition"]] = c('PGRN_KO')
PGRN_KO_2[["Sex"]] = c('F')
rm(PGRN_KO_2.counts)
#vizualize QC metrics and filtering====
PGRN_KO_2[["percent.mt"]] <- PercentageFeatureSet(object = PGRN_KO_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- PGRN_KO_2

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("PGRN_KO_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("PGRN_KO_2_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("PGRN_KO_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "PGRN_KO_2_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8118) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_495", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("PGRN_KO_2_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("PGRN_KO_2_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_495" #visualizing the singlet vs doublet cells
pdf("PGRN_KO_2_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"PGRN_KO_2_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"PGRN_KO_2_singlets.rds")
singlets<-readRDS("PGRN_KO_2_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("PGRN_KO_2_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("PGRN_KO_2_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"PGRN_KO_2_singlets_PCA.rds")
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/cellranger/LG332_228/outs')
sc = autoEstCont(sc)
PGRN_KO_3.counts = adjustCounts(sc)
PGRN_KO_3 <- CreateSeuratObject(counts = PGRN_KO_3.counts, project = "LG332_228", min.cells = 3, min.features = 200)
PGRN_KO_3[["Genotype"]] = c('PGRN KO: -/-')
PGRN_KO_3[["Sample_Name"]] = c('PGRN_KO_3')
PGRN_KO_3[["Condition"]] = c('PGRN_KO')
PGRN_KO_3[["Sex"]] = c('M')
rm(PGRN_KO_3.counts)
#vizualize QC metrics and filtering====
PGRN_KO_3[["percent.mt"]] <- PercentageFeatureSet(object = PGRN_KO_3, pattern = "^mt-") #recognize mitochondrial transcripts
all <- PGRN_KO_3

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("PGRN_KO_3_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("PGRN_KO_3_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("PGRN_KO_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "PGRN_KO_3_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8529) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.29, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.29, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.29_520", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("PGRN_KO_3_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("PGRN_KO_3_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.29_520" #visualizing the singlet vs doublet cells
pdf("PGRN_KO_3_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"PGRN_KO_3_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"PGRN_KO_3_singlets.rds")
singlets<-readRDS("PGRN_KO_3_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("PGRN_KO_3_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("PGRN_KO_3_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"PGRN_KO_3_singlets_PCA.rds")
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/cellranger/LG332_239/outs')
sc = autoEstCont(sc)
PGRN_KO_4.counts = adjustCounts(sc)
PGRN_KO_4 <- CreateSeuratObject(counts = PGRN_KO_4.counts, project = "LG332_239", min.cells = 3, min.features = 200)
PGRN_KO_4[["Genotype"]] = c('PGRN KO: -/-')
PGRN_KO_4[["Sample_Name"]] = c('PGRN_KO_4')
PGRN_KO_4[["Condition"]] = c('PGRN_KO')
PGRN_KO_4[["Sex"]] = c('F')
rm(PGRN_KO_4.counts)
#vizualize QC metrics and filtering====
PGRN_KO_4[["percent.mt"]] <- PercentageFeatureSet(object = PGRN_KO_4, pattern = "^mt-") #recognize mitochondrial transcripts
all <- PGRN_KO_4

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("PGRN_KO_4_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("PGRN_KO_4_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("PGRN_KO_4_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "PGRN_KO_4_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.069*9133) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_630", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("PGRN_KO_4_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("PGRN_KO_4_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_630" #visualizing the singlet vs doublet cells
pdf("PGRN_KO_4_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"PGRN_KO_4_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"PGRN_KO_4_singlets.rds")
singlets<-readRDS("PGRN_KO_4_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("PGRN_KO_4_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("PGRN_KO_4_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"PGRN_KO_4_singlets_PCA.rds")
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/cellranger/LG324_16/outs')
sc = autoEstCont(sc)
MerTK_KO_1.counts = adjustCounts(sc)
MerTK_KO_1 <- CreateSeuratObject(counts = MerTK_KO_1.counts, project = "LG324_16", min.cells = 3, min.features = 200)
MerTK_KO_1[["Genotype"]] = c('MerTK KO: -/-')
MerTK_KO_1[["Sample_Name"]] = c('MerTK_KO_1')
MerTK_KO_1[["Condition"]] = c('MerTK_KO')
MerTK_KO_1[["Sex"]] = c('M')
rm(MerTK_KO_1.counts)
#vizualize QC metrics and filtering====
MerTK_KO_1[["percent.mt"]] <- PercentageFeatureSet(object = MerTK_KO_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- MerTK_KO_1

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("MerTK_KO_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("MerTK_KO_1_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("MerTK_KO_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "MerTK_KO_1_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.076*10195) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.23, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.23, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.23_775", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("MerTK_KO_1_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("MerTK_KO_1_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.23_775" #visualizing the singlet vs doublet cells
pdf("MerTK_KO_1_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"MerTK_KO_1_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"MerTK_KO_1_singlets.rds")
singlets<-readRDS("MerTK_KO_1_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("MerTK_KO_1_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("MerTK_KO_1_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"MerTK_KO_1_singlets_PCA.rds")
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/cellranger/LG324_17/outs')
sc = autoEstCont(sc)
MerTK_KO_2.counts = adjustCounts(sc)
MerTK_KO_2 <- CreateSeuratObject(counts = MerTK_KO_2.counts, project = "LG324_17", min.cells = 3, min.features = 200)
MerTK_KO_2[["Genotype"]] = c('MerTK KO: -/-')
MerTK_KO_2[["Sample_Name"]] = c('MerTK_KO_2')
MerTK_KO_2[["Condition"]] = c('MerTK_KO')
MerTK_KO_2[["Sex"]] = c('M')
rm(MerTK_KO_2.counts)
#vizualize QC metrics and filtering====
MerTK_KO_2[["percent.mt"]] <- PercentageFeatureSet(object = MerTK_KO_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- MerTK_KO_2

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("MerTK_KO_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("MerTK_KO_2_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("MerTK_KO_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "MerTK_KO_2_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7070) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.04, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.04, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.04_382", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("MerTK_KO_2_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("MerTK_KO_2_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.04_382" #visualizing the singlet vs doublet cells
pdf("MerTK_KO_2_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"MerTK_KO_2_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"MerTK_KO_2_singlets.rds")
singlets<-readRDS("MerTK_KO_2_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("MerTK_KO_2_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("MerTK_KO_2_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"MerTK_KO_2_singlets_PCA.rds")
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/cellranger/LG324_18/outs')
sc = autoEstCont(sc)
MerTK_KO_3.counts = adjustCounts(sc)
MerTK_KO_3 <- CreateSeuratObject(counts = MerTK_KO_3.counts, project = "LG324_18", min.cells = 3, min.features = 200)
MerTK_KO_3[["Genotype"]] = c('MerTK KO: -/-')
MerTK_KO_3[["Sample_Name"]] = c('MerTK_KO_3')
MerTK_KO_3[["Condition"]] = c('MerTK_KO')
MerTK_KO_3[["Sex"]] = c('F')
rm(MerTK_KO_3.counts)
#vizualize QC metrics and filtering====
MerTK_KO_3[["percent.mt"]] <- PercentageFeatureSet(object = MerTK_KO_3, pattern = "^mt-") #recognize mitochondrial transcripts
all <- MerTK_KO_3

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("MerTK_KO_3_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("MerTK_KO_3_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("MerTK_KO_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "MerTK_KO_3_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8861) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.22, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.22, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.22_541", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("MerTK_KO_3_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("MerTK_KO_3_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.22_541" #visualizing the singlet vs doublet cells
pdf("MerTK_KO_3_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"MerTK_KO_3_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"MerTK_KO_3_singlets.rds")
singlets<-readRDS("MerTK_KO_3_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("MerTK_KO_3_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("MerTK_KO_3_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"MerTK_KO_3_singlets_PCA.rds")
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/cellranger/LG324_20/outs')
sc = autoEstCont(sc)
MerTK_KO_4.counts = adjustCounts(sc)
MerTK_KO_4 <- CreateSeuratObject(counts = MerTK_KO_4.counts, project = "LG324_20", min.cells = 3, min.features = 200)
MerTK_KO_4[["Genotype"]] = c('MerTK KO: -/-')
MerTK_KO_4[["Sample_Name"]] = c('MerTK_KO_4')
MerTK_KO_4[["Condition"]] = c('MerTK_KO')
MerTK_KO_4[["Sex"]] = c('F')
rm(MerTK_KO_4.counts)
#vizualize QC metrics and filtering====
MerTK_KO_4[["percent.mt"]] <- PercentageFeatureSet(object = MerTK_KO_4, pattern = "^mt-") #recognize mitochondrial transcripts
all <- MerTK_KO_4

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("MerTK_KO_4_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("MerTK_KO_4_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("MerTK_KO_4_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "MerTK_KO_4_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.076*10887) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.2, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.2, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.2_827", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("MerTK_KO_4_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("MerTK_KO_4_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.2_827" #visualizing the singlet vs doublet cells
pdf("MerTK_KO_4_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"MerTK_KO_4_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"MerTK_KO_4_singlets.rds")
singlets<-readRDS("MerTK_KO_4_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("MerTK_KO_4_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("MerTK_KO_4_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"MerTK_KO_4_singlets_PCA.rds")
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/cellranger/LG332_115/outs')
sc = autoEstCont(sc)
Double_KO_1.counts = adjustCounts(sc)
Double_KO_1 <- CreateSeuratObject(counts = Double_KO_1.counts, project = "LG332_115", min.cells = 3, min.features = 200)
Double_KO_1[["Genotype"]] = c('PGRN KO: -/-; MerTK KO: -/-')
Double_KO_1[["Sample_Name"]] = c('Double_KO_1')
Double_KO_1[["Condition"]] = c('Double_KO')
Double_KO_1[["Sex"]] = c('M')
rm(Double_KO_1.counts)
#vizualize QC metrics and filtering====
Double_KO_1[["percent.mt"]] <- PercentageFeatureSet(object = Double_KO_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- Double_KO_1

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Double_KO_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Double_KO_1_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Double_KO_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "Double_KO_1_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*6903) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_318", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Double_KO_1_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Double_KO_1_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_318" #visualizing the singlet vs doublet cells
pdf("Double_KO_1_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Double_KO_1_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Double_KO_1_singlets.rds")
singlets<-readRDS("Double_KO_1_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Double_KO_1_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Double_KO_1_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Double_KO_1_singlets_PCA.rds")
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/cellranger/LG332_266/outs')
sc = autoEstCont(sc)
Double_KO_2.counts = adjustCounts(sc)
Double_KO_2 <- CreateSeuratObject(counts = Double_KO_2.counts, project = "LG332_266", min.cells = 3, min.features = 200)
Double_KO_2[["Genotype"]] = c('PGRN KO: -/-; MerTK KO: -/-')
Double_KO_2[["Sample_Name"]] = c('Double_KO_2')
Double_KO_2[["Condition"]] = c('Double_KO')
Double_KO_2[["Sex"]] = c('M')
rm(Double_KO_2.counts)
#vizualize QC metrics and filtering====
Double_KO_2[["percent.mt"]] <- PercentageFeatureSet(object = Double_KO_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- Double_KO_2

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Double_KO_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Double_KO_2_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Double_KO_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "Double_KO_2_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.069*9882) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.05, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.05, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.05_682", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Double_KO_2_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Double_KO_2_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.05_682" #visualizing the singlet vs doublet cells
pdf("Double_KO_2_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Double_KO_2_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Double_KO_2_singlets.rds")
singlets<-readRDS("Double_KO_2_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Double_KO_2_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Double_KO_2_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Double_KO_2_singlets_PCA.rds")
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/cellranger/LG332_254/outs')
sc = autoEstCont(sc)
Double_KO_3.counts = adjustCounts(sc)
Double_KO_3 <- CreateSeuratObject(counts = Double_KO_3.counts, project = "LG332_254", min.cells = 3, min.features = 200)
Double_KO_3[["Genotype"]] = c('PGRN KO: -/-; MerTK KO: -/-')
Double_KO_3[["Sample_Name"]] = c('Double_KO_3')
Double_KO_3[["Condition"]] = c('Double_KO')
Double_KO_3[["Sex"]] = c('F')
rm(Double_KO_3.counts)
#vizualize QC metrics and filtering====
Double_KO_3[["percent.mt"]] <- PercentageFeatureSet(object = Double_KO_3, pattern = "^mt-") #recognize mitochondrial transcripts
all <- Double_KO_3

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Double_KO_3_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Double_KO_3_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Double_KO_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "Double_KO_3_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8107) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_495", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Double_KO_3_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Double_KO_3_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_495" #visualizing the singlet vs doublet cells
pdf("Double_KO_3_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Double_KO_3_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Double_KO_3_singlets.rds")
singlets<-readRDS("Double_KO_3_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Double_KO_3_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Double_KO_3_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Double_KO_3_singlets_PCA.rds")
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/cellranger/LG332_255/outs')
sc = autoEstCont(sc)
Double_KO_4.counts = adjustCounts(sc)
Double_KO_4 <- CreateSeuratObject(counts = Double_KO_4.counts, project = "LG332_255", min.cells = 3, min.features = 200)
Double_KO_4[["Genotype"]] = c('PGRN KO: -/-; MerTK KO: -/-')
Double_KO_4[["Sample_Name"]] = c('Double_KO_4')
Double_KO_4[["Condition"]] = c('Double_KO')
Double_KO_4[["Sex"]] = c('F')
rm(Double_KO_4.counts)
#vizualize QC metrics and filtering====
Double_KO_4[["percent.mt"]] <- PercentageFeatureSet(object = Double_KO_4, pattern = "^mt-") #recognize mitochondrial transcripts
all <- Double_KO_4

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Double_KO_4_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Double_KO_4_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Double_KO_4_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "Double_KO_4_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7926) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_428", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Double_KO_4_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Double_KO_4_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_428" #visualizing the singlet vs doublet cells
pdf("Double_KO_4_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Double_KO_4_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Double_KO_4_singlets.rds")
singlets<-readRDS("Double_KO_4_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Double_KO_4_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Double_KO_4_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Double_KO_4_singlets_PCA.rds")

#############################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/cellranger/Axl_1/outs')
sc = autoEstCont(sc)
Axl_KO_1.counts = adjustCounts(sc)
Axl_KO_1 <- CreateSeuratObject(counts = Axl_KO_1.counts, project = "Axl_1", min.cells = 3, min.features = 200)
Axl_KO_1[["Genotype"]] = c('Axl KO: -/-')
Axl_KO_1[["Sample_Name"]] = c('Axl_KO_1')
Axl_KO_1[["Condition"]] = c('Axl_KO')
rm(Axl_KO_1.counts)
#vizualize QC metrics and filtering====
Axl_KO_1[["percent.mt"]] <- PercentageFeatureSet(object = Axl_KO_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- Axl_KO_1

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Axl_KO_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Axl_KO_1_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Axl_KO_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "Axl_KO_1_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8296) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.27, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.27, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.27_506", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Axl_KO_1_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Axl_KO_1_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.27_506" #visualizing the singlet vs doublet cells
pdf("Axl_KO_1_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
#saveRDS(all,"Axl_KO_1_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
#saveRDS(singlets,"Axl_KO_1_singlets.rds")
#singlets<-readRDS("Axl_KO_1_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Axl_KO_1_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Axl_KO_1_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Axl_KO_1_singlets_PCA.rds")
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/cellranger/Axl_2/outs')
sc = autoEstCont(sc)
Axl_KO_2.counts = adjustCounts(sc)
Axl_KO_2 <- CreateSeuratObject(counts = Axl_KO_2.counts, project = "Axl_2", min.cells = 3, min.features = 200)
Axl_KO_2[["Genotype"]] = c('Axl KO: -/-')
Axl_KO_2[["Sample_Name"]] = c('Axl_KO_2')
Axl_KO_2[["Condition"]] = c('Axl_KO')
rm(Axl_KO_2.counts)
#vizualize QC metrics and filtering====
Axl_KO_2[["percent.mt"]] <- PercentageFeatureSet(object = Axl_KO_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- Axl_KO_2

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Axl_KO_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Axl_KO_2_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Axl_KO_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "Axl_KO_2_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.069*9031) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.04, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.04, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.04_623", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Axl_KO_2_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Axl_KO_2_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.04_623" #visualizing the singlet vs doublet cells
pdf("Axl_KO_2_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
#saveRDS(all,"Axl_KO_2_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
#saveRDS(singlets,"Axl_KO_2_singlets.rds")
#singlets<-readRDS("Axl_KO_2_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Axl_KO_2_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Axl_KO_2_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Axl_KO_2_singlets_PCA.rds")
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/cellranger/Axl_3/outs')
sc = autoEstCont(sc)
Axl_KO_3.counts = adjustCounts(sc)
Axl_KO_3 <- CreateSeuratObject(counts = Axl_KO_3.counts, project = "Axl_3", min.cells = 3, min.features = 200)
Axl_KO_3[["Genotype"]] = c('Axl KO: -/-')
Axl_KO_3[["Sample_Name"]] = c('Axl_KO_3')
Axl_KO_3[["Condition"]] = c('Axl_KO')
rm(Axl_KO_3.counts)
#vizualize QC metrics and filtering====
Axl_KO_3[["percent.mt"]] <- PercentageFeatureSet(object = Axl_KO_3, pattern = "^mt-") #recognize mitochondrial transcripts
all <- Axl_KO_3

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Axl_KO_3_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Axl_KO_3_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Axl_KO_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "Axl_KO_3_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8704) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.02, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.02, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.02_531", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Axl_KO_3_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Axl_KO_3_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.02_531" #visualizing the singlet vs doublet cells
pdf("Axl_KO_3_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
#saveRDS(all,"Axl_KO_3_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
#saveRDS(singlets,"Axl_KO_3_singlets.rds")
#singlets<-readRDS("Axl_KO_3_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Axl_KO_3_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Axl_KO_3_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Axl_KO_3_singlets_PCA.rds")
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/cellranger/Axl_4/outs')
sc = autoEstCont(sc)
Axl_KO_4.counts = adjustCounts(sc)
Axl_KO_4 <- CreateSeuratObject(counts = Axl_KO_4.counts, project = "Axl_4", min.cells = 3, min.features = 200)
Axl_KO_4[["Genotype"]] = c('Axl KO: -/-')
Axl_KO_4[["Sample_Name"]] = c('Axl_KO_4')
Axl_KO_4[["Condition"]] = c('Axl_KO')
rm(Axl_KO_4.counts)
#vizualize QC metrics and filtering====
Axl_KO_4[["percent.mt"]] <- PercentageFeatureSet(object = Axl_KO_4, pattern = "^mt-") #recognize mitochondrial transcripts
all <- Axl_KO_4

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Axl_KO_4_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Axl_KO_4_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Axl_KO_4_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "Axl_KO_4_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8161) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.21, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.21, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.21_498", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Axl_KO_4_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Axl_KO_4_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.21_498" #visualizing the singlet vs doublet cells
pdf("Axl_KO_4_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
#saveRDS(all,"Axl_KO_4_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
#saveRDS(singlets,"Axl_KO_4_singlets.rds")
#singlets<-readRDS("Axl_KO_4_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Axl_KO_4_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Axl_KO_4_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Axl_KO_4_singlets_PCA.rds")
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/cellranger/LG333_167/outs')
sc = autoEstCont(sc)
PGRN_Axl_KO_1.counts = adjustCounts(sc)
PGRN_Axl_KO_1 <- CreateSeuratObject(counts = PGRN_Axl_KO_1.counts, project = "LG333_167", min.cells = 3, min.features = 200)
PGRN_Axl_KO_1[["Genotype"]] = c('PGRN KO: -/-; Axl KO: -/-')
PGRN_Axl_KO_1[["Sample_Name"]] = c('PGRN_Axl_KO_1')
PGRN_Axl_KO_1[["Condition"]] = c('PGRN_Axl_KO')
rm(PGRN_Axl_KO_1.counts)
#vizualize QC metrics and filtering====
PGRN_Axl_KO_1[["percent.mt"]] <- PercentageFeatureSet(object = PGRN_Axl_KO_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- PGRN_Axl_KO_1

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("PGRN_Axl_KO_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("PGRN_Axl_KO_1_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#PGRN_Axlt finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/PGRN_AxltFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("PGRN_Axl_KO_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "PGRN_Axl_KO_1_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8130) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.04, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.04, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.04_496", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("PGRN_Axl_KO_1_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("PGRN_Axl_KO_1_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.04_496" #visualizing the singlet vs doublet cells
pdf("PGRN_Axl_KO_1_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
#saveRDS(all,"PGRN_Axl_KO_1_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
#saveRDS(singlets,"PGRN_Axl_KO_1_singlets.rds")
#singlets<-readRDS("PGRN_Axl_KO_1_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("PGRN_Axl_KO_1_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("PGRN_Axl_KO_1_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"PGRN_Axl_KO_1_singlets_PCA.rds")
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/cellranger/LG333_173/outs')
sc = autoEstCont(sc)
PGRN_Axl_KO_2.counts = adjustCounts(sc)
PGRN_Axl_KO_2 <- CreateSeuratObject(counts = PGRN_Axl_KO_2.counts, project = "LG333_173", min.cells = 3, min.features = 200)
PGRN_Axl_KO_2[["Genotype"]] = c('PGRN KO: -/-; Axl KO: -/-')
PGRN_Axl_KO_2[["Sample_Name"]] = c('PGRN_Axl_KO_2')
PGRN_Axl_KO_2[["Condition"]] = c('PGRN_Axl_KO')
rm(PGRN_Axl_KO_2.counts)
#vizualize QC metrics and filtering====
PGRN_Axl_KO_2[["percent.mt"]] <- PercentageFeatureSet(object = PGRN_Axl_KO_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- PGRN_Axl_KO_2

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("PGRN_Axl_KO_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("PGRN_Axl_KO_2_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#PGRN_Axlt finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/PGRN_AxltFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("PGRN_Axl_KO_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "PGRN_Axl_KO_2_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7726) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.03, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.03, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.03_417", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("PGRN_Axl_KO_2_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("PGRN_Axl_KO_2_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.03_417" #visualizing the singlet vs doublet cells
pdf("PGRN_Axl_KO_2_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
#saveRDS(all,"PGRN_Axl_KO_2_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
#saveRDS(singlets,"PGRN_Axl_KO_2_singlets.rds")
#singlets<-readRDS("PGRN_Axl_KO_2_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("PGRN_Axl_KO_2_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("PGRN_Axl_KO_2_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"PGRN_Axl_KO_2_singlets_PCA.rds")
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/cellranger/LG333_174/outs')
sc = autoEstCont(sc)
PGRN_Axl_KO_3.counts = adjustCounts(sc)
PGRN_Axl_KO_3 <- CreateSeuratObject(counts = PGRN_Axl_KO_3.counts, project = "LG333_174", min.cells = 3, min.features = 200)
PGRN_Axl_KO_3[["Genotype"]] = c('PGRN KO: -/-; Axl KO: -/-')
PGRN_Axl_KO_3[["Sample_Name"]] = c('PGRN_Axl_KO_3')
PGRN_Axl_KO_3[["Condition"]] = c('PGRN_Axl_KO')
rm(PGRN_Axl_KO_3.counts)
#vizualize QC metrics and filtering====
PGRN_Axl_KO_3[["percent.mt"]] <- PercentageFeatureSet(object = PGRN_Axl_KO_3, pattern = "^mt-") #recognize mitochondrial transcripts
all <- PGRN_Axl_KO_3

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("PGRN_Axl_KO_3_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("PGRN_Axl_KO_3_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#PGRN_Axlt finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/PGRN_AxltFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("PGRN_Axl_KO_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "PGRN_Axl_KO_3_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8099) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.03, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.03, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.03_494", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("PGRN_Axl_KO_3_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("PGRN_Axl_KO_3_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.03_494" #visualizing the singlet vs doublet cells
pdf("PGRN_Axl_KO_3_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
#saveRDS(all,"PGRN_Axl_KO_3_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
#saveRDS(singlets,"PGRN_Axl_KO_3_singlets.rds")
#singlets<-readRDS("PGRN_Axl_KO_3_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("PGRN_Axl_KO_3_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("PGRN_Axl_KO_3_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"PGRN_Axl_KO_3_singlets_PCA.rds")
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Mouse_pgrn_mertk/cellranger/LG333_178/outs')
sc = autoEstCont(sc)
PGRN_Axl_KO_4.counts = adjustCounts(sc)
PGRN_Axl_KO_4 <- CreateSeuratObject(counts = PGRN_Axl_KO_4.counts, project = "LG333_178", min.cells = 3, min.features = 200)
PGRN_Axl_KO_4[["Genotype"]] = c('PGRN KO: -/-; Axl KO: -/-')
PGRN_Axl_KO_4[["Sample_Name"]] = c('PGRN_Axl_KO_4')
PGRN_Axl_KO_4[["Condition"]] = c('PGRN_Axl_KO')
rm(PGRN_Axl_KO_4.counts)
#vizualize QC metrics and filtering====
PGRN_Axl_KO_4[["percent.mt"]] <- PercentageFeatureSet(object = PGRN_Axl_KO_4, pattern = "^mt-") #recognize mitochondrial transcripts
all <- PGRN_Axl_KO_4

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("PGRN_Axl_KO_4_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("PGRN_Axl_KO_4_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#PGRN_Axlt finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/PGRN_AxltFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("PGRN_Axl_KO_4_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "PGRN_Axl_KO_4_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8926) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.03, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.03, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.03_544", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("PGRN_Axl_KO_4_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("PGRN_Axl_KO_4_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.03_544" #visualizing the singlet vs doublet cells
pdf("PGRN_Axl_KO_4_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
#saveRDS(all,"PGRN_Axl_KO_4_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
#saveRDS(singlets,"PGRN_Axl_KO_4_singlets.rds")
#singlets<-readRDS("PGRN_Axl_KO_4_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("PGRN_Axl_KO_4_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("PGRN_Axl_KO_4_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"PGRN_Axl_KO_4_singlets_PCA.rds")






