

library(data.table)
library(dplyr)
library(sqldf)
library(lubridate)
library(RMySQL)
library(tibble)
library(Seurat)  ##version 3.0
library(dplyr)

#download Matrix data https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE125970&format=file&file=GSE125970%5Fraw%5FUMIcounts%2Etxt%2Egz from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125970
time_fread <- system.time(
  test <- fread("./GSE125970_raw_UMIcounts.txt",)
)

paste("数据的大小为：",format(object.size(test),units="auto"))


test= column_to_rownames(test,"GENE")


View(head(mydata_meta@meta.data))


mydata <- CreateSeuratObject(counts = test, min.cells = 3, min.features =300, project = "PLC_scRNAseq")


pbmc =mydata_meta 



#1-QC

pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.2)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used for anything calculated by
# the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))




pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 50)
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


#Scaling the data

all.genes <- rownames(x = pbmc)
pbmc<- ScaleData(object =pbmc, features = all.genes)

#Perform linear dimensional reduction

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(object = pbmc, dims = 1:2, reduction = "pca")


DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)


#Determine the ‘dimensionality’ of the dataset

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)



JackStrawPlot(pbmc, dims = 1:15)

ElbowPlot(pbmc)


#Cluster the cells

pbmc <- FindNeighbors(pbmc, dims = 1:12)
pbmc <- FindClusters(pbmc, resolution = 0.4)

head(Idents(pbmc), 5)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- RunTSNE(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "tsne")






metadata =fread("GSE125970_cell_inf_clulster.txt",)
row.names(metadata) = metadata$UniqueCell_ID

colnames(metadata) [4]= "seurat_clusters"



pbmc <- AddMetaData(object = pbmc, metadata = metadata)


pbmc@meta.data$seurat_clusters =as.factor(pbmc@meta.data$seurat_clusters)
pbmc@meta.data$RNA_snn_res.0.5 =(pbmc@meta.data$seurat_clusters)


pbmc@meta.data$Sample = str_split_fixed(pbmc@meta.data$Sample_ID,"-",2)[,1]


active.ident_meta =  column_to_rownames(metadata[,-c(2:3)],var = "UniqueCell_ID")
Idents(object = pbmc) = active.ident_meta

pdf("PLC in gut.pdf",width = 10,height = 10)

DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "tsne" ,pt.size = 0.5,group.by = "CellType",label = TRUE) 
DimPlot(pbmc, reduction = "umap" ,pt.size = 0.5,group.by = "CellType",label = TRUE) 
DimPlot(pbmc, reduction = "umap" ,pt.size = 0.5,group.by = "Sample_ID",label = TRUE) 
DimPlot(pbmc, reduction = "umap" ,pt.size = 0.5,group.by = "Sample",label = TRUE) 



pdf("PLC hightlight in gut.pdf",width = 10,height = 10)

DimPlot(subset(pbmc, Sample == c("Ileum","Colon")), label=T, reduction = "umap" ,group.by="CellType", 
        cells.highlight=WhichCells(pbmc, idents = c("4")), 
        cols.highlight = c("#54B401"), 
        cols= "lightgrey")


dev.off()


# active.ident_meta$seurat_clusters = as.factor(active.ident_meta$seurat_clusters)
# pbmc@active.ident =  as.factor(active.ident_meta)
# 
# class(mydata_meta@active.ident)
# str(mydata_meta@active.ident)
# 
# view(pbmc@active.ident)


saveRDS(mydata_meta, "mydata_meta.rds")
saveRDS(pbmc, "mydata_meta_pbmc.rds")

save.image("PLC.RData")

dev.off()

# 提取PLC 细胞亚群

pbmc_subPLC_cell <- subset(pbmc, idents = c("4"))

saveRDS(pbmc_subPLC_cell, "pbmc_subPLC_cell.rds")


pbmc_subPLC_cell@meta.data$Sample = str_split_fixed(pbmc_subPLC_cell@meta.data$Sample_ID,"-",2)[,1]

pdf("PLC cell subcell.pdf",width = 10,height = 10)

DimPlot(pbmc_subPLC_cell, reduction = "umap" ,pt.size = 0.5,group.by = "CellType",label = TRUE) 
DimPlot(pbmc_subPLC_cell, reduction = "umap" ,pt.size = 0.5,group.by = "Sample",label = TRUE) 

dev.off()
list =c("Wnt1","Wnt3","Wnt3a","Wnt7a","Wnt8a" , "Wnt8b","Wnt2","Wnt4","Wnt5a","Wnt5b","Wnt6","Wnt7b" ,"Wnt11")

list = str_to_upper(list)

pdf("wnt in PLC1 .pdf",width = 10,height = 20)
a= FeaturePlot(pbmc_subPLC_cell, features = list,split.by   = "Sample",label = TRUE )
print(a)
dev.off()
pdf("wnt in PLC2.pdf",width = 20,height = 10)

b= VlnPlot(pbmc_subPLC_cell, features  = list,group.by = "Sample" )
print(b)
dev.off()



pdf("PLC_colonsubcell.pdf",width = 10,height = 10)

pbmc_subPLC_cell_COLON <- subset(pbmc_subPLC_cell, Sample_ID == c("Colon-2" , "Colon-1"))

DimPlot(pbmc_subPLC_cell_COLON, reduction = "umap" ,pt.size = 0.5,group.by = "CellType",label = TRUE) 
DimPlot(pbmc_subPLC_cell_COLON, reduction = "umap" ,pt.size = 0.5,group.by = "Sample_ID",label = TRUE) 
dev.off()

saveRDS(pbmc_subPLC_cell_COLON, "pbmc_subPLC_cell_COLON.rds")




pdf("colon_ISC_PLC_sub.pdf",width = 10,height = 10)

pbmc_COLON_subPLC_ISCcell <- subset(pbmc, Sample_ID == c("Colon-2" , "Colon-1"),idents = c("4","6"))
pbmc_COLON_subPLC_ISCcell@meta.data$Sample = str_split_fixed(pbmc_COLON_subPLC_ISCcell@meta.data$Sample_ID,"-",2)[,1]

DimPlot(pbmc_COLON_subPLC_ISCcell, reduction = "umap" ,pt.size = 0.5,group.by = "CellType",label = TRUE) 
DimPlot(pbmc_COLON_subPLC_ISCcell, reduction = "umap" ,pt.size = 0.5,group.by = "Sample_ID",label = TRUE) 
DimPlot(pbmc_COLON_subPLC_ISCcell, reduction = "umap" ,pt.size = 0.5,group.by = "Sample",label = TRUE) 

dev.off()

saveRDS(pbmc_COLON_subPLC_ISCcell, "pbmc_COLON_subPLC_ISCcell.rds")


# PLC  GSEA ssGVSA---------
library(GSVA)
library(GSEABase)

genesets <- getGmt("h.all.v7.0.symbols.gmt")




df.data <- GetAssayData(object = pbmc_subPLC_cell, slot = "data")


df.group <- data.frame(umi = names(Idents(pbmc_subPLC_cell)), 
                       cluster = as.character(pbmc_subPLC_cell@meta.data$Sample_ID), 
                       stringsAsFactors = F)
library(stringr)
df.group$sample = str_split_fixed(df.group$cluster,"-",2)[,1]

head(df.group)


gsva
gsvascore <- gsva(data.matrix(df.data), genesets, parallel.sz = 2)


gsvascore[1:5, 1:5]



library(ComplexHeatmap)
pdf("wnt signal in gut pheatmap.pdf",width = 10,height = 16)

ha.t <- HeatmapAnnotation(Cluster = df.group$sample)
Heatmap(as.matrix(gsvascore), 
        show_column_names = F, 
        cluster_rows = T, 
        cluster_columns = T, 
        top_annotation = ha.t, 
        column_split = df.group$sample, 
        row_names_gp = gpar(fontsize = 8), 
        row_names_max_width = max_text_width(rownames(gsvascore), 
                                             gp = gpar(fontsize = 8)))

dev.off()
rownames(sigPathways)[1]


count <- gsvascore["HALLMARK_WNT_BETA_CATENIN_SIGNALING", , drop = FALSE]
count <- as.data.frame(t(count))
colnames(count) <- "geneset"
count$sample <- as.character(df.group$sample)
title.name ="HALLMARK_WNT_BETA_CATENIN_SIGNALING"


count.geneset.group1 <- count$geneset[count$sample == "Ileum" ]
count.geneset.group2 <- count$geneset[count$sample == "Colon"]
count.geneset.group3 <- count$geneset[count$sample == "Rectum"]


ysegment1 <- max(count.geneset.group1)
ysegment2 <- max(count.geneset.group2)
ysegmenz3 <- max(count.geneset.group3)
ysegment.max <- max(ysegment1, ysegment2,ysegment2)



pval <- sigPathways$P.Value[1]

if (pval < 0.001) {
  pval.label = "***"
} else if (pval < 0.005) {
  pval.label = "**"
} else if (pval < 0.05) {
  pval.laben = "*"
} else if (pval >= 0.05) {
  pval.label = "non.sig"
}

blue <- "#619CD6"
green <- "#89C32E"
pdf("wnt signal in gut.pdf",width = 5,height = 8)
library(ggplot2)
p <- ggplot(count, aes(x = sample, y = geneset, fill = sample)) +
  geom_violin() +
  geom_boxplot(width = 0.2,color = "#FFFFFF")+
  scale_fill_manual(values = c("#440054", "#31688E","#35B779")) + # 用自定义颜色填充
  theme_classic() +
  theme(panel.grid = element_blank(), 
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), 
        axis.title.x = element_text(color = "black", size = 20), 
        axis.title.y = element_blank(), 
        axis.text = element_text(color = "black", size = 16), 
        axis.line = element_line(colour = "black", size = 0.6), 
        plot.title = element_text(size = 20, hjust = 0.5)) + 
  ggtitle(title.name) +
  guides(fill = F)
p

dev.off()





