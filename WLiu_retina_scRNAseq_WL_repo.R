setwd("/Users/wei/Genomics/Mouse_CellRangeOut/Analysis_Wei/")
options(stringsAsFactors = FALSE)
Sys.getenv('R_MAX_VSIZE')

library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(clusterProfiler)
library(patchwork)
library(KEGG.db)
library(DOSE)
library(org.Mm.eg.db)
library(monocle3)

#KO
EC1.data <- Read10X(data.dir = "/Users/wei/Genomics/Mouse_CellRangeOut/EC1_filtered_feature_bc_matrix/")
EC1 <- CreateSeuratObject(counts = EC1.data, project = "EC1")
EC1##12672 cells
EC1$condition <- "KO"
EC1[["percent.mt"]] <- PercentageFeatureSet(EC1, pattern = "^mt-")
EC1 <- RenameCells(EC1, add.cell.id = "KO")


##Control
EC7.data <- Read10X(data.dir = "/Users/wei/Genomics/Mouse_CellRangeOut/EC7_filtered_feature_bc_matrix/")
EC7 <- CreateSeuratObject(counts = EC7.data, project = "EC7")
EC7##10315 cells
EC7$condition <- "Control"
EC7[["percent.mt"]] <- PercentageFeatureSet(EC7, pattern = "^mt-")
EC7 <- RenameCells(EC7, add.cell.id = "Control")

##QC
plot1 <- FeatureScatter(EC1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(EC1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("KO_QC.pdf", width = 12, height = 6)
plot1 + plot2
dev.off()


plot1 <- FeatureScatter(EC7, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(EC7, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("Control_QC.pdf", width = 12, height = 6)
plot1 + plot2
dev.off()

##filter the data
EC1 <- subset(EC1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
length(rownames(EC1@meta.data))

EC7 <- subset(EC7, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
length(rownames(EC7@meta.data))

##standard stimulated vs control intergration method
EC1 <- NormalizeData(EC1)
EC1 <- FindVariableFeatures(EC1, selection.method = "vst", nfeatures = 2000)
EC7 <- NormalizeData(EC7)
EC7 <- FindVariableFeatures(EC7, selection.method = "vst", nfeatures = 2000)


std.anchors <- FindIntegrationAnchors(object.list = c(EC1,EC7), dims = 1:20)
std.combined <- IntegrateData(anchorset = std.anchors, dims = 1:20)



# Run the standard workflow for visualization and clustering
std.combined <- ScaleData(std.combined, verbose = FALSE)
std.combined <- RunPCA(std.combined, npcs = 30, verbose = FALSE)
ElbowPlot(std.combined)

# UMAP and Clustering
std.combined <- RunUMAP(std.combined, reduction = "pca", dims = 1:20)
std.combined <- FindNeighbors(std.combined, reduction = "pca", dims = 1:20)
std.combined <- FindClusters(std.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(std.combined, reduction = "umap", group.by = "orig.ident") 
p2 <- DimPlot(std.combined, reduction = "umap", label = TRUE)

p1

pdf("Control_KO_intergrated_clustering.pdf", width = 12, height = 6)
plot_grid(p1, p2)
dev.off()

pdf("Control_KO_intergrated_clustering_seperated_plot.pdf", width = 12, height = 6)
DimPlot(std.combined, reduction = "umap", split.by = "orig.ident")
dev.off()

##cluster analysis 
table(std.combined$orig.ident, std.combined$seurat_clusters)

##fraction  plot
stat.percentage = round(prop.table(table(std.combined$orig.ident, std.combined$seurat_clusters), margin = 1), 3)
stack_cluster = data.frame(t(as.data.frame.matrix(stat.percentage)))

write.csv(stat.percentage,"stat.percentage.csv")

pdf("Control_vs_KO_fraction_in_clusters.pdf", width = 8, height = 8)
print(ggplot(stack_cluster, aes(x=EC1, y=EC7)) +geom_point() + geom_text(label=rownames(stack_cluster),position = position_nudge(x = 0.005,y = 0.005))+ labs(y="Fraction of Control cells per cluster", x = "Fraction of KO cells per cluster"))
dev.off()

## differential list 
std.combined.ALLmarkers <- FindAllMarkers(std.combined, only.pos = F,verbose = TRUE)
write.csv(std.combined.ALLmarkers,"intergrated_cluster_markers.csv")






##Generate cluster associated heatmap
top10 <- std.combined.ALLmarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

pdf("std.combined_heatmap_with_TopClusterMarkers.pdf", width = 12, height = 18)
DoHeatmap(std.combined, features = top10$gene) + NoLegend()
dev.off()




 
  