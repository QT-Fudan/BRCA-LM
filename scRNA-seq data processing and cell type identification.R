library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(RCurl)
library(cowplot)
library(sctransform)


#Quality control
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
filtered_seurat <- subset(seurat_object, subset = nFeature_RNA > 250)
dim(filtered_seurat)

#SCT(sctransform): normalization and variance stabilization
split_seurat <- SplitObject(filtered_seurat, split.by = "sample_name")
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=cc.genes$g2m.genes, s.features=cc.genes$s.genes)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c('percent.mt'))
}
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000)
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
seurat_object <- IntegrateData(anchorset = integ_anchors, 
                               normalization.method = "SCT")

#PCA
seurat_object <- RunPCA(object = seurat_object)
ElbowPlot(seurat_object,ndims = 50)
seurat_object <- FindNeighbors(seurat_object, dims = 1:30)
DefaultAssay(seurat_object) <- 'integrated'
seurat_object <- FindClusters(seurat_object, resolution = 0.5) 

#umap
seurat_object <- RunUMAP(seurat_object, 
                         dims = 1:30,
                         reduction = "pca")

library(paletteer)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggsci)
library(scales)

pal <- c(pal_npg('nrc')(10),pal_jco('default')(10),pal_aaas('default')(10),pal_npg('nrc',alpha = 0.5)(10))

#dot plot (UMAP-sample information)
dimplot <- DimPlot(seurat_object,label = F,
                   reduction = 'umap',group.by='rename_clusters',#split.by = 'rename.clusters.all',
                   pt.size = 1,
                   cols =pal,
                   repel = T)
#dot plot (UMAP-gene expression)
markerplot_list <- list()
for(x in 1:length(genes)){
  markerplot <- FeaturePlot(seurat_object,
                            features = genes[x],
                            pt.size = 0.1,
                            cols = c('lightgrey','navy'),#split.by='sample_location',ncol = 3,
                            coord.fixed = T)
  markerplot_list[[x]] <- markerplot
}

#heatmap & bubble plot (marker gene expression)
##top10 marker genes
seurat_object <- SetIdent(seurat_object,value = seurat_object@meta.data$seurat_clusters)
markers <- FindAllMarkers(seurat_object,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25,return.thresh = 0.01)
top10 <- markers%>%group_by(cluster)%>%slice_max(n=10,order_by=avg_log2FC)
write.table(top10,file = 'top10_markergene_seurat_clusters.xlsx',sep='\t')

##heatmap
library(ComplexHeatmap)
library(gplots)
library(Seurat)
library(dplyr)
library(circlize)
top10 <- top10$V1
information <- data.frame(rename.clusters.all=seurat_object@meta.data$rename.clusters.all,
                          rename_clusters=seurat_object@meta.data$rename_clusters)
information <- information[!duplicated(information$rename.clusters.all),]
rownames(information) <- information$rename.clusters.all

set_color <- pal[1:length(levels(factor(information$rename_clusters)))]
names(set_color) <- levels(factor(information$rename_clusters))
set_color3 <- pal[1:length(levels(factor(information$rename.clusters.all)))]
names(set_color3) <- levels(factor(information$rename.clusters.all))	
column_annotation <- HeatmapAnnotation(df = information, 
                                       col = list(rename.clusters.all=set_color3,
                                                  rename_clusters=set_color ))
seurat_object<- SetIdent(seurat_object,value=seurat_object@meta.data$rename.clusters.all) 
mat_t <- AverageExpression(seurat_object,
                           group.by = 'ident')
mat <- mat_t$RNA
mat <- as.matrix(mat[top10,])
mat <- t(scale(t(mat)))
complexheatmap <- Heatmap(mat, 
                          col=colorRamp2(c(-2,0,5),c('navy','white','firebrick3')), 
                          column_names_side = "top", 
                          cluster_columns=F, 
                          border=T,
                          cluster_rows=F, 
                          show_column_names = FALSE, 
                          show_row_names = T,
                          column_split = information$rename_clusters,
                          row_names_gp = gpar(fontsize = 8), 
                          column_names_gp = gpar(fontsize = 8), 
                          top_annotation = column_annotation)
##bubble plot
marker_gene_DotPlot <- DotPlot(seurat_object, features = top10,cols = c('lightgrey','navy'))+RotatedAxis()

