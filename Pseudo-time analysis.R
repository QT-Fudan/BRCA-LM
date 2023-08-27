library(Seurat)
library(monocle)
library(ggplot2)
library(ggpubr)
library(dplyr)

data = as(as.matrix(seurat_object@assays$RNA@counts),'sparseMatrix')
pd = new('AnnotatedDataFrame',data = seurat_object@meta.data)
fData = data.frame(gene_short_name = row.names(data),row.names = row.names(data))
fd = new('AnnotatedDataFrame',data = fData)
monocle_object = newCellDataSet(data,phenoData = pd,featureData = fd, lowerDetectionLimit = 0, expressionFamily=negbinomial.size())

monocle_object <- estimateSizeFactors(monocle_object)
monocle_object <- estimateDispersions(monocle_object)
disp_table <- dispersionTable(monocle_object)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
monocle_object <- setOrderingFilter(monocle_object, unsup_clustering_genes$gene_id)
monocle_object <- reduceDimension(
  monocle_object,
  max_components = 2,
  method = 'DDRTree')
monocle_object <- orderCells(monocle_object)

monocle_object@phenoData@data$sample_code <- monocle_object@phenoData@data$sample_name
monocle_object@phenoData@data$sample_name <- NULL
monocle_object@phenoData@data$rename_clusters <- seurat_object@meta.data$rename_clusters

#plot_cell_trajectory
monocle_object_State <- plot_cell_trajectory(monocle_object,color_by = 'State', cell_size = 0.2)
monocle_object_Pseudotime <- plot_cell_trajectory(monocle_object,color_by = 'Pseudotime', cell_size = 0.2)
monocle_object_seurat_clusters <- plot_cell_trajectory(monocle_object,color_by = 'rename.clusters.all', cell_size = 0.2)
monocle_object_sample_location <- plot_cell_trajectory(monocle_object,color_by = 'sample_location', cell_size = 0.2)

gene <- 'CLEC9A'
markers_fData <- subset(fData(monocle_object), gene_short_name %in% gene)
markers_exprs <- reshape2::melt(as.matrix(exprs(monocle_object[row.names(markers_fData), ])))
markers_exprs$value_log <- log10(markers_exprs$value+0.1)
midpoint <- (max(markers_exprs$value_log)-min(markers_exprs$value_log))/2+min(markers_exprs$value_log)
CLEC9A<- plot_cell_trajectory(monocle_object,markers = 'CLEC9A',use_color_gradient = T,cell_size = 0.2)+scale_color_gradient2(low = 'navy',mid = 'white',high = 'firebrick3')#+facet_wrap(~sample_location,nrow=1)
