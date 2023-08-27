rm(list = ls())
#signature score~sample type
library(GSVA)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(clusterProfiler)
library(limma)
library(ggrepel)
library(Seurat)


seurat_object_copy <- SetIdent(seurat_object_copy,value = seurat_object_copy$celltype)
selected_celltype <- 'cDC1'
seurat_object <- subset(seurat_object_copy,idents=selected_celltype)
seurat_object <- SetIdent(seurat_object,value = seurat_object$sample_name)
gene_expr <- AverageExpression(seurat_object)
gene_expr <- gene_expr$RNA

target_gene_1 <- read.table('MHC.txt',sep = '\t',quote = '')
target_genelist <- list()
target_genelist[['MHC']] <- target_gene_1$V1

gsva_result <- gsva(gene_expr,target_genelist,method='zscore')

info <- data.frame(sample_name =seurat_object$sample_name,
                   sample_location=seurat_object$sample_location)
info <- info[!duplicated(info$sample_name),]

info$paired_group <- ''
for(x in 1:length(info$sample_name)){
  info$paired_group[x] <- substring(info$sample_name[x],1,7)
}

gsva_result_info <- cbind(info,t(gsva_result))
paired_group <- gsva_result_info %>% group_by(paired_group) %>% summarise(freq = n()) %>% filter(freq>1) %>% select(paired_group)
gsva_result_info <- gsva_result_info[gsva_result_info$paired_group %in% paired_group$paired_group,]


##paired t test
GS <- 'MHC'
boxplot <- ggpaired(gsva_result_info, x="sample_location", y=GS, id='paired_group',
                    fill = "sample_location", line.color = "gray",
              line.size = 0.4)+ stat_compare_means(method = "t.test",paired = TRUE)+
  scale_fill_manual(values = c("firebrick3", "navy"))+
  labs(title=paste(selected_celltype,GS),x='',y="Signature score")
