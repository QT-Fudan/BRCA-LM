#DE gene analysis
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(ggpubr)

seurat_object_copy<- SetIdent(seurat_object_copy,value=seurat_object_copy@meta.data$functional_clusters)
seurat_clusters <- levels(seurat_object_copy@meta.data$functional_clusters)
seurat_object <- subset(seurat_object_copy,idents =seurat_clusters[1])
seurat_object <- SetIdent(seurat_object,value=seurat_object@meta.data$sample_location)
results.classified <- FindMarkers(seurat_object,ident.1 = 'liver-metastasis',ident.2 = 'adjacent tissue')

results.classified$threshold <- as.factor(ifelse(results.classified$p_val_adj<=0.05 & abs(results.classified$avg_log2FC)>=0.1,ifelse(results.classified$avg_log2FC>=0.1,'Up','Down'),'Not'))
results.classified$significant <- as.factor(results.classified$p_val_adj<=0.05 & abs(results.classified$avg_log2FC)>=0.1)
results.classified$label <-''
label=row.names(results.classified[results.classified$p_val_adj < 0.05 & results.classified$significant != "FALSE",])
results.classified$label[results.classified$p_val_adj < 0.05 & results.classified$significant != "FALSE"] = as.character(label) 

write.table(results.classified,file = 'results.classified_diff.genes_Neutrophils.xlsx',sep = '\t',quote = F,row.names = T)

label <- read.table('results.classified_diff.genes_CD4.T.cells.txt',sep = '\t',quote = '',header = F)
label$V2 <- label$V1
label <- label[!duplicated(label$V1),]
results.classified$genes <- rownames(results.classified)
results.classified_label <- merge(results.classified,label,by.x = 'genes',by.y = 'V1',all.x = T)

vocano <- ggplot(results.classified_label,aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold))+
  geom_point()+
  scale_color_manual(values = c('navy','gray','firebrick3'))+
  geom_vline(xintercept = c(-0.1,0.1),lty=4,col='gray')+
  geom_hline(yintercept = -log10(0.05),lty=4,col='gray')+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border = element_blank(),axis.line = element_line(colour = "black"))+
  geom_text_repel(aes(label = V2),size = 2.5,max.overlaps = 40)