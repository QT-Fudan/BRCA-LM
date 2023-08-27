library(Seurat)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(ggsci)
pal <- c(pal_npg('nrc')(10),pal_jco('default')(10),pal_aaas('default')(10))

cell_df <- data.frame(rename_clusters=as.character(seurat_object$rename.clusters.all),
                      sample_name=as.character(seurat_object$sample_name))
cell_df_immune <- cell_df[which(cell_df$rename_clusters!='Epithelial cells'&cell_df$rename_clusters!='Endothelial cells'&cell_df$rename_clusters!='Hepatocytes'&cell_df$rename_clusters!='Mesenchymal stem cells'&cell_df$rename_clusters!='Tumor cells'),]

cell_info <- data.frame(sample_name=as.character(seurat_object$sample_name),
                        sample_location=as.character(seurat_object$sample_location))
cell_info <- cell_info[!duplicated(cell_info$sample_name),]

cell_ratio <- prop.table(table(cell_df_immune$rename_clusters,cell_df_immune$sample_name),margin = 2)
cell_ratio <- as.data.frame(cell_ratio)
cell_ratio_m <- melt(cell_ratio)
colnames(cell_ratio_m) <- c('rename_clusters','sample_name','v1','ratio')
cell_ratio_m <- merge(cell_ratio_m,cell_info,by.x = 'sample_name',by.y = 'sample_name',all.x = T)
for(x in 1:length(cell_ratio_m$sample_name)){
  cell_ratio_m$paired_group[x] <- substring(cell_ratio_m$sample_name[x],1,7)
}

#pair t test
##box plot
boxplot <- ggpaired(cell_ratio_m, x="sample_location", y="ratio", id='paired_group',
              fill = "sample_location", line.color = "gray",
              line.size = 0.4, palette = "jco")+ stat_compare_means(method = "t.test",paired = TRUE)+
  scale_fill_manual(values = c("firebrick3", "navy"))+
  facet_wrap(~rename_clusters,scales = 'free',ncol = 7)+
  ylab('Ratio among total cells')
##bar plot
barplot <- ggplot(cell_ratio_m)+
  geom_bar(aes(x = ratio, y = sample_name, fill = rename_clusters),stat = "identity")+
  scale_fill_manual(values =pal[1:29])+
  #geom_text(aes(x = ratio, y = sample_name,label=rename_clusters),position = position_stack(vjust = 0.5))+
  facet_grid(sample_location~.,scales = 'free')

#propeller
library(speckle)
library(limma)
library(ggplot2)

  ##in total cells
  cell_info <- data.frame(clusters = seurat_object$rename_clusters,
                             samples = seurat_object$sample_name,
                             group = seurat_object$sample_location)
  
  result <- propeller(clusters = cell_info$clusters, sample = cell_info$samples, 
                      group = cell_info$group)
  
  ## in immune cells
  cell_info_immune <- cell_info[which(cell_info$clusters!='Epithelial cells'&cell_info$clusters!='Endothelial cells'&cell_info$clusters!='Hepatocytes'&cell_info$clusters!='Mesenchymal stem cells'&cell_info$clusters!='Tumor cells'),]
  result <- propeller(clusters = cell_info_immune$clusters, sample = cell_info_immune$samples, 
                      group = cell_info_immune$group)
