library(Seurat)
library(CellChat)

DefaultAssay(seurat_object) <- 'RNA'
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object_copy <- seurat_object
seurat_object_copy <- SetIdent(seurat_object_copy,value=seurat_object_copy$sample_location)

value <- 'liver-metastasis'#adjacent tissue,liver-metastasis
seurat_object <- subset(seurat_object_copy,idents=value)
cellchat <- createCellChat(object=seurat_object@assays$RNA@data,meta=seurat_object@meta.data,group.by = 'rename.clusters.all')
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human  
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat,PPI.human)
cellchat <- computeCommunProb(cellchat,raw.use = F,population.size = T)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

#merge cellchat
object.list <- list(Adj = Adj, 
                    LM = LM)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#signaling changes in each cell type
celltype <- levels(LM@idents)
plot <- list()
for(x in 1:length(celltype)){
  gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = celltype[x])
  plot[[x]] <- gg1
}

#cellchat
LM.net <- subsetCommunication(cellchat)
LM.net$sample_location <- rep('LM',nrow(LM.net))
LM.netp <- subsetCommunication(cellchat,slot.name='netP')
LM.netp$sample_location <- rep('LM',nrow(LM.netp))

ADJ.net <- subsetCommunication(cellchat)
ADJ.net$sample_location <- rep('ADJ',nrow(ADJ.net))
ADJ.netp <- subsetCommunication(cellchat,slot.name='netP')
ADJ.netp$sample_location <- rep('ADJ',nrow(ADJ.netp))

net <- merge(LM.net,ADJ.net,all = T)
netp <- merge(LM.netp,ADJ.netp,all = T)
##heatmap(pathway)
library(dplyr)
library(ggplot2)
df <- netp
unique(df$source)
selected_source <- c('Activated DC','cDC1','cDC2','FCN3+ Macrophages','Germinal center B cells',
                     'Monocytes','Naive B cells','TREM2+ Macrophages')
selected_target <- c('CD4+ T cells','CD8+ cytotoxic T cells','CD8+ NK T cells','CD8+MHCII+ T cells')
selected_pathway <- LM@netP$pathways
selected_source <- c('Activated DC','cDC1','FCN3+ Macrophages')
selected_target <- c('CD8+ cytotoxic T cells','CD8+ NK T cells','CD8+MHCII+ T cells')

df_selected <- df[which(df$source %in% selected_source),]
df_selected <- df_selected[which(df_selected$target %in% selected_target),]
df_selected <- df_selected[which(df_selected$pathway_name %in% selected_pathway),]
df_selected$interaction <- paste(df_selected$source,'|',df_selected$target,sep = '')
df_selected$log.prob <- log10(df_selected$prob)
midpoint <- (max(df_selected$log.prob)+min(df_selected$log.prob))/2
heatmap <- ggplot(df_selected,aes(x = interaction,y = reorder(pathway_name,log.prob),fill=log.prob))+
  geom_tile(color = 'white')+
  scale_fill_gradient2(low = 'navy',mid='white',high = 'firebrick3',midpoint = midpoint)+
  facet_wrap(vars(sample_location))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -90,vjust = 0,hjust = 0))
##bubble plot(LR-pair)
df <- net
selected_pathway <- c('MHC-II')
selected_source <- c('Activated DC','cDC1','FCN3+ Macrophages','Germinal center B cells','Naive B cells')
selected_target <- c('CD4+ T cells')

df_selected <- df[which(df$source %in% selected_source),]
df_selected <- df_selected[which(df_selected$target %in% selected_target),]
df_selected <- df_selected[which(df_selected$pathway_name %in% selected_pathway),]
df_selected$interaction <- paste(df_selected$source,'|',df_selected$target,sep = '')
df_selected$log.prob <- log10(df_selected$prob)

midpoint <- (max(df_selected$log.prob)+min(df_selected$log.prob))/2
bubble <- ggplot(df_selected,aes(x=interaction,y=interaction_name))+
  geom_point(aes(size = -log10(pval+0.001),colour = log.prob))+
  facet_wrap(vars(sample_location))+
  scale_size_continuous(range = c(1,3))+
  scale_color_gradient2(high="firebrick3",mid = "white",low ="navy",midpoint = midpoint)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -90,vjust = 0,hjust = 0),
        text = element_text(size = 10))
#heatmap(gene expression)
library(Seurat)
library(reshape2)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggsci)
df <- net
selected_pathway <- c('MHC-II','MHC-I')
selected_source <- c('Activated DC','cDC1','FCN3+ Macrophages','Germinal center B cells','Naive B cells')
selected_target <- c('CD4+ T cells','CD8+ cytotoxic T cells','CD8+ NK T cells','CD8+MHCII+ T cells')
df_selected <- df[which(df$source %in% selected_source),]
df_selected <- df_selected[which(df_selected$target %in% selected_target),]
df_selected <- df_selected[which(df_selected$pathway_name %in% selected_pathway),]
genes <- unique(df_selected$ligand)

seurat_object_copy <- SetIdent(seurat_object_copy,value = seurat_object_copy$rename.clusters.all)

selected_celltype <- 'FCN3+ Macrophages'
seurat_object <- subset(seurat_object_copy,idents=selected_celltype)
seurat_object <- SetIdent(seurat_object,value = seurat_object$sample_name)
gene_expr <- AverageExpression(seurat_object)
gene_expr <- gene_expr$RNA
gene_expr_df <- as.data.frame(gene_expr)
gene_expr_df <- gene_expr_df[genes,]
gene_expr_df <- t(scale(t(gene_expr_df)))
gene_expr_df <- as.data.frame(gene_expr_df)
gene_expr_df$gene <- row.names(gene_expr_df)
gene_expr_df <- melt(gene_expr_df)
colnames(gene_expr_df) <- c('gene','sample_name','value')
gene_expr_df$sample_location <- ''
for(x in 1:length(gene_expr_df$sample_name)){
  gene_expr_df$sample_location[x] <- substring(gene_expr_df$sample_name[x],9,9)
}
gene_expr_df$paired_group <- ''
for(x in 1:length(gene_expr_df$sample_name)){
  gene_expr_df$paired_group[x] <- substring(gene_expr_df$sample_name[x],1,7)
}

heatmap <- ggplot(gene_expr_df,aes(x = sample_name,y = gene,fill=value))+
  geom_tile(color = 'white')+
  scale_fill_gradient2(low = 'navy',mid='white',high = 'firebrick3',midpoint = 0)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -90,vjust = 0,hjust = 0))+
  labs(title = selected_celltype)