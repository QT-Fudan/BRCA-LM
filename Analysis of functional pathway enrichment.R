#Functional pathway enrichment(between different cell clusters)(GO)
  library(Seurat)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  seurat_object<- SetIdent(seurat_object,value=seurat_object@meta.data$rename_clusters)
  markers <- FindAllMarkers(seurat_object,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25,return.thresh = 0.01)
  markers <- markers[which(markers$p_val_adj<=0.05),]
  markers_list <- list()
  for (x in unique(markers$cluster)){
    markers_t <- markers[which(markers$cluster == x),]
    markers_list[[x]] <- markers_t$gene
  }
  
  go_list <- list()
  for (x in 1:length(markers_list)){
    gene.df <- bitr(markers_list[[x]], fromType="SYMBOL",
                    toType="ENTREZID", 
                    OrgDb = "org.Hs.eg.db")
    go <- enrichGO(gene = gene.df$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all")
    #go <- enrichKEGG(gene = gene.df$ENTREZID, organism = "hsa")
    go_list[[x]] <- go
  }
  dotplot_go <- list()
  for (x in 1:length(go_list)){
    dotplot <- dotplot(go_list[[x]])
    dotplot_go[[x]] <- dotplot
  }

#Functional pathway enrichment(between different sample locations)(GSEA)
  library(Seurat)
  library(stringr)
  library(clusterProfiler)
  library(ggplot2)
  library(enrichplot)
  library(ggpubr)
  
  seurat_object_copy <- seurat_object
  DefaultAssay(seurat_object_copy) <- 'SCT'
  seurat_object_copy <- SetIdent(seurat_object_copy,value = seurat_object_copy@meta.data$rename_clusters)
  seurat_clusters <- levels(factor(seurat_object_copy@meta.data$rename_clusters))
  
  for(x in 1:length(seurat_clusters)){
    subcluster <- seurat_clusters[x]
    seurat_object <- subset(seurat_object_copy,idents = subcluster)
    seurat_object <- SetIdent(seurat_object,value = seurat_object@meta.data$sample_location)
    markers <- FindMarkers(seurat_object, ident.1 = "liver-metastasis", ident.2 = "adjacent tissue", verbose = FALSE,logfc.threshold = 0,only.pos =F,min.pct = 0)
    
    #GSEA
    gene_list <- markers$avg_log2FC
    names(gene_list) <- rownames(markers)
    gene_list <- sort(gene_list,decreasing = T)
    
    gsea_result_list <- list()
    dotplot_list <- list()
    sig_geneset_list <- list()
    for(x in 1:length(list.files(swd))){
      geneset <- list.files(swd)[x]
      gs <- read.gmt(geneset)
      gsea_result <- GSEA(geneList = gene_list,
                          TERM2GENE = gs,
                          pvalueCutoff = 0.05)
      gsea_result_list[x] <- gsea_result
      a <- gsea_result@result
      a <- a[which(a$p.adjust<0.05),]
      if(nrow(a)!=0){
        sig_geneset_list[x] <- a$ID
      }
    }
  }
  
  ##barplot
  library(clusterProfiler)
  library(dplyr)
  library(stringr)
  library(Hmisc)
  a <- gsea_result_list[[5]]@result
  a <- a[which(a$p.adjust<0.05),]
  a <- arrange(a,a$NES)
  a$ID <- str_replace(a$ID,'GOBP_','')
  a$ID <- str_replace(a$ID,'GOCC_','')
  a$ID <- str_replace(a$ID,'GOMF_','')
  a$ID <- str_replace_all(a$ID,'_',' ')
  a$ID <- tolower(a$ID)
  a$ID <- capitalize(a$ID)
  rownames(a) <- a$ID
  a$ID <- as.factor(a$ID)
  a$threshold <- as.factor(ifelse(a$NES>=0 ,'Up','Down'))
  up <- nrow(a[which(a$threshold=='Up'),])
  down <- nrow(a[which(a$threshold=='Down'),])
  gsea_barplot <- ggplot(a,aes(x = reorder(ID,NES),y = NES))+
    geom_bar(aes(fill=threshold),stat='identity') + 
    scale_fill_manual(values = c('Up'='firebrick3','Down'='navy'))+
    coord_flip()+
    guides(fill='none')+
    xlab('')+
    ylab('NES, LM vs adjacent tissue')+
    labs(title = "GSEA in C5")+
    geom_text(data = a[1:down,],aes(x=ID,y=0.1,label=ID),hjust=0,color='black',size = 3)+
    geom_text(data = a[(down+1):(up+down),],aes(x=ID,y=-0.1,label=ID),hjust=1,color='black',size = 3)+
    theme_bw()+ #背景变为白色
    theme(panel.border = element_blank(),
          axis.line.x = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
          panel.grid.major = element_blank(),   #不显示网格线
          panel.grid.minor = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  ##curve
  library(enrichplot)
  cell <- 'cDC2'
  b <- gsea_result_list[[5]]
  a <- b@result
  GS <- 'GOBP_LYMPHANGIOGENESIS'
  gseaplot<- gseaplot2(b,
                       which(b@result$ID==GS, arr.ind=TRUE),
                       color="green",pvalue_table = T)
  filename <- paste(GS,'_',cell,'.pdf',sep = '')
  pdf(filename,width = 7,height = 6)
  print(gseaplot)
  dev.off()
  
  ##heatmap
  pathway_list = lapply(list.files(pattern = '*.Rdata$', recursive = T)
                        , function(x){
                          load(x)
                          result <- gsea_result@result
                          pathway <- result[which(result$p.adjust<=0.05),]$ID
                          return(pathway)
                        })
  pathway <- do.call(c,pathway_list)
  pathway <- pathway[!duplicated(pathway)]
  pathway <- read.table('pathway.txt',sep = '\t',quote = '')
  pathway <- pathway$V1
  library(stringr)
  gsea_p = lapply(list.files(pattern = '*.Rdata$', recursive = T), function(x){
    load(x)
    cell <- str_split(x,'_')[[1]][1]
    result <- gsea_result@result
    result <- result[pathway,]
    result <- data.frame(result$ID,result$p.adjust)
    colnames(result) <- c('ID',cell)
    return(result)
  })
  gsea_p <- do.call(cbind,gsea_p)
  head(gsea_p)
  rownames(gsea_p) <- pathway
  gsea_p <- gsea_p[,seq(2,ncol(gsea_p),by=2)]
  gsea_p <- as.matrix(gsea_p)
  gsea_p[is.na(gsea_p)] <- 0.99
  gsea_NES = lapply(list.files(pattern = '*.Rdata$', recursive = T), function(x){
    load(x)
    cell <- str_split(x,'_')[[1]][1]
    result <- gsea_result@result
    result <- result[pathway,]
    result <- data.frame(result$ID,result$NES)
    colnames(result) <- c('ID',cell)
    return(result)
  })
  gsea_NES <- do.call(cbind,gsea_NES)
  head(gsea_NES)
  rownames(gsea_NES) <- pathway
  gsea_NES <- gsea_NES[,seq(2,ncol(gsea_NES),by=2)]
  gsea_NES <- as.matrix(gsea_NES)
  gsea_NES[is.na(gsea_NES)] <- 0.01
  
  gsea_p_d <- as.data.frame(gsea_p)
  gsea_NES_d <- as.data.frame(gsea_NES)
  gsea_p_d_log <- -log10(gsea_p_d)
  gsea_NES_d_abs <- gsea_NES_d/abs(gsea_NES_d)
  gsea_p_d_log_NES <- gsea_p_d_log*gsea_NES_d_abs
  max(gsea_p_d_log_NES)
  min(gsea_p_d_log_NES)
  gsea_p_d_log_NES <- as.matrix(gsea_p_d_log_NES)
  library(ComplexHeatmap)
  library(circlize)
  complexheatmap <- Heatmap(gsea_p_d_log_NES, 
                            col=colorRamp2(c(-7,log10(0.05),-log10(0.05),7),c('navy','white','white','firebrick3')), 
                            column_names_side = "top", 
                            cluster_columns=T, 
                            border=T,
                            cluster_rows=T, 
                            show_column_names = T, 
                            show_row_names = T,
                            row_names_gp = gpar(fontsize = 6), 
                            column_names_gp = gpar(fontsize = 8))