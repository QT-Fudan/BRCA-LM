library(survival)
library(survminer)
library(ggplot2)

dat <- RNA_df_d
dat_d <- as.data.frame(t(dat))

gene1 <- 'MHC'
gene2 <- 'XCR1'
cancer <- 'BRCA'
level <- 'low'

survplot_list <- list()
for(x in 1:length(gene1)){
  gene <- gene1[x]
  dat_t <- dat_d[order(dat_d[,gene2]),]
  nlow <-  round(nrow(dat_t)*0.5,digits = 0)
  nhigh <-nrow(dat_t)-round(nrow(dat_t)*0.5,digits = 0)
  dat_t$gene2_level <- c(rep('Low',nlow),rep('High',nhigh))
  dat_t <- dat_t[order(dat_t$gene2_level,dat_t[,gene]),]
  dat_t$gene1_level <- c(rep('Low',round(nlow*0.25,digits = 0)),rep('medium',nlow-2*round(nlow*0.25,digits = 0)),rep('High',round(nlow*0.25,digits = 0)),
                         rep('Low',round(nhigh*0.25,digits = 0)),rep('medium',nhigh-2*round(nhigh*0.25,digits = 0)),rep('High',round(nhigh*0.25,digits = 0)))
  dat_t$gene_level <- paste(dat_t$gene1_level,'_',dat_t$gene2_level,sep = '')
  dat_t$patient <- row.names(dat_t)
  clinical_df$bcr_patient_uuid <- tolower(clinical_df$bcr_patient_uuid)
  dat_t_clinical <- merge(dat_t,clinical_df,by.x = 'patient',by.y = 'bcr_patient_uuid')
  dat_survival <- data.frame(patient=dat_t_clinical$patient,
                             gene_level=dat_t_clinical$gene_level,
                             days_to_death=dat_t_clinical$days_to_death,
                             days_to_last_followup=dat_t_clinical$days_to_last_followup)
  
  dat_lost <- dat_survival[which((dat_survival$days_to_death=='') & (dat_survival$days_to_last_followup!='')),]
  dat_lost <- dat_lost[,-3]
  dat_lost$status <-rep(1,nrow(dat_lost)) 
  colnames(dat_lost) <- c('patient','gene_level','time','status')
  
  dat_death <- dat_survival[which(dat_survival$days_to_death!=''),]
  dat_death <- dat_death[,-4]
  dat_death$status <-rep(2,nrow(dat_death)) 
  colnames(dat_death) <- c('patient','gene_level','time','status')
  
  survival <- merge(dat_lost,dat_death,all = T)
  survival$time <- as.numeric(survival$time)
  survival_gene2_high <- survival[which(survival$gene_level=='High_High'|survival$gene_level=='Low_High'),]
  survival_gene2_low <- survival[which(survival$gene_level=='High_Low'|survival$gene_level=='Low_Low'),]
  if(level=='low'){
    survival_target <- survival_gene2_low
  }else{
    survival_target <- survival_gene2_high
  }
  fit <- survfit(Surv(time,status) ~ gene_level,  # 创建生存对象 
                 data = survival_target) # 数据集来???
  title <- paste(gene,gene2,level, cancer)
  survplot <- ggsurvplot(fit, data = survival_target,
                         risk.table = TRUE,
                         pval = TRUE,
                         palette=c('red','blue'))+
    labs(title=title)
  survplot_list[[x]] <- survplot
}

