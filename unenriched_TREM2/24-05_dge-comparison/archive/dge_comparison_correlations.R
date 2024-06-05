library(tidyr)
library(dplyr)
library(ggplot2)

setwd("/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/dge/")

# prev, new, comparison factor prev, comparison factor new
DEs<-list(c('/archive/previous_results/de_results_diagnosis_cngeneson_pc_mito_sex_brain_region_apoe_CD33_group_age_PMD/',
            '/DGE/de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group/',
            'diagnosis',
            'NeuropathologicalDiagnosis'),
          c('/archive/previous_results/de_results_TREM2_stratified_PHF1_cngeneson_pc_mito_sex_brain_region_apoe_CD33_group/',
            '/DGE/de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group/',
            'PHF1', 
            'pctPHF1PositiveArea'))

TREM2Variants<-c('CV','R47H','R62H')
celltypes<-c('Micro','Astro')

de_names<-c()
corr_data <- data.frame('DE','logFC','padj')
names(corr_data)<-c('DE','logFC','padj')

for (DE in DEs){
  for (celltype in celltypes){
    for (TREM2Variant in TREM2Variants){  
      
      prev<-read.csv(paste0(getwd(), "/", DE[1],celltype,'/',celltype,'_',DE[3],'_in_',TREM2Variant,'.tsv',sep=''),sep = '\t')
      new<-read.csv(paste0(getwd(), "/", DE[2],celltype,'/',celltype,'_',DE[4],'_in_',TREM2Variant,'.tsv',sep=''),sep = '\t')
      
      prev_int<-na.omit(prev)
      new_int<-na.omit(new)
      
      intersect_prev_new<-intersect(prev_int$ensembl_gene_id,new_int$ensembl_gene_id) 
      
      row.names(prev_int)<-prev_int$ensembl_gene_id
      row.names(new_int)<-new_int$ensembl_gene_id
      
      prev_int<-prev_int[c(intersect_prev_new),]
      new_int<-new_int[c(intersect_prev_new),]
      
      prev_int<-prev_int[ order(row.names(prev_int)), ]
      new_int<-new_int[ order(row.names(new_int)), ]
      
      corr_data<-rbind(corr_data,c(paste(DE[3],celltype,TREM2Variant,sep='_'),
                                   cor(c(new_int$logFC),
                                       c(prev_int$logFC),
                                       method='pearson'),
                                   cor(c(new_int$padj),
                                       c(prev_int$padj),
                                       method='pearson')))
      
      
      dt_old <- prev
      
      x <- dt_old %>%
        mutate(de = case_when(padj <= 0.05 & logFC >= 0.25 ~ "Up",
                              padj <= 0.05 & logFC <= -0.25 ~ "Down",
                              TRUE ~ "NS"
        ))
      table(x$de)
      
      
      de_old <- dt_old %>%
        select(c(gene, logFC, pval, padj)) %>%
        mutate(term = "old",
               de = case_when(padj <= 0.05 & abs(logFC) >= 0.25 ~ "de",
                              TRUE ~ "NS"))
      
      
      
      dt_new <- new
      
      y <- dt_new %>%
        mutate(de = case_when(padj <= 0.05 & logFC >= 0.25 ~ "Up",
                              padj <= 0.05 & logFC <= -0.25 ~ "Down",
                              TRUE ~ "NS"
        ))
      table(y$de)
      
      de_new <- dt_new %>%
        select(c(gene, logFC, pval, padj)) %>%
        mutate(term = "new",
               de = case_when(padj <= 0.05 & abs(logFC) >= 0.25 ~ "de",
                              TRUE ~ "NS")) 
      
      dt <- rbind(de_old, de_new)
      
      
      label1 <- "prev"
      label2 <- "new"
      
      dt <- pivot_wider(dt, names_from = term, values_from = c(2:4,6))
      dt <- dt %>%
        mutate(de = case_when(de_old == "de" & de_new == "de" ~ "Both",
                              de_old == "de" & de_new != "de" ~ label1,
                              de_old != "de" & de_new == "de" ~ label2,
                              de_old != "de" & de_new != "de" ~ "NS",
        )
        ) %>%
        mutate(de = factor(de, levels = c(label1, label2,
                                          "Both", "NS"))) %>%
        arrange(desc(de))
      
      
      palette_choice <- paletteer::paletteer_d("ggsci::default_nejm")
      
      
      plot<-ggplot(dt) +
        geom_point(aes(x = logFC_old, y = logFC_new, color = de)) +
        geom_hline(yintercept = 0) +
        geom_vline(xintercept = 0) +
        xlab(label1)+
        ylab(label2)+
        scale_colour_manual(name = NULL,
                            aesthetics = c("colour", "fill"),
                            values = c(palette_choice[1], palette_choice[2], palette_choice[3], "grey"),
                            breaks = c(label1, label2, "Both", "NS")
        ) +
        scale_alpha_manual(values = c(1, 1, 1, 0.01)) +
        ggtitle(paste0("LogFC of Prev vs New ", TREM2Variant, " in ", celltype))
      
      
      ggsave(paste(DE[3],celltype,TREM2Variant,".png",sep='_'))
      
    }}}

corr_data<-corr_data[-1,]
write.csv(corr_data,file='corr_table.csv',row.names = FALSE)
