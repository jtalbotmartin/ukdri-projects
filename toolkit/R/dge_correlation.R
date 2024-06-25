dge_correlations <- function(celltype_list_of_df_1, celltype_list_of_df_2, dependent_var, dir_path, df1_label, df2_label){

    # de_names<-c()
    corr_data <- data.frame('DE','logFC','padj')
    names(corr_data)<-c('DE','logFC','padj')


    for (i in 1:length(celltype_list_of_df_1)){ 
        
        # extract variant from file
        split_celltype_list_df1 <- strsplit(names(celltype_list_of_df_1)[[i]], "_")
        split_celltype_list_df2 <- strsplit(names(celltype_list_of_df_2)[[i]], "_")
        TREM2Variant_df1 <- tail(split_celltype_list_df1[[1]], n=1)
        TREM2Variant_df2 <- tail(split_celltype_list_df2[[1]], n=1)
        
        # check correspondence of variants
        if (TREM2Variant_df1 == TREM2Variant_df2){
            TREM2Variant <- TREM2Variant_df1
        }
        
        # subset and assign data
        df1 <- celltype_list_of_df_1[[i]]
        df2 <- celltype_list_of_df_2[[i]]
        
        # handle NA data
        df1_int <- na.omit(df1)
        df2_int <- na.omit(df2)
        
        # intersect gene lists from df1 and df2 and rename rows
        intersect_df1_df2 <- intersect(df1_int$ensembl_gene_id, df2_int$ensembl_gene_id) 
        row.names(df1_int) <- df1_int$ensembl_gene_id
        row.names(df2_int) <- df2_int$ensembl_gene_id
        
        # subset and order dfs by intersecting genes
        df1_int <- df1_int[c(intersect_df1_df2),]
        df2_int <- df2_int[c(intersect_df1_df2),]
        
        df1_int <- df1_int[ order(row.names(df1_int)), ]
        df2_int <- df2_int[ order(row.names(df2_int)), ]
        
        # generate pearson correlation data for logFC and padj 

        corr_data <- rbind(corr_data,c(paste(dependent_var,celltype,TREM2Variant,sep='_'),
                                       cor(c(df2_int$logFC),
                                           c(df1_int$logFC),
                                           method='pearson'),
                                       cor(c(df2_int$padj),
                                           c(df1_int$padj),
                                           method='pearson')))
        
            
        x <- df1 %>%
        mutate(de = case_when(padj <= 0.05 & logFC >= 0.25 ~ "Up",
                              padj <= 0.05 & logFC <= -0.25 ~ "Down",
                              TRUE ~ "NS"))
        table(x$de)
    
        de_df1 <- df1 %>%
            select(c(gene, logFC, pval, padj)) %>%
            mutate(term = "df1",
                   de = case_when(padj <= 0.05 & abs(logFC) >= 0.25 ~ "de",
                                  TRUE ~ "NS"))
    
        y <- df2 %>%
        mutate(de = case_when(padj <= 0.05 & logFC >= 0.25 ~ "Up",
                              padj <= 0.05 & logFC <= -0.25 ~ "Down",
                              TRUE ~ "NS"))
        table(y$de)
        
        de_df2 <- df2 %>%
            select(c(gene, logFC, pval, padj)) %>%
            mutate(term = "df2",
                   de = case_when(padj <= 0.05 & abs(logFC) >= 0.25 ~ "de",
                                  TRUE ~ "NS")) 

        dt <- rbind(de_df1, de_df2)
    
        label1 <- df1_label
        label2 <- df2_label
        
        dt <- pivot_wider(dt, names_from = term, values_from = c(2:4,6))
        dt <- dt %>%
            mutate(de = case_when(de_df1 == "de" & de_df2 == "de" ~ "Both",
                                  de_df1 == "de" & de_df2 != "de" ~ label1,
                                  de_df1 != "de" & de_df2 == "de" ~ label2,
                                  de_df1 != "de" & de_df2 != "de" ~ "NS",)) %>%
            mutate(de = factor(de, levels = c(label1, label2,
                                              "Both", "NS"))) %>%
            arrange(desc(de))
        
        # define list that will contain top 5 (or less) up and down regulated genes across each significance group
        top_genes <- list()

        # calculate mean of the logFC values
        dt <- dt |>
            mutate(mean_logFC = rowMeans(cbind(logFC_df1, logFC_df2), na.rm = TRUE))
        
        # loop over significance groups to calculate top DE genes
        for (group in c(df1_label, df2_label, "Both")){

            # filter list, arrange my aggregated logFC values and extract to list
            up <- dt |>
                filter(de == group) |>
                filter(mean_logFC >0) |>
                arrange(desc(mean_logFC)) |> 
                head(5) |>
                select(gene)

            top_genes <- append(top_genes, up)
            
            down <- dt |>
                filter(de == group) |>
                filter(mean_logFC <0) |>
                arrange(mean_logFC) |> 
                head(5) |>
                select(gene)
                
            top_genes <- append(top_genes, down)
        }

        # flatten list of lists, and discard duplicate values
        top_genes <- unique(unlist(top_genes)) 

        # add label for only genes present in the top_genes of each significance group
        dt <- dt |>
            mutate(label = ifelse(gene %in% top_genes, gene, NA))


        palette_choice <- paletteer::paletteer_d("ggsci::default_nejm")
        
        # set options to avoid limiting the number of labels
        options(ggrepel.max.overlaps = Inf)
        
        # plot figure including labels of top DE genes

        plot<-ggplot(dt) +
            geom_point(aes(x = logFC_df1, y = logFC_df2, color = de)) +
            geom_hline(yintercept = 0) +
            geom_vline(xintercept = 0) +
            ggrepel::geom_label_repel(aes(x = logFC_df1, y = logFC_df2, label = label), 
                       color = "black", 
                       min.segment.length = 0, 
                       size = 2.5, 
                       box.padding = 0.5, 
                       point.padding = 0.5,
                       force = 3,
                       show.legend = FALSE) +
            xlab(label1)+
            ylab(label2)+
            scale_colour_manual(name = NULL,
                                aesthetics = c("colour", "fill"),
                                values = c(palette_choice[1], palette_choice[2], palette_choice[3], "grey"),
                                breaks = c(label1, label2, "Both", "NS")
            ) +
            scale_alpha_manual(values = c(1, 1, 1, 0.01)) +
            ggtitle(paste0("LogFC of ", df1_label, " vs ", df2_label, " ", TREM2Variant, " in ", head(split_celltype_list_df1[[1]], n=1)))
            
            
            ggsave(paste(dependent_var, head(split_celltype_list_df1[[1]], n=1), TREM2Variant, ".png",sep='_'), path = dir_path)
        
    }

    corr_data<-corr_data[-1,]
    write.csv(corr_data, file=paste0(dir_path, "/corr_table.csv"),row.names = FALSE)

}