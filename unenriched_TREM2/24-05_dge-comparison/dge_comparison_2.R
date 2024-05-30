library(dplyr)
library(ggplot2)
library(tidyr)

setwd("~/Desktop/Biogen_data/de_results")

outdir <- "bar_plot_DEG"
dir.create(outdir)

input_dir <- "de_results_amyloid_beta_cngeneson_pc_mito_sex_brain_region_apoe_CD33_group"
#celltype <- c("Astro", "Micro", "EN", "IN")
celltype <- dir(input_dir)
#celltype <- celltype[grep("EN-|IN-", celltype)]
#contrast <- c("CV", "R47H", "R62H")

res_all <- list()


for( i in celltype){
  
  file_path <- sprintf("%s/%s", input_dir, i)
  
  file_list <- list.files(path = file_path, pattern =".tsv" , full.names = TRUE)[-4]
  
  dt_l <- lapply(file_list, read.delim)
  names(dt_l) <- basename(file_list) %>% tools::file_path_sans_ext()
  
  res_l <- lapply(names(dt_l), function(x){
    
    de_up <- dt_l[[x]] %>%
      dplyr::filter(padj <= 0.05, logFC >= 0.25) %>%
      dplyr::pull(gene) %>%
      as.character() %>%
      length()
    
    de_down <- dt_l[[x]] %>%
      dplyr::filter(padj <= 0.05, logFC <= -0.25) %>%
      dplyr::pull(gene) %>%
      as.character() %>%
      length()
    
    
    tmp_dt <- data.frame(contrast = x,
                         total_expressed_gene = nrow(dt_l[[x]]),
                         up = de_up,
                         down = de_down,
                         pct_up = de_up*100/nrow(dt_l[[x]]),
                         pct_down = de_down*100/nrow(dt_l[[x]])
    )
  })
  
  res <- do.call(rbind, res_l)
  
  res_all[[i]] <- res
  
}

res <- do.call(rbind, res_all)
res <- res %>%
  mutate(celltype = purrr::map_chr(contrast, ~strsplit(., "_")[[1]][1]),
         trem2 = purrr::map_chr(contrast, ~strsplit(., "_")[[1]][5]),
         pct_down = -1*pct_down)

res_l <- pivot_longer(res, 
                      cols = c(pct_up,pct_down),
                      names_to = "pct_de",
                      values_to = "value"
) %>%
  mutate(label = ifelse(value > 0, paste0(round(value, 2), "%"), paste0(-1* round(value, 2), "%")),
         label_y = ifelse(value > 0, value + 1, value -1)
  )



palette_choice <- paletteer::paletteer_d("ggsci::nrc_npg")

ggplot(res_l, aes(x = trem2, y = value, fill = trem2)) +
  geom_bar(stat = "identity", position = position_dodge(0.8)) +
  facet_wrap(~ celltype) +
  geom_hline(yintercept = 0) +
  scale_color_manual(name = "TREM2", values = palette_choice, 
                     aesthetics = c("colour", "fill")) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(size = 16, colour = "black"),
    axis.title = ggplot2::element_text(size = 18),
    legend.text = ggplot2::element_text(size = 16),
    legend.title = ggplot2::element_text(size = 16),
    strip.text = ggplot2::element_text(size = 16)
  ) +
  geom_text(aes(label = label, y = label_y)) 

####Astro ####
res_l %>%
  filter(celltype %in% c("Astro") )%>%
  ggplot(., aes(x = trem2, y = value, fill = trem2)) +
  geom_bar(stat = "identity", position = position_dodge(0.8)) +
  xlab("") +
  ylab("") +
  facet_wrap(~ celltype) +
  geom_hline(yintercept = 0) +
  scale_color_manual(name = "TREM2", values = palette_choice, 
                     aesthetics = c("colour", "fill")) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(size = 16, colour = "black"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = ggplot2::element_text(size = 18),
    legend.text = ggplot2::element_text(size = 16),
    legend.title = ggplot2::element_text(size = 16),
    strip.text = ggplot2::element_text(size = 16),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm")
  ) +
  geom_text(aes(label = label, y = label_y)) +
  coord_cartesian( clip="off") +
  annotation_custom(grob = linesGrob(arrow=arrow(type="open", 
                                                 ends="first", 
                                                 length=unit(3,"mm")),
                                     gp=gpar(col="black", lwd=2)),
                    xmin = 0.25, xmax = 0.25, 
                    ymax = -1, ymin = -7)+
  annotation_custom(grob = textGrob(label = "% Down",
                                    hjust=0.5, 
                                    rot = 90,
                                    gp=gpar(col="black", cex=1.2)),
                    xmin = 0.1, xmax = 0.1, 
                    ymin = -1, ymax = -6) +
  annotation_custom(grob = linesGrob(arrow=arrow(type="open",
                                                 ends="last",
                                                 length=unit(3,"mm")),
                                     gp=gpar(col="black", lwd=2)),
                    xmin = 0.25, xmax = 0.25, 
                    ymax = 1, ymin = 6) +
  annotation_custom(grob = textGrob(label = "% Up",
                                    hjust=0.5, 
                                    rot = 90,
                                    gp=gpar(col="black", cex=1.2)),
                    xmin = 0.1,xmax = 0.1, 
                    ymin = 1, ymax = 6) 

ggsave("bar_plot_DEG/Astro_amyloid_beta_trem2_stratified.png", 
       width = 5, height = 5, units = "in", dpi = 300)