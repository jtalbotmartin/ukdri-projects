##########################################################################
# function to plot individual TREM2 variant enrichment scores
##########################################################################

library(scales)

plot_facet_celltype_facet_by_trem2 <- function(trem2var_stratified_df, title_label){

    # plot geneset enrichment scores, by each celltype, faceted by each contrast
    # colour by z-score, outline by significance

    ggplot(trem2var_stratified_df, aes(x = contrast_truncated, y = geneset, fill = z.std, stroke = ifelse(significance, 0.75, 0))) +
    geom_point(aes(size = abs(z.std)), shape = 21) +
    xlab("") +
    ylab("") +
    ggtitle(title_label) +
    scale_x_discrete(position = "bottom") +
    facet_grid(cols = vars(trem2)) +
    scale_fill_gradient2(name = "Z-score", 
                       low = "#5164C2", mid = "white", high = "#D71157", 
                       midpoint = 0, space = "Lab", 
                       na.value = "grey50", guide = "colourbar") +
    scale_size(range = c(2, 7)) +
    guides(size = "none", color = "none") +
    theme_bw() +
    theme(axis.text.x = element_text(colour = "black", size = 11, angle = 90, vjust = 0.5, hjust = 1),
            axis.text.y = element_text(colour = "black", size = 11), 
            strip.text.x = element_text(size = 12), 
            strip.background = element_blank(),
            panel.spacing = unit(0.1, "lines")
  )
}