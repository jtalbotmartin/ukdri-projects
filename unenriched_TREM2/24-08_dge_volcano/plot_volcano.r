#' Plot volcano plot for differential expression analysis
#'
#' @param dt Differential expression result table from perform_de() function.
#' @param lfc_threshold Log-fold change threshold for the volcano plot. This will be adjusted and plotted as the log2 fold change.  Default is 0.25.
#' @param fdr_threshold The adjusted p-value cut-off for the volcano plot. Default is 0.05.
#' @param n_label The number of top up & down differentially expressed genes to be labeled. Default is 10.
#' 
#' @return No return. The plot is printed out.
#' 
#' @family differential gene expression
#' 
#' @importFrom ggplot2 ggplot geom_point aes coord_cartesian
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr %>% filter top_n
#'
#' @export
volcano_plot <- function(dt,
                         lfc_threshold = 0.25,
                         fdr_threshold = 0.05,
                         n_label = 10) {
  dt <- dt[!is.na(dt$padj), ]
  dt <- dt[!is.nan(dt$logFC), ]
  dt <- dt[order(dt$padj, decreasing = FALSE), ]
  dt$de <- "Not sig"
  dt$de[dt$padj <= fdr_threshold & dt$logFC >= lfc_threshold] <- "Up"
  dt$de[dt$padj <= fdr_threshold & dt$logFC <= -lfc_threshold] <- "Down"
  dt$de <- factor(dt$de, levels = c("Up", "Down", "Not sig"))
  n_up <- sum(dt$de == "Up" & dt$gene_biotype == "protein_coding")
  n_down <- sum(dt$de == "Down" & dt$gene_biotype == "protein_coding")
  
  dt$label <- NA
  if (n_up > 0) {
    top_up <- dt %>%
      dplyr::filter(de == "Up" & gene_biotype == "protein_coding") %>%
      dplyr::top_n(min(n_up, n_label), wt = -padj)
    top_up_logfc <- dt %>%
      dplyr::filter(de == "Up" & gene_biotype == "protein_coding") %>%
      dplyr::top_n(min(n_up, n_label), wt = abs(logFC))
    top_up <- rbind(top_up, top_up_logfc)
    dt$label[dt$gene %in% unique(top_up$gene)] <- "Yes"
  }
  if (n_down > 0) {
    top_down <- dt %>%
      dplyr::filter(de == "Down" & gene_biotype == "protein_coding") %>%
      dplyr::top_n(min(n_down, n_label), wt = -padj)
    top_down_logfc <- dt %>%
      dplyr::filter(de == "Down" & gene_biotype == "protein_coding") %>%
      dplyr::top_n(min(n_down, n_label), wt = abs(logFC))
    top_down <- rbind(top_down, top_down_logfc)
    dt$label[dt$gene %in% unique(top_down$gene)] <- "Yes"
  }
  
  dt$padj[dt$padj == 0] <- min(dt[dt$padj > 0, "padj"]) # to prevent infinite points
  
  #max_logFC <- max(abs(dt$logFC))
  max_logFC <- 3
  
  ggplot2::ggplot(dt) +
    ggplot2::geom_point(ggplot2::aes(x = logFC, y = -log10(padj), fill = de, colour = de),
                        show.legend = T, alpha = 0.5) +
    ggrepel::geom_text_repel(
      data = dt,
      ggplot2::aes(logFC, y = -log10(padj), label = ifelse(label == "Yes", as.character(.data[["gene"]]), "")),
      max.iter = 1000, size = 4, na.rm = TRUE
    ) +
    ggplot2::xlab(bquote(Log[2]*" (fold-change)")) +
    ggplot2::ylab(bquote("-"*Log[10]*" (adjusted p-value)")) +
    ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
                        linetype = 2, size = 0.5, alpha = 0.5) +
    ggplot2::geom_hline(yintercept = -log10(fdr_threshold),
                        linetype = 2, size = 0.5, alpha = 0.5) +
    ggplot2::scale_colour_manual(name = NULL,
                                 aesthetics = c("colour", "fill"),
                                 values = c("#DC0000FF", "#3C5488FF", "grey"),
                                 label = c("Up-regulated", "Down-regulated"),
                                 breaks = c("Up", "Down")) +
    ggplot2::scale_y_continuous(limits = c(0, 3)) +
    ggplot2::scale_x_continuous(limits = c(-max_logFC, max_logFC)) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3))) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "black", size = 16),
      axis.title = ggplot2::element_text(color = "black", size = 18),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.background = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA),
      legend.text = ggplot2::element_text(size = 12, colour = "black"),
      legend.position = "top",
      legend.direction="horizontal",
      plot.margin = ggplot2::margin(c(1, 1, 1, 1), unit = "cm")
    )
}

