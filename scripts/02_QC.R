# =============================================================================
#
# 02_QC.R: QUALITY CONTROL PLOTTING FUNCTIONS
#
# Description:
# This script contains a function to generate and save a suite of
# standard quality control plots for a DESeq2 analysis.
#
# =============================================================================

#' Generate and Save Quality Control Plots for a DESeq2 Analysis
#'
#' @param dds_processed A `DESeqDataSet` object that has been run through `DESeq()`.
#' @param analysis_name A string used for naming output files.
#' @param gtf_map A data frame for mapping gene IDs to gene names.
#' @param cfg The project configuration list.
#'
#' @return This function does not return a value; it saves plots to files.
#'
generate_qc_plots <- function(dds_processed,
                              analysis_name,
                              gtf_map,
                              cfg) {
  message(paste("\n--- Generating QC plots for:", analysis_name, "---"))

  vsd <- vst(dds_processed, blind = TRUE)
  intgroup <- intersect(cfg$main_vars, names(colData(dds_processed)))
  if (length(intgroup) < 1) stop("None of the main variables are present.")

  # --- PCA Plot ---
  pca_data <- plotPCA(vsd, intgroup = intgroup, returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))

  # Initialize the base plot
  p <- ggplot(pca_data, aes(x = .data$PC1, y = .data$PC2)) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    labs(title = paste("PCA:", analysis_name)) +
    theme_bw(base_size = 14) +
    coord_fixed()

  # <<< FIX: Apply aesthetics and scales based on the number of variables
  # This new logic defines aesthetics directly inside geom_point for robustness.
  if (length(intgroup) >= 3) {
    p <- p + geom_point(aes(
      color = .data[[intgroup[1]]], # Driver
      shape = .data[[intgroup[2]]], # Model
      fill = .data[[intgroup[3]]] # Host
    ), size = 5, alpha = 0.9, stroke = 1.2) +
      ggsci::scale_color_jco(name = intgroup[1]) +
      ggsci::scale_fill_jco(name = intgroup[3]) +
      scale_shape_manual(name = intgroup[2], values = c(21, 22, 23, 24, 25)[1:length(unique(pca_data[[intgroup[2]]]))])
  } else if (length(intgroup) == 2) {
    p <- p + geom_point(aes(color = .data[[intgroup[1]]], shape = .data[[intgroup[2]]]), size = 5, alpha = 0.8) +
      ggsci::scale_color_jco(name = intgroup[1])
  } else if (length(intgroup) == 1) {
    p <- p + geom_point(aes(color = .data[[intgroup[1]]]), size = 5, alpha = 0.8) +
      ggsci::scale_color_jco(name = intgroup[1])
  }

  pca_filename_png <- here::here(cfg$dir_graphs_png, paste0(analysis_name, "_PCA.png"))
  pca_filename_pdf <- here::here(cfg$dir_graphs_pdf, paste0(analysis_name, "_PCA.pdf"))
  ggsave(filename = pca_filename_png, plot = p, width = 10, height = 8, dpi = 300)
  ggsave(filename = pca_filename_pdf, plot = p, width = 10, height = 8, device = "pdf")


  # --- Sample Distance Heatmap ---
  sample_dists <- dist(t(assay(vsd)))
  sample_dist_matrix <- as.matrix(sample_dists)
  annotation_df <- as.data.frame(colData(dds_processed)[, intgroup, drop = FALSE])

  annotation_colors <- lapply(intgroup, function(var) {
    levels <- levels(annotation_df[[var]])
    colors <- ggsci::pal_jco()(length(levels))
    names(colors) <- levels
    return(colors)
  })
  names(annotation_colors) <- intgroup

  dend <- as.dendrogram(hclust(sample_dists))
  dend <- dendextend::rotate(dend, order = rownames(sample_dist_matrix))

  heatmap_filename_png <- here::here(cfg$dir_graphs_png, paste0(analysis_name, "_Heatmap.png"))
  heatmap_filename_pdf <- here::here(cfg$dir_graphs_pdf, paste0(analysis_name, "_Heatmap.pdf"))

  pheatmap::pheatmap(sample_dist_matrix,
    cluster_cols = as.hclust(dend), cluster_rows = as.hclust(dend),
    annotation_col = annotation_df, annotation_colors = annotation_colors,
    main = paste("Sample-to-Sample Distance:", analysis_name), color = cfg$ryb,
    border_color = NA, filename = heatmap_filename_png
  )
  pheatmap::pheatmap(sample_dist_matrix,
    cluster_cols = as.hclust(dend), cluster_rows = as.hclust(dend),
    annotation_col = annotation_df, annotation_colors = annotation_colors,
    main = paste("Sample-to-Sample Distance:", analysis_name), color = cfg$ryb,
    border_color = NA, filename = heatmap_filename_pdf
  )


  # --- Individual Gene Plots ---
  for (gene_name in cfg$qc_plot_settings$project_genes) {
    gene_id <- gtf_map$gene_id[gtf_map$gene_name == gene_name]
    if (length(gene_id) == 1 && gene_id %in% rownames(dds_processed)) {
      plot_data <- plotCounts(dds_processed, gene = gene_id, intgroup = intgroup[1], returnData = TRUE)

      p_gene <- ggplot(plot_data, aes(x = .data[[intgroup[1]]], y = count, fill = .data[[intgroup[1]]])) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.8, size = 3) +
        ggsci::scale_fill_jco() +
        labs(
          title = paste("Normalized Counts for", gene_name),
          subtitle = analysis_name,
          x = intgroup[1]
        ) +
        scale_y_log10() +
        theme_bw(base_size = 14) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

      gene_plot_filename_png <- here::here(cfg$dir_graphs_png, paste0(analysis_name, "_GenePlot_", gene_name, ".png"))
      gene_plot_filename_pdf <- here::here(cfg$dir_graphs_pdf, paste0(analysis_name, "_GenePlot_", gene_name, ".pdf"))

      ggsave(filename = gene_plot_filename_png, plot = p_gene, width = 8, height = 7, dpi = 300)
      ggsave(filename = gene_plot_filename_pdf, plot = p_gene, width = 8, height = 7, device = "pdf")
    }
  }
}
