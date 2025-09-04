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

  # Base plot with color and shape
  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = .data[[intgroup[1]]], shape = .data[[intgroup[2]]])) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    labs(title = paste("PCA:", analysis_name)) +
    ggsci::scale_color_jco() +
    theme_bw(base_size = 14) +
    coord_fixed()

  # Add fill aesthetic if a third variable is present
  if (length(intgroup) > 2) {
    pca_plot <- pca_plot +
      geom_point(aes(fill = .data[[intgroup[3]]]), size = 5, alpha = 0.8) +
      scale_shape_manual(values = c(21, 22, 23, 24, 25)[1:length(unique(pca_data[[intgroup[2]]]))]) +
      ggsci::scale_fill_jco(alpha = 0.5)
  } else {
    pca_plot <- pca_plot + geom_point(size = 4, alpha = 0.8)
  }

  ggsave(
    filename = here::here(cfg$dir_graphs_png, paste0(cfg$date, "_", analysis_name, "_PCA.png")),
    plot = pca_plot, width = 11, height = 8.5, dpi = 300
  )
  ggsave(
    filename = here::here(cfg$dir_graphs_pdf, paste0(cfg$date, "_", analysis_name, "_PCA.pdf")),
    plot = pca_plot, width = 11, height = 8.5, device = "pdf"
  )

  # --- Sample Distance Heatmap ---
  sample_dists <- dist(t(assay(vsd)))
  sample_dist_matrix <- as.matrix(sample_dists)
  annotation_df <- as.data.frame(colData(dds_processed)[, intgroup, drop = FALSE])

  pheatmap::pheatmap(sample_dist_matrix,
    clustering_distance_rows = sample_dists, clustering_distance_cols = sample_dists,
    annotation_col = annotation_df, main = paste("Sample-to-Sample Distance:", analysis_name),
    border_color = NA,
    filename = here::here(cfg$dir_graphs_png, paste0(cfg$date, "_", analysis_name, "_Heatmap.png"))
  )
  pheatmap::pheatmap(sample_dist_matrix,
    clustering_distance_rows = sample_dists, clustering_distance_cols = sample_dists,
    annotation_col = annotation_df, main = paste("Sample-to-Sample Distance:", analysis_name),
    border_color = NA,
    filename = here::here(cfg$dir_graphs_pdf, paste0(cfg$date, "_", analysis_name, "_Heatmap.pdf"))
  )

  # --- Individual Gene Plots ---
  for (gene_name in cfg$qc_plot_settings$project_genes) {
    gene_id <- gtf_map$gene_id[gtf_map$gene_name == gene_name]
    if (length(gene_id) == 1 && gene_id %in% rownames(dds_processed)) {
      plot_data <- plotCounts(dds_processed, gene = gene_id, intgroup = intgroup, returnData = TRUE)

      gene_plot <- ggplot(plot_data, aes(x = .data[[intgroup[1]]], y = count, fill = .data[[intgroup[1]]])) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.8, size = 3) +
        labs(title = paste("Normalized Counts for", gene_name), subtitle = analysis_name, x = intgroup[1]) +
        ggsci::scale_fill_jco() +
        scale_y_log10() +
        theme_bw(base_size = 14) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

      if (length(intgroup) > 1) {
        facet_formula <- as.formula(paste("~", paste(intgroup[-1], collapse = " + ")))
        gene_plot <- gene_plot + facet_wrap(facet_formula)
      }

      ggsave(
        filename = here::here(cfg$dir_graphs_png, paste0(cfg$date, "_", analysis_name, "_GenePlot_", gene_name, ".png")),
        plot = gene_plot, width = 10, height = 8, dpi = 300
      )
      ggsave(
        filename = here::here(cfg$dir_graphs_pdf, paste0(cfg$date, "_", analysis_name, "_GenePlot_", gene_name, ".pdf")),
        plot = gene_plot, width = 10, height = 8, device = "pdf"
      )
    }
  }
  message("...QC plots saved.")
}
