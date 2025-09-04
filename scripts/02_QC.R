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
#' @param analysis_name A string for naming output files.
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

  # Determine which annotation variables are present in the current data subset
  intgroup <- intersect(cfg$main_vars, names(colData(dds_processed)))
  if (length(intgroup) < 1) {
    stop("None of the main variables are present in this data subset.")
  }

  # --- PCA Plot ---
  pca_data <- plotPCA(vsd, intgroup = intgroup, returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))
  pca_plot <- NULL

  if (length(intgroup) >= 3) {
    message("...Generating 3-variable PCA plot (color=Driver, shape=Model, fill=Host).")
    pca_plot <- ggplot(
      pca_data,
      aes(
        x = .data$PC1,
        y = .data$PC2,
        color = .data[[intgroup[1]]], # e.g., Driver
        shape = .data[[intgroup[3]]], # e.g., Model
        fill = .data[[intgroup[2]]]  # e.g., Host
      )
    ) +
      geom_point(size = 5, alpha = 0.9, stroke = 1.2) +
      xlab(paste0("PC1: ", percent_var[1], "% variance")) +
      ylab(paste0("PC2: ", percent_var[2], "% variance")) +
      labs(title = paste("PCA:", analysis_name)) +
      scale_shape_manual(values = c(21, 22, 23, 24, 25)[1:length(levels(pca_data[[intgroup[3]]]))]) +
      ggsci::scale_color_jco() + # Use JCO palette for outline
      ggsci::scale_fill_jco() +  # Use JCO palette for fill
      theme_bw(base_size = 14) +
      coord_fixed()
  } else if (length(intgroup) >= 2) {
    message("...Generating 2-variable PCA plot.")
    pca_plot <- ggplot(
      pca_data,
      aes(
        x = .data$PC1,
        y = .data$PC2,
        color = .data[[intgroup[1]]],
        shape = .data[[intgroup[2]]]
      )
    ) +
      geom_point(size = 4, alpha = 0.8) +
      xlab(paste0("PC1: ", percent_var[1], "% variance")) +
      ylab(paste0("PC2: ", percent_var[2], "% variance")) +
      labs(title = paste("PCA:", analysis_name)) +
      ggsci::scale_color_jco() + # Use JCO palette
      theme_bw(base_size = 14) +
      coord_fixed()
  } else {
    pca_plot <- ggplot(pca_data, aes(x = .data$PC1, y = .data$PC2, color = .data[[intgroup[1]]])) +
      geom_point(size = 4, alpha = 0.8) +
      xlab(paste0("PC1: ", percent_var[1], "% variance")) +
      ylab(paste0("PC2: ", percent_var[2], "% variance")) +
      labs(title = paste("PCA:", analysis_name)) +
      ggsci::scale_color_jco() + # Use JCO palette
      theme_bw(base_size = 14) +
      coord_fixed()
  }

  # Save plots...
  # ... (ggsave calls remain the same) ...

  # --- Sample Distance Heatmap ---
  sample_dists <- dist(t(assay(vsd)))
  sample_dist_matrix <- as.matrix(sample_dists)
  annotation_df <- as.data.frame(colData(dds_processed)[, intgroup, drop = FALSE])

  # --- Generate pheatmap colors from ggsci palettes ---
  annotation_colors_heatmap <- list()
  for (var in intgroup) {
    n_levels <- length(levels(colData(dds_processed)[[var]]))
    # Extract JCO colors for the number of levels in the factor
    palette_colors <- ggsci::pal_jco("default")(n_levels)
    names(palette_colors) <- levels(colData(dds_processed)[[var]])
    annotation_colors_heatmap[[var]] <- palette_colors
  }

  # Use the generated color list for the heatmap
  pheatmap(
    sample_dist_matrix,
    # ... other pheatmap args ...
    annotation_col = annotation_df,
    annotation_colors = annotation_colors_heatmap,
    # ... rest of pheatmap args ...
  )
  # ... (saving logic remains the same) ...

  # --- Individual Gene Plots ---
  for (gene_name in cfg$qc_plot_settings$project_genes) {
    # ... (logic to get plot_data remains the same) ...
    gene_plot <- ggplot(plot_data, aes(
      x = .data[[group_var]],
      y = .data$count,
      fill = .data[[group_var]]
    )) +
      geom_boxplot(outlier.shape = NA, alpha = 0.5) +
      geom_jitter(width = 0.2, alpha = 0.7, size = 3) +
      labs(
        title = paste("Normalized Counts for", gene_name),
        subtitle = analysis_name,
        x = group_var
      ) +
      ggsci::scale_fill_jco() + # Use JCO palette for box plots
      scale_y_log10() +
      theme_bw(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )

    # ... (facet and save logic remains the same) ...
  }
  message("...QC plots saved.")
}

# =============================================================================
# End of 02_QC.R
# =============================================================================