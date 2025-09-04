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
  
  # Determine which annotation variables are present in the current data subset
  intgroup <- intersect(cfg$main_vars, names(colData(dds_processed)))
  if (length(intgroup) < 1) {
    stop("None of the main variables are present in this data subset.")
  }
  if (length(intgroup) < 2) {
    intgroup <- c(intgroup, intgroup[1]) # Duplicate if only one for syntax
  }
  
  # --- FIX: Dynamically create a color list for the CURRENT data subset ---
  # This robustly handles cases where not all factor levels are present.
  plot_colors <- list()
  for (var in intgroup) {
    # Get the factor levels that are actually in the current data
    levels_in_data <- levels(colData(dds_processed)[[var]])
    # From the master color list in cfg, pull the colors for ONLY those levels
    if (var %in% names(cfg$qc_plot_settings$annotation_colors)) {
      plot_colors[[var]] <- cfg$qc_plot_settings$annotation_colors[[var]][levels_in_data]
      # Remove any NA values that might result if a level has no specified color
      plot_colors[[var]] <- plot_colors[[var]][!is.na(plot_colors[[var]])]
    }
  }
  
  # --- PCA Plot ---
  pca_data <- plotPCA(vsd, intgroup = intgroup, returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))
  
  pca_plot <- ggplot(pca_data,
                     aes(
                       x = .data$PC1,
                       y = .data$PC2,
                       color = .data[[intgroup[1]]],
                       shape = .data[[intgroup[2]]]
                     )) +
    geom_point(size = 4, alpha = 0.8) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    labs(title = paste("PCA:", analysis_name)) +
    # Use the dynamically created color list for this plot
    scale_color_manual(values = plot_colors[[intgroup[1]]]) +
    theme_bw(base_size = 14) +
    coord_fixed()
  
  pca_filename <- here::here(cfg$dir_graphs,
                             paste0(cfg$date, "_", analysis_name, "_PCA.pdf"))
  ggsave(
    filename = pca_filename,
    plot = pca_plot,
    width = 10,
    height = 8,
    device = "pdf"
  )
  
  # --- Sample Distance Heatmap ---
  sample_dists <- dist(t(assay(vsd)))
  sample_dist_matrix <- as.matrix(sample_dists)
  annotation_df <- as.data.frame(colData(dds_processed)[, intgroup, drop = FALSE])
  
  heatmap_filename <- here::here(cfg$dir_graphs,
                                 paste0(cfg$date, "_", analysis_name, "_Heatmap.pdf"))
  pdf(heatmap_filename, width = 10, height = 8)
  pheatmap(
    sample_dist_matrix,
    clustering_distance_rows = sample_dists,
    clustering_distance_cols = sample_dists,
    annotation_col = annotation_df,
    # Use the dynamically created color list for the heatmap
    annotation_colors = plot_colors,
    main = paste("Sample-to-Sample Distance:", analysis_name)
  )
  dev.off()
  
  # --- Individual Gene Plots ---
  for (gene_name in cfg$qc_plot_settings$project_genes) {
    gene_id <- gtf_map$gene_id[gtf_map$gene_name == gene_name]
    if (length(gene_id) == 1 &&
        gene_id %in% rownames(dds_processed)) {
      plot_data <- plotCounts(
        dds_processed,
        gene = gene_id,
        intgroup = intgroup,
        returnData = TRUE
      )
      
      gene_plot <- ggplot(plot_data, aes(
        x = .data[[intgroup[1]]],
        y = .data$count,
        fill = .data[[intgroup[1]]]
      )) +
        geom_boxplot(outlier.shape = NA, alpha = 0.5) +
        geom_jitter(width = 0.2,
                    alpha = 0.7,
                    size = 3) +
        labs(
          title = paste("Normalized Counts for", gene_name),
          subtitle = analysis_name,
          x = intgroup[1]
        ) +
        scale_fill_manual(values = plot_colors[[intgroup[1]]]) +
        scale_y_log10() +
        theme_bw(base_size = 14) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none")
      
      if (length(intgroup) > 1) {
        facet_formula <- as.formula(paste("~", paste(intgroup[-1], collapse = " + ")))
        gene_plot <- gene_plot + facet_wrap(facet_formula)
      }
      
      gene_plot_filename <- here::here(
        cfg$dir_graphs,
        paste0(
          cfg$date,
          "_",
          analysis_name,
          "_GenePlot_",
          gene_name,
          ".pdf"
        )
      )
      ggsave(
        filename = gene_plot_filename,
        plot = gene_plot,
        width = 10,
        height = 8,
        device = "pdf"
      )
    }
  }
  message("...QC plots saved.")
}