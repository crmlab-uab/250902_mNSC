# =============================================================================
#
# 02_QC.R: QUALITY CONTROL PLOTTING FUNCTIONS
#
# =============================================================================




#' Generate Quality Control Plots
#'
#' This function generates a variety of QC plots for a given DESeq2 analysis,
#' including a PCA plot, a sample distance heatmap, and boxplots for specified
#' genes of interest.
#'
#' @param dds_processed A processed DESeqDataSet object.
#' @param analysis_name A string used to name the output files.
#' @param gtf_map A data frame mapping gene IDs to gene names.
#' @param cfg A list containing the project configuration.
#'
#' @return This function does not return a value but saves plots to files.
generate_qc_plots <- function(dds_processed,
                              analysis_name,
                              gtf_map,
                              cfg) {
  vsd <- vst(dds_processed, blind = FALSE)
  num_genes <- nrow(dds_processed)

  # PCA Plot
  pca_data <- plotPCA(vsd, intgroup = cfg$main_vars, returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))

  # Create a custom color mapping for the Model variable
  model_levels <- levels(pca_data$Model)
  model_colors <- ggsci::pal_jco("default")(length(model_levels))
  names(model_colors) <- model_levels

  # Create a custom fill variable for the PCA plot
  pca_data$fill_color <- ifelse(pca_data$Host == "BL6", model_colors[as.character(pca_data$Model)], "white")

  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
    geom_point(
      mapping = aes(
        shape = .data[[cfg$main_vars[1]]], # Driver
        color = .data[[cfg$main_vars[2]]], # Model
        fill = I(fill_color)
      ),
      size = 5,
      stroke = 1.5
    ) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    labs(
      title = paste("PCA:", gsub("_", " ", analysis_name)),
      subtitle = paste(num_genes, "genes included"),
      shape = cfg$main_vars[1],
      color = cfg$main_vars[2]
    ) +
    scale_shape_manual(values = c(21, 22, 24)) + # Fillable shapes
    scale_color_manual(values = model_colors) +
    coord_fixed() +
    theme_bw(base_size = 14)

  save_plot_formats(
    plot_object = pca_plot,
    dir_png = cfg$dir_graphs_png,
    dir_pdf = cfg$dir_graphs_pdf,
    filename_base = paste0(cfg$date, "_", analysis_name, "_PCA"),
    width = 11,
    height = 8.5
  )

  # Sample Distance Heatmap
  sample_dists <- dist(t(assay(vsd)))
  sample_dist_matrix <- as.matrix(sample_dists)

  track_order <- intersect(c("Driver", "Model", "Host"), names(colData(dds_processed)))
  annotation_df <- as.data.frame(colData(dds_processed)[, track_order, drop = FALSE]) %>%
    droplevels()

  ann_colors <- list()
  jco_palette <- ggsci::pal_jco("default")(10)
  color_idx <- 1

  for (var in names(annotation_df)) {
    levels_in_data <- levels(annotation_df[[var]])
    if (length(levels_in_data) > 0) {
      var_colors <- jco_palette[color_idx:(color_idx + length(levels_in_data) - 1)]
      names(var_colors) <- levels_in_data
      ann_colors[[var]] <- var_colors
      color_idx <- color_idx + length(levels_in_data)
    }
  }

  dend <- as.dendrogram(hclust(sample_dists))
  dend <- dendextend::rotate(dend, order = rownames(sample_dist_matrix))

  p_heatmap <- pheatmap::pheatmap(
    sample_dist_matrix,
    main = paste("Sample-to-Sample Distance:", gsub("_", " ", analysis_name)),
    annotation_col = annotation_df,
    annotation_colors = ann_colors,
    color = cfg$ryb,
    border_color = NA,
    cluster_rows = as.hclust(dend),
    cluster_cols = as.hclust(dend),
    show_rownames = FALSE,
    show_colnames = FALSE,
    silent = TRUE
  )

  save_plot_formats(
    plot_object = p_heatmap,
    dir_png = cfg$dir_graphs_png,
    dir_pdf = cfg$dir_graphs_pdf,
    filename_base = paste0(cfg$date, "_", analysis_name, "_Heatmap"),
    width = 11,
    height = 8.5
  )

  # GOI Boxplots
  for (gene_name in cfg$qc_plot_settings$project_genes) {
    gene_id <- gtf_map$gene_id[gtf_map$gene_name_upper == toupper(gene_name)]
    if (length(gene_id) == 1 &&
      gene_id %in% rownames(dds_processed)) {
      plot_data <- plotCounts(
        dds_processed,
        gene = gene_id,
        intgroup = cfg$main_vars,
        returnData = TRUE
      )

      p <- ggplot(plot_data, aes(x = .data[[cfg$main_vars[1]]], y = count)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.5) +
        geom_jitter(
          aes(
            shape = .data[[cfg$main_vars[1]]],
            color = .data[[cfg$main_vars[2]]],
            fill = .data[[cfg$main_vars[3]]]
          ),
          width = 0.2,
          size = 4,
          stroke = 1.2
        ) +
        labs(
          title = paste("Normalized Counts for", gene_name),
          subtitle = gsub("_", " ", analysis_name),
          x = cfg$main_vars[1],
          y = "Normalized Counts"
        ) +
        scale_y_log10() +
        scale_shape_manual(values = c(21, 22, 24)) +
        scale_fill_manual(values = c("BL6" = "black", "NSG" = "grey80")) +
        ggsci::scale_color_jco() +
        theme_bw(base_size = 14) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

      save_plot_formats(
        plot_object = p,
        dir_png = cfg$dir_graphs_png,
        dir_pdf = cfg$dir_graphs_pdf,
        filename_base = paste0(cfg$date, "_", analysis_name, "_GenePlot_", gene_name),
        width = 10,
        height = 8
      )
    } else {
      message(paste("Gene of interest not found:", gene_name))
    }
  }
}
