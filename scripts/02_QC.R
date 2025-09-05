# =============================================================================
#
# 02_QC.R: QUALITY CONTROL PLOTTING FUNCTIONS
#
# =============================================================================

generate_qc_plots <- function(dds_processed,
                              analysis_name,
                              gtf_map,
                              cfg) {
  vsd <- vst(dds_processed, blind = FALSE)
  num_genes <- nrow(dds_processed)

  # PCA Plot
  pca_data <- plotPCA(vsd, intgroup = cfg$main_vars, returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))
  pca_plot <- ggplot(pca_data) +
    geom_point(
      mapping = aes(
        x = PC1,
        y = PC2,
        color = .data[[cfg$main_vars[1]]],
        fill = .data[[cfg$main_vars[3]]],
        shape = .data[[cfg$main_vars[2]]]
      ),
      size = 5,
      stroke = 1.5
    ) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    labs(
      title = paste("PCA:", gsub("_", " ", analysis_name)),
      subtitle = paste(num_genes, "genes included"),
      color = cfg$main_vars[1],
      fill = cfg$main_vars[3],
      shape = cfg$main_vars[2]
    ) +
    scale_shape_manual(values = c(21, 22, 24, 25)) +
    scale_fill_manual(values = c("BL6" = "black", "NSG" = "white")) +
    ggsci::scale_color_jco() +
    coord_fixed() +
    theme_bw(base_size = 14)

  ggsave(
    filename = here::here(
      cfg$dir_graphs_png,
      paste0(cfg$date, "_", analysis_name, "_PCA.png")
    ),
    plot = pca_plot,
    width = 11,
    height = 8.5,
    dpi = 300
  )
  ggsave(
    filename = here::here(
      cfg$dir_graphs_pdf,
      paste0(cfg$date, "_", analysis_name, "_PCA.pdf")
    ),
    plot = pca_plot,
    width = 11,
    height = 8.5,
    device = "pdf"
  )

  # Sample Distance Heatmap
  sample_dists <- dist(t(assay(vsd)))
  sample_dist_matrix <- as.matrix(sample_dists)

  # <<< FIX: Explicitly set the track order and droplevels() to prevent color mismatch errors.
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

  pheatmap::pheatmap(
    sample_dist_matrix,
    main = paste("Sample-to-Sample Distance:", gsub("_", " ", analysis_name)),
    subtitle = paste(num_genes, "genes included"),
    annotation_col = annotation_df,
    annotation_colors = ann_colors,
    color = cfg$ryb,
    border_color = NA,
    cluster_rows = dend,
    cluster_cols = dend,
    show_rownames = FALSE,
    show_colnames = FALSE,
    filename = here::here(
      cfg$dir_graphs_png,
      paste0(cfg$date, "_", analysis_name, "_Heatmap.png")
    )
  )
  pheatmap::pheatmap(
    sample_dist_matrix,
    main = paste("Sample-to-Sample Distance:", gsub("_", " ", analysis_name)),
    subtitle = paste(num_genes, "genes included"),
    annotation_col = annotation_df,
    annotation_colors = ann_colors,
    color = cfg$ryb,
    border_color = NA,
    cluster_rows = dend,
    cluster_cols = dend,
    show_rownames = FALSE,
    show_colnames = FALSE,
    filename = here::here(
      cfg$dir_graphs_pdf,
      paste0(cfg$date, "_", analysis_name, "_Heatmap.pdf")
    )
  )


  # GOI Boxplots
  for (gene_name in cfg$qc_plot_settings$project_genes) {
    gene_id <- gtf_map$gene_id[gtf_map$gene_name == gene_name]
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
            color = .data[[cfg$main_vars[1]]],
            fill = .data[[cfg$main_vars[3]]],
            shape = .data[[cfg$main_vars[2]]]
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
        scale_shape_manual(values = c(21, 22, 24, 25)) +
        scale_fill_manual(values = c("BL6" = "black", "NSG" = "white")) +
        ggsci::scale_color_jco() +
        theme_bw(base_size = 14) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

      ggsave(
        filename = here::here(
          cfg$dir_graphs_png,
          paste0(
            cfg$date,
            "_",
            analysis_name,
            "_GenePlot_",
            gene_name,
            ".png"
          )
        ),
        plot = p,
        width = 10,
        height = 8,
        dpi = 300
      )
      ggsave(
        filename = here::here(
          cfg$dir_graphs_pdf,
          paste0(
            cfg$date,
            "_",
            analysis_name,
            "_GenePlot_",
            gene_name,
            ".pdf"
          )
        ),
        plot = p,
        width = 10,
        height = 8,
        device = "pdf"
      )
    }
  }
}
