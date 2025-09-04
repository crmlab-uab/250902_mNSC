# =============================================================================
#
# 02_QC.R: QUALITY CONTROL PLOTTING FUNCTIONS
#
# =============================================================================

#' Generate and Save Quality Control Plots for a DESeq2 Analysis
#'
generate_qc_plots <- function(dds_processed,
                              analysis_name,
                              gtf_map,
                              cfg) {
  message(paste("\n--- Generating QC plots for:", analysis_name, "---"))
  vsd <- vst(dds_processed, blind = TRUE)
  
  # --- PCA Plot ---
  pca_data <- plotPCA(vsd, intgroup = cfg$main_vars, returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))
  
  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = Driver, fill = Host, shape = Model), size = 5, stroke = 1.2) +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    ggsci::scale_color_jco() +
    ggsci::scale_fill_jco(na.value = "white") +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    labs(title = paste("PCA:", gsub("_", " ", analysis_name))) +
    theme_bw(base_size = 14) +
    coord_fixed()
  
  ggsave(here::here(cfg$dir_graphs_png, paste0(analysis_name, "_PCA.png")), pca_plot, width = 11, height = 8)
  ggsave(here::here(cfg$dir_graphs_pdf, paste0(analysis_name, "_PCA.pdf")), pca_plot, width = 11, height = 8)
  
  # --- Sample Distance Heatmap ---
  sample_dists <- dist(t(assay(vsd)))
  dend <- hclust(sample_dists) %>% as.dendrogram() %>% dendextend::rotate_dend()
  
  annotation_df <- as.data.frame(colData(dds_processed)[, c("Driver", "Model", "Host"), drop = FALSE])
  
  annotation_colors <- lapply(names(annotation_df), function(var) {
    levels <- levels(as.factor(annotation_df[[var]]))
    colors <- ggsci::pal_jco()(length(levels)); names(colors) <- levels; return(colors)
  })
  names(annotation_colors) <- names(annotation_df)
  
  pheatmap::pheatmap(as.matrix(sample_dists),
                     main = paste("Sample Distance:", gsub("_", " ", analysis_name)),
                     annotation_col = annotation_df, annotation_colors = annotation_colors,
                     color = cfg$ryb, clustering_distance_rows = sample_dists,
                     cluster_cols = as.hclust(dend), border_color = NA, 
                     filename = here::here(cfg$dir_graphs_png, paste0(analysis_name, "_Heatmap.png"))
  )
  pheatmap::pheatmap(as.matrix(sample_dists),
                     main = paste("Sample Distance:", gsub("_", " ", analysis_name)),
                     annotation_col = annotation_df, annotation_colors = annotation_colors,
                     color = cfg$ryb, clustering_distance_rows = sample_dists,
                     cluster_cols = as.hclust(dend), border_color = NA,
                     filename = here::here(cfg$dir_graphs_pdf, paste0(analysis_name, "_Heatmap.pdf"))
  )
  
  # --- GOI Box Plots ---
  for (gene_name in cfg$qc_plot_settings$project_genes) {
    gene_id <- gtf_map$gene_id[gtf_map$gene_name == gene_name]
    if (length(gene_id) == 1 && gene_id %in% rownames(dds_processed)) {
      plot_data <- plotCounts(dds_processed, gene = gene_id, intgroup = c("Driver", "Model", "Host"), returnData = TRUE)
      
      p <- ggplot(plot_data, aes(x = Driver, y = count, fill = Driver)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.8, size = 3) +
        ggsci::scale_fill_jco() +
        facet_grid(. ~ Model + Host) + scale_y_log10() +
        labs(title = paste("Normalized Counts for", gene_name), subtitle = gsub("_", " ", analysis_name)) +
        theme_bw(base_size = 14) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
      
      ggsave(here::here(cfg$dir_graphs_png, paste0(analysis_name, "_GenePlot_", gene_name, ".png")), p, width = 10, height = 6)
      ggsave(here::here(cfg$dir_graphs_pdf, paste0(analysis_name, "_GenePlot_", gene_name, ".pdf")), p, width = 10, height = 6)
    }
  }
}

