# 03_QC.R
# Contains a function to generate and save a suite of QC plots.

#' Generate and Save Quality Control Plots
#'
#' @param dds_models The list of DESeqDataSet objects.
#' @param gene_map Data frame to map gene IDs to names for plotting.
#' @param project_genes A character vector of specific genes to plot.
#' @param anno_df A data frame with sample annotations for plotting.
#' @param anno_colors A list of color palettes for annotations.
#' @param output_dir Path to the directory where plots will be saved.

generate_qc_plots <- function(dds_models, gene_map, project_genes, anno_df, anno_colors, output_dir) {
  message("--- Running generate_qc_plots function ---\n")
  
  # Define color scales once
  scale_fill_driver <- scale_fill_manual(values = anno_colors$Driver)
  scale_color_driver <- scale_color_manual(values = anno_colors$Driver)
  
  for (model_name in names(dds_models)) {
    for (subset_name in names(dds_models[[model_name]])) {
      message(paste("...Generating QC plots for:", model_name, "-", subset_name))
      dds <- dds_models[[model_name]][[subset_name]]
      vsd <- vst(dds, blind = TRUE)
      
      # --- 1. PCA Plot ---
      pcaData <- plotPCA(vsd, intgroup = c("Driver", "Host", "Model"), returnData = TRUE)
      percentVar <- round(100 * attr(pcaData, "percentVar"))
      
      pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = Driver, shape = Host)) +
        geom_point(size = 4) + facet_wrap(~Model) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        labs(title = paste("PCA:", model_name, "-", subset_name)) +
        scale_color_driver + theme_bw() + coord_fixed()
      
      ggsave(file.path(output_dir, paste0("QC_PCA_", model_name, "_", subset_name, ".png")), pca_plot, width = 10, height = 8)
      
      # --- 2. Sample-to-Sample Distance Heatmap ---
      sample_dists <- dist(t(assay(vsd)))
      heatmap_plot <- pheatmap(as.matrix(sample_dists),
                               main = paste("Sample Distance:", model_name, "-", subset_name),
                               annotation_col = anno_df, annotation_colors = anno_colors,
                               clustering_distance_rows = sample_dists,
                               clustering_distance_cols = sample_dists)
      
      ggsave(file.path(output_dir, paste0("QC_Heatmap_", model_name, "_", subset_name, ".png")), heatmap_plot, width = 10, height = 8)
      
      # --- 3. Project-Specific Gene Boxplots ---
      for (gene_name in project_genes) {
        gene_id <- gene_map$gene_id[gene_map$gene_name == gene_name]
        if (length(gene_id) == 1 && gene_id %in% rownames(dds)) {
          plot_data <- plotCounts(dds, gene = gene_id, intgroup = c("Driver", "Host", "Model"), returnData = TRUE)
          
          p <- ggplot(plot_data, aes(x = Driver, y = count, fill = Driver)) +
            geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
            facet_wrap(~ Model + Host) + scale_y_log10() +
            labs(title = paste("Counts for", gene_name), subtitle = paste(model_name, "|", subset_name)) +
            scale_fill_driver + theme_bw()
          
          ggsave(file.path(output_dir, paste0("QC_Gene_", gene_name, "_", model_name, "_", subset_name, ".png")), p, width = 10, height = 8)
        }
      }
    }
  }
  message("\n--- Completed QC plot generation ---")
}