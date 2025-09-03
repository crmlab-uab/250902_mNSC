# 03_QC.R
# This script performs comprehensive quality control on DESeq2 objects,
# generating and saving a suite of standard and custom plots.

message("--- Running 03_QC.R: Comprehensive Quality Control ---\n")

# Load required libraries
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(DESeq2)
library(rtracklayer)
library(tidyr)
library(tibble)
library(dplyr)
library(RNAseqQC)

# Ensure the main data object from 02_data_prep.R exists
if (!exists("dds_models")) {
  stop("dds_models object not found. Please ensure 02_data_prep.R has been run.")
}

# --- Comprehensive QC Plotting Loop ---
# Iterate through each statistical model and data subset
for (model_name in names(dds_models)) {
  for (subset_name in names(dds_models[[model_name]])) {
    
    message(paste("...Generating QC plots for:", model_name, "-", subset_name))
    
    # Get the specific DESeqDataSet object for the current model and subset
    dds <- dds_models[[model_name]][[subset_name]]
    
    # --- Define Color Scales from Parent RMD ---
    if (exists("anno_colors") && is.list(anno_colors)) {
      scale_fill_driver <- if (!is.null(anno_colors$Driver)) scale_fill_manual(values = anno_colors$Driver) else NULL
      scale_color_driver <- if (!is.null(anno_colors$Driver)) scale_color_manual(values = anno_colors$Driver) else NULL
      project_gene_color_scale <- if (!is.null(anno_colors$Model)) {
        scale_color_manual(values = anno_colors$Model)
      } else {
        scale_color_manual(values = c("black", "grey50")) # Fallback
      }
      heatmap_anno_colors <- anno_colors
    } else {
      warning("'anno_colors' object not found in parent RMD. Using default color schemes.", call. = FALSE)
      scale_fill_driver <- NULL
      scale_color_driver <- NULL
      project_gene_color_scale <- scale_color_manual(values = c("black", "grey50")) # Fallback
      heatmap_anno_colors <- NA # pheatmap uses NA for defaults
    }
    
    # --- PCA and Heatmap Plots ---
    message("    -> Performing VST for PCA and Heatmap...")
    vsd <- vst(dds, blind = TRUE)
    
    # --- PCA Plot ---
    # CORRECTED: Changed 'Type' to 'Model' to match colData
    pcaData <- plotPCA(vsd, intgroup = c("Driver", "Host", "Model"), returnData = TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = Driver, shape = Host)) +
      geom_point(size = 4) +
      # CORRECTED: Changed 'Type' to 'Model' to match colData
      facet_wrap(~ Model) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      coord_fixed() +
      labs(title = paste("PCA on VST", model_name, "-", subset_name)) +
      scale_color_driver + # Use variable color scale
      theme_bw()
    print(pca_plot)
    
    # --- Sample-to-Sample Distance Heatmap ---
    sample_dists <- dist(t(assay(vsd)))
    sample_dist_matrix <- as.matrix(sample_dists)
    
    heatmap_plot <- pheatmap(
      sample_dist_matrix,
      clustering_distance_rows = sample_dists,
      clustering_distance_cols = sample_dists,
      main = paste("Sample-to-Sample Distance:", model_name, "-", subset_name),
      annotation_col = anno,
      annotation_colors = heatmap_anno_colors, # Use variable color scale
      silent = TRUE
    )
    print(heatmap_plot)
    
    # --- Custom Counts Distribution Boxplot ---
    counts_long <- as.data.frame(counts(dds)) %>%
      tibble::rownames_to_column("gene_id") %>%
      tidyr::pivot_longer(-gene_id, names_to = "Sample", values_to = "count") %>%
      mutate(count = count + 1)
    
    counts_long <- merge(counts_long, as.data.frame(colData(dds)), by.x = "Sample", by.y = "row.names")
    
    counts_plot <- ggplot(counts_long, aes(x = Sample, y = count, fill = Driver)) +
      geom_boxplot(outlier.shape = NA) +
      scale_y_log10() +
      scale_fill_driver + # Use variable color scale
      labs(title = paste("Counts Distribution:", model_name, "-", subset_name), y = "Raw Counts (log10 scale)") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    print(counts_plot)
    
    # --- Project-Specific Gene Boxplots ---
    if (exists("project_genes") && is.vector(project_genes) && length(project_genes) > 0) {
      for (gene_id in project_genes) {
        if (!gene_id %in% rownames(dds)) {
          warning(paste("Gene '", gene_id, "' not found in '", subset_name, "' subset. Skipping."), call. = FALSE)
          next
        }
        
        # CORRECTED: Changed 'Type' to 'Model' to match colData
        plot_data <- plotCounts(dds, gene = gene_id, intgroup = c("Driver", "Host", "Model"), returnData = TRUE)
        
        p <- ggplot(plot_data, aes(x = interaction(Driver, Host), y = count, fill = Driver)) +
          geom_boxplot(outlier.shape = NA) +
          # CORRECTED: Changed 'Type' to 'Model' to match colData
          geom_jitter(aes(color = Model), width = 0.2, size = 3) +
          scale_y_log10() +
          scale_fill_driver + # Use variable fill scale
          project_gene_color_scale + # Use variable color scale
          labs(
            title = paste("Normalized Counts for Gene:", gene_id),
            subtitle = paste("Model:", model_name, "| Subset:", subset_name),
            x = "Group (Driver.Host)",
            y = "Normalized Counts (log10 scale)"
          ) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        print(p)
      }
    }
  }
}

message("\n--- Completed 03_QC.R ---")