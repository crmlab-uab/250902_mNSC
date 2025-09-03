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
      scale_fill_driver <- if (!is.null(anno_colors$Driver))
        scale_fill_manual(values = anno_colors$Driver)
      else
        NULL
      scale_color_driver <- if (!is.null(anno_colors$Driver))
        scale_color_manual(values = anno_colors$Driver)
      else
        NULL
      project_gene_color_scale <- if (!is.null(anno_colors$Model)) {
        scale_color_manual(values = anno_colors$Model)
      } else {
        scale_color_manual(values = c("black", "grey50")) # Fallback
      }
      heatmap_anno_colors <- anno_colors
    } else {
      warning(
        "'anno_colors' object not found in parent RMD. Using default color schemes.",
        call. = FALSE
      )
      scale_fill_driver <- NULL
      scale_color_driver <- NULL
      project_gene_color_scale <- scale_color_manual(values = c("black", "grey50")) # Fallback
      heatmap_anno_colors <- NA # pheatmap uses NA for defaults
    }
    
    # --- PCA and Heatmap Plots ---
    message("    -> Performing VST for PCA and Heatmap...")
    vsd <- vst(dds, blind = TRUE)
    
    # --- PCA Plot ---
    pcaData <- plotPCA(vsd,
                       intgroup = c("Driver", "Host", "Model"),
                       returnData = TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    pca_plot <- ggplot(pcaData, aes(
      x = PC1,
      y = PC2,
      color = Driver,
      shape = Host
    )) +
      geom_point(size = 4) +
      facet_wrap( ~ Model) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      coord_fixed() +
      labs(title = paste("PCA on VST", model_name, "-", subset_name)) +
      scale_color_driver +
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
      annotation_colors = heatmap_anno_colors,
      silent = TRUE
    )
    
    grid::grid.newpage()
    grid::grid.draw(heatmap_plot$gtable)
    
    # --- Custom Counts Distribution Boxplot ---
    counts_long <- as.data.frame(counts(dds)) %>%
      tibble::rownames_to_column("gene_id") %>%
      tidyr::pivot_longer(-gene_id, names_to = "Sample", values_to = "count") %>%
      mutate(count = count + 1)
    
    counts_long <- merge(counts_long,
                         as.data.frame(colData(dds)),
                         by.x = "Sample",
                         by.y = "row.names")
    
    counts_plot <- ggplot(counts_long, aes(x = Sample, y = count, fill = Driver)) +
      geom_boxplot(outlier.shape = NA) +
      scale_y_log10() +
      scale_fill_driver +
      labs(
        title = paste("Counts Distribution:", model_name, "-", subset_name),
        y = "Raw Counts (log10 scale)"
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ))
    print(counts_plot)
    
    # --- [UPDATED] Project-Specific Gene Boxplots (Faceted by Model and Host) ---
    if (exists("project_genes") &&
        is.vector(project_genes) && length(project_genes) > 0) {
      message("    -> Generating faceted boxplots for project-specific genes...")
      
      for (gene_name_to_plot in project_genes) {
        # Find the corresponding Ensembl ID from the gene_map
        gene_id_to_plot <- gene_map$gene_id[gene_map$gene_name == gene_name_to_plot]
        
        # Ensure we found a unique, valid ID before plotting
        if (length(gene_id_to_plot) == 1 &&
            gene_id_to_plot %in% rownames(dds)) {
          # Get normalized counts for the specific gene
          plot_data <- plotCounts(
            dds,
            gene = gene_id_to_plot,
            intgroup = c("Driver", "Host", "Model"),
            returnData = TRUE
          )
          
          # Create the faceted plot
          p <- ggplot(plot_data, aes(
            x = Driver,
            y = count,
            fill = Driver
          )) +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter(width = 0.2,
                        alpha = 0.7,
                        shape = 16) +
            # Facet by Model and Host to create subsets
            facet_wrap( ~ Model + Host, scales = "free_y") +
            scale_y_log10() +
            # Use pre-defined colors for consistency
            scale_fill_driver +
            labs(
              # Use the gene NAME for the title
              title = paste("Normalized Counts for Gene:", gene_name_to_plot),
              subtitle = paste(
                "Statistical Model:",
                model_name,
                "| Data Subset:",
                subset_name
              ),
              x = "Driver",
              y = "Normalized Counts (log10 scale)"
            ) +
            theme_bw() +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1),
              strip.background = element_rect(fill = "grey90", color = "black"),
              legend.position = "bottom"
            )
          
          print(p)
          
        } else {
          warning(
            paste(
              "Gene '",
              gene_name_to_plot,
              "' not found or not unique in the current data subset ('",
              subset_name,
              "'). Skipping plot."
            ),
            call. = FALSE
          )
        }
      }
    }
  }
}

message("\n--- Completed 03_QC.R ---")