# 03_QC.R
message("--- Running 03_QC.R: Comprehensive Quality Control ---\n")
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(DESeq2)

# Ensure the dds_models object from 02_data_prep.R exists
if (!exists("dds_models")) {
  stop("dds_models object not found. Please ensure 02_data_prep.R has been run.")
}

# --- Comprehensive QC Plotting Loop ---
# This script now uses nested loops to generate QC plots for:
# 1. Each model design ('model_1', 'model_2')
# 2. Each data subset ('all', 'filt', 'pcg')

# Initialize the qc_plots list to store RNAseqQC outputs
qc_plots <- list()

for (model_name in names(dds_models)) {
  for (subset_name in names(dds_models[[model_name]])) {
    
    message(paste("...Generating QC plots for:", model_name, "-", subset_name))
    
    # Get the specific DESeqDataSet object for the current model and subset
    dds <- dds_models[[model_name]][[subset_name]]
    
    # --- RNAseqQC plots ---
    if (requireNamespace("RNAseqQC", quietly = TRUE)) {
      message("    -> Generating RNAseqQC plots...")
      
      list_entry_name <- paste(model_name, subset_name, sep = "_")
      
      qc_plots[[list_entry_name]] <- RNAseqQC::run_RNAseqQC(
        dds,
        sample_info = samples,
        factors = c("Driver", "Host", "Type")
      )
    }
    
    # --- PCA and Heatmap Plots ---
    # These plots require variance stabilized data.
    message("    -> Performing VST for PCA and Heatmap...")
    vsd <- vst(dds, blind = TRUE)
    
    # PCA Plot
    pca_plot <- plotPCA(vsd, intgroup = c("Driver", "Host", "Type")) +
      labs(title = paste("PCA on VST", model_name, "-", subset_name)) +
      geom_point(size = 4) +
      theme_bw()
    
    # Print plot to the knitted document
    print(pca_plot)
    
    # Save the PCA plot as a dated PDF
    pca_filename <- file.path(dir_qc, paste0(date, "_PCA_", model_name, "_", subset_name, ".pdf"))
    ggsave(filename = pca_filename, plot = pca_plot, width = 10, height = 8)
    message(paste("    -> PCA plot saved to:", pca_filename))
    
    
    # Sample-to-Sample Distance Heatmap
    sample_dists <- dist(t(assay(vsd)))
    sample_dist_matrix <- as.matrix(sample_dists)
    
    heatmap_plot <- pheatmap(
      sample_dist_matrix,
      clustering_distance_rows = sample_dists,
      clustering_distance_cols = sample_dists,
      main = paste("Sample-to-Sample Distance:", model_name, "-", subset_name),
      annotation_col = anno,
      annotation_colors = anno_colors,
      silent = TRUE 
    )
    
    # Print plot to the knitted document
    print(heatmap_plot)
    
    # Save the heatmap plot as a dated PDF
    heatmap_filename <- file.path(dir_qc, paste0(date, "_Heatmap_", model_name, "_", subset_name, ".pdf"))
    pdf(heatmap_filename, width = 10, height = 8)
    print(heatmap_plot)
    dev.off()
    message(paste("    -> Heatmap plot saved to:", heatmap_filename))
    
  }
}

# --- Save RNAseqQC plots from the list to PDF files ---
# This part of the script was already correct.
if (exists("qc_plots") && length(qc_plots) > 0) {
  message("\n...Saving RNAseqQC plots to files...")
  for (name in names(qc_plots)) { # e.g., name is "model_1_all"
    for (plot_name in names(qc_plots[[name]])) {
      pdf_file <- file.path(dir_qc, paste0(date, "_RNAseqQC_", name, "_", plot_name, ".pdf"))
      pdf(pdf_file, width = 10, height = 8)
      print(qc_plots[[name]][[plot_name]])
      dev.off()
    }
  }
  message("...RNAseqQC plots saved successfully.")
}

message("\n--- Completed 03_QC.R ---")