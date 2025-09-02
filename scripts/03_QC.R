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

for (model_name in names(dds_models)) {
  for (subset_name in names(dds_models[[model_name]])) {
    message(paste("...Generating QC plots for:", model_name, "-", subset_name))
    
    # Get the specific DESeqDataSet object for the current model and subset
    dds <- dds_models[[model_name]][[subset_name]]
    
    # --- PCA and Heatmap Plots ---
    # These plots require variance stabilized data.
    message("    -> Performing VST for PCA and Heatmap...")
    vsd <- vst(dds, blind = TRUE)
    
    # --- PCA Plot ---
    # This robust method extracts data from plotPCA and builds the ggplot manually
    # to ensure full control over aesthetics like color.
    
    # 1. Get the data used for plotting instead of the plot object itself
    pcaData <- plotPCA(vsd,
                       intgroup = c("Driver", "Host", "Type"),
                       returnData = TRUE)
    
    # 2. Get the percent variance explained for the axes labels
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    # 3. Build the plot manually using ggplot, giving us full control
    pca_plot <- ggplot(pcaData, aes(
      x = PC1,
      y = PC2,
      color = Driver,
      shape = Host
    )) +
      geom_point(size = 4) +
      # Create subplots for each 'Type'
      facet_wrap(~ Type) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      coord_fixed() +
      labs(title = paste("PCA on VST", model_name, "-", subset_name)) +
      scale_color_manual(values = anno_colors$Driver) +
      theme_bw()
    
    # Print plot to the knitted document
    print(pca_plot)
    
    # Save the PCA plot as a dated PDF
    pca_filename <- file.path(dir_graph,
                              paste0(date, "_PCA_", model_name, "_", subset_name, ".pdf"))
    ggsave(
      filename = pca_filename,
      plot = pca_plot,
      width = 10,
      height = 8
    )
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
    
    # This line ensures the plot is drawn on a fresh canvas
    grid::grid.newpage()
    
    # Print plot to the knitted document
    print(heatmap_plot)
    
    # Save the heatmap plot as a dated PDF
    heatmap_filename <- file.path(dir_graph,
                                  paste0(date, "_Heatmap_", model_name, "_", subset_name, ".pdf"))
    pdf(heatmap_filename, width = 10, height = 8)
    print(heatmap_plot)
    dev.off()
    message(paste("    -> Heatmap plot saved to:", heatmap_filename))
    
    # --- Annotation-based QC (requires GTF file) ---
    if (exists("file_gtf") && file.exists(file_gtf)) {
      message("    -> Found GTF file, generating annotation-based QC plots...")
      
      # 4. Gene Biotypes
      # This plot summarizes the percentage of reads mapping to different gene types.
      message("    -> Plotting gene biotypes...")
      biotype_plot <- RNAseqQC::plot_biotypes(
        dds = dds,
        gtf = file_gtf
      ) +
        labs(title = paste("Gene Biotype Distribution:", model_name, "-", subset_name)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      print(biotype_plot)
      
      # 5. Chromosomal Expression
      # This plot shows the relative expression level for each chromosome.
      message("    -> Plotting chromosomal expression...")
      chr_exp_plot <- RNAseqQC::plot_chromosome_expression(
        dds = dds,
        gtf = file_gtf
      ) +
        labs(title = paste("Chromosomal Expression:", model_name, "-", subset_name))
      
      print(chr_exp_plot)
      
    } else {
      message("    -> GTF file not found or path not set, skipping annotation-based QC.")
    }
  }
}

# --- RNAseqQC plots ---
if (requireNamespace("RNAseqQC", quietly = TRUE)) {
  message("    -> Generating RNAseqQC plots...")
  
  list_entry_name <- paste(model_name, subset_name, sep = "_")
  
  # This line runs all standard RNAseqQC analyses and stores the plots
  qc_plots[[list_entry_name]] <- RNAseqQC::run_RNAseqQC(
    dds,
    sample_info = samples,
    factors = c("Driver", "Host", "Type")
  )
  
  # --- Print requested QC plots to the HTML document ---
  message("    -> Printing selected QC plots to report...")
  
  # 1. Total Sample Counts
  # This plot shows the distribution of raw read counts for each sample.
  print(
    qc_plots[[list_entry_name]][["Counts_Boxplot"]] +
      labs(title = paste("Total Read Counts:", model_name, "-", subset_name)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  
  # 2. Library Complexity
  # This shows the diversity of transcripts captured. More complex libraries have curves that rise quickly.
  print(
    qc_plots[[list_entry_name]][["Library_Complexity"]] +
      labs(title = paste("Library Complexity:", model_name, "-", subset_name))
  )
  
  # 3. Gene Detection (Custom Plot)
  # This bar chart shows the number of genes with more than 10 reads in each sample.
  message("    -> Generating custom gene detection plot...")
  detected_genes <- colSums(counts(dds) > 10) 
  detected_df <- data.frame(Sample = names(detected_genes),
                            Detected_Genes = detected_genes)
  detected_df <- merge(detected_df, as.data.frame(colData(dds)), by.x = "Sample", by.y = "row.names")
  
  gene_detection_plot <- ggplot(detected_df, aes(x = reorder(Sample, -Detected_Genes), y = Detected_Genes, fill = Driver)) +
    geom_col(stat = "identity") +
    scale_fill_manual(values = anno_colors$Driver) +
    labs(title = paste("Gene Detection (>10 Reads):", model_name, "-", subset_name),
         x = "Sample",
         y = "Number of Detected Genes") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  print(gene_detection_plot)
}

message("\n--- Completed 03_QC.R ---")