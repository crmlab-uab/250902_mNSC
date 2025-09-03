# 03_QC.R
# This script performs comprehensive quality control on DESeq2 objects,
# generating and saving a suite of standard and custom plots.
# It is designed to be run as a child script in an R Markdown environment
# and can dynamically use variables like 'anno_colors' and 'project_genes'.

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
    # This block robustly checks for an 'anno_colors' list from the parent RMD.
    # If found, it creates ggplot scale layers; otherwise, it falls back to
    # defaults and issues a single warning for the entire run.
    if (!exists("color_warning_shown")) { # Initialize a flag
      color_warning_shown <- FALSE
    }
    
    if (exists("anno_colors") && is.list(anno_colors)) {
      # --- Use Custom Colors Defined in Parent RMD ---
      scale_fill_driver <- if (!is.null(anno_colors$Driver)) scale_fill_manual(values = anno_colors$Driver) else NULL
      scale_color_driver <- if (!is.null(anno_colors$Driver)) scale_color_manual(values = anno_colors$Driver) else NULL
      
      # Specific color scale for 'Type' aesthetic in the project genes plot
      project_gene_color_scale <- if (!is.null(anno_colors$Type)) {
        scale_color_manual(values = anno_colors$Type)
      } else {
        scale_color_manual(values = c("black", "grey50")) # Fallback
      }
      
      heatmap_anno_colors <- anno_colors
      
    } else {
      # --- Fallback to Default Colors ---
      if (!color_warning_shown) {
        warning("'anno_colors' object not found in parent RMD. Using default color schemes.", call. = FALSE)
        color_warning_shown <- TRUE # Prevent repeated warnings
      }
      scale_fill_driver <- NULL
      scale_color_driver <- NULL
      project_gene_color_scale <- scale_color_manual(values = c("black", "grey50")) # Fallback
      heatmap_anno_colors <- NA # pheatmap uses NA for defaults
    }
    
    # --- PCA and Heatmap Plots ---
    message("    -> Performing VST for PCA and Heatmap...")
    vsd <- vst(dds, blind = TRUE)
    
    # --- PCA Plot ---
    pcaData <- plotPCA(vsd, intgroup = c("Driver", "Host", "Type"), returnData = TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = Driver, shape = Host)) +
      geom_point(size = 4) +
      facet_wrap(~ Type) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      coord_fixed() +
      labs(title = paste("PCA on VST", model_name, "-", subset_name)) +
      scale_color_driver + # Use variable color scale
      theme_bw()
    
    print(pca_plot)
    
    pca_filename <- file.path(dir_graph, paste0(date, "_PCA_", model_name, "_", subset_name, ".pdf"))
    ggsave(filename = pca_filename, plot = pca_plot, width = 12, height = 9)
    message(paste("    -> PCA plot saved to:", pca_filename))
    
    # --- PCA Loadings Plot for PC1 ---
    message("    -> Generating PCA loadings plot for PC1...")
    loadings_plot_pc1 <- RNAseqQC::plot_loadings(
      vsd = vsd,
      pca = pcaData,
      PC = 1,
      annotate_top_n = 5
    ) + labs(title = paste("Top 5 PC1 Loadings:", model_name, "-", subset_name))
    
    print(loadings_plot_pc1)
    
    loadings_filename_pc1 <- file.path(dir_graph, paste0(date, "_PC1_Loadings_", model_name, "_", subset_name, ".pdf"))
    ggsave(filename = loadings_filename_pc1, plot = loadings_plot_pc1, width = 10, height = 8)
    message(paste("    -> PC1 loadings plot saved to:", loadings_filename_pc1))
    
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
    
    grid::grid.newpage()
    print(heatmap_plot)
    
    heatmap_filename <- file.path(dir_graph, paste0(date, "_Heatmap_", model_name, "_", subset_name, ".pdf"))
    pdf(heatmap_filename, width = 10, height = 8)
    print(heatmap_plot)
    dev.off()
    message(paste("    -> Heatmap plot saved to:", heatmap_filename))
    
    # --- Additional QC Plots ---
    message("    -> Generating additional QC plots...")
    
    # 1. Counts Distribution Boxplot
    message("    -> Generating manual counts distribution boxplot...")
    counts_long <- as.data.frame(counts(dds)) %>%
      tibble::rownames_to_column("gene_id") %>%
      tidyr::pivot_longer(-gene_id, names_to = "Sample", values_to = "count")
    counts_long$count <- counts_long$count + 1
    counts_long <- merge(counts_long, as.data.frame(colData(dds)), by.x = "Sample", by.y = "row.names")
    
    counts_plot <- ggplot(counts_long, aes(x = Sample, y = count, fill = Driver)) +
      geom_boxplot(outlier.shape = NA) +
      scale_y_log10() +
      scale_fill_driver + # Use variable color scale
      labs(
        title = paste("Counts Distribution per Sample:", model_name, "-", subset_name),
        x = "Sample",
        y = "Raw Counts (log10 scale)"
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    print(counts_plot)
    
    counts_filename <- file.path(dir_graph, paste0(date, "_CountsBoxplot_", model_name, "_", subset_name, ".pdf"))
    ggsave(filename = counts_filename, plot = counts_plot, width = 12, height = 8)
    message(paste("    -> Counts boxplot saved to:", counts_filename))
    
    # 2. Library Complexity
    lib_comp_plot <- RNAseqQC::plot_library_complexity(dds = dds) +
      labs(title = paste("Library Complexity:", model_name, "-", subset_name))
    print(lib_comp_plot)
    
    lib_comp_filename <- file.path(dir_graph, paste0(date, "_LibComplexity_", model_name, "_", subset_name, ".pdf"))
    ggsave(filename = lib_comp_filename, plot = lib_comp_plot, width = 10, height = 8)
    message(paste("    -> Library complexity plot saved to:", lib_comp_filename))
    
    # 3. Gene Detection
    message("    -> Generating custom gene detection plot...")
    detected_genes <- colSums(counts(dds) > 10) 
    detected_df <- data.frame(Sample = names(detected_genes), Detected_Genes = detected_genes)
    detected_df <- merge(detected_df, as.data.frame(colData(dds)), by.x = "Sample", by.y = "row.names")
    
    gene_detection_plot <- ggplot(detected_df, aes(x = reorder(Sample, -Detected_Genes), y = Detected_Genes, fill = Driver)) +
      geom_col(stat = "identity") +
      scale_fill_driver + # Use variable color scale
      labs(title = paste("Gene Detection (>10 Reads):", model_name, "-", subset_name), x = "Sample", y = "Number of Detected Genes") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    print(gene_detection_plot)
    
    gene_detect_filename <- file.path(dir_graph, paste0(date, "_GeneDetection_", model_name, "_", subset_name, ".pdf"))
    ggsave(filename = gene_detect_filename, plot = gene_detection_plot, width = 12, height = 8)
    message(paste("    -> Gene detection plot saved to:", gene_detect_filename))
    
    # --- Annotation-based QC (requires GTF file) ---
    if (exists("file_gtf") && file.exists(file_gtf)) {
      message("    -> Found GTF file, preparing annotations...")
      
      gtf_data <- rtracklayer::import(file_gtf)
      genes_gtf <- gtf_data[gtf_data$type == "gene"]
      
      # Safely find common genes to prevent 'subscript contains NAs' error
      common_genes <- intersect(rownames(dds), genes_gtf$gene_id)
      message(paste("    -> Found", length(common_genes), "common genes between DDS object and GTF for annotation."))
      
      if (length(common_genes) == 0) {
        message("    -> WARNING: No matching genes found. Skipping annotation-based QC.")
      } else {
        # Create a temporary, annotated DDS object for these plots
        dds_annotated <- dds[common_genes, ]
        genes_gtf_subset <- genes_gtf[genes_gtf$gene_id %in% common_genes, ]
        rowRanges(dds_annotated) <- genes_gtf_subset[match(rownames(dds_annotated), genes_gtf_subset$gene_id), ]
        
        message("    -> Generating annotation-based QC plots...")
        
        # 4. Gene Biotypes
        message("    -> Plotting gene biotypes...")
        biotype_plot <- RNAseqQC::plot_biotype(dds = dds_annotated) +
          labs(title = paste("Gene Biotype Distribution:", model_name, "-", subset_name)) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        print(biotype_plot)
        
        biotype_filename <- file.path(dir_graph, paste0(date, "_Biotypes_", model_name, "_", subset_name, ".pdf"))
        ggsave(filename = biotype_filename, plot = biotype_plot, width = 12, height = 8)
        message(paste("    -> Biotypes plot saved to:", biotype_filename))
        
        # 5. Chromosomal Expression
        message("    -> Plotting chromosomal expression...")
        chr_exp_plot <- RNAseqQC::plot_chromosome_expression(dds = dds_annotated) +
          labs(title = paste("Chromosomal Expression:", model_name, "-", subset_name))
        print(chr_exp_plot)
        
        chrexp_filename <- file.path(dir_graph, paste0(date, "_ChrExpression_", model_name, "_", subset_name, ".pdf"))
        ggsave(filename = chrexp_filename, plot = chr_exp_plot, width = 12, height = 8)
        message(paste("    -> Chromosomal expression plot saved to:", chrexp_filename))
      }
    } else {
      message("    -> GTF file not found or path not set, skipping annotation-based QC.")
    }
    
    # --- Project-Specific Gene Boxplots ---
    if (exists("project_genes") && is.vector(project_genes) && length(project_genes) > 0) {
      message("    -> Generating boxplots for project-specific genes...")
      
      project_genes_filename <- file.path(dir_graph, paste0(date, "_Project_Gene_Boxplots_", model_name, "_", subset_name, ".pdf"))
      
      pdf(project_genes_filename, width = 11, height = 8.5)
      
      for (gene_id in project_genes) {
        if (!gene_id %in% rownames(dds)) {
          warning(paste("Gene '", gene_id, "' not found in '", subset_name, "' subset. Skipping."), call. = FALSE)
          next
        }
        
        plot_data <- plotCounts(dds, gene = gene_id, intgroup = c("Driver", "Host", "Type"), returnData = TRUE)
        
        p <- ggplot(plot_data, aes(x = interaction(Driver, Host), y = count, fill = Driver)) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(aes(color = Type), width = 0.2, size = 3) +
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
      
      dev.off()
      message(paste("    -> Project-specific gene plots saved to:", project_genes_filename))
    }
  }
}

message("\n--- Completed 03_QC.R ---")