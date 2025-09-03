# 04_DESeq_runner.R
# Contains a reusable function to perform a complete DESeq2 analysis.

#' Run a complete DESeq2 differential expression analysis workflow.
#'
#' @param samples_df A data frame with the sample metadata.
#' @param files_sf_vec A named character vector of paths to Salmon 'quant.sf' files.
#' @param design_formula A formula specifying the design of the experiment.
#' @param tx2gene_map A data frame with transcript_id and gene_id columns.
#' @param comparisons_list A list defining the pairwise comparisons to perform.
#' @param analysis_name A string for naming output files.
#' @param gene_map A data frame for mapping gene IDs to names for plot labels.
#' @param dir_graphs The path to the directory where graphs should be saved.
#' @param qval The adjusted p-value (padj) cutoff.
#' @param lfc The log2 fold change cutoff.
#' @return The final DESeqDataSet object after running DESeq().

run_deseq_analysis <- function(samples_df, files_sf_vec, design_formula, tx2gene_map,
                               comparisons_list, analysis_name, gene_map, dir_graphs,
                               qval = 0.05, lfc = 0.58) {
  message(paste("\n--- Starting DESeq2 Analysis for:", analysis_name, "---"))
  
  # 1. Import data and create DESeqDataSet
  txi_data <- tximport(files_sf_vec, type = "salmon", tx2gene = tx2gene_map, ignoreTxVersion = TRUE)
  dds <- DESeqDataSetFromTximport(txi_data, colData = samples_df, design = design_formula)
  
  # 2. Pre-filter and run DESeq
  keep <- rowSums(counts(dds) >= 10) >= min(table(samples_df$Driver))
  dds_filt <- dds[keep, ]
  dds_results <- DESeq(dds_filt)
  
  message("\n--- Performing specified pairwise comparisons ---")
  
  # 3. Loop through comparisons, get results, and plot
  for (comp in comparisons_list) {
    res_name <- paste0(analysis_name, "_", comp$factor, "_", comp$group1, "_vs_", comp$group2)
    message(paste("...Running comparison:", res_name))
    
    res <- results(dds_results, contrast = c(comp$factor, comp$group1, comp$group2), alpha = qval, lfcThreshold = lfc)
    
    # Merge results with gene names for labeling
    res_df <- as.data.frame(res) %>% 
      tibble::rownames_to_column("gene_id") %>% 
      dplyr::left_join(gene_map, by = "gene_id")
    
    # Create Volcano Plot
    volcano_plot <- EnhancedVolcano(
      res_df,
      lab = res_df$gene_name,
      x = 'log2FoldChange',
      y = 'padj',
      pCutoff = qval,
      FCcutoff = lfc,
      title = res_name,
      subtitle = paste("Analysis:", analysis_name)
    )
    
    # Save the plot using the provided directory path
    output_filename <- file.path(dir_graphs, paste0(res_name, "_volcano.png"))
    ggsave(filename = output_filename, plot = volcano_plot, width = 12, height = 10, dpi = 300)
    message(paste("...Saved plot to:", output_filename))
  }
  
  message(paste("\n--- Completed analysis:", analysis_name, "---"))
  return(dds_results)
}