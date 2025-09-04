# 04_DESeq_runner.R

#' 1. Run the core (and slow) DESeq2 analysis.
#' This function takes raw data and returns a processed DESeqDataSet object.
run_deseq_core <- function(samples_df,
                           files_sf_vec,
                           design_formula,
                           tx2gene_map) {
  message(paste0(
    "\n--- Running Core DESeq Analysis for: ",
    deparse(design_formula),
    " ---"
  ))
  
  # Import data
  txi_data <- tximport(
    files_sf_vec,
    type = "salmon",
    tx2gene = tx2gene_map,
    ignoreTxVersion = TRUE
  )
  dds <- DESeqDataSetFromTximport(txi_data, colData = samples_df, design = design_formula)
  
  # Pre-filter for low counts to speed up the analysis
  # Keep genes where at least the number of samples in the smallest group have a count of 10 or more
  smallest_group_size <- min(table(samples_df$Driver))
  keep <- rowSums(counts(dds) >= filter_min_count) >= smallest_group_size
  dds_filt <- dds[keep, ]
  
  # Run the main DESeq2 function
  dds_processed <- DESeq(dds_filt)
  
  message("--- Core analysis complete. ---")
  return(dds_processed)
}

#' 2. Extract results and generate volcano plots for multiple comparisons.
#' This function takes a processed DESeqDataSet object and is very fast.
extract_comparisons <- function(dds_processed,
                                comparisons_list,
                                analysis_name,
                                gene_map,
                                dir_graphs,
                                qval = 0.05,
                                lfc = 1.0) {
  message(paste("\n--- Extracting comparisons for:", analysis_name, "---"))
  
  results_list <- list()
  
  for (comp in comparisons_list) {
    res_name <- paste0(analysis_name,
                       "_",
                       comp$factor,
                       "_",
                       comp$group1,
                       "_vs_",
                       comp$group2)
    message(paste("...Extracting results for:", res_name))
    
    # Extract results from the pre-computed object
    res <- results(
      dds_processed,
      contrast = c(comp$factor, comp$group1, comp$group2),
      alpha = qval,
      lfcThreshold = lfc
    )
    
    # Merge with gene names and format
    res_df <- as.data.frame(res) %>%
      tibble::rownames_to_column("gene_id") %>%
      dplyr::left_join(gene_map, by = "gene_id")
    
    results_list[[res_name]] <- res_df
    
    # Create and save Volcano Plot
    volcano_plot <- EnhancedVolcano(
      res_df,
      lab = res_df$gene_name,
      x = 'log2FoldChange',
      y = 'padj',
      pCutoff = qval,
      FCcutoff = lfc,
      title = res_name
    )
    
    output_filename <- here::here(dir_graphs, paste0(res_name, "_volcano.png"))
    ggsave(
      filename = output_filename,
      plot = volcano_plot,
      width = 12,
      height = 10,
      dpi = 300
    )
    message(paste("...Saved volcano plot to:", output_filename))
  }
  
  message("--- All comparisons extracted. ---")
  return(results_list)
}