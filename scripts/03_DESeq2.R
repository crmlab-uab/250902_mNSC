# =============================================================================
#
# 03_DESeq2.R: DESeq2 ANALYSIS FUNCTIONS
#
# Description:
# This script contains reusable functions for DESeq2 analysis.
# All project-specific variables are passed as arguments via the cfg list.
#
# =============================================================================

#' Run the Core DESeq2 Analysis (Modular)
#'
#' @param samples_df A data frame with sample metadata.
#' @param files_sf_vec A named character vector of paths to 'quant.sf' files.
#' @param design_formula A formula specifying the experimental design.
#' @param tx2gene_map A data frame with `transcript_id` and `gene_id` columns.
#' @param filter_by_factor character. The column name in `samples_df` for the low-count filter.
#' @param filter_min_count numeric. The minimum count for the low-count filter.
#'
#' @return A processed `DESeqDataSet` object after running `DESeq()`.
#'
run_deseq_core <- function(samples_df,
                           files_sf_vec,
                           design_formula,
                           tx2gene_map,
                           filter_by_factor,
                           filter_min_count) {
  message(paste0(
    "\n--- Running Core DESeq Analysis for: ",
    deparse(design_formula),
    " ---"
  ))
  txi_data <- tximport(
    files_sf_vec,
    type = "salmon",
    tx2gene = tx2gene_map,
    ignoreTxVersion = TRUE
  )
  dds <- DESeqDataSetFromTximport(txi_data, colData = samples_df, design = design_formula)
  
  if (!filter_by_factor %in% names(samples_df)) {
    stop(
      paste(
        "The filter_by_factor '",
        filter_by_factor,
        "' is not a column in the samples data frame."
      )
    )
  }
  
  smallest_group_size <- min(table(samples_df[[filter_by_factor]]))
  keep <- rowSums(counts(dds) >= filter_min_count) >= smallest_group_size
  dds_filt <- dds[keep, ]
  
  dds_processed <- DESeq(dds_filt)
  message("--- Core analysis complete. ---")
  return(dds_processed)
}


#' Extract Pairwise Comparisons, Save Results, and Generate Plots
#'
#' @param dds_processed A `DESeqDataSet` object that has been run through `DESeq()`.
#' @param comparisons_list A list defining the pairwise comparisons to perform.
#' @param analysis_name A string for naming output files.
#' @param gene_map A data frame for mapping gene IDs to names for plot labels.
#' @param cfg The project configuration list.
#'
#' @return A named list of data frames containing DESeq2 results.
#'
extract_comparisons <- function(dds_processed,
                                comparisons_list,
                                analysis_name,
                                gene_map,
                                cfg) {
  message(paste(
    "\n--- Extracting comparisons and saving results for:",
    analysis_name,
    "---"
  ))
  results_list <- list()
  
  for (comp in comparisons_list) {
    res_name <- paste0(analysis_name,
                       "_",
                       comp$factor,
                       "_",
                       comp$group1,
                       "_vs_",
                       comp$group2)
    message(paste("...Processing comparison:", res_name))
    
    res <- try(results(
      dds_processed,
      contrast = c(comp$factor, comp$group1, comp$group2),
      alpha = cfg$qval_threshold,
      lfcThreshold = cfg$lfc_threshold
    ),
    silent = TRUE)
    
    if (inherits(res, "try-error")) {
      message(
        paste(
          "......Skipping comparison:",
          res_name,
          "(likely redundant or not estimable)."
        )
      )
      next
    }
    
    res_df <- as.data.frame(res) %>%
      tibble::rownames_to_column("gene_id") %>%
      dplyr::left_join(gene_map, by = "gene_id") %>%
      dplyr::arrange(padj)
    
    results_list[[res_name]] <- res_df
    
    # Save full results table
    csv_filename <- here::here(cfg$dir_results_csv,
                               paste0(cfg$date, "_", res_name, "_results.csv"))
    write.csv(res_df, file = csv_filename, row.names = FALSE)
    
    # Create and save Volcano Plot
    volcano_plot <- EnhancedVolcano(
      res_df,
      lab = res_df$gene_name,
      x = 'log2FoldChange',
      y = 'padj',
      pCutoff = cfg$qval_threshold,
      FCcutoff = cfg$lfc_threshold,
      title = res_name,
      subtitle = analysis_name
    )
    plot_filename <- here::here(cfg$dir_graphs,
                                paste0(cfg$date, "_", res_name, "_volcano.png"))
    ggsave(
      filename = plot_filename,
      plot = volcano_plot,
      width = 12,
      height = 10,
      dpi = 300
    )
  }
  
  message("--- All valid comparisons processed and saved. ---")
  return(results_list)
}