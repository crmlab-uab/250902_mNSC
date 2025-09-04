# 02_data_prep.R
# Contains a function to prepare DESeq2 data objects.

#' Prepare DESeq2 Datasets
#'
#' @param files_sf Named vector of paths to Salmon quant.sf files.
#' @param samples_df Data frame with sample metadata.
#' @param tx2gene_maps A list of tx2gene data frames (e.g., list(all = ..., pcg = ...)).
#' @param design_formulas A named list of design formulas for DESeq2.
#' @param filter_cutoff Minimum count for the low-count filter.
#' @return A list containing 'dds_models' and 'vsd_models'.

prepare_deseq_objects <- function(files_sf,
                                  samples_df,
                                  tx2gene_maps,
                                  design_formulas,
                                  filter_cutoff) {
  message("--- Running prepare_deseq_objects function ---")
  
  dds_models <- list()
  vsd_models <- list()
  
  for (model_name in names(design_formulas)) {
    message(paste("\n... processing model:", model_name))
    design_formula <- design_formulas[[model_name]]
    
    # --- SF File Logic ---
    # The 'files_sf' vector (created in the Rmd) is used here by tximport.
    # The function itself doesn't know the file structure; it just receives the paths.
    txi_all <- tximport(
      files_sf,
      type = "salmon",
      tx2gene = tx2gene_maps$all,
      ignoreTxVersion = TRUE
    ) #
    dds_all <- DESeqDataSetFromTximport(txi_all, colData = samples_df, design = design_formula) #
    
    txi_pcg <- tximport(
      files_sf,
      type = "salmon",
      tx2gene = tx2gene_maps$pcg,
      ignoreTxVersion = TRUE
    ) #
    dds_pcg <- DESeqDataSetFromTximport(txi_pcg, colData = samples_df, design = design_formula) #
    # --- END SF File Logic ---
    
    # Apply a consistent low-count filter
    smallest_group_size <- min(table(samples_df$Driver))
    keep_all <- rowSums(counts(dds_all) >= filter_cutoff) >= smallest_group_size #
    dds_all_filt <- dds_all[keep_all, ] #
    
    keep_pcg <- rowSums(counts(dds_pcg) >= filter_cutoff) >= smallest_group_size #
    dds_pcg_filt <- dds_pcg[keep_pcg, ] #
    
    current_model_dds_list <- list(all = dds_all,
                                   filt = dds_all_filt,
                                   pcg = dds_pcg_filt)
    
    # Estimate size factors and perform VST
    dds_models[[model_name]] <- lapply(current_model_dds_list, estimateSizeFactors) #
    message(paste("... applying VST to filtered data for", model_name))
    vsd_models[[model_name]] <- vst(dds_models[[model_name]]$filt, blind = TRUE) #
  }
  
  message("\n--- Completed DESeq2 object preparation ---")
  return(list(dds_models = dds_models, vsd_models = vsd_models))
}