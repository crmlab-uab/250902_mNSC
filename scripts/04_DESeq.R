# 04_DESeq.R
# This script includes robust checks to ensure all coefficients
# required for a contrast exist within the model results before attempting to
# extract them. This prevents "subscript contains invalid names" errors
# that occur when DESeq2 drops terms after filtering.

message("--- Running 04_DESeq.R: Differential Expression Analysis ---")
library(DESeq2)
library(BiocParallel)
library(tidyverse)

# --- Main DESeq2 Analysis Loop ---
for (model_name in names(dds_models)) {
  for (subset_name in names(dds_models[[model_name]])) {
    
    message(paste(
      "\\n\\n===== Running DESeq2 analysis for:",
      model_name,
      "- Subset:",
      subset_name,
      "=====\\n"
    ))
    
    dds_to_analyze <- dds_models[[model_name]][[subset_name]]
    
    # Run DESeq2
    des_results <- DESeq(dds_to_analyze,
                         parallel = TRUE,
                         BPPARAM = MulticoreParam(nc))
    
    message("DESeq2 analysis complete for ", model_name, " - ", subset_name)
    current_results_names <- resultsNames(des_results)
    message("--> Available model coefficients:")
    print(current_results_names)
    
    # --- CONTRAST ANALYSIS ---
    
    # --- Analysis for model_1 (main effects only) ---
    if (model_name == "model_1") {
      message("\\n--- Extracting main effects for model_1 ---\\n")
      # Check if the coefficient exists before trying to extract it
      if ("Type_Tumor_vs_GBO" %in% current_results_names) {
        res <- results(des_results, name = "Type_Tumor_vs_GBO", alpha = qval, lfcThreshold = lfc)
        res_name <- paste0("model_1_", subset_name, "_Type_Tumor_vs_GBO")
        message(paste("Summary for", res_name))
        print(summary(res))
      } else {
        warning(paste("SKIPPING: Coefficient 'Type_Tumor_vs_GBO' not found in", model_name, "-", subset_name))
      }
    }
    
    # --- Analysis for model_2 (interaction effects) ---
    if (model_name == "model_2") {
      message("\\n--- Extracting interaction effects for model_2 ---\\n")
      
      # Analysis 1: Host effect within each Driver type
      message("--- Host effect analysis by Driver type ---\\n")
      res_host_by_driver <- list()
      
      for (drv in unique(samples$Driver)) {
        message(paste("\\nProcessing Driver:", drv))
        
        # Define the required coefficients for the current driver
        required_terms <- c("Host_NSG_vs_BL6")
        if (drv != "EGFRvIII") {
          required_terms <- c(required_terms, paste0("HostNSG.Driver", drv))
        }
        
        # *** THE DEFINITIVE FIX ***
        # Check if ALL required terms are present in the current model's results.
        if (all(required_terms %in% current_results_names)) {
          
          # Build the contrast based on the number of terms
          if (length(required_terms) == 1) {
            # This is the main effect for the reference driver
            res <- results(des_results, name = required_terms[1], alpha = qval, lfcThreshold = lfc)
          } else {
            # This is the sum of main effect + interaction for non-reference drivers
            res <- results(des_results, contrast = list(required_terms), alpha = qval, lfcThreshold = lfc)
          }
          
          res_name <- paste0("model_2_", subset_name, "_Host_NSG_vs_BL6_in_", drv)
          res_host_by_driver[[res_name]] <- res
          message(paste("Summary for", res_name))
          print(summary(res))
          
        } else {
          # If any term is missing, issue a detailed warning and skip.
          missing_terms <- required_terms[!required_terms %in% current_results_names]
          warning(paste0(
            "SKIPPING contrast for Driver '", drv, "' in subset '", subset_name,
            "' because the following coefficients were dropped by DESeq2: ",
            paste(missing_terms, collapse = ", ")
          ))
          next
        }
      }
      
      # Analysis 2: Driver comparisons within BL6 host
      message("\\n--- Driver comparisons in BL6 host ---\\n")
      
      dds_bl6 <- dds_to_analyze[, dds_to_analyze$Host == "BL6"]
      dds_bl6$Driver <- droplevels(dds_bl6$Driver)
      design(dds_bl6) <- formula(~ SeqBatch + Type + Driver)
      dds_bl6 <- DESeq(dds_bl6, parallel = TRUE, BPPARAM = MulticoreParam(nc))
      
      res_driver_in_bl6 <- list()
      driver_combos <- combn(levels(dds_bl6$Driver), 2)
      
      for (i in 1:ncol(driver_combos)) {
        drv1 <- driver_combos[1, i]
        drv2 <- driver_combos[2, i]
        contrast_vec <- c("Driver", drv1, drv2)
        
        res <- results(
          dds_bl6,
          contrast = contrast_vec,
          alpha = qval,
          lfcThreshold = lfc
        )
        
        res_name <- paste0("model_2_", subset_name, "_Driver_", drv1, "_vs_", drv2, "_in_BL6")
        res_driver_in_bl6[[res_name]] <- res
        
        message(paste("Summary for", res_name))
        print(summary(res))
      }
    }
  }
}

message("--- Completed 04_DESeq.R ---")