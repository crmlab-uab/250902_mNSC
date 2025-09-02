# 02_data_prep.R
message("--- Running 02_data_prep.R: Preparing DESeq2 Objects ---")
library(tximport)
library(DESeq2)

# Import salmon counts
txi <- tximport(
  files_sf,
  type = "salmon",
  tx2gene = tx_gene_symbol,
  ignoreAfterBar = TRUE
)

# Initialize master lists to hold the results for each model
dds_models <- list()
vsd_models <- list()

message("Creating and processing DESeqDataSet for each model...")

# Loop through each design formula provided by the parent RMD
for (model_name in names(designs)) {
  message(paste("... processing model:", model_name))
  design_formula <- designs[[model_name]]
  
  # Create the initial full DESeqDataSet
  # Corrected function call to include `colData = samples`
  dds_all <- DESeqDataSetFromTximport(txi, colData = samples, design = design_formula)
  
  # Pre-filtering: Keep genes with at least 'filt' counts in the smallest group of samples
  smallest_group_size <- 3
  keep <- rowSums(counts(dds_all) >= filt) >= smallest_group_size
  dds_filt <- dds_all[keep, ]
  
  # Filter for protein-coding genes
  dds_pcg <- dds_filt[rownames(dds_filt) %in% pcg, ]
  
  # Create a named list of the dds objects for the current model
  current_model_dds_list <- list(
    all = dds_all,
    filt = dds_filt,
    pcg = dds_pcg
  )
  
  # Estimate size factors for each dataset within the current model
  current_model_dds_list <- lapply(current_model_dds_list, estimateSizeFactors)
  
  # Add the processed list of dds objects to the master list
  dds_models[[model_name]] <- current_model_dds_list
  
  # Perform Variance Stabilizing Transformation (VST) on the filtered data for this model
  message(paste("... applying VST to", model_name))
  vsd_models[[model_name]] <- vst(dds_models[[model_name]]$filt, blind = TRUE)
}

# # --- Save Data Objects ---
# message("Saving processed data objects...")
# saveRDS(dds_models, file = paste0(dir_output, date, "_dds_models_list.rds"))
# saveRDS(vsd_models, file = paste0(dir_output, date, "_vsd_models_list.rds"))

message("Data preparation complete.")

