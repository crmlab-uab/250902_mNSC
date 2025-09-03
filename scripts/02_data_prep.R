# 02_data_prep.R
message("--- Running 02_data_prep.R: Preparing DESeq2 Objects from new tx2gene maps ---")
library(tximport)
library(DESeq2)

# --- Create a list of the tx2gene maps to process ---
tx2gene_maps <- list(
  all = tx2gene_all,
  pcg = tx2gene_pcg
)

# Initialize master lists to hold the results for each model
dds_models <- list()
vsd_models <- list()

message("Creating and processing DESeqDataSet for each model and gene subset...")

# Loop through each of the main DESeq2 statistical models
for (model_name in names(designs)) {
  message(paste("\n... processing model:", model_name))
  design_formula <- designs[[model_name]]
  
  # Create an inner list to hold the dds objects for this model (all, filt, pcg)
  current_model_dds_list <- list()
  
  # --- Step 1: Import counts using the 'all' and 'pcg' maps ---
  # CORRECTED: Added ignoreTxVersion = TRUE to both tximport calls
  txi_all <- tximport(files_sf, type = "salmon", tx2gene = tx2gene_maps$all, ignoreTxVersion = TRUE)
  dds_all <- DESeqDataSetFromTximport(txi_all, colData = samples, design = design_formula)
  
  txi_pcg <- tximport(files_sf, type = "salmon", tx2gene = tx2gene_maps$pcg, ignoreTxVersion = TRUE)
  dds_pcg <- DESeqDataSetFromTximport(txi_pcg, colData = samples, design = design_formula)
  
  # --- Step 2: Apply a low-count filter to the 'all' and 'pcg' sets ---
  # Keep genes with at least 'filt' counts in the smallest group of samples
  smallest_group_size <- 3 
  keep_all <- rowSums(counts(dds_all) >= filt) >= smallest_group_size
  dds_all_filt <- dds_all[keep_all, ]
  
  keep_pcg <- rowSums(counts(dds_pcg) >= filt) >= smallest_group_size
  dds_pcg_filt <- dds_pcg[keep_pcg, ]
  
  # --- Step 3: Assemble the list of datasets for this model ---
  # Note: we use 'dds_all_filt' as the 'filt' subset for QC and analysis
  current_model_dds_list <- list(
    all = dds_all,
    filt = dds_all_filt,
    pcg = dds_pcg_filt # This is now the refined AND filtered pcg set
  )
  
  # Estimate size factors for each dataset
  current_model_dds_list <- lapply(current_model_dds_list, estimateSizeFactors)
  
  # Add the processed list to the master list
  dds_models[[model_name]] <- current_model_dds_list
  
  # Perform VST on the filtered data for this model
  message(paste("... applying VST to", model_name))
  vsd_models[[model_name]] <- vst(dds_models[[model_name]]$filt, blind = TRUE)
}

message("\n--- Completed 02_data_prep.R ---")