# 06_DESeq_BL6_subset.R
# This script performs a DESeq2 analysis on the pre-defined 'samples_bl6' subset.

message("--- Running 06_DESeq_BL6_subset.R: Analysis on 'BL6' Host Subset ---")
library(DESeq2)
library(BiocParallel)
library(EnhancedVolcano)
library(tidyverse)
library(tibble)

# --- The script now expects 'samples_subset' and 'files_sf_subset' to exist ---
samples_bl6 <- samples_subset
files_sf_bl6 <- files_sf_subset

# --- 1. Main Effects Analysis (Model vs. Driver across all BL6 samples) ---
message("--- Part 1: Running Main Effects Analysis for BL6 Subset ---")
message("...Importing counts and creating DESeqDataSet for BL6 subset...")
txi_bl6 <- tximport(
  files_sf_bl6,
  type = "salmon",
  tx2gene = tx2gene_pcg,
  ignoreTxVersion = TRUE
)
dds_bl6 <- DESeqDataSetFromTximport(txi_bl6,
                                    colData = samples_bl6,
                                    design = ~ SeqBatch + Model + Driver)

# Filter for low counts
keep <- rowSums(counts(dds_bl6) >= filt) >= 3
dds_bl6_filt <- dds_bl6[keep, ]

# Run DESeq
message("...Running DESeq()...")
dds_bl6_results <- DESeq(dds_bl6_filt)
message("DESeq2 analysis for BL6 subset complete.")

# --- 2. Extract Pairwise DEG Results for Main Effects ---
# Helper function to run and plot a comparison
run_comparison <- function(dds_results,
                           factor,
                           group1,
                           group2,
                           subset_name) {
  res_name <- paste0(factor, "_", group1, "_vs_", group2, "_in_", subset_name)
  message("\n--- Running comparison: ", res_name, " ---")
  
  res <- results(
    dds_results,
    contrast = c(factor, group1, group2),
    alpha = qval,
    lfcThreshold = lfc
  )
  print(summary(res))
  
  res_df <- as.data.frame(res) %>% rownames_to_column("gene_id") %>% dplyr::left_join(gene_map, by = "gene_id")
  labels <- ifelse(!is.na(res_df$gene_name), res_df$gene_name, res_df$gene_id)
  
  volcano_plot <- EnhancedVolcano(
    res,
    lab = labels,
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = qval,
    FCcutoff = lfc,
    title = res_name,
    subtitle = paste("Comparison in", subset_name, "subset")
  )
  print(volcano_plot)
}

# --- Programmatically loop through all Model comparisons ---
message("\n--- Processing Main Effect: Model Comparisons in BL6 Host ---")
model_combos <- combn(levels(samples_bl6$Model), 2)
for (i in 1:ncol(model_combos)) {
  run_comparison(dds_bl6_results,
                 "Model",
                 model_combos[1, i],
                 model_combos[2, i],
                 "BL6_Host")
}

# --- Programmatically loop through all Driver comparisons ---
message("\n--- Processing Main Effect: Driver Comparisons in BL6 Host ---")
driver_combos <- combn(levels(samples_bl6$Driver), 2)
for (i in 1:ncol(driver_combos)) {
  run_comparison(dds_bl6_results,
                 "Driver",
                 driver_combos[1, i],
                 driver_combos[2, i],
                 "BL6_Host")
}


# --- [NEW] Part 2: Stratified Analysis (Driver effects within each Model) ---
message("\n\n--- Part 2: Running Stratified Analysis: Driver Effects Within Each Model ---")

# Loop through each Model type present in the BL6 subset
for (model_level in levels(samples_bl6$Model)) {
  message(paste("\n--- Stratifying by Model:", model_level, "---"))
  
  # Create a further subset for the current model
  samples_bl6_model_strat <- samples_bl6[samples_bl6$Model == model_level, ]
  
  # Check if there are multiple drivers to compare within this model
  if (nlevels(factor(samples_bl6_model_strat$Driver)) < 2) {
    message(
      paste(
        "...Skipping Model '",
        model_level,
        "' as it does not contain multiple Driver types for comparison."
      )
    )
    next
  }
  
  # Prepare data for DESeq2
  files_sf_bl6_model_strat <- files_sf_bl6[rownames(samples_bl6_model_strat)]
  
  txi_strat <- tximport(
    files_sf_bl6_model_strat,
    type = "salmon",
    tx2gene = tx2gene_pcg,
    ignoreTxVersion = TRUE
  )
  
  # The design is simpler now as we are only looking at Driver effects
  dds_strat <- DESeqDataSetFromTximport(txi_strat,
                                        colData = samples_bl6_model_strat,
                                        design = ~ SeqBatch + Driver)
  
  # Filter and run DESeq
  keep_strat <- rowSums(counts(dds_strat) >= filt) >= 3
  dds_strat_filt <- dds_strat[keep_strat, ]
  dds_strat_results <- DESeq(dds_strat_filt)
  
  # Perform pairwise Driver comparisons for this model stratum
  driver_combos_strat <- combn(levels(factor(samples_bl6_model_strat$Driver)), 2)
  for (i in 1:ncol(driver_combos_strat)) {
    run_comparison(
      dds_strat_results,
      "Driver",
      driver_combos_strat[1, i],
      driver_combos_strat[2, i],
      paste0("BL6_Host-", model_level)
    )
  }
}

message("\n--- Completed 06_DESeq_BL6_subset.R ---")