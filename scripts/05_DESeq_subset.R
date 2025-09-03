# 05_DESeq_generic_subset.R
# This is a generic script to perform a DESeq2 analysis on any pre-defined subset.
# It uses a list named 'analysis_params' that must be defined in the parent environment.

message(
  paste(
    "\n--- Running Subset Analysis for:",
    analysis_params$subset_name,
    "---"
  )
)
library(DESeq2)
library(BiocParallel)
library(EnhancedVolcano)
library(tidyverse)
library(tibble)

# --- 1. Unpack Parameters from Parent Environment ---
samples_subset <- analysis_params$samples_df
files_sf_subset <- analysis_params$files_sf_vec
design_formula <- analysis_params$design_formula
factors_to_compare <- analysis_params$factors_to_compare
subset_name <- analysis_params$subset_name

# --- 2. Run DESeq2 Analysis ---
message(paste("...Importing counts for", subset_name, "subset..."))
txi_subset <- tximport(
  files_sf_subset,
  type = "salmon",
  tx2gene = tx2gene_pcg,
  ignoreTxVersion = TRUE
)
dds_subset <- DESeqDataSetFromTximport(txi_subset, colData = samples_subset, design = design_formula)

# Filter for low counts
keep <- rowSums(counts(dds_subset) >= filt) >= 3 # Using a sensible default
dds_subset_filt <- dds_subset[keep, ]

# Run DESeq
message("...Running DESeq()...")
dds_subset_results <- DESeq(dds_subset_filt)
message("DESeq2 analysis complete for this subset.")

# --- 3. Extract Pairwise DEG Results and Generate Volcano Plots ---
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
    title = res_name
  )
  print(volcano_plot)
}

# --- Loop through all factors specified for comparison ---
for (factor_to_test in factors_to_compare) {
  message(paste(
    "\n--- Processing comparisons for factor:",
    factor_to_test,
    "---"
  ))
  
  factor_levels <- levels(samples_subset[[factor_to_test]])
  
  if (length(factor_levels) > 1) {
    level_combos <- combn(factor_levels, 2)
    for (i in 1:ncol(level_combos)) {
      run_comparison(
        dds_subset_results,
        factor_to_test,
        level_combos[1, i],
        level_combos[2, i],
        subset_name
      )
    }
  } else {
    message(
      paste(
        "...Skipping factor '",
        factor_to_test,
        "' as it has less than two levels in this subset."
      )
    )
  }
}

message(
  paste(
    "\n--- Completed Generic Subset Analysis for:",
    analysis_params$subset_name,
    "---"
  )
)