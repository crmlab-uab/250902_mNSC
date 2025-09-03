# 06_DESeq_BL6_subset.R
# This script performs a DESeq2 analysis on the subset of samples
# where Host == "BL6".

message("--- Running 06_DESeq_BL6_subset.R: Subset Analysis for Host == 'BL6' ---")
library(DESeq2)
library(BiocParallel)
library(EnhancedVolcano)
library(tidyverse)

# --- 1. Create Data Subset ---
message("--- Creating subset for Host == 'BL6' ---")
samples_bl6 <- samples[samples$Host == "BL6", ]

# Drop unused factor levels to ensure clean analysis
samples_bl6$Host <- droplevels(samples_bl6$Host)
samples_bl6$Model <- factor(samples_bl6$Model)
samples_bl6$Driver <- factor(samples_bl6$Driver)

# Subset the quantification file paths
files_sf_bl6 <- files_sf[rownames(samples_bl6)]
message(paste("Subset contains", nrow(samples_bl6), "samples."))

# --- 2. Run DESeq2 Analysis ---
# Since Host is constant, the design formula includes Model and Driver effects.
message("...Importing counts and creating DESeqDataSet...")
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
smallest_group_size <- 3 # Based on experimental design
keep <- rowSums(counts(dds_bl6) >= filt) >= smallest_group_size
dds_bl6_filt <- dds_bl6[keep, ]

# Run DESeq
message("...Running DESeq()...")
dds_bl6_results <- DESeq(dds_bl6_filt)
message("DESeq2 analysis for BL6 subset complete.")

# --- 3. Extract Pairwise DEG Results and Generate Volcano Plots ---
message("--- Extracting pairwise DEG results for BL6 Host subset ---")

# A helper function to run and plot a comparison
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
  
  # Volcano Plot
  volcano_plot <- EnhancedVolcano(
    res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = qval,
    FCcutoff = lfc,
    title = res_name,
    subtitle = paste("Comparison in", subset_name, "subset")
  )
  print(volcano_plot)
  ggsave(
    file.path(dir_graph, paste0(date, "_Volcano_", res_name, ".pdf")),
    plot = volcano_plot,
    width = 10,
    height = 10
  )
  
  return(res)
}

# --- Programmatically loop through all Model comparisons in BL6 Host ---
message("\n--- Processing Model Comparisons ---")
model_levels <- levels(samples_bl6$Model)
if (length(model_levels) > 1) {
  model_combos <- combn(model_levels, 2)
  for (i in 1:ncol(model_combos)) {
    run_comparison(dds_bl6_results,
                   "Model",
                   model_combos[1, i],
                   model_combos[2, i],
                   "BL6_Host")
  }
}

# --- Programmatically loop through all Driver comparisons in BL6 Host ---
message("\n--- Processing Driver Comparisons ---")
driver_levels <- levels(samples_bl6$Driver)
if (length(driver_levels) > 1) {
  driver_combos <- combn(driver_levels, 2)
  for (i in 1:ncol(driver_combos)) {
    run_comparison(dds_bl6_results,
                   "Driver",
                   driver_combos[1, i],
                   driver_combos[2, i],
                   "BL6_Host")
  }
}

message("\n--- Completed 06_DESeq_BL6_subset.R ---")