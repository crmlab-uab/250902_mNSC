# 05_DESeq_NS1_subset.R
# This script performs a DESeq2 analysis on the subset of samples
# where Model == "NS1".

message("--- Running 05_DESeq_NS1_subset.R: Subset Analysis ---")
library(DESeq2)
library(BiocParallel)
library(EnhancedVolcano)

# --- 1. Create Data Subset ---
message("--- Creating subset for Model == 'NS1' ---")
samples_ns1 <- samples[samples$Model == "NS1", ]

# Drop unused factor levels to avoid errors in the model
samples_ns1$Model <- droplevels(samples_ns1$Model)
samples_ns1$Driver <- factor(samples_ns1$Driver)
samples_ns1$Host <- factor(samples_ns1$Host)

# Subset the quantification file paths
files_sf_ns1 <- files_sf[rownames(samples_ns1)]
message(paste("Subset contains", nrow(samples_ns1), "samples."))

# --- 2. Run DESeq2 Analysis ---
# Design formula includes Host and Driver as sources of variation.
message("...Importing counts and creating DESeqDataSet...")
txi_ns1 <- tximport(
  files_sf_ns1,
  type = "salmon",
  tx2gene = tx2gene_pcg,
  ignoreTxVersion = TRUE
)
dds_ns1 <- DESeqDataSetFromTximport(txi_ns1,
                                    colData = samples_ns1,
                                    design = ~ SeqBatch + Host + Driver)

# Filter for low counts
smallest_group_size <- min(table(samples_ns1$Driver))
keep <- rowSums(counts(dds_ns1) >= filt) >= smallest_group_size
dds_ns1_filt <- dds_ns1[keep, ]

# Run DESeq
message("...Running DESeq()...")
dds_ns1_results <- DESeq(dds_ns1_filt)
message("DESeq2 analysis for NS1 subset complete.")

# --- 3. Extract Pairwise DEG Results and Generate Volcano Plots ---
message("--- Looping through all pairwise Driver comparisons in NS1 subset ---")

# Get all unique pairwise combinations of Driver levels
driver_levels <- levels(samples_ns1$Driver)
driver_combos <- combn(driver_levels, 2)

# Loop through each combination
for (i in 1:ncol(driver_combos)) {
  group1 <- driver_combos[1, i]
  group2 <- driver_combos[2, i]
  
  res_name <- paste0("Driver_", group1, "_vs_", group2, "_in_NS1")
  message("\n--- Running comparison: ", res_name, " ---")
  
  # Extract results
  res <- results(
    dds_ns1_results,
    contrast = c("Driver", group1, group2),
    alpha = qval,
    lfcThreshold = lfc
  )
  
  print(summary(res))
  
  # Generate and save Volcano Plot
  volcano_plot <- EnhancedVolcano(
    res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = qval,
    FCcutoff = lfc,
    title = res_name
  )
  print(volcano_plot)
  ggsave(
    file.path(dir_graph, paste0(date, "_Volcano_", res_name, ".pdf")),
    plot = volcano_plot,
    width = 10,
    height = 10
  )
}

message("\n--- Completed 05_DESeq_NS1_subset.R ---")