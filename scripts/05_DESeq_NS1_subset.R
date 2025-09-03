# 05_DESeq_NS1_subset.R
# This script performs a DESeq2 analysis on the pre-defined 'samples_ns1' subset.

message("--- Running 05_DESeq_NS1_subset.R: Analysis on 'NS1' Subset ---")
library(DESeq2)
library(BiocParallel)
library(EnhancedVolcano)
library(tibble)

# --- The script now expects 'samples_subset' and 'files_sf_subset' to exist ---
samples_ns1 <- samples_subset
files_sf_ns1 <- files_sf_subset

# --- Part 1: Main Effects Analysis (Driver effects across all Hosts in NS1) ---
message("--- Part 1: Running Main Effects Analysis for NS1 Subset ---")
message("...Importing counts and creating DESeqDataSet for the NS1 subset...")
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
keep <- rowSums(counts(dds_ns1) >= filt) >= min(table(samples_ns1$Driver))
dds_ns1_filt <- dds_ns1[keep, ]

# Run DESeq
message("...Running DESeq()...")
dds_ns1_results <- DESeq(dds_ns1_filt)
message("DESeq2 analysis for the main NS1 subset complete.")

# --- Loop through all pairwise Driver comparisons ---
message("--- Processing Main Effect: Driver Comparisons in NS1 Subset ---")
driver_combos <- combn(levels(samples_ns1$Driver), 2)

for (i in 1:ncol(driver_combos)) {
  group1 <- driver_combos[1, i]
  group2 <- driver_combos[2, i]
  res_name <- paste0("Driver_", group1, "_vs_", group2, "_in_NS1_MainEffect")
  message("\n--- Running comparison: ", res_name, " ---")
  
  res <- results(
    dds_ns1_results,
    contrast = c("Driver", group1, group2),
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

# --- [NEW] Part 2: Stratified Analysis (Driver effects within each Host) ---
message("\n\n--- Part 2: Running Stratified Analysis: Driver Effects Within Each Host ---")

# Loop through each Host type present in the NS1 subset
for (host_level in levels(samples_ns1$Host)) {
  message(paste("\n--- Stratifying by Host:", host_level, "---"))
  
  # Create a further subset for the current host
  samples_ns1_host_strat <- samples_ns1[samples_ns1$Host == host_level, ]
  
  # Check if there are multiple drivers to compare within this host
  if (nlevels(factor(samples_ns1_host_strat$Driver)) < 2) {
    message(
      paste(
        "...Skipping Host '",
        host_level,
        "' as it does not contain multiple Driver types for comparison."
      )
    )
    next
  }
  
  # Prepare data for this specific stratum
  files_sf_ns1_host_strat <- files_sf_ns1[rownames(samples_ns1_host_strat)]
  txi_strat <- tximport(
    files_sf_ns1_host_strat,
    type = "salmon",
    tx2gene = tx2gene_pcg,
    ignoreTxVersion = TRUE
  )
  
  # Design is simpler as Host is now constant
  dds_strat <- DESeqDataSetFromTximport(txi_strat,
                                        colData = samples_ns1_host_strat,
                                        design = ~ SeqBatch + Driver)
  
  # Filter and run DESeq
  keep_strat <- rowSums(counts(dds_strat) >= filt) >= min(table(samples_ns1_host_strat$Driver))
  dds_strat_filt <- dds_strat[keep_strat, ]
  dds_strat_results <- DESeq(dds_strat_filt)
  
  # Perform pairwise Driver comparisons for this host stratum
  driver_combos_strat <- combn(levels(factor(samples_ns1_host_strat$Driver)), 2)
  for (i in 1:ncol(driver_combos_strat)) {
    group1 <- driver_combos_strat[1, i]
    group2 <- driver_combos_strat[2, i]
    res_name <- paste0("Driver_", group1, "_vs_", group2, "_in_NS1-", host_level)
    message("\n--- Running comparison: ", res_name, " ---")
    
    res <- results(
      dds_strat_results,
      contrast = c("Driver", group1, group2),
      alpha = qval,
      lfcThreshold = lfc
    )
    print(summary(res))
    
    res_df <- as.data.frame(res) %>% rownames_to_column("gene_id") %>% dplyr::left_join(gene_map, by = "gene_id")
    labels <- ifelse(!is.na(res_df$gene_name),
                     res_df$gene_name,
                     res_df$gene_id)
    
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
}

message("\n--- Completed 05_DESeq_NS1_subset.R ---")