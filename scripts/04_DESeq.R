# 04_DESeq.R
message("--- Running 04_DESeq.R: Differential Expression Analysis ---")
library(DESeq2)
library(BiocParallel)
library(tidyverse)
library(EnhancedVolcano)
library(tibble)

# --- Main DESeq2 Analysis Loop ---
for (model_name in names(dds_models)) {
  for (subset_name in names(dds_models[[model_name]])) {
    message(
      paste(
        "\n\n===== Running DESeq2 analysis for:",
        model_name,
        "- Subset:",
        subset_name,
        "=====\n"
      )
    )
    
    dds_to_analyze <- dds_models[[model_name]][[subset_name]]
    des_results <- DESeq(dds_to_analyze,
                         parallel = TRUE,
                         BPPARAM = MulticoreParam(nc))
    current_results_names <- resultsNames(des_results)
    
    # --- Analysis for model_1 (main effects only) ---
    if (model_name == "model_1") {
      if ("Model_Tumor_vs_GBO" %in% current_results_names) {
        res <- results(
          des_results,
          name = "Model_Tumor_vs_GBO",
          alpha = qval,
          lfcThreshold = lfc
        )
        res_name <- paste0("model_1_", subset_name, "_Model_Tumor_vs_GBO")
        print(summary(res))
        
        # Create labels with gene names
        res_df <- as.data.frame(res) %>% rownames_to_column("gene_id") %>% left_join(gene_map, by = "gene_id")
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
    
    # --- Analysis for model_2 (interaction effects) ---
    if (model_name == "model_2") {
      # Host effect within each Driver type
      for (drv in unique(samples$Driver)) {
        required_terms <- c("Host_NSG_vs_BL6")
        if (drv != "EGFRvIII") {
          required_terms <- c(required_terms, paste0("HostNSG.Driver", drv))
        }
        
        if (all(required_terms %in% current_results_names)) {
          res <- if (length(required_terms) == 1) {
            results(
              des_results,
              name = required_terms[1],
              alpha = qval,
              lfcThreshold = lfc
            )
          } else {
            results(
              des_results,
              contrast = list(required_terms),
              alpha = qval,
              lfcThreshold = lfc
            )
          }
          res_name <- paste0("model_2_",
                             subset_name,
                             "_Host_NSG_vs_BL6_in_",
                             drv)
          print(summary(res))
          
          res_df <- as.data.frame(res) %>% rownames_to_column("gene_id") %>% left_join(gene_map, by = "gene_id")
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
      
      # Driver comparisons within BL6 host
      dds_bl6 <- dds_to_analyze[, dds_to_analyze$Host == "BL6"]
      dds_bl6$Driver <- droplevels(dds_bl6$Driver)
      design(dds_bl6) <- formula( ~ SeqBatch + Model + Driver)
      dds_bl6 <- DESeq(dds_bl6, parallel = TRUE, BPPARAM = MulticoreParam(nc))
      
      driver_combos <- combn(levels(dds_bl6$Driver), 2)
      for (i in 1:ncol(driver_combos)) {
        drv1 <- driver_combos[1, i]
        drv2 <- driver_combos[2, i]
        res <- results(
          dds_bl6,
          contrast = c("Driver", drv1, drv2),
          alpha = qval,
          lfcThreshold = lfc
        )
        res_name <- paste0("model_2_",
                           subset_name,
                           "_Driver_",
                           drv1,
                           "_vs_",
                           drv2,
                           "_in_BL6")
        print(summary(res))
        
        res_df <- as.data.frame(res) %>% rownames_to_column("gene_id") %>% left_join(gene_map, by = "gene_id")
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
  }
}
message("--- Completed 04_DESeq.R ---")