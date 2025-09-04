# =============================================================================
#
# 03_DESeq2.R: DESeq2 ANALYSIS FUNCTIONS
#
# =============================================================================

library(tximport)

#' Run the Core DESeq2 Analysis
#'
run_deseq_core <- function(samples_df,
                           files_sf_vec,
                           design_formula,
                           tx2gene_map,
                           filter_by_factor,
                           filter_min_count) {
  message(paste0("\n--- Running Core DESeq Analysis for: ", deparse(design_formula), " ---"))
  txi_data <- tximport(files_sf_vec, type = "salmon", tx2gene = tx2gene_map, ignoreTxVersion = TRUE)
  dds <- DESeqDataSetFromTximport(txi_data, colData = samples_df, design = design_formula)
  
  smallest_group_size <- min(table(samples_df[[filter_by_factor]]))
  keep <- rowSums(counts(dds) >= filter_min_count) >= smallest_group_size
  dds <- DESeq(dds[keep, ])
  message("--- Core analysis complete. ---")
  return(dds)
}

#' Extract Pairwise Comparisons, Save Results, and Generate Plots
#'
extract_comparisons <- function(dds_processed,
                                comparisons_list,
                                analysis_name,
                                gene_map,
                                cfg,
                                kinase_genes = NULL) {
  message(paste("\n--- Extracting comparisons for:", analysis_name, "---"))
  results_list <- list()
  
  for (comp in comparisons_list) {
    res_name <- paste0(analysis_name, "_", comp$factor, "_", comp$group1, "_vs_", comp$group2)
    message(paste("...Processing comparison:", res_name))
    
    res <- try(results(dds_processed, contrast = c(comp$factor, comp$group1, comp$group2), alpha = cfg$qval_threshold), silent = TRUE)
    if (inherits(res, "try-error")) {
      message("......Skipping comparison (redundant/not estimable).")
      next
    }
    
    res_df <- as.data.frame(res) %>%
      tibble::rownames_to_column("gene_id") %>%
      dplyr::left_join(gene_map, by = "gene_id") %>%
      dplyr::arrange(padj)
    results_list[[res_name]] <- res_df
    
    write.csv(res_df, here::here(cfg$dir_results_csv, paste0(res_name, "_results.csv")), row.names = FALSE)
    
    sig_genes_padj <- res_df %>% dplyr::filter(padj < cfg$qval_threshold & !is.na(padj))
    sig_genes_df <- sig_genes_padj %>% dplyr::filter(abs(log2FoldChange) > cfg$lfc_threshold)
    message(paste("......Found", nrow(sig_genes_df), "DEGs (padj <", cfg$qval_threshold, " & LFC >", cfg$lfc_threshold, ")"))
    
    if (nrow(sig_genes_df) > 0) {
      top_pcg_labs <- sig_genes_df %>% dplyr::arrange(padj, desc(abs(log2FoldChange))) %>% head(10) %>% pull(gene_name)
      top_kinase_labs <- sig_genes_df %>% dplyr::filter(toupper(gene_name) %in% kinase_genes | toupper(human_ortholog) %in% kinase_genes) %>% dplyr::arrange(padj, desc(abs(log2FoldChange))) %>% head(10) %>% pull(gene_name)
      labels_to_plot <- unique(c(top_pcg_labs, top_kinase_labs))
      
      jco_palette <- ggsci::pal_jco()(4)
      volcano_plot <- EnhancedVolcano(
        res_df, lab = res_df$gene_name, x = 'log2FoldChange', y = 'padj',
        pCutoff = cfg$qval_threshold, FCcutoff = cfg$lfc_threshold,
        col = c("grey30", jco_palette[3], jco_palette[2], jco_palette[1]),
        title = gsub("_", " ", res_name), selectLab = labels_to_plot,
        drawConnectors = TRUE, widthConnectors = 0.5
      )
      ggsave(here::here(cfg$dir_graphs_png, paste0(res_name, "_volcano.png")), volcano_plot, width = 12, height = 10)
      ggsave(here::here(cfg$dir_graphs_pdf, paste0(res_name, "_volcano.pdf")), volcano_plot, width = 12, height = 10)
    }
    
    vsd <- vst(dds_processed, blind = FALSE)
    samples_in_comparison <- rownames(colData(dds_processed)[colData(dds_processed)[[comp$factor]] %in% c(comp$group1, comp$group2), ])
    
    if (nrow(sig_genes_df) > 1) {
      top_50_degs <- sig_genes_df %>% dplyr::arrange(padj, desc(abs(log2FoldChange))) %>% head(50)
      generate_deg_heatmap(top_50_degs, vsd, dds_processed, samples_in_comparison, cfg, res_name, "DEG_heatmap_top50", paste("Top 50 DEGs:", gsub("_", " ", res_name)))
      
      if (nrow(sig_genes_df) > 50) {
        top_250_degs <- sig_genes_df %>% dplyr::arrange(padj, desc(abs(log2FoldChange))) %>% head(250)
        generate_deg_heatmap(top_250_degs, vsd, dds_processed, samples_in_comparison, cfg, res_name, "DEG_heatmap_top250", paste("Top 250 DEGs:", gsub("_", " ", res_name)))
      }
    }
    
    if (!is.null(kinase_genes)) {
      sig_kinases_df <- sig_genes_df %>% dplyr::filter(toupper(gene_name) %in% kinase_genes | toupper(human_ortholog) %in% kinase_genes)
      message(paste("......Found", nrow(sig_kinases_df), "differentially expressed kinases."))
      
      if (nrow(sig_kinases_df) > 1) {
        generate_deg_heatmap(sig_kinases_df, vsd, dds_processed, samples_in_comparison, cfg, res_name, "Kinase_DEG_heatmap", paste("Significant Kinase DEGs:", gsub("_", " ", res_name)))
      }
      if (nrow(sig_kinases_df) > 0) {
        generate_kinase_boxplots(sig_kinases_df, dds_processed, comp, cfg, res_name)
      }
    }
  }
  return(results_list)
}

#' Generate and Save a Heatmap of DEGs
#'
generate_deg_heatmap <- function(heatmap_df, vsd, dds, samples_in_comparison, cfg, res_name, plot_suffix, title) {
  heatmap_mat <- assay(vsd)[heatmap_df$gene_id, samples_in_comparison, drop = FALSE]
  rownames(heatmap_mat) <- heatmap_df$gene_name
  annotation_df <- as.data.frame(colData(dds)[samples_in_comparison, c("Driver", "Model", "Host"), drop = FALSE])
  
  annotation_colors <- lapply(names(annotation_df), function(var) {
    levels <- levels(as.factor(annotation_df[[var]]))
    colors <- ggsci::pal_jco()(length(levels)); names(colors) <- levels; return(colors)
  })
  names(annotation_colors) <- names(annotation_df)
  
  pheatmap::pheatmap(heatmap_mat,
                     main = title, annotation_col = annotation_df, annotation_colors = annotation_colors,
                     color = cfg$ryb, scale = "row", show_colnames = FALSE,
                     fontsize_row = 8, border_color = NA, 
                     filename = here::here(cfg$dir_graphs_png, paste0(res_name, "_", plot_suffix, ".png"))
  )
  pheatmap::pheatmap(heatmap_mat,
                     main = title, annotation_col = annotation_df, annotation_colors = annotation_colors,
                     color = cfg$ryb, scale = "row", show_colnames = FALSE,
                     fontsize_row = 8, border_color = NA,
                     filename = here::here(cfg$dir_graphs_pdf, paste0(res_name, "_", plot_suffix, ".pdf"))
  )
}

#' Generate and Save Box Plots for Individual DE Kinases
#'
generate_kinase_boxplots <- function(sig_kinases_df, dds_processed, comp, cfg, res_name) {
  for (i in 1:nrow(sig_kinases_df)) {
    gene_id <- sig_kinases_df$gene_id[i]; gene_name <- sig_kinases_df$gene_name[i]
    plot_data <- plotCounts(dds_processed, gene = gene_id, intgroup = comp$factor, returnData = TRUE)
    
    p <- ggplot(plot_data, aes(x = .data[[comp$factor]], y = count, fill = .data[[comp$factor]])) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.8, size = 3) +
      ggsci::scale_fill_jco() + scale_y_log10() +
      labs(title = paste("Normalized Counts for", gene_name), subtitle = gsub("_", " ", res_name)) +
      theme_bw(base_size = 14) + theme(legend.position = "none")
    
    ggsave(here::here(cfg$dir_graphs_png_kinase_boxplots, paste0(res_name, "_Boxplot_", gene_name, ".png")), p, width = 7, height = 6)
    ggsave(here::here(cfg$dir_graphs_pdf_kinase_boxplots, paste0(res_name, "_Boxplot_", gene_name, ".pdf")), p, width = 7, height = 6)
  }
}

