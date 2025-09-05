# =============================================================================
#
# 03_DESeq2.R: DESeq2 ANALYSIS AND PLOTTING FUNCTIONS
#
# =============================================================================

run_deseq_core <- function(samples_df,
                           files_sf_vec,
                           design_formula,
                           tx2gene_map,
                           filter_by_factor,
                           filter_min_count) {
  txi_data <- tximport::tximport(
    files_sf_vec,
    type = "salmon",
    tx2gene = tx2gene_map,
    ignoreTxVersion = TRUE
  )
  dds <- DESeq2::DESeqDataSetFromTximport(txi_data, colData = samples_df, design = design_formula)
  keep <- rowSums(counts(dds) >= filter_min_count) >= min(table(samples_df[[filter_by_factor]]))
  DESeq2::DESeq(dds[keep, ])
}

extract_comparisons <- function(dds_processed,
                                comparisons_list,
                                analysis_name,
                                gene_map,
                                cfg,
                                kinase_genes = NULL) {
  results_list <- list()
  for (comp in comparisons_list) {
    comp_name <- paste(analysis_name,
      comp$factor,
      comp$group1,
      "vs",
      comp$group2,
      sep = "_"
    )
    res <- DESeq2::results(
      dds_processed,
      contrast = c(comp$factor, comp$group1, comp$group2),
      alpha = cfg$qval_threshold
    )
    res_df <- as.data.frame(res) %>%
      tibble::rownames_to_column("gene_id") %>%
      dplyr::left_join(gene_map, by = "gene_id") %>%
      dplyr::filter(!is.na(padj))

    results_list[[comp_name]] <- res_df

    sig_genes <- res_df %>% dplyr::filter(padj < cfg$qval_threshold &
      abs(log2FoldChange) > cfg$lfc_threshold)
    sig_kinases <- sig_genes %>% dplyr::filter(gene_name_upper %in% kinase_genes |
      toupper(human_ortholog) %in% kinase_genes)

    # Volcano Plot
    genes_to_label <- (sig_genes %>% dplyr::arrange(padj) %>% head(10))$gene_name
    if (nrow(sig_kinases) > 0) {
      genes_to_label <- unique(c(
        genes_to_label,
        (sig_kinases %>% arrange(padj) %>% head(10))$gene_name
      ))
    }
    EnhancedVolcano::EnhancedVolcano(
      res_df,
      lab = res_df$gene_name,
      x = "log2FoldChange",
      y = "padj",
      selectLab = genes_to_label,
      pCutoff = cfg$qval_threshold,
      FCcutoff = cfg$lfc_threshold,
      title = gsub("_", " ", comp_name),
      subtitle = NULL,
      col = ggsci::pal_jco()(4),
      drawConnectors = TRUE,
      widthConnectors = 0.5
    )
    ggsave(
      here::here(
        cfg$dir_graphs_png,
        paste0(cfg$date, "_", comp_name, "_volcano.png")
      ),
      width = 12,
      height = 10
    )

    # Heatmaps
    vsd <- vst(dds_processed, blind = FALSE)
    if (nrow(sig_genes) > 1) {
      generate_deg_heatmap(
        sig_genes %>% head(50),
        vsd,
        dds_processed,
        comp_name,
        cfg,
        "_DEG_heatmap_top50"
      )
      if (nrow(sig_genes) > 50) {
        generate_deg_heatmap(
          sig_genes %>% head(250),
          vsd,
          dds_processed,
          comp_name,
          cfg,
          "_DEG_heatmap_top250"
        )
      }
    }
    if (nrow(sig_kinases) > 1) {
      generate_deg_heatmap(
        sig_kinases,
        vsd,
        dds_processed,
        comp_name,
        cfg,
        "_Kinase_DEG_heatmap"
      )
    }

    # Kinase Boxplots
    if (nrow(sig_kinases) > 0) {
      for (i in 1:nrow(sig_kinases)) {
        gene_id <- sig_kinases$gene_id[i]
        gene_name <- sig_kinases$gene_name[i]
        plot_data <- plotCounts(
          dds_processed,
          gene = gene_id,
          intgroup = cfg$main_vars,
          returnData = TRUE
        )
        p <- ggplot(plot_data, aes(
          x = Driver,
          y = count,
          fill = Host
        )) +
          geom_boxplot(outlier.shape = NA, alpha = 0.7) +
          geom_jitter(
            aes(shape = Model),
            width = 0.2,
            alpha = 0.9,
            size = 3
          ) +
          scale_fill_jco() +
          scale_y_log10() +
          labs(title = gene_name, y = "Normalized Counts") +
          theme_bw(base_size = 14)
        ggsave(
          here::here(
            cfg$dir_graphs_png_kinases,
            paste0(
              cfg$date,
              "_",
              comp_name,
              "_Boxplot_",
              gene_name,
              ".png"
            )
          ),
          p,
          width = 8,
          height = 6
        )
      }
    }
  }
  return(results_list)
}

generate_deg_heatmap <- function(df, vsd, dds, comp_name, cfg, suffix) {
  mat <- assay(vsd)[df$gene_id, ]
  rownames(mat) <- df$gene_name
  dend_rows <- as.dendrogram(hclust(dist(mat))) %>% dendextend::rotate(order = rev(labels(.)))

  pheatmap::pheatmap(
    mat,
    main = paste(
      gsub("_", " ", comp_name),
      "\nTop",
      nrow(df),
      gsub("_", " ", suffix)
    ),
    scale = "row",
    show_colnames = FALSE,
    annotation_col = as.data.frame(colData(dds)[, c("Driver", "Model", "Host")]),
    annotation_colors = list(Driver = pal_jco()(length(
      levels(colData(dds)$Driver)
    )), Model = pal_jco()(length(
      levels(colData(dds)$Model)
    ))),
    color = cfg$ryb,
    cluster_rows = as.hclust(dend_rows),
    filename = here::here(
      cfg$dir_graphs_png,
      paste0(cfg$date, "_", comp_name, suffix, ".png")
    )
  )
}
