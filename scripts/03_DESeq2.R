# =============================================================================
#
# 03_DESeq2.R: DESeq2 ANALYSIS FUNCTIONS
#
# Description:
# This script contains reusable functions for DESeq2 analysis.
# All project-specific variables are passed as arguments via the cfg list.
#
# =============================================================================

library(tximport)

#' Run the Core DESeq2 Analysis (Modular)
#'
#' @param samples_df A data frame with sample metadata.
#' @param files_sf_vec A named character vector of paths to 'quant.sf' files.
#' @param design_formula A formula specifying the experimental design.
#' @param tx2gene_map A data frame with `transcript_id` and `gene_id` columns.
#' @param filter_by_factor character. The column name in `samples_df` for the low-count filter.
#' @param filter_min_count numeric. The minimum count for the low-count filter.
#'
#' @return A processed `DESeqDataSet` object after running `DESeq()`.
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

  if (!filter_by_factor %in% names(samples_df)) {
    stop(paste("The filter_by_factor '", filter_by_factor, "' is not a column in the samples data frame."))
  }

  smallest_group_size <- min(table(samples_df[[filter_by_factor]]))
  keep <- rowSums(counts(dds) >= filter_min_count) >= smallest_group_size
  dds_filt <- dds[keep, ]

  dds_processed <- DESeq(dds_filt)
  message("--- Core analysis complete. ---")
  return(dds_processed)
}


#' Extract Pairwise Comparisons, Save Results, and Generate Plots
#'
#' @param dds_processed A `DESeqDataSet` object that has been run through `DESeq()`.
#' @param comparisons_list A list defining the pairwise comparisons to perform.
#' @param analysis_name A string for naming output files.
#' @param gene_map A data frame for mapping gene IDs to names for plot labels.
#' @param cfg The project configuration list.
#' @param kinase_genes A character vector of kinase gene names.
#'
#' @return A named list of data frames containing DESeq2 results.
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

    # <<< FIX: Removed lfcThreshold from the initial results call. >>>
    # This is the standard, more sensitive method. We test for significance
    # first (padj), then filter by effect size (log2FoldChange).
    res <- try(results(
      dds_processed,
      contrast = c(comp$factor, comp$group1, comp$group2),
      alpha = cfg$qval_threshold
    ), silent = TRUE)

    if (inherits(res, "try-error")) {
      message(paste("......Skipping comparison:", res_name, "(redundant/not estimable)."))
      next
    }

    res_df <- as.data.frame(res) %>%
      tibble::rownames_to_column("gene_id") %>%
      dplyr::left_join(gene_map, by = "gene_id") %>%
      dplyr::arrange(padj)

    results_list[[res_name]] <- res_df

    csv_filename <- here::here(cfg$dir_results_csv, paste0(res_name, "_results.csv"))
    write.csv(res_df, file = csv_filename, row.names = FALSE)

    # Now define significant genes using a post-hoc filter
    sig_genes_padj <- res_df %>% dplyr::filter(padj < cfg$qval_threshold & !is.na(padj))
    sig_genes_df <- sig_genes_padj %>% dplyr::filter(abs(log2FoldChange) > cfg$lfc_threshold)

    message(paste("......Found", nrow(sig_genes_padj), "genes with padj <", cfg$qval_threshold))
    message(paste("......Found", nrow(sig_genes_df), "genes after applying LFC >", cfg$lfc_threshold, "filter."))

    # --- Volcano Plot ---
    jco_palette <- ggsci::pal_jco()(4)
    volcano_plot <- EnhancedVolcano(
      res_df,
      lab = res_df$gene_name, x = "log2FoldChange", y = "padj",
      pCutoff = cfg$qval_threshold, FCcutoff = cfg$lfc_threshold,
      col = c("grey30", jco_palette[3], jco_palette[2], jco_palette[1]),
      title = gsub("_", " ", res_name)
    )
    volcano_png <- here::here(cfg$dir_graphs_png, paste0(res_name, "_volcano.png"))
    volcano_pdf <- here::here(cfg$dir_graphs_pdf, paste0(res_name, "_volcano.pdf"))
    ggsave(volcano_png, volcano_plot, width = 12, height = 10, dpi = 300)
    ggsave(volcano_pdf, volcano_plot, width = 12, height = 10, device = "pdf")

    vsd <- vst(dds_processed, blind = FALSE)
    samples_in_comparison <- colData(dds_processed) %>%
      as.data.frame() %>%
      dplyr::filter(.data[[comp$factor]] %in% c(comp$group1, comp$group2)) %>%
      rownames()

    if (nrow(sig_genes_df) > 1) {
      top_50_degs <- sig_genes_df %>%
        dplyr::arrange(padj, desc(abs(log2FoldChange))) %>%
        head(50)
      generate_deg_heatmap(
        heatmap_df = top_50_degs, vsd = vsd, dds = dds_processed,
        samples_in_comparison = samples_in_comparison, cfg = cfg, res_name = res_name,
        plot_suffix = "DEG_heatmap", title = paste("Top 50 DEGs:", gsub("_", " ", res_name))
      )
    }

    if (!is.null(kinase_genes)) {
      # Check for significant kinases before and after LFC filter for debugging
      sig_kinases_padj <- sig_genes_padj %>% dplyr::filter(gene_name %in% kinase_genes)
      sig_kinases_df <- sig_genes_df %>% dplyr::filter(gene_name %in% kinase_genes)

      message(paste("......Found", nrow(sig_kinases_padj), "kinases with padj <", cfg$qval_threshold))
      message(paste("......Found", nrow(sig_kinases_df), "kinases after applying LFC filter."))

      if (nrow(sig_kinases_df) > 1) {
        generate_deg_heatmap(
          heatmap_df = sig_kinases_df, vsd = vsd, dds = dds_processed,
          samples_in_comparison = samples_in_comparison, cfg = cfg, res_name = res_name,
          plot_suffix = "Kinase_DEG_heatmap", title = paste("Significant Kinase DEGs:", gsub("_", " ", res_name))
        )
      }

      if (nrow(sig_kinases_df) > 0) {
        generate_kinase_boxplots(
          sig_kinases_df = sig_kinases_df, dds_processed = dds_processed, comp = comp,
          cfg = cfg, res_name = res_name
        )
      }
    }
  }
  return(results_list)
}


#' Generate and Save a Heatmap of Differentially Expressed Genes
#'
#' @param heatmap_df A data frame of genes to include in the heatmap.
#' @param vsd A `DESeqTransform` object.
#' @param dds A `DESeqDataSet` object.
#' @param samples_in_comparison A character vector of sample names.
#' @param cfg The project configuration list.
#' @param res_name The name of the comparison for file naming.
#' @param plot_suffix A suffix for the output filenames.
#' @param title The main title for the heatmap plot.
#'
generate_deg_heatmap <- function(heatmap_df, vsd, dds, samples_in_comparison, cfg, res_name, plot_suffix, title) {
  heatmap_mat <- assay(vsd)[heatmap_df$gene_id, samples_in_comparison, drop = FALSE]
  rownames(heatmap_mat) <- heatmap_df$gene_name
  annotation_df <- as.data.frame(colData(dds)[samples_in_comparison, cfg$main_vars, drop = FALSE])

  annotation_colors <- lapply(cfg$main_vars, function(var) {
    if (var %in% names(annotation_df)) {
      levels <- levels(annotation_df[[var]])
      colors <- ggsci::pal_jco()(length(levels))
      names(colors) <- levels
      return(colors)
    }
  })
  names(annotation_colors) <- cfg$main_vars

  png_filename <- here::here(cfg$dir_graphs_png, paste0(res_name, "_", plot_suffix, ".png"))
  pdf_filename <- here::here(cfg$dir_graphs_pdf, paste0(res_name, "_", plot_suffix, ".pdf"))

  pheatmap::pheatmap(heatmap_mat,
    main = title, annotation_col = annotation_df, annotation_colors = annotation_colors,
    color = cfg$ryb, scale = "row", show_colnames = FALSE,
    fontsize_row = 8, border_color = NA, filename = png_filename
  )
  pheatmap::pheatmap(heatmap_mat,
    main = title, annotation_col = annotation_df, annotation_colors = annotation_colors,
    color = cfg$ryb, scale = "row", show_colnames = FALSE,
    fontsize_row = 8, border_color = NA, filename = pdf_filename
  )
}

#' Generate and Save Box Plots for Individual DE Kinases
#'
#' @param sig_kinases_df A data frame of significant kinase genes.
#' @param dds_processed A `DESeqDataSet` object.
#' @param comp The current comparison list object.
#' @param cfg The project configuration list.
#' @param res_name The name of the comparison for file naming.
#'
generate_kinase_boxplots <- function(sig_kinases_df, dds_processed, comp, cfg, res_name) {
  comparison_factor <- comp$factor
  for (i in 1:nrow(sig_kinases_df)) {
    gene_data <- sig_kinases_df[i, ]
    gene_id <- gene_data$gene_id
    gene_name <- gene_data$gene_name

    if (is.na(gene_name)) gene_name <- gene_id

    plot_data <- plotCounts(dds_processed, gene = gene_id, intgroup = comparison_factor, returnData = TRUE)

    p <- ggplot(plot_data, aes(x = .data[[comparison_factor]], y = count, fill = .data[[comparison_factor]])) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.8, size = 3) +
      ggsci::scale_fill_jco() +
      scale_y_log10() +
      labs(
        title = paste("Normalized Counts for", gene_name),
        subtitle = gsub("_", " ", res_name),
        x = comparison_factor, y = "Normalized Counts (log10 scale)"
      ) +
      theme_bw(base_size = 14) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

    plot_filename_png <- here::here(cfg$dir_graphs_png_kinase_boxplots, paste0(res_name, "_Boxplot_", gene_name, ".png"))
    plot_filename_pdf <- here::here(cfg$dir_graphs_pdf_kinase_boxplots, paste0(res_name, "_Boxplot_", gene_name, ".pdf"))

    ggsave(filename = plot_filename_png, plot = p, width = 7, height = 6, dpi = 300)
    ggsave(filename = plot_filename_pdf, plot = p, width = 7, height = 6, device = "pdf")
  }
}
