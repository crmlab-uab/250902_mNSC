# =============================================================================
#
# 03_DESeq2.R: CORE DIFFERENTIAL EXPRESSION ANALYSIS FUNCTIONS
#
# =============================================================================

#' Run Core DESeq2 Workflow
#'
#' This function takes sample metadata and quantification files to run the
#' standard DESeq2 workflow, including pre-filtering and normalization.
#'
#' @param samples_df Data frame of sample metadata.
#' @param files_sf Named vector of salmon quantification file paths.
#' @param design_formula Formula for the DESeq2 model.
#' @param tx2gene_map Transcript-to-gene mapping data frame.
#' @param filter_factor The column name in colData to use for determining the
#'   smallest group size for filtering.
#' @param min_count The minimum read count for a gene to be kept.
#' @param parallel A logical value indicating whether to use parallel processing.
#' @return A processed DESeqDataSet object after running DESeq().
run_deseq_core <- function(samples_df,
                           files_sf,
                           design_formula,
                           tx2gene_map,
                           filter_factor,
                           min_count,
                           parallel = FALSE) {
  txi <- tximport::tximport(
    files_sf,
    type = "salmon",
    tx2gene = tx2gene_map,
    ignoreTxVersion = TRUE
  )
  dds <- DESeq2::DESeqDataSetFromTximport(txi, colData = samples_df, design = design_formula)

  # Pre-filtering based on the smallest group size
  smallest_group_size <- min(table(dds[[filter_factor]]))
  keep <- rowSums(counts(dds) >= min_count) >= smallest_group_size
  dds <- dds[keep, ]

  # Run DESeq2
  if (parallel) {
    dds <- DESeq2::DESeq(dds, parallel = TRUE)
  } else {
    dds <- DESeq2::DESeq(dds)
  }

  return(dds)
}


#' Perform Differential Expression Analysis and Generate Outputs
#'
#' This function takes a processed DESeqDataSet object and performs differential
#' expression analysis for specified comparisons. It generates volcano plots,
#' heatmaps, and result tables.
#'
#' @param dds A processed DESeqDataSet object.
#' @param comparisons A list of comparisons to be made.
#' @param gtf_map A data frame mapping gene IDs to gene names.
#' @param p_adj_threshold The adjusted p-value threshold for significance.
#' @param lfc_threshold The log2 fold change threshold for significance.
#' @param gene_sets A list of custom gene sets for heatmap analysis.
#' @param analysis_name A string used to name output files.
#' @param cfg A list containing the project configuration.
#' @return A list of data frames, each containing the results for a comparison.
perform_differential_analysis <- function(dds,
                                          comparisons,
                                          gtf_map,
                                          p_adj_threshold,
                                          lfc_threshold,
                                          gene_sets,
                                          analysis_name,
                                          cfg) {
  all_results <- list()

  for (comp in comparisons) {
    comp_name <- paste0(analysis_name, "_", comp[1], "_", comp[2], "_vs_", comp[3])
    message(paste("Performing analysis for comparison:", comp_name))

    res <- DESeq2::results(
      dds,
      contrast = comp,
      alpha = p_adj_threshold
    )
    res_df <- as.data.frame(res) %>%
      tibble::rownames_to_column("gene_id") %>%
      dplyr::left_join(gtf_map, by = "gene_id") %>%
      dplyr::mutate(gene_name_upper = toupper(gene_name))

    all_results[[comp_name]] <- res_df

    res_df_sig <- res_df %>% dplyr::filter(padj < p_adj_threshold, abs(log2FoldChange) > lfc_threshold)
    
    # MODIFICATION: Added detailed logging to show how many significant genes were found.
    message(paste("Comparison:", comp_name, "- Found", nrow(res_df_sig), "significant genes with padj <", p_adj_threshold, "and |LFC| >", lfc_threshold))

    # MODIFICATION: Added an 'else' block to explicitly state when plots are being skipped.
    if (nrow(res_df_sig) > 0) {
      # Save significant results to CSV
      write.csv(
        res_df_sig,
        file = here::here(
          cfg$dir_data_csv,
          paste0(cfg$date, "_", comp_name, "_DEG_sig.csv")
        )
      )
      
      # Save all results to CSV
      write.csv(
        res_df,
        file = here::here(
          cfg$dir_data_csv,
          paste0(cfg$date, "_", comp_name, "_DEG_all.csv")
        )
      )

      # Volcano Plot
      p_volcano <- EnhancedVolcano::EnhancedVolcano(
        res_df,
        lab = res_df$gene_name,
        x = "log2FoldChange",
        y = "padj",
        pCutoff = p_adj_threshold,
        FCcutoff = lfc_threshold,
        title = comp_name
      )
      save_plot_formats(
        plot_object = p_volcano,
        dir_png = cfg$dir_graphs_png,
        dir_pdf = cfg$dir_graphs_pdf,
        filename_base = paste0(cfg$date, "_", comp_name, "_volcano"),
        width = 12,
        height = 10
      )

      # Heatmaps of DEGs
      vsd <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = FALSE))
      
      # Top 50 DEG Heatmap
      top_degs_50 <- res_df_sig %>% dplyr::arrange(padj) %>% head(50)
      mat_50 <- vsd[top_degs_50$gene_id, ]
      rownames(mat_50) <- top_degs_50$gene_name
      
      p_heatmap_top50 <- pheatmap::pheatmap(
          mat_50,
          scale = "row",
          annotation_col = as.data.frame(colData(dds)[, cfg$main_vars]),
          main = paste("Top 50 DEGs:", comp_name)
      )
      save_plot_formats(
          plot_object = p_heatmap_top50,
          dir_png = cfg$dir_graphs_png,
          dir_pdf = cfg$dir_graphs_pdf,
          filename_base = paste0(cfg$date, "_", comp_name, "_DEG_heatmap_top50"),
          width = 8,
          height = 10
      )

      # Top 250 DEG Heatmap
      if (nrow(res_df_sig) > 50) {
        top_degs_250 <- res_df_sig %>% dplyr::arrange(padj) %>% head(250)
        mat_250 <- vsd[top_degs_250$gene_id, ]
        rownames(mat_250) <- top_degs_250$gene_name
        
        p_heatmap_top250 <- pheatmap::pheatmap(
            mat_250,
            scale = "row",
            annotation_col = as.data.frame(colData(dds)[, cfg$main_vars]),
            main = paste("Top 250 DEGs:", comp_name)
        )
        save_plot_formats(
            plot_object = p_heatmap_top250,
            dir_png = cfg$dir_graphs_png,
            dir_pdf = cfg$dir_graphs_pdf,
            filename_base = paste0(cfg$date, "_", comp_name, "_DEG_heatmap_top250"),
            width = 8,
            height = 15
        )
      }


      # Gene Set Analysis
      if (length(gene_sets) > 0) {
        for (gs_name in names(gene_sets)) {
          gs_genes <- gene_sets[[gs_name]]
          sig_gs <- res_df_sig %>% dplyr::filter(gene_name_upper %in% gs_genes | toupper(human_ortholog) %in% gs_genes)

          if (nrow(sig_gs) > 0) {
            mat_gs <- vsd[sig_gs$gene_id, ]
            rownames(mat_gs) <- sig_gs$gene_name
            p_heatmap_gs <- pheatmap::pheatmap(
              mat_gs,
              scale = "row",
              annotation_col = as.data.frame(colData(dds)[, cfg$main_vars]),
              main = paste("DEGs in", gs_name, ":", comp_name),
              cluster_cols = TRUE
            )
            save_plot_formats(
              plot_object = p_heatmap_gs,
              dir_png = cfg$dir_graphs_png,
              dir_pdf = cfg$dir_graphs_pdf,
              filename_base = paste0(cfg$date, "_", comp_name, "_", gs_name, "_DEG_heatmap"),
              width = 8,
              height = 10
            )

            # Box plots for individual DE genes in the gene set
            for (i in 1:nrow(sig_gs)) {
              gene_id <- sig_gs$gene_id[i]
              gene_name <- sig_gs$gene_name[i]
              plot_data <- DESeq2::plotCounts(
                dds,
                gene = gene_id,
                intgroup = cfg$main_vars,
                returnData = TRUE
              )
              p_box <- ggplot2::ggplot(
                plot_data,
                ggplot2::aes_string(
                  x = "Driver",
                  y = "count",
                  fill = "Model"
                )
              ) +
                ggplot2::geom_boxplot(outlier.shape = NA) +
                ggplot2::geom_jitter(width = 0.2) +
                ggplot2::facet_wrap(~Host) +
                ggplot2::ggtitle(paste("Expression of", gene_name)) +
                ggplot2::theme_bw()
              save_plot_formats(
                plot_object = p_box,
                dir_png = cfg$dir_graphs_png,
                dir_pdf = cfg$dir_graphs_pdf,
                filename_base = paste0(cfg$date, "_", comp_name, "_", gs_name, "_BoxPlot_", gene_name),
                width = 8,
                height = 6
              )
            }
          }
        }
      }
    } else {
        message(paste("--> Skipping plot and table generation for", comp_name, "due to no significant results."))
    }
  }
  return(all_results)
}
