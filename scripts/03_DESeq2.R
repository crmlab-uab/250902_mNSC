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
  smallest_group_size <- min(table(colData(dds)[[filter_factor]]))
  keep <- rowSums(counts(dds) >= min_count) >= smallest_group_size
  dds <- dds[keep, ]

  # Run the main DESeq2 function
  dds <- DESeq2::DESeq(dds, parallel = parallel)
  return(dds)
}


#' Extract and Save DESeq2 Comparisons
#'
#' Iterates through a list of comparisons, extracts results from the DESeq2
#' object, saves them to CSV, and calls plotting functions.
#'
#' @param dds_processed The processed DESeqDataSet object from run_deseq_core.
#' @param comparisons_list A list of comparison definitions. Each element of the
#'   list should be a list with three elements: `factor`, `group1`, and `group2`.
#' @param analysis_name A string identifying the current analysis run.
#' @param gene_map A data frame mapping gene IDs to gene names.
#' @param cfg The project's central configuration list.
#' @param gene_sets A named list of gene sets, where each element is a character
#'   vector of gene symbols.
#' @param filt_name The name of the gene filter used (e.g., "protein_coding"),
#'   used for creating output subdirectories.
#' @return A list of data frames, with each data frame containing the results
#'   for one comparison.
extract_comparisons <- function(dds_processed,
                                comparisons_list,
                                analysis_name,
                                gene_map,
                                cfg,
                                gene_sets,
                                filt_name) {
  results_list <- list()

  for (comp in comparisons_list) {
    factor <- comp$factor
    group1 <- comp$group1
    group2 <- comp$group2

    if (is.null(factor) || is.null(group1) || is.null(group2) || !is.character(factor) || !is.character(group1) || !is.character(group2)) {
      warning(paste("Skipping invalid comparison:", paste(comp, collapse = " ")))
      next
    }

    contrast_vector <- c(factor, group2, group1)
    if (length(contrast_vector) != 3) {
      warning(paste("Skipping comparison with invalid contrast vector:", paste(contrast_vector, collapse = " ")))
      next
    }

    comp_name <- paste(analysis_name, factor, group2, "vs", group1, sep = "_")

    message(paste("...EXTRACTING results for", comp_name))

    res <- DESeq2::results(
      dds_processed,
      contrast = contrast_vector,
      alpha = cfg$qval_threshold
    )
    res_df <- as.data.frame(res) %>%
      tibble::rownames_to_column(var = "gene_id") %>%
      dplyr::left_join(gene_map, by = "gene_id") %>%
      dplyr::mutate(gene_name_upper = toupper(gene_name)) %>%
      dplyr::arrange(padj)

    # Create a subdirectory for the specific gene filter type and save CSV
    csv_dir <- here::here(cfg$dir_output, "csv", filt_name)
    dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
    file_path_csv <- here::here(
      csv_dir,
      paste0(cfg$date, "_", comp_name, "_DEG_results.csv")
    )
    write.csv(res_df, file_path_csv, row.names = FALSE)
    message(paste("...Saving results to", file_path_csv))

    # Generate all plots for this specific comparison
    generate_comparison_plots(
      res_df,
      dds_processed,
      comp_name,
      cfg,
      gene_sets,
      gene_map
    )

    results_list[[comp_name]] <- res_df
  }
  return(results_list)
}


#' Generate Plots for a Single Comparison
#'
#' Creates volcano plots, heatmaps, and gene set-specific plots for a given
#' differential expression result.
#'
#' @param res_df The results data frame for the comparison.
#' @param dds The DESeqDataSet object.
#' @param comp_name The unique name of the comparison.
#' @param cfg The project configuration list.
#' @param gene_sets A named list of gene sets.
#' @param gene_map A data frame mapping gene IDs to gene names.
generate_comparison_plots <- function(res_df,
                                      dds,
                                      comp_name,
                                      cfg,
                                      gene_sets,
                                      gene_map) {
  # Volcano Plot
  p_volcano <- EnhancedVolcano::EnhancedVolcano(
    res_df,
    lab = res_df$gene_name,
    x = "log2FoldChange",
    y = "padj",
    pCutoff = cfg$qval_threshold,
    FCcutoff = cfg$lfc_threshold,
    title = gsub("_", " ", comp_name)
  )
  save_plot_formats(
    plot_object = p_volcano,
    dir_png = cfg$dir_graphs_png,
    dir_pdf = cfg$dir_graphs_pdf,
    filename_base = paste0(cfg$date, "_", comp_name, "_volcano"),
    width = 10,
    height = 10
  )

  sig_genes <- res_df %>% dplyr::filter(padj < cfg$qval_threshold &
    abs(log2FoldChange) > cfg$lfc_threshold)
  if (nrow(sig_genes) < 2) {
    message(paste(
      "...Skipping heatmaps for",
      comp_name,
      "due to < 2 significant genes."
    ))
    return() # Not enough significant genes to plot heatmaps
  }

  # Use variance-stabilized counts for heatmaps
  norm_counts <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = FALSE))

  # Top 50 DEG Heatmap
  top50_genes <- head(sig_genes$gene_id, 50)
  p_heatmap50 <- pheatmap::pheatmap(
    norm_counts[top50_genes, ],
    annotation_col = as.data.frame(colData(dds)[, cfg$main_vars, drop = FALSE]),
    labels_row = gene_map$gene_name[match(top50_genes, gene_map$gene_id)],
    main = paste("Top 50 DEGs:", gsub("_", " ", comp_name)),
    show_colnames = FALSE,
    scale = "row",
    silent = TRUE
  )
  save_plot_formats(
    plot_object = p_heatmap50,
    dir_png = cfg$dir_graphs_png,
    dir_pdf = cfg$dir_graphs_pdf,
    filename_base = paste0(cfg$date, "_", comp_name, "_DEG_heatmap_top50"),
    width = 8,
    height = 12
  )

  # Top 250 DEG Heatmap (if applicable)
  if (nrow(sig_genes) > 50) {
    top250_genes <- head(sig_genes$gene_id, 250)
    p_heatmap250 <- pheatmap::pheatmap(
      norm_counts[top250_genes, ],
      annotation_col = as.data.frame(colData(dds)[, cfg$main_vars, drop = FALSE]),
      labels_row = gene_map$gene_name[match(top250_genes, gene_map$gene_id)],
      main = paste("Top 250 DEGs:", gsub("_", " ", comp_name)),
      show_colnames = FALSE,
      scale = "row",
      show_rownames = (length(top250_genes) <= 100), # Only show names if not too cluttered
      silent = TRUE
    )
    save_plot_formats(
      plot_object = p_heatmap250,
      dir_png = cfg$dir_graphs_png,
      dir_pdf = cfg$dir_graphs_pdf,
      filename_base = paste0(cfg$date, "_", comp_name, "_DEG_heatmap_top250"),
      width = 8,
      height = 18
    )
  }

  # Gene set-specific analysis
  if (length(gene_sets) > 0) {
    for (gs_name in names(gene_sets)) {
      gs_genes <- gene_sets[[gs_name]]
      sig_gs <- sig_genes %>% dplyr::filter(gene_name_upper %in% gs_genes |
        toupper(human_ortholog) %in% gs_genes)
      if (nrow(sig_gs) > 1) {
        p_heatmap_gs <- pheatmap::pheatmap(
          norm_counts[sig_gs$gene_id, ],
          annotation_col = as.data.frame(colData(dds)[, cfg$main_vars, drop = FALSE]),
          labels_row = sig_gs$gene_name,
          main = paste("DE", gs_name, ":", gsub("_", " ", comp_name)),
          show_colnames = FALSE,
          scale = "row",
          silent = TRUE
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
            filename_base = paste0(
              cfg$date,
              "_",
              comp_name,
              "_Boxplot_",
              gene_name
            ),
            width = 8,
            height = 6
          )
        }
      }
    }
  }
}