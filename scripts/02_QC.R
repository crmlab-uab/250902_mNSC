# =============================================================================
#
# 02_QC.R: QUALITY CONTROL AND DIFFERENTIAL EXPRESSION ANALYSIS FUNCTIONS
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

#' Helper to generate and save a DEG heatmap
#'
#' @param mat Matrix of expression data.
#' @param annotation_df Data frame for column annotations.
#' @param main_title The main title for the heatmap.
#' @param cfg The project configuration list.
#' @param filename_base The base name for the output file.
#' @param width The width of the plot.
#' @param height The height of the plot.
generate_and_save_heatmap <- function(mat, annotation_df, main_title, cfg, filename_base, width, height) {
  p_heatmap <- pheatmap::pheatmap(
    mat,
    scale = "row",
    annotation_col = annotation_df,
    main = main_title
  )
  save_plot_formats(
    plot_object = p_heatmap,
    dir_png = cfg$dir_graphs_png,
    dir_pdf = cfg$dir_graphs_pdf,
    filename_base = filename_base,
    width = width,
    height = height
  )
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
    # --- START: MODIFIED CONTRAST HANDLING ---
    res_args <- list(
      object = dds,
      alpha = p_adj_threshold,
      cooksCutoff = FALSE,
      independentFiltering = TRUE
    )
    comp_label <- ""

    if (!is.null(comp$name)) { # by name
      res_args$name <- comp$name
      comp_label <- comp$name
      message("...Running comparison (name): ", comp_label)
    } else if (!is.null(comp$contrast)) { # by list
      res_args$contrast <- comp$contrast
      comp_label <- comp$label %||% paste(comp$contrast, collapse = "_vs_")
      message("...Running comparison (contrast list): ", comp_label)
    } else { # by pair
      f <- comp$factor; g1 <- comp$group1; g2 <- comp$group2
      res_args$contrast <- c(f, g1, g2)
      comp_label <- paste(f, g1, "vs", g2, sep = "_")
      message("...Running comparison (pair): ", comp_label)
    }

    comp_name <- paste(analysis_name, comp_label, sep = "_")

    res_obj <- try(do.call(DESeq2::results, res_args), silent = TRUE)

    if (inherits(res_obj, "try-error")) {
      warning("......FAILED contrast (skipping): ", comp_name,
              " | Error: ", conditionMessage(attr(res_obj, "condition")))
      next
    }
    # --- END: MODIFIED CONTRAST HANDLING ---

    res_df <- as.data.frame(res_obj) %>%
      tibble::rownames_to_column("gene_id") %>%
      dplyr::left_join(gtf_map, by = "gene_id") %>%
      dplyr::mutate(gene_name_upper = toupper(gene_name))

    all_results[[comp_name]] <- res_df

    res_df_sig <- res_df %>% dplyr::filter(padj < p_adj_threshold, abs(log2FoldChange) > lfc_threshold)

    message(paste("...", comp_name, "found", nrow(res_df_sig), "significant genes."))

    if (nrow(res_df_sig) > 0) {
      # Create analysis-specific output directories
      safe_analysis_name <- sanitize_analysis_name(analysis_name)
      csv_dir <- file.path(cfg$dir_data_csv, safe_analysis_name)
      png_dir <- file.path(cfg$dir_graphs_png, safe_analysis_name)
      pdf_dir <- file.path(cfg$dir_graphs_pdf, safe_analysis_name)
      dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

      # Save significant results to CSV
      write.csv(
        res_df_sig,
        file = file.path(csv_dir, paste0(cfg$date, "_", comp_label, "_DEG_sig.csv")),
        row.names = FALSE
      )

      # Save all results to CSV
      write.csv(
        res_df,
        file = file.path(csv_dir, paste0(cfg$date, "_", comp_label, "_DEG_all.csv")),
        row.names = FALSE
      )

      # Volcano Plot
      p_volcano <- EnhancedVolcano::EnhancedVolcano(
        res_df,
        lab = res_df$gene_name,
        x = "log2FoldChange",
        y = "padj",
        pCutoff = p_adj_threshold,
        FCcutoff = lfc_threshold,
        title = comp_label
      )
      save_plot_formats(
        plot_object = p_volcano,
        dir_png = png_dir,
        dir_pdf = pdf_dir,
        filename_base = paste0(cfg$date, "_", comp_label, "_volcano"),
        width = 12,
        height = 10
      )

      # Heatmaps of DEGs
      vsd <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = FALSE))
      annotation_df <- as.data.frame(colData(dds)[, cfg$main_vars, drop = FALSE])
      
      # Top 50 DEG Heatmap
      top_degs_50 <- res_df_sig %>% dplyr::arrange(padj) %>% head(50)
      mat_50 <- vsd[top_degs_50$gene_id, ]
      rownames(mat_50) <- ifelse(is.na(top_degs_50$gene_name), top_degs_50$gene_id, top_degs_50$gene_name)
      generate_and_save_heatmap(
        mat_50, annotation_df, paste("Top 50 DEGs:", comp_label), cfg,
        paste0(cfg$date, "_", comp_label, "_DEG_heatmap_top50"), 8, 10
      )

      # Top 250 DEG Heatmap
      if (nrow(res_df_sig) > 50) {
        top_degs_250 <- res_df_sig %>% dplyr::arrange(padj) %>% head(250)
        mat_250 <- vsd[top_degs_250$gene_id, ]
        rownames(mat_250) <- ifelse(is.na(top_degs_250$gene_name), top_degs_250$gene_id, top_degs_250$gene_name)
        generate_and_save_heatmap(
          mat_250, annotation_df, paste("Top 250 DEGs:", comp_label), cfg,
          paste0(cfg$date, "_", comp_label, "_DEG_heatmap_top250"), 8, 15
        )
      }


      # Gene Set Analysis
      if (length(gene_sets) > 0) {
        for (gs_name in names(gene_sets)) {
          gs_genes <- gene_sets[[gs_name]]
          sig_gs <- res_df_sig %>% dplyr::filter(gene_name_upper %in% gs_genes | toupper(human_ortholog) %in% gs_genes)

          if (nrow(sig_gs) > 0) {
            mat_gs <- vsd[sig_gs$gene_id, ]
            rownames(mat_gs) <- ifelse(is.na(sig_gs$gene_name), sig_gs$gene_id, sig_gs$gene_name)
            generate_and_save_heatmap(
              mat_gs, annotation_df, paste("DEGs in", gs_name, ":", comp_label), cfg,
              paste0(cfg$date, "_", comp_label, "_", gs_name, "_DEG_heatmap"), 8, 10
            )

            # Box plots for individual DE genes in the gene set
            # Use the main plotting utility for consistency
            # This requires a temporary modification of the config
            temp_cfg <- cfg
            temp_cfg$qc_plot_settings$project_genes <- sig_gs$gene_name
            
            # Suppress the other QC plots (PCA, heatmap) since we only want boxplots
            # A more elegant solution would be a dedicated boxplot function,
            # but this reuses existing code effectively.
            suppressMessages(generate_qc_plots(
              dds = dds,
              analysis_name = file.path(analysis_name, paste0(comp_label, "_", gs_name, "_boxplots")),
              gene_name_map = gtf_map,
              cfg = temp_cfg
            ))
          }
        }
      }
    } else {
        message(paste("...Skipping plot and table generation for", comp_name, "due to no significant results."))
    }
  }
  return(all_results)
}
