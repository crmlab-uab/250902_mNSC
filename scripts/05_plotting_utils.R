# =============================================================================
#
# 05_plotting_utils.R: PLOTTING UTILITY FUNCTIONS
#
# =============================================================================

#' Save a plot in multiple formats (PNG and PDF)
#'
#' @param plot_object The plot object to save (ggplot or pheatmap).
#' @param dir_png Directory for PNG output.
#' @param dir_pdf Directory for PDF output.
#' @param filename_base The base name for the output file, without extension.
#' @param width The width of the plot.
#' @param height The height of the plot.
#' Save a ggplot (or pheatmap object) to PNG and PDF with consistent handling.
save_plot_formats <- function(plot_object,
                              dir_png,
                              dir_pdf,
                              filename_base,
                              width,
                              height) {
  if (!dir.exists(dir_png)) dir.create(dir_png, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(dir_pdf)) dir.create(dir_pdf, recursive = TRUE, showWarnings = FALSE)

  # Convert pheatmap to ggplot if needed
  if (inherits(plot_object, "pheatmap")) {
    plot_object <- ggplotify::as.ggplot(plot_object)
  }

  png_path <- file.path(dir_png, paste0(filename_base, ".png"))
  pdf_path <- file.path(dir_pdf, paste0(filename_base, ".pdf"))

  ggplot2::ggsave(
    filename = png_path,
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300
  )
  ggplot2::ggsave(
    filename = pdf_path,
    plot = plot_object,
    width = width,
    height = height,
    device = "pdf"
  )
}

# (continuation of plotting utility functions)

# Helper to produce a filesystem-safe directory name
sanitize_analysis_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]", "_", x)
}

#
## Additional plotting utilities follow

#' Generate PCA, sample distance heatmap, and GOI boxplots for a DESeqDataSet.
#' Saves into per-analysis subdirectories under graphs/png and graphs/pdf.
generate_qc_plots <- function(dds, analysis_name, gene_name_map, cfg) {
  stopifnot(methods::is(dds, "DESeqDataSet"))

  safe_name <- sanitize_analysis_name(analysis_name)
  png_dir <- file.path(cfg$dir_graphs_png, safe_name)
  pdf_dir <- file.path(cfg$dir_graphs_pdf, safe_name)
  if (!dir.exists(png_dir)) dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(pdf_dir)) dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

  vst_mat <- DESeq2::vst(dds, blind = TRUE)
  vst_assay <- SummarizedExperiment::assay(vst_mat)

  pca_df <- DESeq2::plotPCA(
    vst_mat,
    intgroup = intersect(cfg$main_vars, colnames(SummarizedExperiment::colData(dds))),
    returnData = TRUE
  )
  percent_var <- round(100 * attr(pca_df, "percentVar"))

  color_factor <- intersect(cfg$main_vars, colnames(pca_df))
  color_factor <- if (length(color_factor)) color_factor[1] else NULL

  pca_plot <- ggplot2::ggplot(
    pca_df,
    ggplot2::aes(x = PC1, y = PC2,
                 color = if (!is.null(color_factor)) .data[[color_factor]] else NULL,
                 label = name)
  ) +
    ggplot2::geom_point(size = 3, alpha = 0.85) +
    ggplot2::labs(
      title = paste("PCA -", gsub("_", " ", analysis_name)),
      x = paste0("PC1 (", percent_var[1], "%)"),
      y = paste0("PC2 (", percent_var[2], "%)")
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "right")

  save_plot_formats(
    plot_object = pca_plot,
    dir_png = png_dir,
    dir_pdf = pdf_dir,
    filename_base = paste0(cfg$date, "_", analysis_name, "_PCA"),
    width = 6,
    height = 5
  )

  sample_dist <- dist(t(vst_assay))
  sample_dist_mat <- as.matrix(sample_dist)
  rownames(sample_dist_mat) <- colnames(sample_dist_mat) <- colnames(vst_assay)

  annotation_df <- as.data.frame(SummarizedExperiment::colData(dds)) %>%
    dplyr::select(any_of(c(cfg$batch_vars, cfg$main_vars))) %>%
    dplyr::mutate(dplyr::across(where(is.character), as.factor))

  heatmap_obj <- pheatmap::pheatmap(
    sample_dist_mat,
    clustering_distance_rows = sample_dist,
    clustering_distance_cols = sample_dist,
    annotation_col = annotation_df,
    main = paste("Sample Distance -", gsub("_", " ", analysis_name)),
    silent = TRUE,
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdYlBu")))(100)
  )

  save_plot_formats(
    plot_object = heatmap_obj,
    dir_png = png_dir,
    dir_pdf = pdf_dir,
    filename_base = paste0(cfg$date, "_", analysis_name, "_Heatmap"),
    width = 6,
    height = 5
  )

  # GOI boxplots
  if (!is.null(cfg$qc_plot_settings$project_genes) &&
      length(cfg$qc_plot_settings$project_genes) > 0) {

    if (!is.null(gene_name_map) &&
        all(c("gene_id", "gene_name") %in% colnames(gene_name_map))) {
      goi_map <- gene_name_map %>%
        dplyr::filter(gene_name %in% cfg$qc_plot_settings$project_genes)
    } else {
      goi_map <- data.frame(gene_id = character(), gene_name = character())
    }

    norm_counts <- DESeq2::counts(dds, normalized = TRUE) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("gene_id")

    if (nrow(goi_map) > 0) {
      plot_group_var <- intersect(cfg$main_vars, colnames(SummarizedExperiment::colData(dds)))
      plot_group <- if (length(plot_group_var)) plot_group_var[1] else NULL
      meta_df <- as.data.frame(SummarizedExperiment::colData(dds)) %>%
        tibble::rownames_to_column("sample_id")

      for (gene_symbol in cfg$qc_plot_settings$project_genes) {
        gene_row <- goi_map %>% dplyr::filter(gene_name == gene_symbol)
        if (nrow(gene_row) == 0) next
        gid <- gene_row$gene_id[1]

        cnt_row <- norm_counts %>% dplyr::filter(gene_id == gid)
        if (nrow(cnt_row) == 0) next

        long_df <- cnt_row %>%
          tidyr::pivot_longer(-gene_id,
                              names_to = "sample_id",
                              values_to = "norm_count") %>%
          dplyr::left_join(meta_df, by = "sample_id")

        bp <- ggplot2::ggplot(
          long_df,
          ggplot2::aes_string(
            x = ifelse(is.null(plot_group), "sample_id", plot_group),
            y = "norm_count",
            fill = ifelse(is.null(plot_group), "sample_id", plot_group)
          )
        ) +
          ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
          ggplot2::geom_jitter(width = 0.15, size = 1.6, alpha = 0.8) +
          ggplot2::scale_y_log10() +
          ggplot2::labs(
            title = paste(gene_symbol, "-", gsub("_", " ", analysis_name)),
            x = ifelse(is.null(plot_group), "Sample", plot_group),
            y = "Normalized Count (log10)"
          ) +
          ggplot2::theme_minimal(base_size = 11) +
          ggplot2::theme(legend.position = "none")

        save_plot_formats(
          plot_object = bp,
          dir_png = png_dir,
          dir_pdf = pdf_dir,
          filename_base = paste0(cfg$date, "_", analysis_name, "_GenePlot_", gene_symbol),
          width = 4,
          height = 4
        )
      }
    }
  }

  invisible(TRUE)
}
