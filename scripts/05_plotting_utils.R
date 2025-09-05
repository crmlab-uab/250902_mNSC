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
save_plot_formats <- function(plot_object,
                              dir_png,
                              dir_pdf,
                              filename_base,
                              width,
                              height) {
  # For pheatmap objects, convert to ggplot-compatible grob
  if (inherits(plot_object, "pheatmap")) {
    plot_object <- ggplotify::as.ggplot(plot_object)
  }

  ggsave(
    filename = here::here(dir_png, paste0(filename_base, ".png")),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300
  )
  ggsave(
    filename = here::here(dir_pdf, paste0(filename_base, ".pdf")),
    plot = plot_object,
    width = width,
    height = height,
    device = "pdf"
  )
}
