# =============================================================================
#
# 00_setup.R: PROJECT PACKAGE LOADER AND VALIDATION
#
# Description:
# This script loads all R packages required for the analysis pipeline and
# contains a validation function to check inputs before running the analysis.
#
# =============================================================================

message("Loading required packages...")

packages_to_load <- c(
  "tidyverse", "here", "rmarkdown", "reactable", "BiocManager",
  "DESeq2", "tximport", "pheatmap", "RColorBrewer", "ggplotify",
  "BiocParallel", "rtracklayer", "EnhancedVolcano", "knitr",
  "ggsci", "styler", "dendextend"
)

suppressPackageStartupMessages({
  for (pkg in packages_to_load) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warning_msg <- paste("Package", pkg, "not found. Installing...")
      warning(warning_msg)
      BiocManager::install(pkg, update = FALSE)

      # After installing, check again and issue a version warning
      if (requireNamespace(pkg, quietly = TRUE)) {
        library(pkg, character.only = TRUE)
        pkg_version <- as.character(packageVersion(pkg))
        version_warning <- paste("Package", pkg, "version", pkg_version, "was installed and loaded.")
        warning(version_warning)
      } else {
        stop(paste("Failed to install package:", pkg))
      }
    } else {
      library(pkg, character.only = TRUE)
    }
  }
})

message("All packages loaded successfully.")


#' Validate Configuration and Input Files
#'
#' This function checks that the `cfg` list and the `samples` data frame are
#' correctly configured before the main analysis begins. It stops with an
#' informative error if a required column is missing.
#'
#' @param cfg The project configuration list.
#' @param samples A data frame with sample metadata, read from samples.csv.
#'
#' @return This function does not return a value; it stops on error.
#'
validate_config_and_inputs <- function(cfg, samples) {
  message("--- Validating configuration and input files ---")

  required_cols <- unique(c(
    cfg$batch_vars,
    cfg$main_vars,
    cfg$deg_design_vars,
    cfg$interaction_vars,
    cfg$filter_by_factor
  ))

  missing_cols <- setdiff(required_cols, names(samples))

  if (length(missing_cols) > 0) {
    stop(paste(
      "\n\nERROR: The samples.csv file is missing the following required columns defined in your cfg list:",
      paste(missing_cols, collapse = ", "),
      "\nPlease check your column names and the configuration variables in the RMD.",
      "\n"
    ))
  }

  message("...Validation successful.")
}
