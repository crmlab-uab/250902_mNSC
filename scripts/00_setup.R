# =============================================================================
#
# 00_setup.R: PROJECT PACKAGE LOADER
#
# Description:
# This script loads all R packages required for the entire analysis pipeline.
# It also contains a function to validate that the configuration list and
# the samples metadata file are compatible.
#
# =============================================================================

message("Loading required packages...")

packages_to_load <- c(
  "tidyverse",
  "here",
  "rmarkdown",
  "reactable",
  "BiocManager",
  "DESeq2",
  "tximport",
  "pheatmap",
  "RColorBrewer",
  "ggplotify",
  "BiocParallel",
  "rtracklayer",
  "EnhancedVolcano",
  "knitr",
  "styler",
  "ggsci" # Added for publication-quality color palettes
)

# Load each package, suppressing startup messages for cleaner output
suppressPackageStartupMessages({
  for (pkg in packages_to_load) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      library(pkg, character.only = TRUE)
    } else {
      warning(paste("Package", pkg, "not found and could not be loaded."))
    }
  }
})

message("All packages loaded successfully.")


#' Validate Project Configuration Against Input Files
#'
#' This function checks that the columns specified in the `cfg` list exist
#' in the samples data frame. It provides a single, early-fail point to
#' catch common configuration errors.
#'
#' @param cfg The project configuration list.
#' @param samples_df The data frame loaded from the samples metadata file.
#'
#' @return This function does not return a value. It will stop with an error
#'   if validation fails.
#'
validate_config_and_inputs <- function(cfg, samples_df) {
  message("--- Validating configuration against sample metadata ---")

  # Gather all column names required by the configuration
  required_cols <- unique(c(
    cfg$batch_vars,
    cfg$main_vars,
    cfg$interaction_vars,
    cfg$filter_by_factor
  ))

  # Check which columns are missing from the samples data frame
  missing_cols <- setdiff(required_cols, names(samples_df))

  if (length(missing_cols) > 0) {
    stop(paste(
      "\n\nERROR: The following required columns were not found in your samples file:\n  -",
      paste(missing_cols, collapse = "\n  - "),
      "\nPlease check the 'file_samples' path and column names in both your config and CSV file."
    ))
  }

  message("...Validation successful.")
}
# =============================================================================
# End of 00_setup.R
# =============================================================================