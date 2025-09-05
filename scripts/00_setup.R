# =============================================================================
#
# 00_setup.R: PROJECT PACKAGE LOADER AND VALIDATION
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
  "ggsci",
  "styler",
  "biomaRt",
  "dendextend",
  "GSEABase"
)

suppressPackageStartupMessages({
  for (pkg in packages_to_load) {
    library(pkg, character.only = TRUE)
  }
})

message("All packages loaded successfully.")

#' Validate Configuration and Input Files
#'
#' This function checks for the existence of critical input files and directories,
#' and ensures that the sample metadata file contains all required columns
#' defined in the configuration list.
#'
#' @param cfg A list containing the project configuration.
#' @param samples A data frame containing the sample metadata.
#'
#' @return This function does not return a value but will stop execution if
#'   critical files are missing or if required columns are not found in the
#'   sample metadata. It will issue a warning if the genesets directory is
#'   missing.
validate_config_and_inputs <- function(cfg, samples) {
  message("Validating configuration and inputs...")

  # Check for critical file existence
  if (!file.exists(cfg$file_samples)) {
    stop("Sample file not found: ", cfg$file_samples)
  }
  if (!file.exists(cfg$file_gtf)) {
    stop("GTF file not found: ", cfg$file_gtf)
  }
  if (!dir.exists(cfg$dir_genesets)) {
    warning("Genesets directory not found: ", cfg$dir_genesets)
  }

  # Check for required columns in samples.csv
  required_cols <- unique(c(
    cfg$batch_vars,
    cfg$main_vars,
    cfg$interaction_vars,
    cfg$filter_by_factor
  ))
  missing_cols <- setdiff(required_cols, names(samples))
  if (length(missing_cols) > 0) {
    stop(
      paste(
        "\nERROR: samples.csv is missing required columns:",
        paste(missing_cols, collapse = ", "),
        "\nPlease check column names and cfg variables (batch_vars, main_vars, etc.)."
      )
    )
  }

  message("...Validation successful.")
}
