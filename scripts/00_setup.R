# =============================================================================
#
# 00_setup.R: PROJECT PACKAGE LOADER AND VALIDATION
#
# Description:
# This script loads all R packages required for the analysis pipeline and
# contains a validation function to check inputs before running the analysis.
# It will automatically attempt to install any missing packages.
#
# =============================================================================

message("Loading required packages...")

packages_to_load <- c(
  "tidyverse", "here", "rmarkdown", "reactable", "BiocManager",
  "DESeq2", "tximport", "pheatmap", "RColorBrewer", "ggplotify",
  "BiocParallel", "rtracklayer", "EnhancedVolcano", "knitr",
  "ggsci", "styler"
)

# --- Automated Package Installation ---
# This loop checks for each package. If it's not installed, it attempts
# to install it using BiocManager, then issues a warning with the version.
for (pkg in packages_to_load) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Package", pkg, "not found. Attempting installation..."))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
    
    # After installation, issue a warning to the user with the version.
    if (requireNamespace(pkg, quietly = TRUE)) {
      pkg_version <- as.character(utils::packageVersion(pkg))
      warning(
        paste0(
          "SETUP WARNING: Package '", pkg, "' was not found and has been installed. ",
          "Version: ", pkg_version
        ),
        call. = FALSE
      )
    } else {
      stop(paste("CRITICAL ERROR: Failed to install package:", pkg, ". Please install it manually and restart."), call. = FALSE)
    }
  }
}

# --- Load all packages ---
suppressPackageStartupMessages({
  for (pkg in packages_to_load) {
    library(pkg, character.only = TRUE)
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
    cfg$interaction_vars,
    cfg$filter_by_factor
  ))
  
  missing_cols <- setdiff(required_cols, names(samples))
  
  if (length(missing_cols) > 0) {
    stop(paste(
      "\n\nERROR: The samples.csv file is missing the following required columns defined in your cfg list:",
      paste(missing_cols, collapse = ", "),
      "\nPlease check your column names and the 'batch_vars', 'main_vars', etc., in the RMD.",
      "\n"
    ))
  }
  
  message("...Validation successful.")
}