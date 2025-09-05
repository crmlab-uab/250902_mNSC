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
  "dendextend"
)

suppressPackageStartupMessages({
  for (pkg in packages_to_load) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Package", pkg, "not found. Attempting to install..."))
      BiocManager::install(pkg, update = FALSE)
    }
    if (requireNamespace(pkg, quietly = TRUE)) {
      # Issue a warning if a package was just installed
      if (!pkg %in% loadedNamespaces()) {
        warning(
          paste0(
            "Package '",
            pkg,
            "' (version ",
            packageVersion(pkg),
            ") was not found and has been installed. ",
            "Please review your environment for reproducibility."
          )
        )
      }
      library(pkg, character.only = TRUE)
    } else {
      stop(paste("Failed to install and load package:", pkg))
    }
  }
})

message("All packages loaded successfully.")

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
      "\nERROR: samples.csv is missing required columns:",
      paste(missing_cols, collapse = ", ")
    ))
  }
  message("...Validation successful.")
}
