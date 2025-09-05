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
      if (requireNamespace(pkg, quietly = TRUE)) {
        library(pkg, character.only = TRUE)
        version <- packageVersion(pkg)
        warning(paste(
          "PACKAGE WAS MISSING:",
          pkg,
          "(v",
          version,
          ") was installed."
        ),
        call. = FALSE)
      } else {
        stop(paste("Failed to install package:", pkg), call. = FALSE)
      }
    } else {
      library(pkg, character.only = TRUE)
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
