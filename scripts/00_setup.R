# =============================================================================
#
# 00_setup.R: PROJECT PACKAGE LOADER
#
# Description:
# This script loads all R packages required for the entire analysis pipeline.
# It is sourced by the parent R Markdown file at the beginning of a run.
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
  "knitr"
)

# Suppress startup messages for a cleaner console output.
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