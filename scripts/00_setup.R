# 00_setup.R
# In a containerized environment, this script simply loads all required packages.
message("Loading required packages...")

packages_to_load <- c(
  "tidyverse", "rmarkdown", "reactable", "BiocManager", "DESeq2", 
  "tximport", "biomaRt", "pheatmap", "RColorBrewer", "ggplotify", 
  "BiocParallel", "RNAseqQC"
)

# Load all libraries, suppressing startup messages
suppressPackageStartupMessages({
  for (pkg in packages_to_load) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      library(pkg, character.only = TRUE)
    } else {
      warning(paste("Package", pkg, "not found and could not be loaded."))
    }
  }
})

message("All packages loaded.")
