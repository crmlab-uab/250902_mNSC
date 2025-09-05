# =============================================================================
#
# 04_genesets.R: GENESET LOADING FUNCTIONS
#
# =============================================================================

#' Load Gene Sets from a Directory
#'
#' This function loads all gene sets from a specified directory, supporting
#' .csv, .tsv, and .gmt file formats.
#'
#' @param genesets_dir Path to the directory containing gene set files.
#'
#' @return A named list of gene sets. Each element of the list is a character
#'   vector of gene symbols, and the name of the element is the name of the
#'   gene set (derived from the file name).
load_gene_sets <- function(genesets_dir) {
  gene_sets <- list()
  files <- list.files(genesets_dir, full.names = TRUE)

  for (file in files) {
    geneset_name <- tools::file_path_sans_ext(basename(file))

    if (endsWith(file, ".csv")) {
      df <- read.csv(file, header = TRUE)
      gene_sets[[geneset_name]] <- toupper(unique(df[, 1]))
    } else if (endsWith(file, ".tsv")) {
      df <- read.delim(file, header = TRUE)
      gene_sets[[geneset_name]] <- toupper(unique(df[, 1]))
    } else if (endsWith(file, ".gmt")) {
      gmt <- GSEABase::getGmt(file)
      gene_sets[[geneset_name]] <- toupper(unique(GSEABase::geneIds(gmt)))
    }
  }

  return(gene_sets)
}
