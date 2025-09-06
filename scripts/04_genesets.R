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
  if (!dir.exists(genesets_dir)) {
    warning("Genesets directory not found: ", genesets_dir)
    return(gene_sets)
  }

  files <- list.files(genesets_dir, full.names = TRUE)

  for (file in files) {
    geneset_name <- tools::file_path_sans_ext(basename(file))
    ext <- tolower(tools::file_ext(file))

    genes <- tryCatch({
      dplyr::case_when(
        ext == "csv" ~ toupper(unique(read.csv(file, header = TRUE)[[1]])),
        ext %in% c("tsv", "txt") ~ toupper(unique(read.delim(file, header = TRUE)[[1]])),
        ext == "gmt" ~ toupper(unique(GSEABase::geneIds(GSEABase::getGmt(file)))),
        TRUE ~ NA_character_
      )
    }, error = function(e) {
      warning("Failed to read geneset file: ", file, "\nError: ", e$message)
      return(NULL)
    })
    if (!is.null(genes) && !all(is.na(genes))) gene_sets[[geneset_name]] <- genes
  }
  return(gene_sets)
}
