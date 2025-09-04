# =============================================================================
#
# 01_genes_mouse.R: GENE ANNOTATION FUNCTIONS
#
# Description:
# This script contains functions to process a GTF annotation file.
# The filtering criteria are parameterized to make the function reusable
# across different projects with different analytical goals.
#
# =============================================================================

#' Create Transcript-to-Gene and Gene ID-to-Name Maps from a GTF File (Modular)
#'
#' This function parses a GTF file to generate mapping tables used by `tximport`.
#' It implements a caching mechanism to speed up subsequent runs.
#'
#' @param gtf_path character. Path to the gzipped GTF annotation file.
#' @param cache_dir character. Directory to save the cached .rds file.
#' @param current_date character. Date string for naming the cache file.
#' @param filter_gene_types character vector. Gene biotypes to keep.
#' @param exclude_patterns character vector. Patterns in gene description to filter out.
#' @param exclude_chromosomes character vector. Chromosomes to exclude.
#'
#' @return A list containing three data frames: 'tx2gene_all', 'tx2gene_pcg', and 'gene_name_map'.
#'
create_tx2gene_maps <- function(gtf_path,
                                cache_dir,
                                current_date,
                                filter_gene_types,
                                exclude_patterns,
                                exclude_chromosomes) {
  message("--- Running create_tx2gene_maps function ---")
  
  file_cache_rds <- file.path(cache_dir, paste0(current_date, "_gtf_data_cache.rds"))
  if (!file.exists(gtf_path))
    stop(paste("GTF file not found:", gtf_path))
  
  if (file.exists(file_cache_rds)) {
    message(paste("...Loading cached GTF data from:", file_cache_rds))
    gtf_data <- readRDS(file_cache_rds)
  } else {
    message("...Importing GTF. This may take a moment...")
    gtf_data <- rtracklayer::import(gtf_path)
    message(paste("...Saving imported GTF data to cache:", file_cache_rds))
    saveRDS(gtf_data, file = file_cache_rds)
  }
  
  message("...Extracting annotations...")
  full_map <- as.data.frame(mcols(gtf_data))
  full_map$chromosome_name <- as.character(seqnames(gtf_data))
  
  base_map <- full_map %>%
    dplyr::select(any_of(
      c(
        "transcript_id",
        "gene_id",
        "gene_name",
        "gene_type",
        "chromosome_name",
        "description"
      )
    )) %>%
    dplyr::filter(!is.na(.data$transcript_id) &
                    !is.na(.data$gene_id))
  
  # --- Create initial maps ---
  tx2gene_all <- base_map %>%
    dplyr::select(.data$transcript_id, .data$gene_id) %>%
    dplyr::distinct()
  
  filtered_map <- base_map
  if (!is.null(filter_gene_types)) {
    filtered_map <- filtered_map %>% dplyr::filter(.data$gene_type %in% filter_gene_types)
  }
  
  if (!is.null(exclude_patterns) &&
      "description" %in% names(filtered_map)) {
    message("...Filtering based on 'description' column.")
    pattern_regex <- paste(exclude_patterns, collapse = "|")
    filtered_map <- filtered_map %>% dplyr::filter(!grepl(pattern_regex, .data$description, ignore.case = TRUE))
  } else if (!is.null(exclude_patterns)) {
    message("...NOTE: 'description' column not found in GTF, skipping pattern exclusion filter.")
  }
  
  if (!is.null(exclude_chromosomes)) {
    filtered_map <- filtered_map %>% dplyr::filter(!(.data$chromosome_name %in% exclude_chromosomes))
  }
  
  tx2gene_filtered <- filtered_map %>%
    dplyr::select(.data$transcript_id, .data$gene_id) %>%
    dplyr::distinct()
  
  gene_name_map <- as.data.frame(mcols(gtf_data)) %>%
    dplyr::filter(.data$type == "gene") %>%
    dplyr::select(.data$gene_id, .data$gene_name) %>%
    dplyr::distinct()
  
  # --- Remove transcript version numbers to ensure matching ---
  # This makes the matching robust against annotation/quantification mismatches.
  message("...Removing transcript version numbers from mapping files.")
  tx2gene_all$transcript_id <- gsub("\\..*$", "", tx2gene_all$transcript_id)
  tx2gene_filtered$transcript_id <- gsub("\\..*$", "", tx2gene_filtered$transcript_id)
  
  message(paste(
    "...tx2gene_all map created with",
    nrow(tx2gene_all),
    "entries."
  ))
  message(paste(
    "...tx2gene_filtered map created with",
    nrow(tx2gene_filtered),
    "entries."
  ))
  
  return(
    list(
      tx2gene_all = tx2gene_all,
      tx2gene_pcg = tx2gene_filtered,
      gene_name_map = gene_name_map
    )
  )
}