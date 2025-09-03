# 01_genes_mouse.R
# Contains a function to create transcript-to-gene mapping files.

#' Create Transcript-to-Gene and Gene ID-to-Name Maps
#'
#' @param gtf_path Path to the GTF annotation file.
#' @param cache_dir Directory to save the cached RDS file.
#' @param current_date A string (e.g., "250903") for naming the cache file.
#' @return A list containing three data frames: 'tx2gene_all', 'tx2gene_pcg',
#'         and 'gene_name_map'.

create_tx2gene_maps <- function(gtf_path, cache_dir, current_date) {
  message("--- Running create_tx2gene_maps function ---")
  
  file_cache_rds <- file.path(cache_dir, paste0(current_date, "_gtf_data_cache.rds"))
  
  if (!file.exists(gtf_path)) {
    stop(paste("GTF file not found at:", gtf_path))
  }
  
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
  
  # Create a robust base map
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
    dplyr::filter(!is.na(transcript_id) & !is.na(gene_id))
  
  # --- 1. Create tx2gene map for ALL genes ---
  tx2gene_all <- base_map %>%
    dplyr::select(transcript_id, gene_id) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      transcript_id = gsub("\\..*", "", transcript_id),
      gene_id = gsub("\\..*", "", gene_id)
    ) %>%
    dplyr::distinct()
  
  # --- 2. Create refined PROTEIN-CODING (pcg) map ---
  tx2gene_pcg <- base_map %>%
    dplyr::filter(gene_type == "protein_coding") %>%
    dplyr::filter(
      !grepl(
        "RIKEN|cDNA sequence|DNA segment|predicted gene",
        description,
        ignore.case = TRUE
      )
    ) %>%
    dplyr::filter(!chromosome_name %in% c("MT", "X", "Y")) %>%
    dplyr::select(transcript_id, gene_id) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      transcript_id = gsub("\\..*", "", transcript_id),
      gene_id = gsub("\\..*", "", gene_id)
    ) %>%
    dplyr::distinct()
  
  # --- 3. Create Gene ID to Gene Name Map for Plotting ---
  gene_name_map <- as.data.frame(mcols(gtf_data)) %>%
    dplyr::filter(type == "gene") %>%
    dplyr::select(gene_id, gene_name) %>%
    dplyr::distinct() %>%
    dplyr::mutate(gene_id = gsub("\\..*", "", gene_id))
  
  message(paste(
    "...tx2gene_all map created with",
    nrow(tx2gene_all),
    "entries."
  ))
  message(paste(
    "...tx2gene_pcg map created with",
    nrow(tx2gene_pcg),
    "entries."
  ))
  message(paste(
    "...gene_name_map created with",
    nrow(gene_name_map),
    "entries."
  ))
  
  return(
    list(
      tx2gene_all = tx2gene_all,
      tx2gene_pcg = tx2gene_pcg,
      gene_name_map = gene_name_map
    )
  )
}