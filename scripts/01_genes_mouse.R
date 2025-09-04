# =============================================================================
#
# 01_genes_mouse.R: GENE ANNOTATION FUNCTIONS
#
# =============================================================================

#' Create Transcript-to-Gene and Gene ID-to-Name Maps from a GTF File
#'
create_tx2gene_maps <- function(gtf_path,
                                cache_dir,
                                filter_gene_types,
                                exclude_patterns,
                                exclude_chromosomes) {
  message("--- Running create_tx2gene_maps function ---")
  gtf_cache_rds <- file.path(cache_dir, "gtf_data_cache.rds")
  ortholog_cache_rds <- file.path(cache_dir, "ortholog_map_cache.rds")
  
  if (file.exists(gtf_cache_rds)) {
    gtf_data <- readRDS(gtf_cache_rds)
  } else {
    message("...Importing GTF. This may take a moment...")
    gtf_data <- rtracklayer::import(gtf_path)
    saveRDS(gtf_data, file = gtf_cache_rds)
  }
  
  gene_name_map_mouse <- as.data.frame(mcols(gtf_data)) %>%
    dplyr::filter(type == "gene") %>%
    dplyr::select(gene_id, gene_name) %>%
    dplyr::distinct() %>%
    dplyr::filter(!is.na(gene_name) & gene_name != "") %>%
    dplyr::filter(!startsWith(gene_name, "Gm")) %>%
    dplyr::filter(!endsWith(gene_name, "Rik")) %>%
    dplyr::mutate(gene_name_upper = toupper(gene_name))
  
  if (file.exists(ortholog_cache_rds)) {
    message("...Loading cached human ortholog map.")
    ortholog_map <- readRDS(ortholog_cache_rds)
  } else {
    message("...Querying Ensembl for human orthologs. This may take a moment...")
    ensembl <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    ortholog_map <- biomaRt::getBM(
      attributes = c("mgi_symbol", "hgnc_symbol"),
      filters = "mgi_symbol",
      values = gene_name_map_mouse$gene_name,
      mart = ensembl
    )
    saveRDS(ortholog_map, file = ortholog_cache_rds)
  }
  
  ortholog_map <- ortholog_map %>%
    dplyr::rename(gene_name = mgi_symbol, human_ortholog = hgnc_symbol) %>%
    dplyr::filter(human_ortholog != "") %>%
    dplyr::distinct(gene_name, .keep_all = TRUE)
  
  gene_name_map <- dplyr::left_join(gene_name_map_mouse, ortholog_map, by = "gene_name")
  
  message(paste("...Found human orthologs for", sum(!is.na(gene_name_map$human_ortholog)), "genes."))
  
  tx2gene_filtered <- as.data.frame(mcols(gtf_data)) %>%
    dplyr::filter(type == "transcript") %>%
    dplyr::select(transcript_id, gene_id) %>%
    dplyr::distinct() %>%
    dplyr::mutate(transcript_id = gsub("\\..*$", "", transcript_id))
  
  return(
    list(
      tx2gene_pcg = tx2gene_filtered,
      gene_name_map = gene_name_map
    )
  )
}

