# =============================================================================
#
# 01_genes_mouse.R: GENE ANNOTATION FUNCTIONS
#
# =============================================================================

create_tx2gene_maps <- function(gtf_path,
                                cache_dir,
                                filter_gene_types,
                                use_cache = TRUE) {
  file_cache_rds <- file.path(cache_dir, "gtf_data_cache.rds")
  ortholog_cache_rds <- file.path(cache_dir, "ortholog_map_cache.rds")

  if (use_cache && file.exists(file_cache_rds)) {
    message("...Loading cached GTF data.")
    gtf_data <- readRDS(file_cache_rds)
  } else {
    message("...Importing GTF data from file.")
    gtf_data <- rtracklayer::import(gtf_path)
    saveRDS(gtf_data, file = file_cache_rds)
  }

  if (use_cache && file.exists(ortholog_cache_rds)) {
    message("...Loading cached ortholog map.")
    ortholog_map <- readRDS(ortholog_cache_rds)
  } else {
    message("...Querying Ensembl for human orthologs. This may take a moment...")
    ortholog_map <- tryCatch(
      {
        ensembl <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")

        map1 <- biomaRt::getBM(
          attributes = c("ensembl_gene_id", "mgi_symbol"),
          mart = ensembl
        )
        map2 <- biomaRt::getBM(
          attributes = c(
            "ensembl_gene_id",
            "hsapiens_homolog_associated_gene_name"
          ),
          mart = ensembl
        )

        ortholog_map_result <- dplyr::full_join(map1, map2, by = "ensembl_gene_id")
        saveRDS(ortholog_map_result, file = ortholog_cache_rds)
        ortholog_map_result
      },
      error = function(e) {
        warning(
          "Failed to query Ensembl for human orthologs. Proceeding without ortholog information.\nError: ",
          e$message
        )
        return(
          data.frame(
            ensembl_gene_id = character(),
            mgi_symbol = character(),
            hsapiens_homolog_associated_gene_name = character()
          )
        )
      }
    )
  }

  full_map <- as.data.frame(mcols(gtf_data))
  full_map$chromosome_name <- as.character(seqnames(gtf_data))

  base_map <- full_map %>%
    dplyr::filter(
      !is.na(transcript_id) & !is.na(gene_id) &
        !grepl("^Gm", gene_name) &
        !grepl("Rik$", gene_name) & !is.na(gene_name)
    )

  if (!is.null(filter_gene_types)) {
    base_map <- base_map %>% dplyr::filter(gene_type %in% filter_gene_types)
  }

  tx2gene_pcg <- base_map %>%
    dplyr::select(transcript_id, gene_id) %>%
    dplyr::distinct() %>%
    dplyr::mutate(transcript_id = gsub("[.].*$", "", transcript_id))

  gene_name_map <- base_map %>%
    dplyr::filter(type == "gene") %>%
    dplyr::select(gene_id, gene_name) %>%
    dplyr::distinct()

  gene_name_map <- gene_name_map %>%
    dplyr::left_join(ortholog_map, by = c("gene_id" = "ensembl_gene_id"))

  if ("hsapiens_homolog_associated_gene_name" %in% names(gene_name_map)) {
    gene_name_map <- gene_name_map %>%
      dplyr::rename(human_ortholog = hsapiens_homolog_associated_gene_name)
  } else {
    warning("Could not retrieve human orthologs. `human_ortholog` column will be empty.")
    gene_name_map$human_ortholog <- NA_character_
  }

  gene_name_map <- gene_name_map %>%
    dplyr::mutate(gene_name_upper = toupper(gene_name))

  return(list(tx2gene_pcg = tx2gene_pcg, gene_name_map = gene_name_map))
}
