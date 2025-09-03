# 01_genes_mouse.R
message("--- Running 01_genes_mouse.R: Creating multiple tx2gene maps with Caching ---")
library(rtracklayer)
library(dplyr)

# --- Gene Annotation with Caching Logic ---
# Define paths for the input GTF and the RDS cache file
file_gtf <- "./input/gencode.vM37.primary_assembly.annotation.gtf.gz"
file_cache_rds <- paste0("./output/", date, "_gtf_data_cache.rds")

if (!file.exists(file_gtf)) {
  stop(paste("GTF file not found at:", file_gtf))
}

# Caching logic: Import the GTF and save as an RDS file if the cache doesn't exist.
if (file.exists(file_cache_rds)) {
  message(paste("...Loading cached GTF data from:", file_cache_rds))
  gtf_data <- readRDS(file_cache_rds)
} else {
  message("...Cached data not found. Importing GTF. This may take a moment...")
  gtf_data <- rtracklayer::import(file_gtf)
  message(paste("...Saving imported GTF data to cache:", file_cache_rds))
  saveRDS(gtf_data, file = file_cache_rds)
}

# --- Extract and Map ---
message("...Extracting transcript, gene, and other annotation information...")
# Combine metadata with chromosome names from the GRanges object
full_map <- as.data.frame(mcols(gtf_data))
full_map$chromosome_name <- as.character(seqnames(gtf_data))

# Ensure required columns exist for filtering
required_cols <- c("transcript_id", "gene_id", "gene_type", "chromosome_name")
if (!"description" %in% colnames(full_map)) {
  full_map$description <- "" # Add dummy column if description is missing
  warning("GTF file does not contain a 'description' column. Some filters will be skipped.")
}

full_map <- full_map %>%
  dplyr::select(any_of(c(required_cols, "description"))) %>%
  dplyr::filter(!is.na(transcript_id) & !is.na(gene_id))

# --- 1. Create tx2gene map for ALL genes ---
message("...Creating tx2gene_all map...")
tx2gene_all <- full_map %>%
  dplyr::select(transcript_id, gene_id) %>%
  dplyr::distinct()

# Remove version numbers for robustness
tx2gene_all$transcript_id <- gsub("\\..*","", tx2gene_all$transcript_id)
tx2gene_all$gene_id <- gsub("\\..*","", tx2gene_all$gene_id)
tx2gene_all <- dplyr::distinct(tx2gene_all)

message(paste("...tx2gene_all map created with", nrow(tx2gene_all), "entries."))
print(head(tx2gene_all))


# --- 2. Gene Subsetting to create a refined PROTEIN-CODING (pcg) map ---
message("\n...Creating tx2gene_pcg map using your specific Gene Subsetting Logic...")

# Apply all filters as requested
tx2gene_pcg_filtered <- full_map %>%
  # Filter for protein-coding biotype
  dplyr::filter(gene_type == "protein_coding") %>%
  # Filter based on description content
  dplyr::filter(!grepl("RIKEN|cDNA sequence|DNA segment|predicted gene", description, ignore.case = TRUE)) %>%
  # Filter out mitochondrial and sex chromosomes
  dplyr::filter(!chromosome_name %in% c("MT", "X", "Y"))

# Finalize the map by selecting the two required columns
tx2gene_pcg <- tx2gene_pcg_filtered %>%
  dplyr::select(transcript_id, gene_id) %>%
  dplyr::distinct()

# Remove version numbers for robustness
tx2gene_pcg$transcript_id <- gsub("\\..*","", tx2gene_pcg$transcript_id)
tx2gene_pcg$gene_id <- gsub("\\..*","", tx2gene_pcg$gene_id)
tx2gene_pcg <- dplyr::distinct(tx2gene_pcg)

message(paste("...Refined tx2gene_pcg map created with", nrow(tx2gene_pcg), "entries."))
print(head(tx2gene_pcg))

message("\n--- Completed 01_genes_mouse.R ---")