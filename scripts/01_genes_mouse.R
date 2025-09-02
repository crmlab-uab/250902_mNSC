# 01_genes_mouse.R
message("--- Running 01_genes_mouse.R ---")

# Load libraries needed for this script
library(biomaRt)

# Load tx2gene mapping file
tx_to_gene_file <- "./input/tx_to_gene_vM37.tsv"
tx_geneID_genename <- read.delim2(file = tx_to_gene_file, header = TRUE, sep = "\t")
tx_gene_symbol <- tx_geneID_genename[, c(1, 3)] # Use gene_name for summarization
colnames(tx_gene_symbol) <- c("TXNAME", "GENEID")

# --- Gene Annotation with Caching Logic ---
# Define the path for the cached gene annotation file
genes_mouse_cache_file <- paste0(
  dir_output,
  date,
  "_genes_mouse_vM37.csv")
  
  # Check if the cached file already exists
  if (file.exists(genes_mouse_cache_file)) {
    # If it exists, load the data from the CSV file
    message("...Loading cached gene annotations from ",
            genes_mouse_cache_file)
    genes_mouse <- read.csv(genes_mouse_cache_file)
    
  } else {
    # If it does not exist, query biomaRt to get the data
    message(
      "...Cached file not found. Fetching gene annotations from biomaRt (this may take a moment)..."
    )
    mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    genes_mouse <- getBM(
      attributes = c(
        "ensembl_gene_id",
        "external_gene_name",
        "description",
        "gene_biotype",
        "chromosome_name"
      ),
      mart = mart
    )
    
    # Save the fetched data to the cache file for future runs
    message("...Saving annotations to ",
            genes_mouse_cache_file,
            " for future use.")
    write.csv(genes_mouse, file = genes_mouse_cache_file, row.names = FALSE)
  }
  
  # Gene Subsetting (LOGIC PRESERVED AS REQUESTED)
  pcg_df <- subset(genes_mouse, genes_mouse$gene_biotype == "protein_coding")
  pcg_df <- pcg_df[!grepl("RIKEN", pcg_df$description), ]
  pcg_df <- pcg_df[!grepl("cDNA sequence", pcg_df$description), ]
  pcg_df <- pcg_df[!grepl("DNA segment", pcg_df$description), ]
  pcg_df <- pcg_df[!grepl("predicted gene", pcg_df$description), ]
  pcg_df <- subset(pcg_df, pcg_df$chromosome_name != "MT")
  pcg_df <- subset(pcg_df, !chromosome_name %in% c("X", "Y"))
  
  # Create a final vector of protein-coding gene names
  pcg <- pcg_df$external_gene_name
  
  message("Gene and transcript inputs loaded and processed.")
  
  