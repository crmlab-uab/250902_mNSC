# main.R - Master Analysis Script (Corrected and Final)

# 1. SETUP ----
# Use here::here() to ensure paths are always relative to the project root (/data in container)
source(here::here("00_setup.R")) 

# 2. CONFIGURATION ----
# Central location for all parameters
cfg <- list(
  # Input file paths
  file_samples = here::here("input", "samples.csv"),
  file_gtf = here::here("input", "gencode.vM37.primary_assembly.annotation.gft.gz"),
  dir_sf = here::here("sf"),
  
  # Output directories
  dir_output = here::here("output"),
  dir_graphs = here::here("graphs"),
  
  # DESeq2 parameters
  pval = 0.05,
  qval = 0.05,
  lfc  = 1.0,
  
  # Project-specific genes
  project_genes = c("Egfr", "Pdgfra", "Pten", "Cdkn2a")
)

# Create output directories
dir.create(cfg$dir_output, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$dir_graphs, recursive = TRUE, showWarnings = FALSE)

# Source all function definitions
source(here::here("01_genes_mouse.R"))
source(here::here("03_QC.R")) 
source(here::here("04_DESeq_runner.R"))


# 3. DATA LOADING AND PREPARATION ----
message("--- Loading and preparing data ---")

# Load and process sample metadata
samples <- read.csv(cfg$file_samples, row.names = 1, header = TRUE)
potential_factors <- names(which(sapply(samples, function(col) {
  n_unique <- length(unique(col))
  n_unique > 1 && n_unique < nrow(samples) && (!is.numeric(col) || n_unique < 15)
})))
samples[, potential_factors] <- lapply(samples[, potential_factors], factor)

# Generate full paths to quantification files
files_sf <- file.path(cfg$dir_sf, paste0(rownames(samples), ".sf"))
names(files_sf) <- rownames(samples)
if (!all(file.exists(files_sf))) {
  stop("One or more Salmon quantification files are missing.")
}

# Create transcript-to-gene maps
gtf_maps <- create_tx2gene_maps(
  gtf_path = cfg$file_gtf,
  cache_dir = cfg$dir_output,
  current_date = format(Sys.Date(), "%y%m%d")
)

# 4. ANALYSIS: FULL DATASET ----
message("\n--- Starting Analysis: Full Dataset ---")

# Define the design formula
design_full <- as.formula("~ SeqBatch + Host + Model + Driver")

# Run the core DESeq analysis ONCE
dds_full_processed <- run_deseq_core(
  samples_df = samples,
  files_sf_vec = files_sf,
  design_formula = design_full,
  tx2gene_map = gtf_maps$tx2gene_pcg
)

# Programmatically generate all pairwise comparisons for key factors
factors_to_compare_full <- c("Driver", "Model", "Host")
comparisons_full <- list()
for (f in factors_to_compare_full) {
  level_pairs <- combn(levels(samples[[f]]), 2, simplify = FALSE)
  for (pair in level_pairs) {
    comparisons_full <- append(comparisons_full, list(list(factor = f, group1 = pair[1], group2 = pair[2])))
  }
}
message(paste("...Generated", length(comparisons_full), "total pairwise comparisons for the full dataset."))

# Extract all comparisons and generate plots
results_full <- extract_comparisons(
  dds_processed = dds_full_processed,
  comparisons_list = comparisons_full,
  analysis_name = "Full_Dataset",
  gene_map = gtf_maps$gene_name_map,
  dir_graphs = cfg$dir_graphs,
  qval = cfg$qval,
  lfc = cfg$lfc
)

# 5. ANALYSIS: NS1 MODEL SUBSET (RESTORED) ----
message("\n--- Starting Analysis: NS1 Model Subset ---")

# Subset the data
samples_ns1 <- droplevels(samples[samples$Model == "NS1", ])
files_sf_ns1 <- files_sf[rownames(samples_ns1)]

# Define the design for this subset
design_ns1 <- as.formula("~ SeqBatch + Host + Driver")

# Run the core analysis ONCE for this subset
dds_ns1_processed <- run_deseq_core(
  samples_df = samples_ns1,
  files_sf_vec = files_sf_ns1,
  design_formula = design_ns1,
  tx2gene_map = gtf_maps$tx2gene_pcg
)

# Programmatically generate all comparisons for varying factors in this subset
factors_to_compare_ns1 <- c("Driver", "Host")
comparisons_ns1 <- list()
for (f in factors_to_compare_ns1) {
  level_pairs <- combn(levels(samples_ns1[[f]]), 2, simplify = FALSE)
  for (pair in level_pairs) {
    comparisons_ns1 <- append(comparisons_ns1, list(list(factor = f, group1 = pair[1], group2 = pair[2])))
  }
}
message(paste("...Generated", length(comparisons_ns1), "total pairwise comparisons for the NS1 subset."))

# Extract results
results_ns1 <- extract_comparisons(
  dds_processed = dds_ns1_processed,
  comparisons_list = comparisons_ns1,
  analysis_name = "NS1_Model_Subset",
  gene_map = gtf_maps$gene_name_map,
  dir_graphs = cfg$dir_graphs
)

# 6. ANALYSIS: BL6 HOST SUBSET (CORRECTED) ----
message("\n--- Starting Analysis: BL6 Host Subset ---")

# Subset the data
samples_bl6 <- droplevels(samples[samples$Host == "BL6", ])
files_sf_bl6 <- files_sf[rownames(samples_bl6)]

# Define the design for this subset
design_bl6 <- as.formula("~ SeqBatch + Model + Driver")

# Run the core analysis ONCE for this subset
dds_bl6_processed <- run_deseq_core(
  samples_df = samples_bl6,
  files_sf_vec = files_sf_bl6,
  design_formula = design_bl6,
  tx2gene_map = gtf_maps$tx2gene_pcg
)

# Programmatically generate all comparisons for varying factors in this subset
factors_to_compare_bl6 <- c("Model", "Driver")
comparisons_bl6 <- list()
for (f in factors_to_compare_bl6) {
  level_pairs <- combn(levels(samples_bl6[[f]]), 2, simplify = FALSE)
  for (pair in level_pairs) {
    comparisons_bl6 <- append(comparisons_bl6, list(list(factor = f, group1 = pair[1], group2 = pair[2])))
  }
}
message(paste("...Generated", length(comparisons_bl6), "total pairwise comparisons for the BL6 subset."))

# Extract results
results_bl6 <- extract_comparisons(
  dds_processed = dds_bl6_processed,
  comparisons_list = comparisons_bl6,
  analysis_name = "BL6_Host_Subset",
  gene_map = gtf_maps$gene_name_map,
  dir_graphs = cfg$dir_graphs
)

# 7. SAVE RESULTS ----
# Save key objects for the report. This allows the .Rmd to load them instantly.
message("\n--- Saving analysis objects to file ---")
save(
  samples,
  gtf_maps,
  dds_full_processed,
  results_full,
  dds_ns1_processed,
  results_ns1,
  dds_bl6_processed,
  results_bl6,
  file = here::here(cfg$dir_output, "deseq_analysis_results.RData")
)

message("--- MASTER SCRIPT COMPLETE ---")