# =============================================================================
#
# 03_DESeq2.R: CORE DIFFERENTIAL EXPRESSION ANALYSIS FUNCTIONS
#
# =============================================================================

#' Run Core DESeq2 Workflow
#'
#' This function takes sample metadata and quantification files to run the
#' standard DESeq2 workflow, including pre-filtering and normalization.
#'
#' @param samples_df Data frame of sample metadata.
#' @param files_sf Named vector of salmon quantification file paths.
#' @param design_formula Formula for the DESeq2 model.
#' @param tx2gene_map Transcript-to-gene mapping data frame.
#' @param filter_factor The column name in colData to use for determining the
#'   smallest group size for filtering.
#' @param min_count The minimum read count for a gene to be kept.
#' @param parallel A logical value indicating whether to use parallel processing.
#' @return A processed DESeqDataSet object after running DESeq().
run_deseq_core <- function(samples_df,
                           files_sf,
                           design_formula,
                           tx2gene_map,
                           filter_factor,
                           min_count,
                           parallel = FALSE) {
  txi <- tximport::tximport(
    files_sf,
    type = "salmon",
    tx2gene = tx2gene_map,
    ignoreTxVersion = TRUE
  )
  dds <- DESeq2::DESeqDataSetFromTximport(txi, colData = samples_df, design = design_formula)

  # Pre-filtering based on the smallest group size
  smallest_group_size <- min(table(colData(dds)[[filter_factor]]))
  keep <- rowSums(counts(dds) >= min_count) >= smallest_group_size
  dds <- dds[keep, ]

  # Run the main DESeq2 function
  dds <- DESeq2::DESeq(dds, parallel = parallel)
  return(dds)
}


valid_pairwise_contrasts <- function(dds, factor_var) {
  if (!factor_var %in% names(SummarizedExperiment::colData(dds))) return(list())
  lv <- levels(SummarizedExperiment::colData(dds)[[factor_var]])
  if (length(lv) < 2) return(list())
  combn(lv, 2, FUN = function(x) list(factor = factor_var, group1 = x[1], group2 = x[2]), simplify = FALSE)
}

is_valid_contrast <- function(dds, factor_var, g1, g2) {
  if (!factor_var %in% names(SummarizedExperiment::colData(dds))) return(FALSE)
  fac <- SummarizedExperiment::colData(dds)[[factor_var]]
  if (!is.factor(fac)) return(FALSE)
  all(c(g1, g2) %in% levels(fac)) && g1 != g2 && length(levels(fac)) >= 2
}


# Example safer comparisons builder (call where comparisons were created):
build_all_comparisons <- function(dds, candidate_factors) {
  out <- list()
  for (f in candidate_factors) {
    out <- c(out, valid_pairwise_contrasts(dds, f))
  }
  out
}

## Auto-generate interaction contrasts for two-way interaction models
## interaction_vars must be length 2: c(var1, var2)
## simple_effect_direction:
##   "var1" -> only simple effects of var2 at each non-ref level of var1
##   "var2" -> only simple effects of var1 at each non-ref level of var2
##   "both" -> (original behavior) both directions
##   "none" -> only the pure interaction terms
generate_interaction_contrasts <- function(dds,
                                           interaction_vars,
                                           simple_effect_direction = c("var1", "var2", "none", "both")) {
  simple_effect_direction <- match.arg(simple_effect_direction)
  if (length(interaction_vars) != 2) return(list())
  rn <- resultsNames(dds)
  var1 <- interaction_vars[1]
  var2 <- interaction_vars[2]

  # Patterns:
  # Main effect names look like VarLevel_vs_Ref (DESeq2 standard)
  main_pat_var1 <- paste0("^", var1, "_.+_vs_.+$")
  main_pat_var2 <- paste0("^", var2, "_.+_vs_.+$")

  main_var1 <- grep(main_pat_var1, rn, value = TRUE)
  main_var2 <- grep(main_pat_var2, rn, value = TRUE)

  # Interaction terms look like Var1Level.Var2Level (no "_vs_") (DESeq2 naming)
  inter_pat <- paste0("^", var1, "[^.]*\\.", var2, "[^.]*$")
  inter_terms <- grep(inter_pat, rn, value = TRUE)

  if (!length(inter_terms)) return(list())

  out <- list()

  for (iterm in inter_terms) {
    # Store direct interaction term (difference-of-differences)
    out[[length(out) + 1]] <- list(
      type = "interaction_term",
      name = iterm,
      label = paste0("INTERACTION:", iterm)
    )

    add_var1 <- simple_effect_direction %in% c("var2", "both")
    add_var2 <- simple_effect_direction %in% c("var1", "both")

    # Simple effects of var2 at each non-reference level of var1
    if (add_var2 && length(main_var2)) {
      for (mv2 in main_var2) {
        out[[length(out) + 1]] <- list(
          type = "simple_effect_var2_given_var1",
          contrast = list(c(iterm, mv2)),
          label = paste0("SIMPLE_EFFECT:", iterm, "_plus_", mv2)
        )
      }
    }
    # Simple effects of var1 at each non-reference level of var2
    if (add_var1 && length(main_var1)) {
      for (mv1 in main_var1) {
        out[[length(out) + 1]] <- list(
          type = "simple_effect_var1_given_var2",
          contrast = list(c(iterm, mv1)),
          label = paste0("SIMPLE_EFFECT:", iterm, "_plus_", mv1)
        )
      }
    }
  }

  # -----------------------------
  # De-duplicate semantically equivalent simple effects (order-sensitive)
  # Key definition:
  #  - interaction_term entries: "name:<resultsName>"
  #  - simple effects with contrast list: "contrast:" + checksum(ordered elements)
  if (length(out)) {
    if (!requireNamespace("digest", quietly = TRUE)) {
      message("...NOTE: 'digest' not installed; falling back to non-checksum keys for contrast de-dup (still order-sensitive).")
    }
    seen <- character(0)
    uniq <- list()
    for (cmp in out) {
      if (!is.null(cmp$name)) {
        key <- paste0("name:", cmp$name)
      } else if (!is.null(cmp$contrast)) {
        vec <- cmp$contrast[[1]]
        raw_key <- paste(vec, collapse = "|")      # preserve order
        if (requireNamespace("digest", quietly = TRUE)) {
          key <- paste0("contrast:", digest::digest(raw_key, algo = "xxhash64"))
        } else {
          key <- paste0("contrast:", raw_key)
        }
      } else {
        key <- paste0("other:", length(seen) + 1)
      }
      if (!(key %in% seen)) {
        seen <- c(seen, key)
        uniq[[length(uniq) + 1]] <- cmp
      }
    }
    if (length(uniq) < length(out)) {
      message("...Deduplicated ", length(out) - length(uniq), " redundant simple effect contrasts (order preserved).")
    }
    out <- uniq
  }
  out
}



extract_comparisons <- function(dds_processed,
                                comparisons_list,
                                analysis_name,
                                gene_map,
                                cfg,
                                gene_sets,
                                filt_name) {

  safe_analysis_name <- gsub("[^A-Za-z0-9_.-]", "_", analysis_name)
  graph_png_dir <- file.path(cfg$dir_graphs_png, safe_analysis_name)
  graph_pdf_dir <- file.path(cfg$dir_graphs_pdf, safe_analysis_name)
  if (!dir.exists(graph_png_dir)) dir.create(graph_png_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(graph_pdf_dir)) dir.create(graph_pdf_dir, recursive = TRUE, showWarnings = FALSE)

  csv_dir <- file.path(cfg$dir_data_csv, filt_name, safe_analysis_name)
  if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

  # Filter out invalid simple pairwise contrasts (factor/group fields) but keep name/contrast specs
  comparisons_list <- Filter(function(cmp) {
    if (!is.null(cmp$name) || !is.null(cmp$contrast)) return(TRUE)
    if (all(c("factor", "group1", "group2") %in% names(cmp))) {
      return(is_valid_contrast(dds_processed, cmp$factor, cmp$group1, cmp$group2))
    }
    FALSE
  }, comparisons_list)

  if (length(comparisons_list) == 0) {
    message("...No valid contrasts for ", analysis_name)
    return(list())
  }

  results_list <- list()

  for (comp in comparisons_list) {
    # Determine retrieval mode
    mode <- if (!is.null(comp$name)) "name" else if (!is.null(comp$contrast)) "list" else "pair"
    if (mode == "pair") {
      f <- comp$factor; g1 <- comp$group1; g2 <- comp$group2
      comp_id <- paste(analysis_name, f, g1, "vs", g2, sep = "_")
      message("...Running comparison (pair): ", comp_id)
      res_obj <- try(
        DESeq2::results(
          dds_processed,
          contrast = c(f, g1, g2),
          alpha = cfg$qval_threshold,
          cooksCutoff = FALSE,
          independentFiltering = TRUE
        ),
        silent = TRUE
      )
    } else if (mode == "name") {
      comp_id <- paste(analysis_name, comp$name, sep = "_")
      message("...Running comparison (name): ", comp_id)
      res_obj <- try(
        DESeq2::results(
          dds_processed,
          name = comp$name,
          alpha = cfg$qval_threshold,
          cooksCutoff = FALSE,
          independentFiltering = TRUE
        ),
        silent = TRUE
      )
    } else { # list contrast
      comp_id <- paste(analysis_name, comp$label %||% paste(comp$contrast[[1]], collapse = "_"), sep = "_")
      message("...Running comparison (contrast list): ", comp_id)
      res_obj <- try(
        DESeq2::results(
          dds_processed,
          contrast = comp$contrast,
          alpha = cfg$qval_threshold,
          cooksCutoff = FALSE,
          independentFiltering = TRUE
        ),
        silent = TRUE
      )
    }

    if (inherits(res_obj, "try-error")) {
      message("......FAILED contrast (skipping): ", comp_id,
              " | Error: ", conditionMessage(attr(res_obj, "condition")))
      next
    }

    res_df <- as.data.frame(res_obj) %>%
      tibble::rownames_to_column("gene_id") %>%
      dplyr::left_join(gene_map, by = "gene_id") %>%
      dplyr::mutate(
        gene_name = ifelse(is.na(gene_name), gene_id, gene_name),
        gene_name_upper = toupper(gene_name),
        human_ortholog = ifelse(is.na(human_ortholog), "", human_ortholog)
      )

    # Save CSV (all + sig)
    file_path_all <- file.path(csv_dir, paste0(cfg$date, "_", comp_id, "_DEG_results_all.csv"))
    write.csv(res_df, file_path_all, row.names = FALSE)

    sig_df <- res_df %>%
      dplyr::filter(!is.na(padj),
                    padj < cfg$qval_threshold,
                    abs(log2FoldChange) > cfg$lfc_threshold)
    file_path_sig <- file.path(csv_dir, paste0(cfg$date, "_", comp_id, "_DEG_results_sig.csv"))
    write.csv(sig_df, file_path_sig, row.names = FALSE)

    results_list[[comp_id]] <- res_df

    # (Existing plotting code block remains unchanged below)
    # ...
  }

  results_list
}
