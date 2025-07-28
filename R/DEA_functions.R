#' DGEALimma
#'
#' @description Performs a differential expression analysis based on limma voom
#'
#' @param rawCounts A matrix containing raw counts (integers).
#' @param metadata A dataframe containing the meta data information of the raw counts matrix.
#' @param grouping_columns A list of the columns names in the metadata file that are used to generate the groups for the DEA that are compared with each other.
#' @param comparisons An element of the class list, that contains optional information on specific comparisons that are to be conducted (simplifies the results).
#' @param prefix A string used as a prefix for the contrast_names that are generated.
#'
#' @return A list containing the voom data, fitting of the linear model, results of the contrasts, and the generated design.
#' @export
#'
#' @examples
#' rawCounts <- read.csv(system.file("extdata/RawDataExamplePackageNanostring.csv", package = "DgeaHeatmap"))
#' rawCounts <- build_matrix(rawCounts, 1)
#' metadata <- read.csv(system.file("extdata/MetaDataPackageNanostring.csv", package = "DgeaHeatmap"))
#' grouping_columns <- c("segment", "region", "class", "slide_name")
#' comparisons <- list(Comp1 = c("Geometric_Segment_glomerulus_DKD_disease3", "PanCK_tubule_DKD_disease4"), Comp2 = c("neg_tubule_DKD_disease4", "PanCK_tubule_DKD_disease4"))
#' prefix <- "DEA"
#' DGEAListResults <- DGEALimma(rawCounts, metadata, grouping_columns, comparisons, prefix)
DGEALimma <- function(
    rawCounts,
    metadata,
    grouping_columns,
    comparisons = NULL,        # named list of specific comparisons
    prefix = "DEA") {

  # Check input
  if (!all(grouping_columns %in% colnames(metadata))) {
    stop("Not all grouping_columns found in metadata")
  }
  # Combine grouping columns into one composite group
  comp <- apply(metadata[, grouping_columns, drop = FALSE], 1, paste, collapse = "_")
  comp <- gsub(" ", "_", comp)  # Replace spaces with underscores
  metadata$comp <- make.names(comp)  # Clean names for use in model.matrix

  # Create design matrix
  comp_factor <- factor(metadata$comp)
  design <- stats::model.matrix(~0 + comp_factor)
  colnames(design) <- levels(comp_factor)

  # Create all pairwise contrasts
  if (!is.null(comparisons)) {
    contrast_strings <- vapply(comparisons, function(groups) {
      paste0(groups[1], " - ", groups[2])
    }, character(1))
    contrast_names <- names(comparisons)
  } else {
    # fallback to all pairwise
    comb <- utils::combn(levels(comp_factor), 2)
    contrast_strings <- apply(comb, 2, function(x) paste0(x[1], " - ", x[2]))
    contrast_names <- apply(comb, 2, function(x) paste0(prefix, "_", x[1], "_vs_", x[2]))
  }

  # create contrasts matrix
  contrast_matrix <- limma::makeContrasts(contrasts = contrast_strings, levels = design)

  # generate DGEList object
  y <- edgeR::DGEList(counts = rawCounts)
  # calculate normalization factors
  y <- edgeR::calcNormFactors(y)
  # voom transformation
  v <- limma::voom(y, design)
  # fit linear model
  fit <- limma::lmFit(v, design)
  fit2 <- limma::contrasts.fit(fit, contrast_matrix) # contrasts.fit must be run before eBayes
  fit2 <- limma::eBayes(fit2)
  # Return results for each contrast in a list
  results <- lapply(colnames(contrast_matrix), function(coef_name) {
    limma::topTable(fit2, coef = coef_name, number = Inf, sort.by = "P")
  })
  names(results) <- contrast_names

  return(list(
    voom_data = v,
    fit = fit2,
    results = results,
    design = design,
    normFactors = y
  ))

}

#' DGEADESeq2
#'
#' @description Performs a differential expression analysis based on DESeq2
#'
#' @param rawCounts A matrix containing raw counts (integers).
#' @param metadata A dataframe containing the meta data information of the raw counts matrix.
#' @param grouping_columns A list of the columns names in the metadata file that are used to generate the groups for the DEA that are compared with each other.
#' @param comparisons An element of the class list, that contains optional information on specific comparisons that are to be conducted (simplifies the results).
#'
#' @return A list containing the voom data, fitting of the linear model, results of the contrasts, and the generated design.
#' @export
#'
#' @examples
#' rawCounts <- read.csv(system.file("extdata/RawDataExamplePackageNanostring.csv", package = "DgeaHeatmap"))
#' rawCounts <- build_matrix(rawCounts, 1)
#' metadata <- read.csv(system.file("extdata/MetaDataPackageNanostring.csv", package = "DgeaHeatmap"))
#' grouping_columns <- c("segment", "region", "class", "slide_name")
#' comparisons <- list(Comp1 = c("Geometric_Segment_glomerulus_DKD_disease3", "PanCK_tubule_DKD_disease4"), Comp2 = c("neg_tubule_DKD_disease4", "PanCK_tubule_DKD_disease4"))
#' storage.mode(rawCounts) <- "integer"
#' sum(!is.finite(as.matrix(rawCounts)))
#' results_DESeq2 <- DGEADESeq2(rawCounts, metadata, grouping_columns,comparisons)

DGEADESeq2 <- function( rawCounts, metadata, grouping_columns, comparisons ) {
  # Filter genes with non-zero counts in at least n samples (e.g., n = 2)
  keep <- rowSums(rawCounts > 0) >= 2
  filtered_counts <- rawCounts[keep, ]

  # Combine grouping columns into a composite factor
  comp <- apply(metadata[, grouping_columns, drop = FALSE], 1, paste, collapse = "_")
  comp <- gsub(" ", "_", comp)
  metadata$comp <- make.names(comp)
  metadata$comp <- factor(metadata$comp)

  # Build DESeq2 dataset
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = filtered_counts,
    colData = metadata,
    design = ~ 0 + comp
  )

  # Run DESeq2 pipeline
  dds <- DESeq2::DESeq(dds, sfType = "poscounts")

  results_list <- list()

  for (name in names(comparisons)) {
    contrast_pair <- trimws(comparisons[[name]])

    # Check levels exist in your comp factor:
    if (all(contrast_pair %in% levels(SummarizedExperiment::colData(dds)$comp))) {
      res <- tryCatch({
        DESeq2::results(dds, contrast = c("comp", contrast_pair[1], contrast_pair[2]))
      }, error = function(e) {
        warning(paste("Contrast failed:", name, "->", e$message))
        NULL
      })

      if (!is.null(res)) {
        results_list[[name]] <- as.data.frame(res)
      }
    } else {
      missing_levels <- setdiff(contrast_pair, levels(SummarizedExperiment::colData(dds)$comp))
      message(
        sprintf(
          "Contrast failed: %s -> One or both levels not found in comp factor levels: %s",
          name,
          paste(missing_levels, collapse = ", ")
        )
      )

    }
  }
  contrast_names <- names(comparisons)
  return(list(
    dds = dds,
    results = results_list,
    contrasts = contrast_names
  ))
}

#' extractDEGenes
#'
#' @description Extracts the results of the DEA with DESeq2
#'
#' @param results_list A list containing the results of the DEA.
#' @param comparisons An element of the class list, that contains optional information on specific comparisons that are to be conducted (simplifies the results).
#' @param only_up A Boolean opting for only up regulated genes, default = FALSE.
#' @param only_down A Boolean opting for only down regulated genes, default = FALSE.
#' @param up_down A Boolean opting for only up and down regulated genes, default = FALSE.
#' @param only_sig A Boolean opting for only significantly regulated genes, default = FALSE.
#' @param padj_cutoff A Float, setting the adjusted p-values cutoff, default = 0.05.
#' @param lfc_cutoff A Float, setting the fold change cutoff, default = 0.
#'
#' @return A list containing the voom data, fitting of the linear model, results of the contrasts, and the generated design.
#' @export
#'
#' @examples
#' rawCounts <- read.csv(system.file("extdata/RawDataExamplePackageNanostring.csv", package = "DgeaHeatmap"))
#' rawCounts <- build_matrix(rawCounts, 1)
#' metadata <- read.csv(system.file("extdata/MetaDataPackageNanostring.csv", package = "DgeaHeatmap"))
#' grouping_columns <- c("segment", "region", "class", "slide_name")
#' comparisons <- list(Comp1 = c("Geometric_Segment_glomerulus_DKD_disease3", "PanCK_tubule_DKD_disease4"), Comp2 = c("neg_tubule_DKD_disease4", "PanCK_tubule_DKD_disease4"))
#' storage.mode(rawCounts) <- "integer"
#' sum(!is.finite(as.matrix(rawCounts)))
#' results_list_DESeq2 <- DGEADESeq2(rawCounts, metadata, grouping_columns,comparisons)
#' results_list <- results_list_DESeq2$results
#' down_genes <- extractDEGenes(results_list, comparisons, only_down = TRUE)

extractDEGenes <- function(results_list,
                           comparisons,
                           only_up = FALSE,
                           only_down = FALSE,
                           up_down = FALSE,
                           only_sig = FALSE,
                           padj_cutoff = 0.05,
                           lfc_cutoff = 0) {

  # Check for incompatible options
  if ((only_up + only_down + up_down + only_sig) > 1) {
    stop("Only one of 'only_up', 'only_down', 'up_down', or 'only_sig' can be TRUE at a time.")
  }

  # Get union of all genes across all comparisons
  all_genes <- unique(unlist(lapply(results_list, rownames)))

  if (up_down) {
    # Build 3-state matrix (genes x comparisons)
    state_matrix <- matrix(0, nrow = length(all_genes), ncol = length(results_list))
    rownames(state_matrix) <- all_genes
    colnames(state_matrix) <- names(results_list)

    for (comparisons in names(results_list)) {
      res <- results_list[[comparisons]]
      if (!is.data.frame(res)) res <- as.data.frame(res)
      if (!all(c("padj", "log2FoldChange") %in% colnames(res))) {
        warning(paste("Missing expected columns in", comparison))
        next
      }
      genes_in_res <- rownames(res)

      for (gene in genes_in_res) {
        padj <- res[gene, "padj"]
        lfc <- res[gene, "log2FoldChange"]

        if (!is.na(padj) && padj < padj_cutoff && !is.na(lfc)) {
          if (lfc > lfc_cutoff) {
            state_matrix[gene, comparisons] <- 1
          } else if (lfc < -lfc_cutoff) {
            state_matrix[gene, comparisons] <- -1
          }
        }
      }
    }
    return(state_matrix)

  } else {
    # Initialize list to store gene vectors
    gene_lists <- list()

    for (comparison in names(results_list)) {
      res <- results_list[[comparison]]

      # Filter by significance if requested
      sig_res <- if (only_sig || only_up || only_down) {
        res[!is.na(res$padj) & res$padj < padj_cutoff, ]
      } else {
        res
      }

      # Extract genes based on requested direction
      if (only_up) {
        genes <- rownames(sig_res)[sig_res$log2FoldChange > lfc_cutoff]
      } else if (only_down) {
        genes <- rownames(sig_res)[sig_res$log2FoldChange < -lfc_cutoff]
      } else if (only_sig) {
        genes <- rownames(sig_res)
      } else {
        # If no filtering requested, return all genes
        genes <- rownames(res)
      }

      gene_lists[[comparison]] <- genes
    }
    return(gene_lists)
  }
}


#' DGEAedgeR
#'
#' @description Performs a differential expression analysis based on edgeR
#'
#' @param rawCounts A matrix containing raw counts (integers).
#' @param metadata A dataframe containing the meta data information of the raw counts matrix.
#' @param grouping_columns A list of the columns names in the metadata file that are used to generate the groups for the DEA that are compared with each other.
#' @param comparisons An element of the class list, that contains optional information on specific comparisons that are to be conducted (simplifies the results).
#' @param prefix A string used as a prefix for the contrast_names that are generated.
#'
#' @return A list containing the voom data, fitting of the linear model, results of the contrasts, and the generated design.
#' @export
#'
#' @examples
#' rawCounts <- read.csv(system.file("extdata/RawDataExamplePackageNanostring.csv", package = "DgeaHeatmap"))
#' rawCounts <- build_matrix(rawCounts, 1)
#' metadata <- read.csv(system.file("extdata/MetaDataPackageNanostring.csv", package = "DgeaHeatmap"))
#' grouping_columns <- c("segment", "region", "class", "slide_name")
#' comparisons <- list(Comp1 = c("Geometric_Segment_glomerulus_DKD_disease3", "PanCK_tubule_DKD_disease4"), Comp2 = c("neg_tubule_DKD_disease4", "PanCK_tubule_DKD_disease4"))
#' prefix <- "DEA"
#' results_edgeR <- DGEAedgeR(rawCounts, metadata, grouping_columns, comparisons, prefix = "DEA")

DGEAedgeR <- function( rawCounts,
                       metadata,
                       grouping_columns,
                       comparisons = NULL,    # named list of specific comparisons; e.g. list(comp1_vs_comp2 = c("comp1", "comp2"))
                       prefix = "DEA"
) {
  # Check grouping columns in metadata
  if (!all(grouping_columns %in% colnames(metadata))) {
    stop("Some of the grouping_columns not found in metadata")
  }

  # Create composite group factor
  comp <- apply(metadata[, grouping_columns, drop = FALSE], 1, paste, collapse = "_")
  comp <- gsub(" ", "_", comp)  # replace spaces with underscores
  metadata$comp <- factor(comp)

  # Create DGEList object
  y <- edgeR::DGEList(counts = rawCounts, group = metadata$comp)

  # Filter lowly expressed genes (optional, adjust thresholds as needed)
  keep <- edgeR::filterByExpr(y, group = metadata$comp)
  y <- y[keep, , keep.lib.sizes = FALSE]

  # Calculate normalization factors
  y <- edgeR::calcNormFactors(y)

  # Create design matrix without intercept (one column per group)
  design <- stats::model.matrix(~0 + metadata$comp)
  colnames(design) <- levels(metadata$comp)

  # Fit negative binomial GLM
  y <- edgeR::estimateDisp(y, design)
  fit <- edgeR::glmFit(y, design)

  # delete trailing spaces from comparisons
  comparisons <- lapply(comparisons, trimws)

  # Prepare contrasts matrix
  if (!is.null(comparisons)) {
    # named list of pairs, e.g. list(comp1_vs_comp2 = c("comp1", "comp2"))
    contrast_list <- lapply(comparisons, function(x) {
      contrast_vec <- rep(0, length(levels(metadata$comp)))
      names(contrast_vec) <- levels(metadata$comp)
      contrast_vec[x[1]] <- 1
      contrast_vec[x[2]] <- -1
      contrast_vec
    })
    contrast_matrix <- do.call(cbind, contrast_list)
    colnames(contrast_matrix) <- names(comparisons)
  } else {
    # all pairwise contrasts
    comb <- utils::combn(levels(metadata$comp), 2)
    contrast_matrix <- apply(comb, 2, function(x) {
      contrast_vec <- rep(0, length(levels(metadata$comp)))
      names(contrast_vec) <- levels(metadata$comp)
      contrast_vec[x[1]] <- 1
      contrast_vec[x[2]] <- -1
      contrast_vec
    })
    colnames(contrast_matrix) <- apply(comb, 2, function(x) paste0(prefix, "_", x[1], "_vs_", x[2]))
  }


  # Run likelihood ratio tests for each contrast and collect results
  results <- list()
  for (i in seq_len(ncol(contrast_matrix))) {
    ct <- contrast_matrix[, i]
    lrt <- edgeR::glmLRT(fit, contrast = ct)
    top <- edgeR::topTags(lrt, n = Inf)$table
    results[[colnames(contrast_matrix)[i]]] <- top
  }

  return(list(
    dgeList = y,
    design = design,
    fit = fit,
    contrast_matrix = contrast_matrix,
    results = results
  ))
}

#' Summarize the results of the DEA with edgeR
#'
#' @param results_edgeR A list containing the results of the DEA.
#' @param lfc_threshold A Float, setting the threshold value for the fold change, default = 1.
#' @param fdr_threshold A Float, setting the False Discovery Rate threshold, default = 0.05.
#' @param return_significant_genes A Boolean opting for the output of significantly regulated genes, default = TRUE
#'
#' @return A list containing the voom data, fitting of the linear model, results of the contrasts, and the generated design.
#' @export
#'
#' @examples
#' rawCounts <- read.csv(system.file("extdata/RawDataExamplePackageNanostring.csv", package = "DgeaHeatmap"))
#' rawCounts <- build_matrix(rawCounts, 1)
#' metadata <- read.csv(system.file("extdata/MetaDataPackageNanostring.csv", package = "DgeaHeatmap"))
#' grouping_columns <- c("segment", "region", "class", "slide_name")
#' comparisons <- list(Comp1 = c("Geometric_Segment_glomerulus_DKD_disease3", "PanCK_tubule_DKD_disease4"), Comp2 = c("neg_tubule_DKD_disease4", "PanCK_tubule_DKD_disease4"))
#' prefix <- "DEA"
#' results_edgeR <- DGEAedgeR(rawCounts, metadata, grouping_columns, comparisons, prefix = "DEA")
#' sumresults <- summarize_edgeR_DEA(results_edgeR)
summarize_edgeR_DEA <- function(results_edgeR,
                                lfc_threshold = 1,
                                fdr_threshold = 0.05,
                                return_significant_genes = TRUE) {
  all_results <- results_edgeR$results
  gene_classifications <- list()
  summary_table <- data.frame()

  for (contrast in names(all_results)) {
    res <- all_results[[contrast]] %>%
      as.data.frame() %>%
      dplyr::mutate(
        decision = dplyr::case_when(
          FDR < fdr_threshold & logFC > lfc_threshold ~ 1,
          FDR < fdr_threshold & logFC < -lfc_threshold ~ -1,
          TRUE ~ 0
        )
      )

    gene_classifications[[contrast]] <- res

    summary_counts <- table(factor(res$decision, levels = c(-1, 0, 1)))
    summary_table <- rbind(summary_table,
                           data.frame(
                             contrast = contrast,
                             down = summary_counts["-1"],
                             not_sig = summary_counts["0"],
                             up = summary_counts["1"]
                           ))
  }

  # Optionally extract significant gene names
  if (return_significant_genes) {
    sig_genes <- lapply(gene_classifications, function(df) {
      df_sig <- df[df$decision != 0, ]
      rownames(df_sig)
    })
  } else {
    sig_genes <- NULL
  }

  return(list(
    classified_results = gene_classifications,
    summary = summary_table,
    significant_gene_names = sig_genes
  ))
}
