#' Performs a differential expression analysis based on limma voom
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
#' rawCounts <- matrix
#' metadata <- data.frame
#' grouping_columns <- c("segment", "region", "class", "slide_name")
#' comparisons <- list("Comp1 = c("Group_a", "Group_b"), "Comp2" = c("Group_c", "Group_d"))
#' prefix = "DEA"
#'
DGEALimma <- function(
    rawCounts,
    metadata,
    grouping_columns,
    comparisons = NULL,        # named list of specific comparisons
    prefix = "DEA") {

  # Check input
  if (!all(grouping_columns %in% colnames(metadata))) {
    stop("Some grouping_columns not found in metadata")
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
    contrast_strings <- sapply(comparisons, function(groups) {
      paste0(groups[1], " - ", groups[2])
    })
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
  y = edgeR::DGEList(counts = rawCounts)
  # calculate normalization factors
  y = edgeR::calcNormFactors(y)
  # voom transformation
  v = limma::voom(y, design)
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

#' Performs a differential expression analysis based on DESeq2
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
#' rawCounts <- matrix
#' metadata <- data.frame
#' grouping_columns <- c("segment", "region", "class", "slide_name")
#' comparisons <- list("Comp1 = c("Group_a", "Group_b"), "Comp2" = c("Group_c", "Group_d"))
#'
DGEA_DESeq2 <- function (rawCounts, metadata, grouping_columns, comparisons){
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
      warning(paste("Contrast failed:", name, "-> One or both levels not found in comp factor levels"))
      print(setdiff(contrast_pair, levels(colData(dds)$comp)))
    }
  }
  contrast_names <- names(comparisons)
  return(list(
    dds = dds,
    results = results_list,
    contrasts = contrast_names
  ))
}
