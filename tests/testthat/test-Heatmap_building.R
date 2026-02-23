library(testthat)
library(stringr)
library(DgeaHeatmap)
library(tinysnapshot)
library(visualTest)

test_that("I can use the 3rd edition", {
  local_edition(3)
  expect_true(TRUE)
})

# basic check
test_that("basic test", {
  expect_true(TRUE)
})

test_that("Matrix is build from an input file and one column is put as rownames", {
  input_data <- data.frame(
    Gene = paste0("Gene", 1:3),
    Sample1 = c(1, 2, 3),
    Sample2 = c(4, 5, 6)
  )
  matrix_class <- c("matrix", "array")
  actual_outcome <- build_matrix(input_data, 1)

  # testing the outcome of the function to be "matrix" "array"
  expect_equal(class(actual_outcome), matrix_class)
  # testing the outcome to have strings as rownames
  expect_equal(class(rownames(actual_outcome)), "character")
})

test_that("An factor dependent individual matrix is build from a bigger matrix", {

  counts_data <- matrix(
    1:12,
    nrow = 3,
    dimnames = list(
      paste0("Gene", 1:3),
      c(
        "DKD_glomerulus_1",
        "DKD_glomerulus_2",
        "DKD_glomerulus_3",
        "DKD_glomerulus_4"
      )
    )
  )

  factors_for_matrix <- list("DKD", "glomerulus")

  actual_outcome <- individual_matrix(factors_for_matrix, counts_data)

  # testing the outcome of the function to be "matrix"
  expect_true(is.matrix(actual_outcome))

  # testing the outcome to have strings as rownames
  expect_equal(rownames(actual_outcome), rownames(counts_data))

  # Check correct number of selected columns
  expect_equal(ncol(actual_outcome), 4)
})


test_that("Variance of each gene is estimated and only the most variable genes are filtered from matrix", {
  set.seed(1)
  counts_data <- matrix(
    rnorm(50),
    nrow = 10,
    dimnames = list(paste0("Gene", 1:10), paste0("Sample", 1:5))
  )

  top_number_of_genes <- 3

  highly_variable_genes <- filtering_for_top_exprGenes(
    counts_data,
    top_number_of_genes
  )

  # Correct number of rows
  expect_equal(nrow(highly_variable_genes), top_number_of_genes)

  # Object type
  expect_true(is.matrix(highly_variable_genes))
})

test_that("Counts are scaled through Z-score scaling", {
  set.seed(1)
  counts_data <- matrix(
    rnorm(50, mean = 10, sd = 5),
    nrow = 10,
    dimnames = list(paste0("Gene", 1:10), paste0("Sample", 1:5))
  )

  scaled_counts <- scale_counts(counts_data)

  # Mean approximately 0 per gene
  expect_true(all(abs(rowMeans(scaled_counts)) < 1e-10))

  # SD approximately 1 per gene
  expect_true(all(abs(apply(scaled_counts, 1, sd) - 1) < 1e-10))
})


test_that("summarizing the biological replicates works", {
  counts_data <- matrix(
    rnorm(30),
    nrow = 5,
    dimnames = list(
      paste0("Gene", 1:5),
      c(
        "DKD_glomerulus_1",
        "DKD_glomerulus_2",
        "DKD_tubule_1",
        "DKD_tubule_2",
        "Other_1",
        "Other_2"
      )
    )
  )

  probes <- list("DKD_glomerulus", "DKD_tubule")

  actual_outcome <- summarise_bio_replicates(counts_data, probes)

  expect_true(is.matrix(actual_outcome))
  expect_equal(ncol(actual_outcome), 2)
  expect_equal(rownames(actual_outcome), rownames(counts_data))
})

test_that("K-mean generation works", {
  SumTable <- matrix(
    rnorm(20),
    nrow = 5,
    dimnames = list(
      paste0("Gene", 1:5),
      paste0("Group", 1:4)
    )
  )
  K_meanTable <- Kmean_generation(SumTable, 3)
  colNumSum <- ncol(SumTable)
  expected_outcome <- colNumSum + 1
  actual_outcome <- ncol(K_meanTable)
  # testing that a new column is added
  expect_equal(actual_outcome, expected_outcome)
  # testing the outcome of the function to be "matrix"
  actual_outcome <- Kmean_generation(SumTable, 1)
  expect_true(is.matrix(K_meanTable))
})

test_that("function makes list of most variable genes (rows) of each cluster", {
  set.seed(1)
  K_meanTable <- matrix(
    rnorm(20),
    nrow = 5,
    dimnames = list(
      paste0("Gene", 1:5),
      paste0("Group", 1:4)
    )
  )
  K_meanTable <- cbind(K_meanTable, cluster = c(1, 1, 1, 1, 1))
  actualOutput <- most_variable_genes(K_meanTable, 1, 2)
  expect_type(actualOutput, "list")
})

test_that("performing k_mean clustering outside of the heatmap works", {
  set.seed(1)
    SumTable <- matrix(
    rnorm(20),
    nrow = 5,
    dimnames = list(
      paste0("Gene", 1:5),
      paste0("Group", 1:4)
    )
  )
  actualOutput <- performing_kMeans(SumTable, 1)
  expect_type(actualOutput, "character")
})

test_that("basic test runs", {
  expect_true(TRUE)
})

test_that("function_complexHeatmap_var runs without error on valid input", {
  # Setup: Load sample data
  input_data <- read.csv(system.file("extdata/testfile_counts.csv",
                                     package = "DgeaHeatmap"
  ))

  x <- 1
  matrixCounts <- build_matrix(input_data, x)
  topGenes_matrix <- scale_counts(matrixCounts)

  # Define other arguments
  probes <- list("disease3", "disease4", "disease1B")
  k <- 1
  seed <- 1
  number_of_annotations_per_cluster <- 5
  Title <- "Heatmap of Data"
  fontsize_rowAnnotation <- 8
  fontsize_columnNames <- 6
  fontsize_rowNames <- 4
  title_heatmapLegend <- "Expression"
  WidthNum <- 4.5
  HeightNum <- 3
  UnitSize <- "cm"
  color_Palette <- "RdBu"

  set.seed(seed)

  # Test: run function silently
  expect_silent({
    suppressMessages({
      suppressWarnings({
        invisible(function_complexHeatmap_var(
          topGenes_matrix = topGenes_matrix,
          probes = probes,
          number_of_annotations_per_cluster = number_of_annotations_per_cluster,
          k = k,
          Title = Title,
          fontsize_rowAnnotation = fontsize_rowAnnotation,
          fontsize_columnNames = fontsize_columnNames,
          fontsize_rowNames = fontsize_rowNames,
          title_heatmapLegend = title_heatmapLegend,
          WidthNum = WidthNum,
          HeightNum = HeightNum,
          UnitSize = UnitSize,
          color_Palette = color_Palette
        ))
      })
    })
  })
})

test_that("get_heatmap_colors returns NULL when colorPalette is NULL", {
  result <- get_heatmap_colors(NULL)
  expect_null(result)
})

test_that("get_heatmap_colors returns a color vector for valid palette", {
  result <- get_heatmap_colors("RdBu")
  expect_type(result, "character") # Should be a character vector of colors
  expect_true(all(grepl("^#", result))) # Should be hex color codes
})

test_that("get_heatmap_colors fails gracefully on invalid palette", {
  expect_error(get_heatmap_colors("NotARealPalette"))
})

test_that("adv_Heatmap works with clustering disabled", {
  input_data <- read.csv(system.file("extdata/testfile_counts.csv",
    package = "DgeaHeatmap"
  ))
  matrixCounts <- build_matrix(input_data, 1)
  ncounts_matrix <- scale_counts(matrixCounts)

  result <- adv_Heatmap(ncounts_matrix, cluster_rows = FALSE, cluster_columns = FALSE)

  expect_s4_class(result, "HeatmapList")
})


test_that("row_clustering works with hierarchical clustering and no k_row", {
  # Setup
  input_data <- read.csv(system.file("extdata/testfile_counts.csv", package = "DgeaHeatmap"))
  matrixCounts <- build_matrix(input_data, 1)
  ncounts_matrix <- scale_counts(matrixCounts)

  # Run clustering
  result <- row_clustering(ncounts_matrix,
                           cluster_rows = TRUE,
                           cluster_method = "hierarchical", distance_method = "euclidean",
                           k_row = NULL, row_split = NULL, row_dend = TRUE
  )

  # Expectations
  expect_type(result, "list")
  expect_true("row_split" %in% names(result))
  expect_true("dend" %in% names(result))
  expect_null(result$row_split)
  expect_s3_class(result$dend, "dendrogram")
})

test_that("row_clustering works with hierarchical clustering and k_row", {
  input_data <- read.csv(system.file("extdata/testfile_counts.csv", package = "DgeaHeatmap"))
  matrixCounts <- build_matrix(input_data, 1)
  ncounts_matrix <- scale_counts(matrixCounts)

  result <- row_clustering(ncounts_matrix,
    cluster_rows = TRUE,
    cluster_method = "hierarchical", distance_method = "euclidean",
    k_row = 3, row_split = NULL, row_dend = TRUE
  )

  expect_type(result, "list")
  expect_s3_class(result$row_split, "factor")
  expect_equal(length(result$row_split), nrow(ncounts_matrix))
  expect_true(!is.null(result$dend))
})


test_that("get_dist returns correct type for euclidean distance", {
  input_data <- read.csv(system.file("extdata/testfile_counts.csv", package = "DgeaHeatmap"))
  matrixCounts <- build_matrix(input_data, 1)
  matrix <- scale_counts(matrixCounts)

  dist_matrix <- get_dist(matrix, "euclidean")

  expect_s3_class(dist_matrix, "dist")
  expect_equal(attr(dist_matrix, "method"), "euclidean")
})

test_that("get_dist returns correct type for manhattan distance", {
  input_data <- read.csv(system.file("extdata/testfile_counts.csv", package = "DgeaHeatmap"))
  matrixCounts <- build_matrix(input_data, 1)
  matrix <- scale_counts(matrixCounts)

  dist_matrix <- get_dist(matrix, "manhattan")

  expect_s3_class(dist_matrix, "dist")
  expect_equal(attr(dist_matrix, "method"), "manhattan")
})

test_that("get_dist returns correct structure for correlation distance", {
  input_data <- read.csv(system.file("extdata/testfile_counts.csv", package = "DgeaHeatmap"))
  matrixCounts <- build_matrix(input_data, 1)
  matrix <- scale_counts(matrixCounts)

  dist_matrix <- get_dist(matrix, "correlation")

  expect_s3_class(dist_matrix, "dist")
  expect_null(attr(dist_matrix, "method")) # correlation does not assign a method attr
  expect_equal(length(dist_matrix), nrow(matrix) * (nrow(matrix) - 1) / 2)
})

#' Setting the column clustering
#'
#' @param ncounts_matrix An input matrix to create the clustering.
#' @param cluster_columns A Boolean switching the optional clustering of the
#' rows on and off, default <- TRUE.
#' @param cluster_method A string setting the cluster method for the heatmap.
#' @param distance_method A string setting the distance method for clustering.
#' @param k_col An integer used to set number of clusters for row clustering.
#' @param col_split An integer seeting the column clustering by kmeans.
#' @param col_dend A Boolean switching column clustering by hierarchical
#' clustering on (default).
#'
#' @return col_split where the rows are split into clusters
#' @export
#'
#' @examples
#' x <- 1
#' seed <- 1
#' input_data <- read.csv(system.file("extdata/testfile_counts.csv",
#'   package = "DgeaHeatmap"
#' ))
#' matrixCounts <- build_matrix(input_data, x)
#' ncounts_matrix <- scale_counts(matrixCounts)
#' set.seed
#' col_split <- NULL
#' col_dend <- TRUE
#' columns_clustered <- column_clustering(ncounts_matrix,
#'   cluster_columns = TRUE, cluster_method = "hierarchical",
#'   distance_method = "euclidean", k_col = NULL, col_split = col_split,
#'   col_dend = col_dend
#' )
column_clustering <- function(ncounts_matrix, cluster_columns = TRUE,
                              cluster_method = "hierarchical",
                              distance_method = "euclidean", k_col = NULL,
                              col_split = NULL, col_dend = TRUE) {
  col_split <- col_split
  col_dend <- col_dend
  if (cluster_columns) {
    col_data <- t(ncounts_matrix)

    if (cluster_method == "hierarchical") {
      col_dist <- get_dist(col_data, distance_method)
      col_clust <- stats::hclust(col_dist, method = "average")
      col_dend <- stats::as.dendrogram(col_clust)
      if (!is.null(k_col)) {
        col_split <- stats::cutree(col_clust, k = k_col)
        col_split <- as.factor(col_split)
        col_dend <- TRUE
      }
    } else if (cluster_method == "kmeans") {
      col_split <- performing_kMeans(col_data, k_col)
    }
  }
  return(list(col_split = col_split, dend = col_dend))
}

test_that("set_sample_annotation returns NULL if sample_metadata is NULL", {
  result <- set_sample_annotation()
  expect_null(result)
})

test_that("set_sample_annotation returns a HeatmapAnnotation object with valid metadata", {
  # Build example input matrix
  input_data <- read.csv(system.file("extdata/testfile_counts.csv", package = "DgeaHeatmap"))
  matrixCounts <- build_matrix(input_data, 1)
  ncounts_matrix <- scale_counts(matrixCounts)

  # Create metadata
  groups <- c("3_DKD_glomerulus_Geometric_S", "1B_DKD_glomerulus_Geometric_S", "2B_DKD_glomerulus_WT")
  sample_names <- colnames(ncounts_matrix)
  group_assignment <- sapply(sample_names, function(sample) {
    matched <- groups[sapply(groups, function(g) grepl(g, sample))]
    if (length(matched) > 0) matched[1] else NA
  })
  sample_metadata <- data.frame(Group = group_assignment, row.names = sample_names)

  # Create color mapping
  group_colors <- list(Group = c(
    "3_DKD_glomerulus_Geometric_S" = "#1b9e77",
    "1B_DKD_glomerulus_Geometric_S" = "#7570b3",
    "2B_DKD_glomerulus_WT" = "#e7298a"
  ))

  # Run function
  annotation <- set_sample_annotation(
    sample_metadata = sample_metadata,
    annotation_colors = group_colors
  )

  # Check result
  expect_s4_class(annotation, "HeatmapAnnotation")
})
