library(testthat)
library(DgeaHeatmap)
library(tinysnapshot)
library(visualTest)

test_that("I can use the 3rd edition", {
  local_edition(3)
  expect_true(TRUE)
})

test_that("Matrix is build from an input file and one column is put as rownames", {

  input_data <- read.csv(test_path("1_Counts_All_Regions_All_Ages_Pos_Neg.csv"))
  x <- 1
  matrix_class <- c("matrix","array")
  actual_outcome <- build_matrix(input_data, x)

  # testing the outcome of the function to be "matrix" "array"
  expect_equal(class(actual_outcome), matrix_class)
  # testing the outcome to have strings as rownames
  expect_equal(class(rownames(actual_outcome)), "character")
})

test_that("An factor dependent individual matrix is build from a bigger matrix", {
  input_data <- read.csv(test_path("1_Counts_All_Regions_All_Ages_Pos_Neg.csv"))
  counts_data <- build_matrix(input_data, 1)
  factors_for_matrix <- list("hippo", "pos")
  matrix_class <- c("matrix","array")
  actual_outcome <- individual_matrix(factors_for_matrix, counts_data)
  # testing the outcome of the function to be "matrix" "array"
  expect_equal(class(actual_outcome), matrix_class)
  # testing the outcome to have strings as rownames
  expect_equal(class(rownames(actual_outcome)), "character")
  # testing if number of columns that contain string in original matrix equals number of column in output matrix
  actual_output <- ncol(actual_outcome)
  expected_output <- 9
  expect_equal(actual_output, expected_output)

})


test_that("Variance of each gene is estimated and only the most variable genes are filtered from matrix", {
  input_data <- read.csv(test_path("1_Counts_All_Regions_All_Ages_Pos_Neg.csv"))
  counts_data <- build_matrix(input_data, 1)
  top_number_of_genes <- 500
  highly_variable_genes <- filtering_for_top_exprGenes(counts_data, top_number_of_genes)

  # testing for correct number of rows
  actual_outcome <- nrow(highly_variable_genes)
  expect_equal(actual_outcome, top_number_of_genes)

  # testing for correct class of outcome
  actual_outcome_class <- class(highly_variable_genes)
  matrix_class <- c("matrix","array")
  expect_equal(actual_outcome_class, matrix_class)


})

test_that("Counts are scaled through Z-score scaling", {
  input_data <- read.csv(test_path("1_Counts_All_Regions_All_Ages_Pos_Neg.csv"))
  counts_data <- build_matrix(input_data, 1)
  scaled_counts <- scale_counts(counts_data)
  if(any(scaled_counts > 2)){
    if(any(scaled_counts < (-2))){
      actual_outcome <- FALSE
    }
    else {
      actual_outcome <- TRUE
    }
  }
  else {
    actual_outcome <- TRUE
  }

  expected_outcome <- FALSE

  expect_equal(actual_outcome, expected_outcome)
})

test_that("data ist distributed equally", {
  #library("tinysnapshot")
  require("tinysnapshot")
  require("visualTest")
  local_edition(3)
  input_data <- read.csv(test_path("1_Counts_All_Regions_All_Ages_Pos_Neg.csv"))
  counts_data <- build_matrix(input_data, 1)
  scaled_counts <- scale_counts(counts_data)
  png(filename = "testShowDataDistribution.png")
  p1 <- show_data_distribution(scaled_counts)
  dev.off()
  tinysnapshot::expect_snapshot_plot(p1 ,label = "test7.pdf")
  #compare_file_text("test7.pdf", "test6.pdf")
  #getFingerprint(file = "testShowDataDistribution.png")
  #getFingerprint(file = "testShowDataDistributionControl.png")
  #isSimilar(file = "testShowDataDistributionControl.png",
   #         fingerprint = getFingerprint("testShowDataDistribution.png"),
    #        threshold = 0.1)

})

#test_that("elbow plot works",{
#  input_data <- read.csv(test_path("1_Counts_All_Regions_All_Ages_Pos_Neg.csv"))
#  counts_data <- build_matrix(input_data, 1)
#  scaled_counts <- scale_counts(counts_data)
#  p1 <- elbow_plot(2, scaled_counts)
#  tinysnapshot::expect_snapshot_plot(p1, label = "test8.pdf")
#})

test_that("summarizing the biological replicates works", {
  input_data <- read.csv(test_path("1_Counts_All_Regions_All_Ages_Pos_Neg.csv"))
  counts_data <- build_matrix(input_data, 1)
  factors_for_matrix <- list("hippo", "pos")
  matrix_class <- c("matrix","array")
  indi_matrix <- individual_matrix(factors_for_matrix, counts_data)
  probes <- list("P0_hippo_Iba1_pos", "P5_hippo_Iba1_pos", "P15_hippo_Iba1_pos","P0_hippo_Iba1_neg", "P5_hippo_Iba1_neg","P15_hippo_Iba1_neg")
  actual_outcome <- summarise_bio_replicates(indi_matrix, probes)
  # testing the outcome of the function to be "matrix" "array"
  expect_equal(class(actual_outcome), matrix_class)
  # testing the outcome to have strings as rownames
  expect_equal(class(rownames(actual_outcome)), "character")
  # testing if number of columns that contain string in original matrix equals number of column in output matrix
  actual_output <- ncol(actual_outcome)
  expected_output <- 6
  expect_equal(actual_output, expected_output)
})

test_that("K-mean generation works",{
  input_data <- read.csv(test_path("1_Counts_All_Regions_All_Ages_Pos_Neg.csv"))
  counts_data <- build_matrix(input_data, 1)
  factors_for_matrix <- list("hippo", "pos")
  matrix_class <- c("matrix","array")
  indi_matrix <- individual_matrix(factors_for_matrix, counts_data)
  probes <- list("P0_hippo_Iba1_pos", "P5_hippo_Iba1_pos", "P15_hippo_Iba1_pos")
  SumTable <- summarise_bio_replicates(indi_matrix, probes)
  K_meanTable <- Kmean_generation(SumTable, 1, 3)
  colNumSum <- ncol(SumTable)
  expected_outcome <- colNumSum +1
  actual_outcome <- ncol(K_meanTable)
  # testing that a new column is added
  expect_equal(actual_outcome, expected_outcome)
  # testing the outcome of the function to be "matrix" "array"
  actual_outcome <- Kmean_generation(SumTable, 1, 3)
  expect_equal(class(actual_outcome), matrix_class)

})

test_that("function makes list of most variable genes (rows) of each cluster", {
  input_data <- read.csv(test_path("1_Counts_All_Regions_All_Ages_Pos_Neg.csv"))
  counts_data <- build_matrix(input_data, 1)
  factors_for_matrix <- list("hippo", "pos")
  matrix_class <- c("matrix","array")
  indi_matrix <- individual_matrix(factors_for_matrix, counts_data)
  probes <- list("P0_hippo_Iba1_pos", "P5_hippo_Iba1_pos", "P15_hippo_Iba1_pos")
  SumTable <- summarise_bio_replicates(indi_matrix, probes)
  K_meanTable <- Kmean_generation(SumTable, 1, 3)
  actualOutput <- most_variable_genes(K_meanTable, 5, 3)
  # testing the outcome is a list
  expected_outcome <- "list"
  expect_equal(class(actualOutput), expected_outcome)

})

test_that("performing k_mean clustering outside of the heatmap works", {
  input_data <- read.csv(test_path("1_Counts_All_Regions_All_Ages_Pos_Neg.csv"))
  counts_data <- build_matrix(input_data, 1)
  factors_for_matrix <- list("hippo", "pos")
  matrix_class <- c("matrix","array")
  indi_matrix <- individual_matrix(factors_for_matrix, counts_data)
  probes <- list("P0_hippo_Iba1_pos", "P5_hippo_Iba1_pos", "P15_hippo_Iba1_pos")
  SumTable <- summarise_bio_replicates(indi_matrix, probes)
  K_meanTable <- Kmean_generation(SumTable, 1, 3)
  actualOutput <- performing_kMeans(SumTable, 3)
  # testing that outcomes is right class
  expect_outcome <- "character"
  expect_equal(class(actualOutput), expect_outcome)


})


