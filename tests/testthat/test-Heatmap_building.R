library(testthat)
library(DgeaHeatmap)
library(tinysnapshot)

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

test_that("elbow plot works",{
  input_data <- read.csv(test_path("1_Counts_All_Regions_All_Ages_Pos_Neg.csv"))
  counts_data <- build_matrix(input_data, 1)
  scaled_counts <- scale_counts(counts_data)
    #elbow_plot(2, scaled_counts)
  #dev.off
  #expected_output <- getFingerprint(file = "testElbowPlot.pdf")
  #actual_output <- getFingerprint(file = "testElbowPlotTrial.pdf")
  #isSimilar(file = "testElbowPlot.pdf",
   #         fingerprint = getFingerprint(file = "testElbowPlotTrial.pdf"))
  #expect_equal(actual_output, expected_output)
})
#input_data <- read.csv(test_path("1_Counts_All_Regions_All_Ages_Pos_Neg.csv"))
#counts_data <- build_matrix(input_data, 1)
#del_layer6 <- grep("layer", colnames(counts_data))
#counts_data <- counts_data[, -del_layer6]
#factors_for_matrix <- list("hippo", "pos")
#indi_matrix <- individual_matrix(factors_for_matrix, counts_data)
#top_number_of_genes <- 500
#highly_variable_genes <- filtering_for_top_exprGenes(indi_matrix, top_number_of_genes)
#scaled_counts <- scale_counts(highly_variable_genes)
#png("testElbowPlotControl.png")
#elbow_plot(1, scaled_counts)
#dev.off()
input_data <- read.csv(test_path("1_Counts_All_Regions_All_Ages_Pos_Neg.csv"))
factors_for_matrix <- list("hippo", "pos")
indi_matrix <- individual_matrix(factors_for_matrix, input_data)
counts_data <- build_matrix(indi_matrix, 1)
scaled_counts <- scale_counts(counts_data)
elbow_plot(2, scaled_counts)
