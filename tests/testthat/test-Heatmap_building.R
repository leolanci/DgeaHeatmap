
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
  expect_snapshot(individual_matrix(factors_for_matrix, counts_data))

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
