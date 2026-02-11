library(testthat)
library(stringr)
library(DgeaHeatmap)
library(tinysnapshot)
library(visualTest)
library(Biobase)

test_that("I can use the 3rd edition", {
  local_edition(3)
  expect_true(TRUE)
})

# basic check
test_that("basic test", {
  expect_true(TRUE)
})


test_that("Build a matrix from an input csv file, with one column as rownames", {
  skip_if_not_installed("GeoMxWorkflows")
  skip_if_not_installed("GeomxTools")

  datadir <- system.file("extdata", "WTA_NGS_Example", package = "GeoMxWorkflows")
  DCCFiles <- dir(file.path(datadir, "dccs"), pattern = ".dcc$", full.names = TRUE, recursive = TRUE)
  PKCFiles <- unzip(zipfile = dir(file.path(datadir, "pkcs"), pattern = ".zip$", full.names = TRUE, recursive = TRUE),
                    exdir = tempdir()
                    )
  SampleAnnotationFile <- dir(file.path(datadir, "annotation"), pattern = ".xlsx$", full.names = TRUE, recursive = TRUE)

  demoData <- suppressWarnings(
    GeomxTools::readNanoStringGeoMxSet(
      dccFiles = DCCFiles,
      pkcFiles = PKCFiles,
      phenoDataFile = SampleAnnotationFile,
      phenoDataSheet = "Template",
      phenoDataDccColName = "Sample_ID",
      protocolDataColNames = c("aoi", "roi"),
      configFile = NULL,
      analyte = "RNA",
      phenoDataColPrefix = "",
      experimentDataColNames = NULL
    )
  )

  # Apply your function
  updated_data <- add_demoElem(demoData)

  # Check that "demoElem" was added
  expect_true("demoElem" %in% Biobase::assayDataElementNames(updated_data))

  # Check that demoElem is log10 of exprs (at a few positions)
  exprs_values <- Biobase::assayDataElement(demoData, "exprs")
  demoElem_values <- Biobase::assayDataElement(updated_data, "demoElem")

  # Compare first 3 rows and 2 columns
  expect_equal(
    demoElem_values[1:3, 1:2],
    log10(exprs_values[1:3, 1:2])
  )
})

test_that("Splitting data by group column with feature, pheno or protocol data to then get the mean.", {
  skip_if_not_installed("GeoMxWorkflows")
  skip_if_not_installed("GeomxTools")
  skip_if_not_installed("NanoStringNCTools")

  # Load example data
  datadir <- system.file("extdata", "WTA_NGS_Example",
                         package = "GeoMxWorkflows")
  DCCFiles <- dir(file.path(datadir, "dccs"),
                  pattern = ".dcc$", full.names = TRUE, recursive = TRUE)
  PKCFiles <- unzip(zipfile = dir(file.path(datadir, "pkcs"),
                                  pattern = ".zip$", full.names = TRUE, recursive = TRUE),
                    exdir = tempdir()
                    )
  SampleAnnotationFile <- dir(file.path(datadir, "annotation"),
                              pattern = ".xlsx$", full.names = TRUE, recursive = TRUE)

  rawDataObject <- suppressWarnings(
    GeomxTools::readNanoStringGeoMxSet(dccFiles = DCCFiles,
                                       pkcFiles = PKCFiles,
                                       phenoDataFile = SampleAnnotationFile,
                                       phenoDataSheet = "Template",
                                       phenoDataDccColName = "Sample_ID",
                                       protocolDataColNames = c("aoi", "roi"),
                                       configFile = NULL,
                                       analyte = "RNA",
                                       phenoDataColPrefix = "",
                                       experimentDataColNames = NULL))

  # Add a dummy expression matrix for testing
  demoData <- rawDataObject
  assayDataElement <- matrix(rnorm(nrow(demoData) * ncol(demoData)),
                             nrow = nrow(demoData),
                             ncol = ncol(demoData))
  dimnames(assayDataElement) <- dimnames(Biobase::exprs(demoData))
  assayDataElement(demoData, "demoElt") <- assayDataElement

  # Set test parameters
  vGroup <- "aoi"
  vElt <- "demoElt"

  # Run the function
  result <- split_data_by_column(demoData, vGroup, vElt)

  # Check that original object is returned
  expect_identical(result, demoData)

  # Optionally, capture the output of esBy to test side effect
  grouped_result <- NanoStringNCTools::esBy(demoData,
                                            GROUP = vGroup,
                                            FUN = function(x) {
                                              NanoStringNCTools::assayDataApply(x, MARGIN = 1, FUN = mean, elt = vElt)
                                            })

})

test_that("Function for automatized quality control.", {
  # Load example data
  datadir <- system.file("extdata", "WTA_NGS_Example", package = "GeoMxWorkflows")
  dcc_files <- dir(file.path(datadir, "dccs"), pattern = ".dcc$", full.names = TRUE)
  pkc_file <- unzip(zipfile = dir(file.path(datadir, "pkcs"), pattern = ".zip$", full.names = TRUE),
                    exdir = tempdir()
                    )
  ann_file <- dir(file.path(datadir, "annotation"), pattern = ".xlsx$", full.names = TRUE)

  demoData <- suppressWarnings(GeomxTools::readNanoStringGeoMxSet(
    dccFiles = dcc_files,
    pkcFiles = pkc_file,
    phenoDataFile = ann_file,
    phenoDataSheet = "Template",
    phenoDataDccColName = "Sample_ID",
    protocolDataColNames = c("aoi", "roi"),
    analyte = "RNA"
  ))

  # Run QC function
  result <- aExprsDataQC(demoData, vFlags = "QCFlags")

  # Tests
  expect_s4_class(result, "NanoStringGeoMxSet")
  expect_true(ncol(result) <= ncol(demoData)) # QC filtering may remove samples
})


test_that("genRawReadCountTable returns a data frame of raw counts", {
  # Load example data
  datadir <- system.file("extdata", "WTA_NGS_Example", package = "GeoMxWorkflows")
  DCCFiles <- dir(file.path(datadir, "dccs"), pattern = ".dcc$", full.names = TRUE)
  PKCFiles <- unzip(zipfile = dir(file.path(datadir, "pkcs"), pattern = ".zip$", full.names = TRUE),
                    exdir = tempdir()
                    )
  SampleAnnotationFile <- dir(file.path(datadir, "annotation"),
                              pattern = ".xlsx$", full.names = TRUE)

  # Read and preprocess data
  rawDataObject <- suppressWarnings(GeomxTools::readNanoStringGeoMxSet(
    dccFiles = DCCFiles,
    pkcFiles = PKCFiles,
    phenoDataFile = SampleAnnotationFile,
    phenoDataSheet = "Template",
    phenoDataDccColName = "Sample_ID",
    protocolDataColNames = c("aoi", "roi"),
    analyte = "RNA"
  ))

  demoData <- add_demoElem(rawDataObject)

  # Optionally do QC before aggregation
  demoData <- aExprsDataQC(demoData, vFlags = "QCFlags")

  # Generate count table
  count_table <- genRawReadCountTable(demoData)

  # Unit test expectations
  expect_s3_class(count_table, "data.frame")
  expect_true(nrow(count_table) > 0)
  expect_true(ncol(count_table) > 0)
})
