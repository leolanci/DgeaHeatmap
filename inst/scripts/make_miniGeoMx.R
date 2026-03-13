## This script creates a minimal NanoStringGeoMxSet
## for use in examples.
## It should be run manually and is NOT executed during package check.

library(GeomxTools)

# Load example data from GeoMxWorkflows
datadir <- system.file(
  "extdata",
  "WTA_NGS_Example",
  package = "GeoMxWorkflows"
)

DCCFiles <- dir(
  file.path(datadir, "dccs"),
  pattern = ".dcc$",
  full.names = TRUE,
  recursive = TRUE
)

PKCFiles <- unzip(
  zipfile = dir(
    file.path(datadir, "pkcs"),
    pattern = ".zip$",
    full.names = TRUE,
    recursive = TRUE
  ),
  exdir = tempdir()
)

SampleAnnotationFile <- dir(
  file.path(datadir, "annotation"),
  pattern = ".xlsx$",
  full.names = TRUE,
  recursive = TRUE
)

demoData <- suppressWarnings(
  GeomxTools::readNanoStringGeoMxSet(
    dccFiles = DCCFiles,
    pkcFiles = PKCFiles,
    phenoDataFile = SampleAnnotationFile,
    phenoDataSheet = "Template",
    phenoDataDccColName = "Sample_ID",
    protocolDataColNames = c("aoi", "roi"),
    analyte = "RNA"
  )
)

# Subset to tiny object (adjust if needed)
miniGeoMx <- demoData[1:5, 1:3]

# Save into inst/extdata
saveRDS(
  miniGeoMx,
  file = file.path("inst", "extdata", "miniGeoMx.rds")
)

message("miniGeoMx.rds successfully created.")
