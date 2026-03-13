## This script generates a minimal edgeR result object
## for use in examples.
## It is run manually and NOT during R CMD check.

library(DgeaHeatmap)

# Load example data from your package
rawCounts <- read.csv(
  system.file(
    "extdata",
    "RawDataExamplePackageNanostring.csv",
    package = "DgeaHeatmap"
  )
)

rawCounts <- build_matrix(rawCounts, 1)

metadata <- read.csv(
  system.file(
    "extdata",
    "MetaDataPackageNanostring.csv",
    package = "DgeaHeatmap"
  )
)

grouping_columns <- c("segment", "region", "class", "slide_name")

comparisons <- list(
  Comp1 = c(
    "Geometric_Segment_glomerulus_DKD_disease3",
    "PanCK_tubule_DKD_disease4"
  )
)

results_edgeR <- DGEAedgeR(
  rawCounts,
  metadata,
  grouping_columns,
  comparisons,
  prefix = "DEA"
)

# Save into inst/extdata
saveRDS(
  results_edgeR,
  file = file.path("inst", "extdata", "mini_edgeR_results.rds")
)

message("mini_edgeR_results.rds successfully created.")
