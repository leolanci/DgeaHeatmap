#' Build a matrix from an input csv file, with one column as rownames
#'
#' @param demoData A NanostringGeoMxSet object, storing all expression, annotation, and probe information linked together.
#'
#' @return A NanostringGeoMxSet object.
#' @export
#'
#' @examples
#' datadir <- system.file("extdata", "WTA_NGS_Example", package="GeoMxWorkflows")
#' DCCFiles <- dir(file.path(datadir, "dccs"), pattern = ".dcc$", full.names = TRUE, recursive = TRUE)
#' PKCFiles <- unzip(zipfile = dir(file.path(datadir, "pkcs"), pattern = ".zip$", full.names = TRUE, recursive = TRUE))
#' SampleAnnotationFile <- dir(file.path(datadir, "annotation"), pattern = ".xlsx$", full.names = TRUE, recursive = TRUE)
#' demoData <- suppressWarnings(GeomxTools::readNanoStringGeoMxSet(dccFiles = DCCFiles, pkcFiles = PKCFiles, phenoDataFile = SampleAnnotationFile, phenoDataSheet = "Template", phenoDataDccColName = "Sample_ID", protocolDataColNames = c("aoi","roi"), configFile = NULL, analyte = "RNA", phenoDataColPrefix = "", experimentDataColNames = NULL))
#' demoData <- add_demoElem(demoData)
add_demoElem <- function(demoData){
  Biobase::assayDataElement(demoData, elt = "exprs")
  Biobase::assayDataElement(demoData, "demoElem") <-
    NanoStringNCTools::assayDataApply(demoData, MARGIN=2, FUN=log, base=10, elt="exprs")
  Biobase::assayDataElement(demoData, "demoElem")[seq_len(3), seq_len(2)]
  return(demoData)
}

#' Splitting data by group column with feature, pheno or protocol data to then get the mean.
#'
#' @param demoData A NanostringGeoMxSet object, storing all expression, annotation, and probe information linked together.
#' @param vGroup Group of features in column.
#' @param vElt Character string of the expression matrix name.
#'
#' @return A NanostringGeoMxSet object.
#' @export
#'
#' @examples
#' datadir <- system.file("extdata", "WTA_NGS_Example", package="GeoMxWorkflows")
#' DCCFiles <- dir(file.path(datadir, "dccs"), pattern = ".dcc$", full.names = TRUE, recursive = TRUE)
#' PKCFiles <- unzip(zipfile = dir(file.path(datadir, "pkcs"), pattern = ".zip$", full.names = TRUE, recursive = TRUE))
#' SampleAnnotationFile <- dir(file.path(datadir, "annotation"), pattern = ".xlsx$", full.names = TRUE, recursive = TRUE)
#' rawDataObject <- suppressWarnings(GeomxTools::readNanoStringGeoMxSet(dccFiles = DCCFiles, pkcFiles = PKCFiles, phenoDataFile = SampleAnnotationFile, phenoDataSheet = "Template", phenoDataDccColName = "Sample_ID", protocolDataColNames = c("aoi","roi"), configFile = NULL, analyte = "RNA", phenoDataColPrefix = "", experimentDataColNames = NULL))
#' demoData <- add_demoElem(rawDataObject)
#' NanoStringNCTools::assayDataApply(demoData, MARGIN=1, FUN=mean, elt="demoElem")[seq_len(5)]
#' VGroup <-"aoi"
#' vElt <- "demoElem"
#' demoData <- split_data_by_column(demoData, VGroup, vElt)
split_data_by_column <- function(demoData, vGroup, vElt){
  utils::head(NanoStringNCTools::esBy(demoData,
            GROUP = vGroup,
            FUN <- function(x) {
              NanoStringNCTools::assayDataApply(x, MARGIN = 1, FUN=mean, elt=vElt)
            }))
  return(demoData)
}

#' Function for automatized quality control.
#'
#' @param demoData A NanostringGeoMxSet object, storing all expression, annotation, and probe information linked together.
#' @param vFlags Character string to set the flags for quality control.
#'
#' @return A matrix.
#' @export
#'
#' @examples
#' demoData <- file
#' vFlags <-"QCFlags"

aExprsDataQC <- function(demoData, vFlags){
  demoData <- GeomxTools::shiftCountsOne(demoData, useDALogic = TRUE)
  demoData <- GeomxTools::setSegmentQCFlags(demoData)
  utils::head(Biobase::protocolData(demoData)[[vFlags]])
  demoData <- GeomxTools::setBioProbeQCFlags(demoData)
  Biobase::featureData(demoData)[[vFlags]][seq_len(5), seq_len(4)]
  QCResultsIndex <- which(apply(Biobase::protocolData(demoData)[[vFlags]],
                                1L , function(x) sum(x) == 0L))
  QCPassed <- demoData[, QCResultsIndex]
  dim(QCPassed)
  return(QCPassed)
}

#' Function to generating a raw read count table.
#'
#' @param demoData A NanostringGeoMxSet object, storing all expression, annotation, and probe information linked together.
#'
#' @return A dataframe containing only the raw read counts that have passed quality control.
#' @export
#'
#' @examples
#' datadir <- system.file("extdata", "WTA_NGS_Example", package="GeoMxWorkflows")
#' DCCFiles <- dir(file.path(datadir, "dccs"), pattern = ".dcc$", full.names = TRUE, recursive = TRUE)
#' PKCFiles <- unzip(zipfile = dir(file.path(datadir, "pkcs"), pattern = ".zip$", full.names = TRUE, recursive = TRUE))
#' SampleAnnotationFile <- dir(file.path(datadir, "annotation"), pattern = ".xlsx$", full.names = TRUE, recursive = TRUE)
#' rawDataObject <- suppressWarnings(GeomxTools::readNanoStringGeoMxSet(dccFiles = DCCFiles, pkcFiles = PKCFiles, phenoDataFile = SampleAnnotationFile, phenoDataSheet = "Template", phenoDataDccColName = "Sample_ID", protocolDataColNames = c("aoi","roi"), configFile = NULL, analyte = "RNA", phenoDataColPrefix = "", experimentDataColNames = NULL))
#' demoData <- add_demoElem(rawDataObject)
#' NanoStringNCTools::assayDataApply(demoData, MARGIN=1, FUN=mean, elt="demoElem")[seq_len(5)]
#' VGroup <-"aoi"
#' vElt <- "demoElem"
#' demoData <- split_data_by_column(demoData, VGroup, vElt)
#' QCPassed <- aExprsDataQC(demoData, "QCFlags")
#' df_Exp <- genRawReadCountTable(demoData)

#' demoData <- file
genRawReadCountTable <- function(demoData){
  target_demoData <- GeomxTools::aggregateCounts(demoData)
  dim(target_demoData)

  GeomxTools::featureType(demoData)
  Biobase::exprs(target_demoData)[seq_len(5), seq_len(5)]

  matrix_Exp <- Biobase::exprs(target_demoData)
  class(matrix_Exp)
  df_Exp <- BiocGenerics::as.data.frame(matrix_Exp)
  return(df_Exp)
}

