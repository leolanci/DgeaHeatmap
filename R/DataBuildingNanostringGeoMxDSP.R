#' Build a matrix from an input csv file, with one column as rownames
#'
#' @param demoData A NanostringGeoMxSet object, storing all expression, annotation, and probe information linked together.
#'
#' @return A NanostringGeoMxSet object.
#' @export
#'
#' @examples
#' demoData <- file
add_demoElem <- function(demoData){
  Biobase::assayDataElement(demoData, elt = "exprs")
  Biobase::assayDataElement(demoData, "demoElem") <-
    NanoStringNCTools::assayDataApply(demoData, MARGIN=2, FUN=log, base=10, elt="exprs")
  Biobase::assayDataElement(demoData, "demoElem")[1:3, 1:2]
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
#' demoData <- file
#' VGroup <-"aoi"
#' elt <- "demoElem"
split_data_by_column <- function(demoData, vGroup, vElt){
  utils::head(NanoStringNCTools::esBy(demoData,
            GROUP = vGroup,
            FUN = function(x) {
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
  Biobase::featureData(demoData)[[vFlags]][1:5, 1:4]
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
#' demoData <- file
genRawReadCountTable <- function(demoData){
  target_demoData <- GeomxTools::aggregateCounts(demoData)
  dim(target_demoData)

  GeomxTools::featureType(demoData)
  Biobase::exprs(target_demoData)[1:5, 1:5]

  matrix_Exp <- Biobase::exprs(target_demoData)
  class(matrix_Exp)
  df_Exp <- BiocGenerics::as.data.frame(matrix_Exp)
  return(df_Exp)
}

