#' Build a matrix from an input csv file, with one column as rownames
#'
#' @param demoData A NanostringGeoMxSet object.
#'
#' @return A NanostringGeoMxSet object.
#' @export
#'
#' @examples
#' demoData <- readRDS(
#'     system.file(
#'         "extdata",
#'         "miniGeoMx.rds",
#'         package = "DgeaHeatmap"
#'     )
#' )
#'
#' demoData <- add_demoElem(demoData)
add_demoElem <- function(demoData) {
    Biobase::assayDataElement(demoData, elt = "exprs")
    Biobase::assayDataElement(demoData, "demoElem") <-
        NanoStringNCTools::assayDataApply(demoData,
            MARGIN = 2,
            FUN = log, base = 10, elt = "exprs"
        )
    Biobase::assayDataElement(
        demoData,
        "demoElem"
    )[seq_len(3), seq_len(2)]
    demoData
}

#' Splitting data by group column with feature,
#' pheno or protocol data to then get the mean.
#'
#' @param demoData A NanostringGeoMxSet object,
#' storing all expression, annotation, and probe information linked together.
#' @param vGroup Group of features in column.
#' @param vElt Character string of the expression matrix name.
#'
#' @return A NanostringGeoMxSet object.
#' @export
#'
#' @examples
#' demoData <- readRDS(
#'     system.file(
#'         "extdata",
#'         "miniGeoMx.rds",
#'         package = "DgeaHeatmap"
#'     )
#' )
#'
#' demoData <- add_demoElem(demoData)
#'
#' NanoStringNCTools::assayDataApply(
#'     demoData,
#'     MARGIN = 1,
#'     FUN = mean,
#'     elt = "demoElem"
#' )[seq_len(5)]
#'
#' demoData <- split_data_by_column(
#'     demoData,
#'     vGroup = "aoi",
#'     vElt = "demoElem"
#' )
split_data_by_column <- function(demoData, vGroup, vElt) {
    utils::head(NanoStringNCTools::esBy(demoData,
        GROUP = vGroup,
        FUN = function(x) {
            NanoStringNCTools::assayDataApply(x,
                MARGIN = 1,
                FUN = mean, elt = vElt
            )
        }
    ))
    demoData
}

#' Function for automatized quality control.
#'
#' @param demoData A NanostringGeoMxSet object,
#' storing all expression, annotation, and probe information linked together.
#' @param vFlags Character string to set the flags for quality control.
#'
#' @return A matrix.
#' @export
#'
#' @examples
#' demoData <- file
#' vFlags <- "QCFlags"
aExprsDataQC <- function(demoData, vFlags) {
    demoData <- GeomxTools::shiftCountsOne(demoData, useDALogic = TRUE)
    demoData <- GeomxTools::setSegmentQCFlags(demoData)
    utils::head(Biobase::protocolData(demoData)[[vFlags]])
    demoData <- GeomxTools::setBioProbeQCFlags(demoData)
    Biobase::featureData(demoData)[[vFlags]][seq_len(5), seq_len(4)]
    QCResultsIndex <- which(apply(
        Biobase::protocolData(demoData)[[vFlags]],
        1L, function(x) sum(x) == 0L
    ))
    QCPassed <- demoData[, QCResultsIndex]
    dim(QCPassed)
    QCPassed
}

#' Function to generating a raw read count table.
#'
#' @param demoData A NanostringGeoMxSet object,
#' storing all expression, annotation, and probe information linked together.
#'
#' @return A dataframe containing only the raw read counts
#' that have passed quality control.
#' @export
#'
#' @examples
#' demoData <- readRDS(
#'     system.file(
#'         "extdata",
#'         "miniGeoMx.rds",
#'         package = "DgeaHeatmap"
#'     )
#' )
#'
#' demoData <- add_demoElem(demoData)
#'
#' NanoStringNCTools::assayDataApply(
#'     demoData,
#'     MARGIN = 1,
#'     FUN = mean,
#'     elt = "demoElem"
#' )[seq_len(5)]
#'
#' demoData <- split_data_by_column(
#'     demoData,
#'     vGroup = "aoi",
#'     vElt = "demoElem"
#' )
#' df_Exp <- genRawReadCountTable(demoData)
genRawReadCountTable <- function(demoData) {
    target_demoData <- GeomxTools::aggregateCounts(demoData)
    dim(target_demoData)

    GeomxTools::featureType(demoData)
    nr <- min(5, nrow(Biobase::exprs(target_demoData)))
    nc <- min(5, ncol(Biobase::exprs(target_demoData)))

    Biobase::exprs(target_demoData)[seq_len(nr), seq_len(nc)]

    matrix_Exp <- Biobase::exprs(target_demoData)
    class(matrix_Exp)
    df_Exp <- BiocGenerics::as.data.frame(matrix_Exp)
    df_Exp
}
