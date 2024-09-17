#' Build a matrix from an input csv file, with one column as rownames
#'
#' @param demoData A demoData object, storing all expression, annotation, and probe information linked together.
#'
#' @return A matrix.
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

