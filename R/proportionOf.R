#' Derive proportion of variance explained by variance components
#' 
#' @param object \code{varExp} object created by \code{\link[varianceExplainedPack:varianceExplained]{varianceExplainedPack::varianceExplained()}}.
#' @export
proportionOf <- function(object){
  structure(list(
    Rx = object$Rx / object$var.y,
    Rxpart = rowSums(object$Rxpart) / object$var.y,
    Rz.1 = object$Rz.1 / object$var.y,
    Rz.2 = object$Rz.2 / object$var.y,
    Rz.pairs = object$Rz.pairs / object$var.y,
    Rxz = object$Rxz / object$var.y,
    se2 = object$se2 / object$var.y,
    error = object$var.y- object$se2 - 
      object$Rx - 
      sum(object$Rz.1) - sum(object$Rz.2) -  
      2 * sum(object$Rz.pairs) - 
      sum(object$Rxz)
  ), 
  class = "varExpProp")
}
