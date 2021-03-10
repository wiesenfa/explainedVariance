#' Derive proportion of variance explained by variance components
#' 
#' @param object \code{varExp} object created by \code{\link[varianceExplainedPack:varianceExplained]{varianceExplainedPack::varianceExplained()}}.
#' @export
proportionOf <- function(object){
  structure(list(
    Rx = object$Rx / object$var.y,
    Rxpart = object$Rxpart / object$var.y,
    Rz.1 = object$Rz.1 / object$var.y,
    Rz.2 = object$Rz.2 / object$var.y,
    Rz.pairs = object$Rz.pairs / object$var.y,
    Rxz = object$Rxz / object$var.y,
    se2 = object$se2 / object$var.y,
    error = object$error / object$var.y,
    var.y = object$var.y
  ), 
  class = "VarExpProp")
}
