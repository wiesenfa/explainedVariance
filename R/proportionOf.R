#' Derive proportion of variance explained by variance components
#' 
#' @param object \code{varExp} object created by \code{\link[varianceExplainedPack:varianceExplained]{varianceExplainedPack::varianceExplained()}}.
#' @export
proportionOf <- function(object, type = "dataset-specific"){
  
  if (type != "population") {
    denom <- object$var.y
    res <- list(
      Rz.2 = object$Rz.2 / denom,
      Rz.pairs = object$Rz.pairs / denom,
      Rxz = object$Rxz / denom,
      error = object$error / denom
    )
  }  else {
    denom <- object$Rx + sum(object$Rz.1) + object$se2
    res <- NULL
  }
  structure(c(res,
               list(
                 Rx = object$Rx / denom,
                 Rxpart = object$Rxpart / denom,
                 Rz.1 = object$Rz.1 / denom,
                 se2 = object$se2 / denom,
                 var.y = denom,
                 type = type
               )), 
            class = "VarExpProp")
}
