#' Derive proportion of variance explained by variance components
#' 
#' @param object \code{varExp} object created by \code{\link[varianceExplainedPack:varianceExplained]{varianceExplainedPack::varianceExplained()}}.
#' @param type specify "dataset-specific" or "population" variation in denominator 
#' @export
#' @rdname proportionOf
proportionOf <- function(object,...) UseMethod("proportionOf")

#' @rdname proportionOf
#' @export
proportionOf.default <- function(object, ...) stop("not implemented for this class")



#' @rdname proportionOf
#' @export
proportionOf.VarExp <- function(object, 
                                type = "dataset-specific",
                                ...){
  
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


#' @rdname proportionOf
#' @export
proportionOf.VarExp.boot <- function(object, 
                                     type = "dataset-specific",
                                     ...){
  
  if (type != "population") {
    object$t0 <-   object$t0 / object$t0["var.y"]
    object$t <- t(apply(object$t,
                        1,
                        function(samp.i) samp.i / samp.i["var.y"]
                        ))
  }  else {
    object$t0 <-   object$t0 / object$t0["Rx"] + sum(object$t0[grep("Rz.1", names(object$t0))]) + object$t0["se2"]
    object$t <- t(apply(object$t,
                        1,
                        function(samp.i) samp.i / (samp.i["Rx"] + sum(samp.i[grep("Rz.1", names(samp.i))]) + samp.i["se2"])
                        ))
  }
  
  attr(object, "type") <- type
  structure(
    object,
    class = c("VarExpProp.boot", "bootMer", "boot" ) 
    )
}


