#' Extract explained covariance matrix by random effects
#' 
#' @param object a \code{VarExp} object 
#' @param ... Currently no additional arguments
#' @return a matrix
#' @export
#' @rdname explainedCovariance
explainedCovariance <- function(object,...) UseMethod("explainedCovariance")

explainedCovariance.default <- function(object, ...) stop("not implemented for this class")

#' @rdname explainedCovariance
explainedCovariance.VarExp <- function(object, ...) {
  diag(object$Rz.pairs) <- object$Rz.1 + object$Rz.2
  object$Rz.pairs 
}
  

