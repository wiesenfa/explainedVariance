#' Summary for (proportion of) variance explained by variance components
#' 
#' @param object \code{varExp} object created by \code{\link[varianceExplainedPack:varianceExplained]{varianceExplainedPack::varianceExplained()}}.or \code{varExpProp} object created by \code{\link[varianceExplainedPack:proportionOf]{varianceExplainedPack::proportionOf()}}
#' @export
#' @rdname summaries
summary.varExp <-function(x,...)  summary(x, ...)


#' @export
#' @rdname summaries
summary.varExpProp <-function(x,...)  summary(x, ...)
