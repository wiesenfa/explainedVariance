#' Summary for (proportion of) variance explained by variance components
#' 
#' @param object \code{varExp} object created by \code{\link[varianceExplainedPack:varianceExplained]{varianceExplainedPack::varianceExplained()}}.or \code{varExpProp} object created by \code{\link[varianceExplainedPack:proportionOf]{varianceExplainedPack::proportionOf()}}
#' @param ... currently not supported
#' @export
#' @rdname summaries
summary.VarExp <-function(object,...){
  fixed <- c(X=object$Rx, X=object$Rxz,Sum=object$Rx+sum(object$Rxz))
  fixedPartial <- rowSums(object$Rxpart)
  
  random <- cbind("weighted ICC"= object$Rz.1,"data-specific deviation"=object$Rz.2, "sum" = object$Rz.1+object$Rz.2,
                  "correlation with fixed"= object$Rxz,
                  "correlation with random"=object$Rz.pairs,
                  "total sum"= object$Rz.1+object$Rz.2+object$Rxz+rowSums(object$Rz.pairs, na.rm = T)
                  )
  unexplained <- object$se2
  total <- object$var.y
  error <- object$error
  return(structure(list(fixed = fixed,
                     fixedPartial = fixedPartial,
                     random = random,
                     unexplained = unexplained,
                     total = total,
                     error = error
  ), 
  class = "summary.VarExp"))
}





#' @export
#' @rdname summaries
summary.VarExpProp <-function(object,...)  {
  return(structure(summary.VarExp(object, ...), 
                   class = "summary.VarExpProp"))
}
  
