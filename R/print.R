#' Print 
#' 
#' @param x \code{VarExp} object created by \code{\link[varianceExplainedPack:varianceExplained]{varianceExplainedPack::varianceExplained()}}.or \code{varExpProp} object created by \code{\link[varianceExplainedPack:proportionOf]{varianceExplainedPack::proportionOf()}}
#' @param ... arguments passed to print()
#' @export
#' @rdname print.varExp
print.VarExp <-function(x,...)  cat("Variance decomposition object. Sum of variance components deviates from var(y) by ", x$error,".",sep = "")

#' @export
#' @rdname print.varExp
print.VarExpProp <-function(x,...)   cat("Proportion of variance explained object. Sum of variance components deviates from var(y) by ", x$error/x$var.y*100,"%.",sep = "")




#' Print  summary
#' @param x object created by summary()
#' @param ... arguments passed to print()
#' @export
#' @rdname print.summary.varExp
print.summary.VarExp <-function(x, ...){
  cat("Fixed effects:\n")
  cat("   Total explained variation by fixed effects (and by correlations with random effects):\n")
  print(x$fixed, ...)
  cat("   Partial explained variation by fixed effects (excluding by correlations with random effects):\n")
  print(x$fixedPartial, ...)
  
  cat("\nExplained variation by random effects:\n")
  print(x$random)
  
  cat("\nUnexplained variation (residual):\n")
  print(x$unexplained, ...)
  cat("\nTotal variation:\n")
  print(x$total, ...)
  cat("\nMismatch:\n")
  print(x$error, ...)
}


#' @export
#' @rdname print.summary.varExp
print.summary.VarExpProp <-function(x,...){
  cat("Fixed effects:\n")
  cat("   Proportion of total explained variation by fixed effects (and by correlations with random effects):\n")
  print(x$fixed, ...)
  cat("   Proportion of partial explained variation by fixed effects (excluding by correlations with random effects):\n")
  print(x$fixedPartial, ...)
  
  cat("\nProportion of explained variation by random effects:\n")
  print(x$random, ...)
  
  cat("\nProportion of unexplained variation (residual):\n")
  print(x$unexplained, ...)
}


