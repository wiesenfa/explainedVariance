#' Coefficient of Determination
#' 
#' @param x a \code{lmerMod} or \code{lmerModLmerTest} object created by \code{\link[lme4:lmer]{lme4::lmer()}} or \code{\link[lmerTest:lmer]{lmerTest::lmer()}}, respectively, or a \code{mmer} object created by  \code{\link[sommer:mmer]{sommer::mmer()}}.
#' @export
#' @rdname R2
R2 <- function(x) UseMethod("R2")

#' @rdname R2
R2.default <- function(x) stop("not implemented for this class")

#' @rdname R2
R2.lmerMod<-  R2.lmerModLmerTest <- function(x){
  #x=varianceExplained(x)
  #sd1=summary(x)
  #a=(sd1$fixed[nrow(sd1$fixed),ncol(sd1$fixed)]+sd1$random[nrow(sd1$random),ncol(sd1$random)]-sum(d1$Rxz))
  #a/(a+x$se2)
  #1-(x$se2/x$var.y)
  var.y <- var(getME(x, "y"))
  se2 <-sigma(x)^2
  1-(se2/var.y)
}
  