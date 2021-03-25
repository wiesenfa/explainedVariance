# auxillary functions such that returned object of varianceExplained is a vector
#' @importFrom reshape2 melt
#' @importFrom stats setNames
matrixToVector <- function(x){
  xx=  melt(x)
  xx$name=paste0(xx$Var1,"___",xx$Var2)
  setNames(xx$value,xx$name)
}

varianceExplainedToVector <- function(x){
  vv <- varianceExplained(x)
  vv$Rx.Sum = vv$Rx+sum(vv$Rxz)
  vv$RxpartRowSums = rowSums(vv$Rxpart)
  vv$Rz.pairsRowSums <- rowSums(vv$Rz.pairs, na.rm = T)
  vv$Rz.sum = vv$Rz.1+vv$Rz.2
  vv$Rz.total =  vv$Rz.sum +vv$Rxz+vv$Rz.pairsRowSums
  vv$Rz.pairs= matrixToVector(vv$Rz.pairs)
  vv$Rxpart= matrixToVector(vv$Rxpart)
  unlist(vv)
}



# auxillary functions to transform vector valued varianceExplained object to usual VarExp object

#' @importFrom tidyr separate
#' @importFrom reshape2 acast
#' @importFrom dplyr %>%
vectorToMatrix <- function(x){
  data.frame(names=names(x),XX=x)%>%
    separate(names, into=c("Var","row"),sep="___")%>% 
    acast(Var~row,value.var="XX")
}


extractMatrix <- function(x, id){
  xx = x[grep(id,names(x))]
  names(xx)= gsub(paste0(id,"."),"",names(xx))
  xx
}


vectorToVarExp <- function(x){
  structure(list(
    var.y=x["var.y"],
    se2=x["se2"],
    Rxpart = vectorToMatrix(extractMatrix(x, "Rxpart")),
    Rx=x["Rx"],
    Rz.1=x[grep("Rz.1",names(x))],
    Rz.2=x[grep("Rz.2",names(x))],
    Rz.pairs = vectorToMatrix(extractMatrix(x, "Rz.pairs")),
    Rxz=x[grep("Rxz",names(x))],
    error=x["error"]
  ), 
  class = "VarExp")
  
}



#' Derive bootstrap confidence intervals for variance decomposition
#' 
#' @param object a \code{lmerMod} or \code{lmerModLmerTest} object created by \code{\link[lme4:lmer]{lme4::lmer()}} or \code{\link[lmerTest:lmer]{lmerTest::lmer()}}, respectively, or a \code{mmer} object created by  \code{\link[sommer:mmer]{sommer::mmer()}}.
#' @param ... arguments passed to  \code{\link[lme4:bootMer]{lme4::bootMer()}}, in particular \code{nsim} for the number of simulations, thy type of bootstrap and arguments for parallel computing
#' @export
#' @rdname bootstrap
bootstrap <- function(object,...) UseMethod("bootstrap")

#' @rdname bootstrap
bootstrap.default <- function(object, ...) stop("not implemented for this class")


#' @importFrom lme4 bootMer
#' @importFrom stats quantile
#' @export
#' @rdname bootstrap
bootstrap.lmerMod <- bootstrap.lmerModLmerTest <- function(object, ...){
  bootobj = bootMer(object,
                    varianceExplainedToVector,
                    ...)
  
  bt=bootobj$t[,-grep("Rxpart.", colnames(bootobj$t), fixed=T)]
  bt=apply(bt, 2, quantile, probs=c(.025,.975), na.rm=TRUE)
  bt0=bootobj$t0[-grep("Rxpart.", names(bootobj$t0), fixed=T)]
  bootobj=as.data.frame(t(rbind(bt0, bt)))
  structure(bootobj,
            class = c("VarExp.boot", "data.frame" ))
}


 