# auxillary functions such that returned object of varianceExplained is a vector
#' @importFrom reshape2 melt
#' @importFrom stats setNames
matrixToVector <- function(x){
  xx=  melt(x)
  xx$name=paste0(xx$Var1,"___",xx$Var2)
  setNames(xx$value,xx$name)
}

varianceExplainedToVector <- function(x, X=NULL,Z=NULL){
  if (inherits(x, c("lmerMod", "lmerModLmerTest"))) vv <- varianceExplained(x)
  else vv <- varianceExplained(x, X=X, Z=Z)
  if(!is.null(vv$Rx)){
    vv$RxzSum= sum(vv$Rxz)/2 
    vv$RxSum = vv$Rx+sum(vv$Rxz)/2
    vv$Rxpart= matrixToVector(vv$Rx.part)
    vv$RxpartRowSums = rowSums(vv$Rx.part)
    vv$RxzpartRowSums = rowSums(vv$Rxz.part)
    
    if (inherits(x, c("lmerMod", "lmerModLmerTest"))){
      vv$RxpartRowSums = reduceFactors.numeric(fixed=vv$RxpartRowSums, object=vv)
      vv$RxzpartRowSums = reduceFactors.numeric(fixed=vv$RxzpartRowSums, object=vv)
    }

    vv$RxpartTotal = vv$RxpartRowSums + vv$RxzpartRowSums
  }
  vv$Rz.sum = vv$Rz.1+vv$Rz.2
  vv$Rz1.combined= sum(vv$Rz.1)
  vv$Rz2.combined= sum(vv$Rz.2)
  vv$Rzsum.combined= sum(vv$Rz.sum)
  if (!is.null(vv$Rz.pairs)) {
    vv$Rz.pairsRowSums <- rowSums(vv$Rz.pairs, na.rm = T)
    vv$Rz.withoutCrossterm.total =  vv$Rz.sum + vv$Rz.pairsRowSums
    vv$Rz.total = vv$Rz.withoutCrossterm.total + ifPresent(vv$Rxz)/2 #vv$Rz.sum  + vv$Rz.pairsRowSums + ifPresent(vv$Rxz)/2
    vv$Rz.pairs= matrixToVector(vv$Rz.pairs)
    vv$RzpairsRowSums.combined = sum(vv$Rz.pairsRowSums)
    vv$Rztotal.withoutCrossterm.combined = sum(vv$Rz.withoutCrossterm.total)
    vv$Rztotal.combined = sum(vv$Rz.total)
  } else {
    vv$Rz.total = vv$Rz.sum + ifPresent(vv$Rxz)/2 
    vv$Rztotal.combined = sum(vv$Rz.total)
  }

  vv$model <- vv$Rx.part <- vv$Rxz.part <- NULL
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
    Rx.part = vectorToMatrix(extractMatrix(x, "Rx.part")),
    Rx=x["Rx"],
    Rz.1=x[grep("Rz.1",names(x))],
    Rz.2=x[grep("Rz.2",names(x))],
    Rz.pairs = vectorToMatrix(extractMatrix(x, "Rz.pairs")),
    Rxz=x[grep("Rxz",names(x))],
    error=x["error"]
  ), 
  class = "VarExp")
  
}



#' Derive bootstrap confidence intervals for variance decomposition (based on percentile method)
#' 
#' @param object a \code{lmerMod} or \code{lmerModLmerTest} object created by \code{\link[lme4:lmer]{lme4::lmer()}} or \code{\link[lmerTest:lmer]{lmerTest::lmer()}}, respectively, or a \code{mmer} object created by  \code{\link[sommer:mmer]{sommer::mmer()}}.
#' @param ... arguments passed to  \code{\link[lme4:bootMer]{lme4::bootMer()}}, in particular \code{nsim} for the number of simulations, the type of bootstrap and arguments for parallel computing.  In case of mmer() objects arguments  X, Z need to be provided
#' @param parallel only mmer objects. TRUE for parallelization. Initiallize parallization e.g. with \code{library(doParallel); registerDoParallel(cores = 8)} before 
#' @param progress only mmer objects. passed to plyr::llply 
#' @export
#' @rdname bootVarianceExplained
bootVarianceExplained <- function(object,...) UseMethod("bootVarianceExplained")


#' @export
#' @rdname bootVarianceExplained
bootVarianceExplained.default <- function(object, ...) stop("not implemented for this class")

#' Derive bootstrap confidence intervals for variance decomposition (based on percentile method)
#' 
#' @param object a \code{VarExp} object from varianceExplained()
#' @param ... arguments passed to  \code{\link[lme4:bootMer]{lme4::bootMer()}}, in particular \code{nsim} for the number of simulations, the type of bootstrap and arguments for parallel computing.
#' @param progress only mmer objects. passed to plyr::llply 
#' @export
#' @rdname boot
boot <- function(object,...) UseMethod("boot")

#' @export
#' @rdname boot
boot.default <- function(object, ...) stop("not implemented for this class")


#' @export
#' @rdname boot
boot.VarExp <- function(object, ...){
  if (!is.null(object$X) & !is.null(object$Z))  bootVarianceExplained(object$model, X = object$X, Z=object$Z, ...) # mmer object
  else bootVarianceExplained(object$model, ...)
}


#' @importFrom lme4 bootMer
#' @importFrom stats quantile
#' @export
#' @rdname bootVarianceExplained
bootVarianceExplained.lmerMod <- bootVarianceExplained.lmerModLmerTest <- function(object, ...){
  bootobj = bootMer(object,
                    varianceExplainedToVector,
                    ...)
  
  bootobj$model <- object
  structure(bootobj,
            class = c("VarExp.boot", "bootMer", "boot"))
}

#' @export
#' @importFrom sommer mmer
#' @importFrom plyr llply 
#' @importFrom stats as.formula rnorm
#' @rdname bootVarianceExplained
bootVarianceExplained.mmer = function(object, X, Z, nsim=1000, parallel=FALSE,
                                      progress="time", ...){
  reform= as.formula(paste("~",paste(names(Z), collapse ="+")))
  se2 <- as.numeric( object$sigma$units )
  su<-  sqrt(unlist(object$sigma[-which(names(object$sigma)=="units")],
                    recursive = F))
  tstar =t(simplify2array(llply(1:nsim,
                                function(.){
                                  object$data$.ystar=
                                    object$fitted+
                                    rnorm(nrow(object$fitted), 0, sqrt(se2))+
                                    rowSums(sapply(names(su), 
                                                   function(i) Z[[i]]%*%rnorm(ncol(Z[[i]]),0, su[i])))
                                  
                                  
                                  varianceExplainedToVector(mmer(.ystar ~ X, random= reform, 
                                                                 verbose=FALSE, 
                                                                 date.warning=FALSE, 
                                                                 data = object$data), 
                                                            X=X, Z=Z)
                                  
                                },
                                .parallel = parallel,
                                .progress = progress
  )))
  t0=varianceExplainedToVector(object, X=X, Z=Z)
  bootobj = list(t=tstar,  t0=t0)
  bootobj$model <- object
  structure(bootobj,
            class = c("VarExp.boot"))#, "boot"))
}

 