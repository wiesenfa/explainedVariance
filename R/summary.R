#' @title
#' Summary for (proportion of) variance explained by variance components
#' 
#' @description
#' Provide description here
#' 
#' @details
#' Provide details here
#' 
#' @param object \code{varExp} object created by \code{\link[varianceExplainedPack:varianceExplained]{varianceExplainedPack::varianceExplained()}}.or \code{varExpProp} object created by \code{\link[varianceExplainedPack:proportionOf]{varianceExplainedPack::proportionOf()}}
#' @param ... currently not supported
#' 
#' @examples 
#' # provide some example code
#' 
#' @export
#' @rdname summaries
summary.VarExp <-function(object,...){
  if (length(object$Rx)>0){
    fixedTotal <- c(X=object$Rx, 
                    XZ=sum(object$Rxz),
                    total=object$Rx+sum(object$Rxz))
    
    fixedPartial <- rowSums(object$Rx.part)
    fixedPartial.XZ <- 2 * rowSums(object$Rxz.part)
    
    fixed <- cbind(" " = fixedPartial, 
                   "Z.*"= fixedPartial.XZ, 
                   "sum"= fixedPartial + fixedPartial.XZ
    )
    
    if (nrow(fixed)>1) fixed <- rbind(fixed,
                                      "total"= fixedTotal
    )
    if (inherits(object$model,c("lmerMod", "lmerModLmerTest"))){
      # combine explained variance of factor levels
        fixed = reduceFactors.matrix(fixed, object)
    }
  } else fixed <- NULL
  
  if (length(object$Rz.1)>0){
    random <- cbind("population"= object$Rz.1,
                    "data-specific deviation"=object$Rz.2, 
                    "data-specific" = object$Rz.1+object$Rz.2,
                    "X.*"= object$Rxz
    )
    if (length(object$Rz.pairs)>0) {
      random <- cbind(random,
                      "explainedCovariance"=rowSums(object$Rz.pairs, na.rm = TRUE),
                      "sum"= object$Rz.1+object$Rz.2 + 
                        ifPresent(object$Rxz) + 
                        rowSums(object$Rz.pairs, na.rm = T)
      )
    }
    if (nrow(random)>1) random <- rbind(random,
                                        total = colSums(random, 
                                                        na.rm = TRUE))
  } else random <- NULL
  unexplained <- object$se2
  total <- object$var.y
  error <- object$error
  return(structure(list(fixed = fixed,
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
  if (object$type != "population"){
    return(structure(summary.VarExp(object, ...), 
                     class = "summary.VarExpProp"))
    
  } else {
    if (length(object$Rx)>0){
      fixedTotal <- c(object$Rx)
      fixedPartial <- rowSums(object$Rx.part)
      fixed <- cbind(" " = fixedPartial   
      )
      
      if (nrow(fixed)>1) fixed <- rbind(fixed,
                                        "total"= fixedTotal
      )
      if (inherits(object$model,c("lmerMod", "lmerModLmerTest"))){
        # combine explained variance of factor levels
          fixed = reduceFactors.matrix(fixed, object)
      }
    } else fixed <- fixedPartial <- NULL
    
    if (length(object$Rz.1)>0){
      random <- cbind("population"= object$Rz.1 )
      if (nrow(random)>1) random <- rbind(random,
                                          total = colSums(random, 
                                                          na.rm = TRUE))
    } else random <- NULL
    unexplained <- object$se2
    total <- object$var.y
    return(structure(list(fixed = fixed,
                          fixedPartial = fixedPartial,
                          random = random,
                          unexplained = unexplained,
                          total = total), 
                     class = "summary.VarExpProp"))
  }
}



#' @export
#' @param probs quantiles
#' @rdname summaries
summary.VarExp.boot <-function(object, 
                               probs = c(.025,.975),
                               ...){
  # uses bootsrap percentile confidence intervals always. could be made more flexible using boot::boot.ci().
  bt=apply(object$t[,!grepl("Rxpart.", colnames(object$t), fixed=T)], 
           2, 
           quantile, 
           probs = probs, 
           na.rm=TRUE)
  bt0=object$t0[!grepl("Rxpart.", names(object$t0), fixed=T)]
  object = as.data.frame(t(rbind(bt0, bt)))
  colnames(object)[1]=""
  fixedPartial <-  object[grepl("RxpartRowSums.",rownames(object), fixed=T),]
  rownames(fixedPartial) <- gsub("RxpartRowSums.", "", rownames(fixedPartial), fixed=F)
  fixedPartial.XZ <-  object[grepl("RxzpartRowSums.",rownames(object), fixed=T),]
  rownames(fixedPartial.XZ) <- gsub("RxzpartRowSums.", "", rownames(fixedPartial.XZ), fixed=F)
  fixedPartialTotal <-  object[grepl("RxpartTotal.",rownames(object), fixed=T),]
  rownames(fixedPartialTotal) <- gsub("RxpartTotal.", "", rownames(fixedPartialTotal), fixed=F)
  
  fixedTotal <- cbind(" "=object["Rx",], 
                      "Z.*"=object["RxzSum",],
                      "sum"=object["RxSum",])
  
  fixed <- cbind(" " = fixedPartial, 
                 "Z.*"= fixedPartial.XZ, 
                 "sum"= fixedPartialTotal
  )
  
  
  if (nrow(fixed)>1) fixed <- rbind(fixed,
                                    "total"= fixedTotal
                                    )
  
  random <- cbind("population" = object[grepl("Rz.1", rownames(object)),],
                  "data-specific deviation" = object[grepl("Rz.2", rownames(object)),],
                  "data-specific" = object[grepl("Rz.sum", rownames(object)),]
                  )
  if (nrow(fixed)>0)  random <- cbind(random,
                                      "X.*" = object[grepl("Rxz.", rownames(object),fixed=T),]
                                      )
  
  if (nrow(object[grep("Rz.pairsRowSums", rownames(object)),])>0) {
    random <- cbind(random,
                    "explainedCovariance" = object[grep("Rz.pairsRowSums", rownames(object)),],
                    "sum"= object[grep("Rz.total", rownames(object)),]
                    )
  }
  random <- t(random)
  colnames(random) <- gsub("Rz.1.","",colnames(random), fixed=F)
  
  randomTotal = cbind("population" = object[grep("Rz1.combined", rownames(object)),],
                       "data-specific deviation" =object[grep("Rz2.combined", rownames(object)),],
                       "data-specific" = object[grep("Rzsum.combined", rownames(object)),]
                      )
  if (nrow(fixed)>0)  randomTotal <- cbind(randomTotal,
                                           "X.*" = object[grep("RxzSum", rownames(object),fixed=T),]  
                                           )
  
  if (nrow(object[grep("Rz.pairsRowSums", rownames(object)),])>0) {
    randomTotal <- cbind(randomTotal,
                         "explainedCovariance" = object[grep("RzpairsRowSums.combined", rownames(object)),],
                         "sum"= object[grep("Rztotal.combined", rownames(object)),]
                         )
  }
  
  randomTotal <- t(randomTotal)
  colnames(randomTotal) <- "total"
  
  if (ncol(random)>1) random <- cbind(random,
                                      total = randomTotal)
  
  # random <- rbind("population: " = object[grep("Rz.1", rownames(object)),],
  #                 "data-specific deviation: " = object[grep("Rz.2", rownames(object)),],
  #                 "data-specific: " = object[grep("Rz.sum", rownames(object)),],
  #                 "X." = object[grep("Rxz.", rownames(object),fixed=T),],
  #                 "explainedCovariance: " = object[grep("Rz.pairsRowSums", rownames(object)),],
  #                 "total: "= object[grep("Rz.total", rownames(object)),]
  # )
  # rownames(random) <- gsub(".Rz.1.","",rownames(random), fixed=F)
  # rownames(random) <- gsub(".Rz.2.","",rownames(random), fixed=F)
  # rownames(random) <- gsub(".Rxz.","",rownames(random), fixed=F)
  # rownames(random) <- gsub(".Rz.pairsRowSums.","",rownames(random), fixed=F)
  # rownames(random) <- gsub(".Rz.total.","",rownames(random), fixed=F)
  # rownames(random) <- gsub(".Rz.sum.","",rownames(random), fixed=F)
  # 
  unexplained <- object[grep("se2",rownames(object)),] 
  total <- object[grep("var.y",rownames(object)),] 
  error <- object[grep("error",rownames(object)),] 
  
  return(
    structure(list(fixed = fixed,
                   fixedPartial = fixedPartial,
                   random = random,
                   unexplained = unexplained,
                   total = total,
                   error = error
    ),
    class = "summary.VarExp")
  )
}


#' @export
#' @rdname summaries
summary.VarExpProp.boot <-function(object, 
                                   probs = c(.025,.975),
                                   ...)  {
  if (attr(object, "type") != "population"){
    return(structure(summary.VarExp.boot(object, ...), 
                     class = "summary.VarExpProp"))
    
  } else {
    bt=apply(object$t[,-grep("Rxpart.", colnames(object$t), fixed=T)], 
             2, 
             quantile, 
             probs = probs, 
             na.rm=TRUE)
    bt0=object$t0[-grep("Rxpart.", names(object$t0), fixed=T)]
    object = as.data.frame(t(rbind(bt0, bt)))
    colnames(object)[1]=""

    fixedTotal <- rbind(X = object["Rx",])
    fixedPartial <- object[grep("RxpartRowSums", rownames(object)),]
    rownames(fixedPartial) <- gsub("RxpartRowSums.", "", rownames(fixedPartial), fixed=F)

    fixed <- cbind(" " =fixedPartial )
    
    if (nrow(fixed)>1) fixed <- rbind(fixed,
                                      "total"= cbind(" "=fixedTotal)
                                      )
    
    random <- cbind("population" = object[grep("Rz.1", rownames(object)),])
    rownames(random) <- gsub("Rz.1.","",rownames(random), fixed=F)
    
    randomTotal =cbind("population" = object[grep("Rz1.combined", rownames(object)),])
                         
    if (nrow(random)>1) random <- rbind(random,
                                        total = randomTotal)
    
    unexplained <- object[grep("se2",rownames(object)),] 
    total <- object[grep("var.y",rownames(object)),] 
    
    return(structure(list(fixed = fixed,
                          fixedPartial = fixedPartial,
                          random = random,
                          unexplained = unexplained,
                          total = total), 
                     class = "summary.VarExpProp"))
  }
}

