#' Summary for (proportion of) variance explained by variance components
#' 
#' @param object \code{varExp} object created by \code{\link[varianceExplainedPack:varianceExplained]{varianceExplainedPack::varianceExplained()}}.or \code{varExpProp} object created by \code{\link[varianceExplainedPack:proportionOf]{varianceExplainedPack::proportionOf()}}
#' @param ... currently not supported
#' @export
#' @rdname summaries
summary.VarExp <-function(object,...){
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
  
  random <- cbind("population"= object$Rz.1,
                  "data-specific deviation"=object$Rz.2, 
                  "data-specific" = object$Rz.1+object$Rz.2,
                  "X.*"= object$Rxz,
                  "explainedCovariance"=rowSums(object$Rz.pairs, na.rm = TRUE),
                  "sum"= object$Rz.1+object$Rz.2+object$Rxz+rowSums(object$Rz.pairs, na.rm = T)
                  )
  if (nrow(random)>1) random <- rbind(random,
                                      total = colSums(random, 
                                                      na.rm = TRUE))
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
    fixedTotal <- c(object$Rx)
    fixedPartial <- rowSums(object$Rx.part)
    fixed <- cbind(" " = fixedPartial   
                   )
    
    if (nrow(fixed)>1) fixed <- rbind(fixed,
                                      "total"= fixedTotal
    )
    
    random <- cbind("population"= object$Rz.1 )
    if (nrow(random)>1) random <- rbind(random,
                                        total = colSums(random, 
                                                        na.rm = TRUE))
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
  bt=apply(object$t[,-grep("Rxpart.", colnames(object$t), fixed=T)], 
           2, 
           quantile, 
           probs = probs, 
           na.rm=TRUE)
  bt0=object$t0[-grep("Rxpart.", names(object$t0), fixed=T)]
  object = as.data.frame(t(rbind(bt0, bt)))
  colnames(object)[1]=""
  fixedPartial <-  object[grep("RxpartRowSums.",rownames(object), fixed=T),]
  rownames(fixedPartial) <- gsub("RxpartRowSums.", "", rownames(fixedPartial), fixed=F)
  fixedPartial.XZ <-  object[grep("RxzpartRowSums.",rownames(object), fixed=T),]
  rownames(fixedPartial.XZ) <- gsub("RxzpartRowSums.", "", rownames(fixedPartial.XZ), fixed=F)
  fixedPartialTotal <-  object[grep("RxpartTotal.",rownames(object), fixed=T),]
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
  

  
  random <- t(cbind("population" = object[grep("Rz.1", rownames(object)),],
                  "data-specific deviation" = object[grep("Rz.2", rownames(object)),],
                  "data-specific" = object[grep("Rz.sum", rownames(object)),],
                  "X.*" = object[grep("Rxz.", rownames(object),fixed=T),],
                  "explainedCovariance" = object[grep("Rz.pairsRowSums", rownames(object)),],
                  "sum"= object[grep("Rz.total", rownames(object)),]
  ))
  colnames(random) <- gsub("Rz.1.","",colnames(random), fixed=F)

  randomTotal =t(cbind("population" = object[grep("Rz1.combined", rownames(object)),],
          "data-specific deviation" =object[grep("Rz2.combined", rownames(object)),],
          "data-specific" = object[grep("Rzsum.combined", rownames(object)),],
          "X.*" = object[grep("RxzSum", rownames(object),fixed=T),],
          "explainedCovariance" = object[grep("RzpairsRowSums.combined", rownames(object)),],
          "sum"= object[grep("Rztotal.combined", rownames(object)),]
  ))
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
    bt=apply(object$t[,-grep("Rx.part.", colnames(object$t), fixed=T)], 
             2, 
             quantile, 
             probs = probs, 
             na.rm=TRUE)
    bt0=object$t0[-grep("Rx.part.", names(object$t0), fixed=T)]
    object = as.data.frame(t(rbind(bt0, bt)))
    
    fixed <- rbind(X = object["Rx",])
    fixedPartial <- object[grep("Rx.partRowSums", rownames(object)),]
    rownames(fixedPartial) <- gsub("Rx.partRowSums.", "", rownames(fixedPartial), fixed=F)
    
    random <- rbind("population: " = object[grep("Rz.1", rownames(object)),])
    rownames(random) <- gsub(".Rz.1.","",rownames(random), fixed=F)
    # if (nrow(random)>1) random <- rbind(random,
    #                                     total = colSums(random, 
    #                                                     na.rm = TRUE))
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

