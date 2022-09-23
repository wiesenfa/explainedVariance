#' @title
#' Add explained variance to model summary
#' 
#' @description
#' Provide description here
#' 
#' @details
#' Results provide proportion of variance explained in percentages
#' 
#' @param object a \code{summary.VarExpProp} object after boostrapping
#' 
#' @param ... Currently no additional arguments
#' 
#' @examples 
#' # provide some example code
#' 
#' @export
#' @importFrom dplyr rename select mutate across everything bind_cols bind_rows
#' @importFrom tibble rownames_to_column
#' @rdname expandResults
expandResults <- function(object, digits, ...) UseMethod("expandResults")

#' @rdname varianceExplained
expandResults.default <- function(object, digits, ...) stop("not implemented for this class. Only applicable to VarExp.boot (i.e. after bootstrapping)")

#' @export
#' @rdname expandResults
expandResults.summary.VarExpProp = function(object,
                                            digits=1, ...){
  lmeObject <- object$model
  toPercentage = function(x) format(round(x*100, digits = digits), nsmall=digits)        
  
  fixed=merge(as.data.frame(summary(lmeObject)$coefficients),
              as.data.frame(toPercentage(object$fixed[,1:3,drop=F])),
              by="row.names",all=T)%>% 
#    select(-"df",-"t value")%>%
    rename(VaExp=` .`) 
  
  random=as.data.frame(object$random[grep("combined.",rownames(object$random),fixed=T),-3,drop=F])%>%
    mutate(across(everything(), toPercentage))
  random= as.data.frame(t(random))%>% 
    rename(VaExp=`combined.`,
           "2.5%"=`combined.2.5%`,
           "97.5%"=`combined.97.5%`) %>%
    rownames_to_column(var = "Row.names")
  
  R_X.Z = 2* t(object$random[grep("X.*",rownames(object$random),fixed=T),"total",drop=F] )
  colnames(R_X.Z)=gsub("X.*.","",colnames(R_X.Z), fixed=T)
  colnames(R_X.Z)[1]="VaExp"
  rownames(R_X.Z)= "$S_{XxZ}$"
  
  random.estimates = unlist(lapply(VarCorr(lmeObject), function(x) attr(x, "stddev")^2),       recursive = F) / var(getME(lmeObject, "y"))
  random2 = bind_cols(random,
                      Estimate = random.estimates)
  random.tab = bind_rows(random2,
                         as.data.frame(  R_X.Z)%>%
                           mutate(across(everything(), toPercentage))%>%
                           rownames_to_column(var = "Row.names"),
                         as.data.frame(  object$unexplained)%>%
                           rename(VaExp="")%>%
                           mutate(across(everything(), toPercentage))%>%
                           mutate("Estimate" = sigma(lmeObject)^2 / var(getME(lmeObject, "y"))) %>%
                           rownames_to_column(var = "Row.names")
  )
  bind_rows(fixed%>% 
              rename(
                "2.5%"=` .2.5%`,
                "97.5%"=` .97.5%`),
            random.tab
  )
}
