
reconstructFactors= function(objcect){
  aa=getME(objcect,"X")
  colIDs=sapply(unique(attr(aa,"assign")), function(x) which(attr(aa,"assign")==x))
  lapply(colIDs, function(x) colnames(aa)[x])
  
}


reduceFactors.numeric = function(fixed, object){
  factor_ids= reconstructFactors(object$model)[-1]
  fixed = sapply(factor_ids, function(fac) ifelse(length(fac)>1, 
                                                  yes=sum(fixed[fac]),
                                                  no=fixed[fac]))
  
  names(fixed)=sapply(factor_ids, function(fac) ifelse(length(fac)>1, 
                                                       yes=gsub(gsub(paste0("[",fac[2],"]"),"",fac[1]), "",fac[1]),
                                                       no=fac))
  fixed
}

reduceFactors.matrix = function(fixed, object){
  factor_ids= reconstructFactors(object$model)[-1]
  for (fac in factor_ids){
    if (length(fac)>1){
      fixed[fac,][1,] <- colSums(fixed[fac,])
      fixed[fac,][-1,] <- NA
      
      rownames(fixed)[which( rownames(fixed) %in% fac)][1] <- gsub(gsub(paste0("[",fac[2],"]"),"",fac[1]), "",fac[1])
      fixed <- fixed[apply(fixed,1, function(x) !all(is.na(x))), ]
    }
  }
  fixed
}


# agrep(pattern = attr(terms(lmm.brain), "term.labels")[1],colnames(aa)
#       )
# aa=getME(lmm.brain,"X")
# 
# colnames(aa)
# reconstructFactors(lmm.brain)
# terms(lmm.brain)
# attr(terms(lmm.brain), "term.labels")[1]
# 
# 
# ve =varianceExplained(lmm.brain)
# ve$Rx.part
# varianceExplainedPack:::summary.VarExp
# 



# factor_ids= reconstructFactors(lmm.brain)[-1]
# new.RXpart = object$Rx.part
# new.RXZpart = object$Rxz.part
# for (fac in factor_ids){
#   if (length(fac)>1){
#     #rows
#     cs=colSums(new.RXpart[fac,])
#     new.RXpart[fac,] <- NA
#     new.RXpart[fac,][1,] <- cs
#     #columns
#     cs=rowSums(new.RXpart[,fac])
#     new.RXpart[,fac] <- NA
#     new.RXpart[,fac][,1] <- cs
#     
#     colnames(new.RXpart)[which( colnames(new.RXpart) %in% fac)][1] <-
#       rownames(new.RXpart)[which( rownames(new.RXpart) %in% fac)][1] <- gsub(gsub(paste0("[",fac[2],"]"),"",fac[1]), "",fac[1])
#     new.RXpart <- new.RXpart[apply(new.RXpart,1, function(x) !all(is.na(x))), apply(new.RXpart,1, function(x) !all(is.na(x)))]
# 
#     cs=colSums(new.RXZpart[fac,])
#     new.RXZpart[fac,] <- NA
#     new.RXZpart[fac,][1,] <- cs
#     
#     rownames(new.RXZpart)[which( rownames(new.RXZpart) %in% fac)][1] <- gsub(gsub(paste0("[",fac[2],"]"),"",fac[1]), "",fac[1])
#     new.RXZpart <- new.RXZpart[apply(new.RXZpart,1, function(x) !all(is.na(x))), ]
#     }
# }
# object$Rx.part=new.RXpart
# object$Rxz.part=new.RXZpart
# 
# fixedTotal <- c(X=object$Rx, 
#                 XZ=sum(object$Rxz),
#                 total=object$Rx+sum(object$Rxz))
# 
# fixedPartial <- rowSums(object$Rx.part)
# fixedPartial.XZ <- 2 * rowSums(object$Rxz.part)
# 
# fixed <- cbind(" " = fixedPartial, 
#                "Z.*"= fixedPartial.XZ, 
#                "sum"= fixedPartial + fixedPartial.XZ
# )
# 
#
