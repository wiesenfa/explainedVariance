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
  vv$Rz.pairs= matrixToVector(vv$Rz.pairs)
  vv$RxpartRowSums = rowSums(vv$Rxpart)
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

# b=bootMer(lmm.class,varianceExplainedToVector,nsim = 3)
# bt=b$t[,-grep("Rxpart.",colnames(b$t), fixed=T)]
# colnames(bt) <- gsub("RxpartRowSums","Rxpart",colnames(bt), fixed=F)
# bt
