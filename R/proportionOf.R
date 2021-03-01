proportionOf <- function(object){
  object$Rx / object$var.y
  rowSums(object$Rxpart) / object$var.y
  object$Rz.1 / object$var.y
  object$Rz.2 / object$var.y
  object$Rz.pairs / object$var.y
  object$Rxz / object$var.y
  object$se2 / object$var.y
  object$var.y- object$se2 - object$Rx - sum(object$Rz.1) - sum(object$Rz.2) -
     2 * sum(object$Rz.pairs) - sum(object$Rxz)

}
