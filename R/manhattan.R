
#' @export
getProfile <- function(object){
  X <-getME(object,"X")[,-1, drop = FALSE]
  #SCALE X?????
  X <- scale(X, center = TRUE, scale = FALSE)
  Z <-lapply(getME(object,"Ztlist"), Matrix::t)
  Z <- lapply(Z, function(x) scale( x, center = TRUE, scale = FALSE))
  
  b.hat <- fixef(object)[-1]
  u.tilde <- unlist(ranef(object,condVar = FALSE), recursive = F)
  
  S.b.hat <- vcov(object, full = FALSE)[-1, -1]
  su <- unlist(lapply(VarCorr(object), function(x) attr(x, "stddev")^2),
               recursive = F)
  
  se2 <-sigma(object)^2
  h1 <- explainedVariance:::compute_h1(Xc = X, Z = Z, su, se2)
  var.u <- sapply(names(su),
                  function(id) su[id]^2 * t(Z[[id]]) %*% h1 %*% Z[[id]] / se2,
                  simplify = FALSE
  )
  
  var.x <- var.part(b = b.hat, 
                    Sx = var(X), 
                    Sb = S.b.hat)
  
  var.xz <- var.z.1 <- var.z.2 <- list()
  for (id in names(su)){
    Sz = var(Z[[id]])
    
    var.xz[[id]] <- var.part_XZ(b = b.hat,
                                u = u.tilde[[id]], 
                                Sxz = cov(X, Z[[id]]))
    var.z.1[[id]] <- var.part.Z1(su = su[id], 
                                 Sz = Sz)
    
    var.z.2[[id]] <- var.part.Z2(u = u.tilde[[id]], 
                                 Sz = Sz,
                                 Su = var.u[[id]] )
  }
  
  # proportion
  var.y <- var(getME(object, "y"))
  var.x <- lapply(var.x, function(x)  x/ var.y)
  var.xz <- lapply(var.xz, function(x)  x/ var.y)
  var.z.2 <- lapply(var.z.2, function(x)  x/ var.y)
  var.z.1 <- lapply(var.z.1, function(x)  x/ var.y)
  
  return(list(
    var.y = var.y,
    var.x = var.x,
    var.xz = var.xz,
    var.z.1 = var.z.1,
    var.z.2 = var.z.2))
}



#' @export
plotProfile <- function(object){
  library(latex2exp)
  for (id in names(object$var.z.1)){
    par(mfrow = c(2, 2), mai = c(1, 1, 1, 1), mar = c(4, 3.5, 2, 1), mgp = c(1.5, 0.5, 0) )
    plot( (object$var.z.1[[id]] + object$var.z.2[[id]] + sum(object$var.xz[[id]])) *100, xlab = id, 
          main = "Parameter-wise Dispersion Relevance Profile", 
          ylab = TeX( r'( $\dot{R}_{z_{\cdot}}^2$ in \%)' )#, 
          #  ylim = c(-0.013, 0.027) 
    )
    abline( h = 0, col = "grey", lwd = 3, lty = 2 )
    plot(object$var.z.1[[id]]*100, xlab = id, 
         main = "Parameter-wise Population Explained Variation", 
         ylab = TeX( r'( ${\dot{R}_{z_{\cdot}}^2}{}^p$ in \%)' )#, 
         #  ylim = c(0, 0.005)  
    )
    plot( (object$var.z.2[[id]] + colSums(object$var.xz[[id]]) ) *100, xlab = id, 
          main = "Parameter-wise Data-set Specific Explained Variation"#,
          #   ylim = c(-0.013, 0.027)
    )
    abline( h = 0, col = "grey", lwd = 3, lty = 2 )
  }
}

