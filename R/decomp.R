
compute_h1 = function(Xc, Z,
                      su, se2,
                      cholesky = TRUE){
  Z <- Z[names(su)]
  q.vec= sapply(Z, ncol)
  Zc <- do.call(cbind, Z)
  n = nrow(Zc)

  if (cholesky) SOLVE <- function(x) chol2inv(chol(x))
  else SOLVE <- function(x) solve(x)
    
  Hcm <- diag(1, n) - Zc %*%
    SOLVE( diag( rep(se2/su,q.vec) ) + crossprod(Zc) ) %*% t(Zc)

  if (ncol(Xc)>0){ # only if fixed effects present
    HcmX <- Hcm %*% Xc
    HcmXm <- SOLVE( t(Xc) %*% HcmX )
    h1 <- Hcm - HcmX %*% HcmXm %*% t(HcmX)
    
  } else h1 <- Hcm
  h1
}

matrix_trace <- function(x) sum(diag(x))


#' @importFrom utils combn
decomp <- function(X, Z,
                   se2, su,
                   b.hat, S.b.hat,
                   u.tilde,
                   var.u, h1 ){

  ##############################################################################
  # Input: Model matrices, some values from the lmer fit
  # Output: Rx (including further decomp), Rz, Rxz
  # Aim: Estimate the decomposition of quadratic form of observations
  ##############################################################################
  n <- max(nrow(X), nrow(Z)) # allows X or Z to be NULL

  # Explained variation by fixed effects
    if (ncol(X)>0){
      Rx.part <- var.part(b = b.hat,
                          Sx = crossprod(X) / (n - 1),
                          Sb = S.b.hat)
      colnames(Rx.part) = rownames(Rx.part) = colnames(X)
      Rx <- sum(Rx.part)
    } else Rx <- Rx.part <- NULL


  # Explained variation by random effects
    if (!is.null(Z)){
      ZtZ <- lapply(names(su),
                    function(id) crossprod(Z[[id]])) 
      names(ZtZ) = names(su)
      Rz.1 <- sapply(names(su),
                     function(id) unname(su[id]) * matrix_trace( ZtZ[[id]] )  / (n - 1)
                     )
      
      Rz.2 <- sapply(names(su),
                     function(id)  {
                       sum(  (Z[[id]] %*% u.tilde[[id]]) ^ 2 ) / (n - 1) -
                         matrix_trace( var.u[[id]] %*% ZtZ[[id]])  / (n - 1)
                     })

      
      
      # all pairwise
        if (length(su)>1){
          combis <- combn(names(su), 2)
  
          Rz.pairs <- matrix(NA, 
                             nrow = length(su),
                             ncol = length(su),
                             dimnames = list(names(su), names(su)))
          for (i in 1:ncol(combis)){
            id1 <- combis[1,i]
            id2 <- combis[2,i]
            
            var.u12 <- su[id1] * su[id2] * t(Z[[id1]]) %*% h1 %*% Z[[id2]] / se2
            Rz.pairs[id1, id2] <- Rz.pairs[id2, id1] <- t(u.tilde[[id1]]) %*% t(Z[[id1]]) %*% Z[[id2]] %*% u.tilde[[id2]] / (n - 1) -
              matrix_trace( Z[[id1]] %*% var.u12 %*% t(Z[[id2]]) )   / (n - 1)
          }
        } else {
          Rz.pairs = NULL
        }
    } else Rz.1 <- Rz.2 <- Rz.pairs <- NULL

  # Explained variation by correlation of fixed and random effects DGP
    if (ncol(X)>0 & !is.null(Z)){
      Rxz <- sapply(names(su),
                    function(id)  2 *  t(b.hat) %*% t(X) %*% Z[[id]] %*% u.tilde[[id]] / (n - 1)
                    )
      Rxz.part <- sapply(names(su),
                        function(id) rowSums(var.part_XZ(b=b.hat, u.tilde[[id]], Sxz= ( t(X) %*% Z[[id]] / (n - 1)   ) )) 
                        )
      if (!is.matrix(Rxz.part)){ # for cases where there is only one fixed effect
        Rxz.part <- t(as.matrix(Rxz.part))
      }
      rownames(Rxz.part) <- colnames(X)
    } else Rxz <- Rxz.part <- NULL

  return( list(se2=se2,
               Rx.part = Rx.part, Rx = Rx,
               Rz.1 = Rz.1, Rz.2 = Rz.2,
               Rz.pairs = Rz.pairs,
               Rxz = Rxz,
               Rxz.part = Rxz.part) )
}


