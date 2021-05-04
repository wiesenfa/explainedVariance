var.part <- function(b, Sx, Sb){
  ##############################################################################
  # Input: estimates FEs, covariance matrix of covariates, estimated covariance
  #        matrix of FEs
  # Output: Matrix including the explained variances due to FEs
  # Aim: Calculate the
  ##############################################################################
  n <- length(b)
  rxb <- matrix(NA, n, n)
  Sb = as.matrix(Sb) # assures that correct dimensions if only single covariate
  Sx = as.matrix(Sx)
  for(i in 1 : n){
    for(j in 1 : n){
      rxb[i, j] <- (b[i] * b[j] - Sb[i, j] ) * Sx[i, j]
    }
  }
  rxb
}


# parital variance including correlations with random effects
var.part_XZ <- function(b, u, Sxz){
  k <- length(b)
  p <- length(u)
  rxz <- matrix(NA, k, p)
  for(i in 1 : k){
    for(j in 1 : p){
      rxz[i, j] <- ( b[i] * u[j] ) * Sxz[i, j] 
    }
  }
  rxz
}


compute_h1 = function(Xc, Z,
                      su, se2){
  Z <- Z[names(su)]
  q.vec= sapply(Z, ncol)
  Zc <- do.call(cbind, Z)
  n = nrow(Zc)

  Hcm <- diag(1, n) - Zc %*%
    solve( diag( rep(se2/su,q.vec) ) + crossprod(Zc) ) %*% t(Zc)

  HcmX <- Hcm %*% Xc
  HcmXm <- solve( t(Xc) %*% HcmX )
  h1 <- Hcm - HcmX %*% HcmXm %*% t(HcmX)
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
  n <- nrow(X)

  # Explained variation by fixed effects
    Rxpart <- var.part(b = b.hat,
                       Sx = crossprod(X) / (n - 1),
                       Sb = S.b.hat)
    colnames(Rxpart) = rownames(Rxpart) = colnames(X)
    Rx <- sum(Rxpart)


  # Explained variation by random effects
    Rz.1 <- sapply(names(su),
                   function(id) unname(su[id]) * matrix_trace( crossprod(Z[[id]]) )  / (n - 1)
                   )
    
    Rz.2 <- sapply(names(su),
                   function(id)  {
                     sum(  (Z[[id]] %*% u.tilde[[id]]) ^ 2 ) / (n - 1) -
                       matrix_trace( var.u[[id]] %*% crossprod(Z[[id]]))  / (n - 1)
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

  # Explained variation by correlation of fixed and random effects DGP
    Rxz <- sapply(names(su),
                  function(id)  2 *  t(b.hat) %*% t(X) %*% Z[[id]] %*% u.tilde[[id]] / (n - 1)
    )


  return( list(se2=se2,
               Rxpart = Rxpart, Rx = Rx,
               Rz.1 = Rz.1, Rz.2 = Rz.2,
               Rz.pairs = Rz.pairs,
               Rxz = Rxz) )
}


