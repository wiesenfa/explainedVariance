
# ------------------------------------------------------------------------------
# Variance Partitioning of the Covariates associated with FE
var.part <- function(b, Sx, Sb){
  ##############################################################################
  # Input: estimates FEs, covariance matrix of covariates, estimated covariance
  #        matrix of FEs
  # Output: Matrix including the explained variances due to FEs
  # Aim: Calculate individual (total) contributions of FE to observed variance; 
  #      direct contributions are on the diagonal, and joint contributions on
  #      off diagonal; see Section 2.3 (components of EV in the LM)
  # Use: Calculate EV of fixed effects and individual total contributions for
  #      tables amending the standard model output
  ##############################################################################
  Sb = as.matrix(Sb) # assures that correct dimensions if only single covariate
  Sx = as.matrix(Sx)
  k <- length(b) # number of FE
  rxb <- matrix(NA, k, k)
  for(i in 1 : k){
    for(j in 1 : k){
      rxb[i, j] <- (b[i] * b[j] - Sb[i, j] ) * Sx[i, j] # see formula (3.12)
    }
  }
  return(rxb) # rowsums or colsums give total contribution (\dot{S}^2_{X_j})
  # sum of all matrix elements gives \dot{S}^2_X
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Variance Partitioning of the Data-set EV between fixed and random covariates
var.part_XZ <- function(b, u, Sxz){
  ##############################################################################
  # Input: estimates FEs, covariance matrix between covariates associated with
  #        fixed and random effects; estimated realizations of REs
  # Output: Matrix including the explained covariances due to FE & RE
  # Aim: Calculate Data-set specific explained covariance between covariates
  #      associated with fixed and random effects
  # Use: Amend the EV by covariates associated with FE with data-set
  #      specific contributions and amend the EV by covariates associated with
  #      RE by data-set specific contributions
  ##############################################################################
  k <- length(b)
  p <- length(u)
  rxz <- matrix(NA, k, p)
  for(i in 1 : k){
    for(j in 1 : p){
      rxz[i, j] <- ( b[i] * u[j] ) * Sxz[i, j] 
    }
  }
  return(rxz)
}
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Variance Partitioning of population EV associated with RE covariates
var.part.Z1 <- function(su, Sz){
  ##############################################################################
  # Input: estimated variance components (REML), covariance matrix of covariates
  #        associated with random effects
  # Output: vector with population contribution of each covariate associated
  #         with the random effect (vector)
  # Aim: Calculate pop contributions of covariates of RE to observed variance; 
  #      see section 3.4 and in particular the paragraph above equation (3.14)
  # Use: Calculate population EV of individual RE covariates eg for manhattan
  #      plots; mainly depends on the variability of each covariate
  ##############################################################################
  d <- diag(Sz)
  sz1 <- su * d
  return(sz1)
}
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Variance Partitioning of data-set specific EV associated with RE covariates
var.part.Z2 <- function(u, Sz, Su){
  ##############################################################################
  # Input: estimated realizations of REs; estimated covariance matrix of 
  #        covariates associated with RE; estimated covariance matrix of eBLUPs
  # Output: vector with data-set specific contribution of each covariate 
  #         associated with the random effect (vector)
  # Aim: Calculate data-set specific contribution of covariates of RE to 
  #       observed variance; see first part of equation (3.14)
  # Use: Calculate data-set specific EV of individual RE covariates eg for 
  #     manhattan plots
  # Warning: This is primarily valid in the case of ONE random effect (vector); 
  #          In the case of several RE (several hierarchies etc) the full 
  #          equation (3.14) has to be implemented in order to also account for
  #          the covariance BETWEEN the random effects 
  ##############################################################################
  M <- diag( Sz %*% ( tcrossprod(u) - Su) ) 
  return(M)
}
# -----------------------------------






