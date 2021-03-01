#' Derive variance decomposition
#' 
#' @param object a \code{lmerMod} or \code{lmerModLmerTest} object created by \code{\link[lme4:lmer]{lme4::lmer()}} or \code{\link[lmerTest:lmer]{lmerTest::lmer()}}, respectively, or a \code{mmer} object created by  \code{\link[sommer:mmer]{sommer::mmer()}}.
#' @param ... Currently no additional arguments
#' @return a varExp object
#' @export
#' @rdname varianceExplained
varianceExplained <- function(object,...) UseMethod("varianceExplained")

#' @rdname varianceExplained
varianceExplained.default <- function(object, ...) stop("not implemented for this class")


#' @importFrom lme4 lFormula ranef fixef VarCorr
#' @importFrom Matrix t
#' @importFrom stats sigma var vcov
#' @export
#' @rdname varianceExplained
varianceExplained.lmerMod <- varianceExplained.lmerModLmerTest <- function(object, ...){

  # get matrices
    mc=object@call
    mc[[1]] <- quote(lme4::lFormula)
    lmod <- eval(mc, parent.frame(1L))
    X <-lmod$X[,-1]
    X <- scale(X, center = TRUE, scale = FALSE)
    Z0 <-lapply(lmod$reTrms$Ztlist,Matrix::t)
    # check whether both random intercept and random slope for same grouping variable
      multipleREs <- lapply(object@cnms, function(x) length(x))

    # lme4 combines all random effects of one grouping variable into a single Z, following code separates this
      Z=list()
      j=1
      for (i in 1:length(Z0)){
        if (multipleREs[[i]]==1) {
          Z[[j]] <- Z0[[i]]
          names(Z)[j] <- paste(names(object@cnms)[i],object@cnms[[i]][1], sep=".")
          j=j+1
        }  else {
          for (k in 1: multipleREs[[i]]){
            Z[[j]] <- Z0[[i]][, seq(k,
                                              by=multipleREs[[i]],
                                              length.out=ncol(Z0[[i]])/multipleREs[[i]]) ]
            names(Z)[j] <- paste(names(object@cnms)[i],object@cnms[[i]][k], sep=".")
            j = j +1

          }
        }
      }
      Z <- lapply(Z, function(x) scale( x, center = TRUE, scale = FALSE))
      rm(Z0)

  # get variance components
    se2 <-sigma(object)^2
    su<-unlist(lapply(VarCorr(object), function(x) attr(x, "stddev")^2),
                   recursive = F)

  # get estimates
    b.hat <- fixef(object)[-1]
    S.b.hat <- vcov(object, full = FALSE)[-1, -1]
    u.tilde <- unlist(ranef(object), recursive = F)



  h1 <- compute_h1(Xc = X, Z = Z, su, se2)

  # variance-covariance matrix for  BLUPs
    var.u <- lapply(1:length(su),
                    function(i) su[i]^2 * t(Z[[i]]) %*% h1 %*% Z[[i]] / se2
                    )

  # decomposition works with centered matrices!
    deco <- decomp(X = X, Z = Z,
                   se2 = se2, su = su,
                   b.hat = b.hat, S.b.hat = S.b.hat,  u.tilde = u.tilde,
                   var.u =var.u, h1 = h1
                   )
  deco <- structure(c(var.y = var(object@frame[,1]),
                      deco), 
                    class = "varExp")
  return(deco)
}






#' @param X Design matrix for fixed effects
#' @param Z List of desgin matrices for random effects (one list element for each effect)
#' @importFrom stats model.frame model.response
#' @export
#' @rdname varianceExplained
varianceExplained.mmer <- function(object, X, Z, ...){   
  # center matrices
    X <- scale(X, center = TRUE, scale = FALSE)
    Z <- lapply(Z, function(x) scale( x, center = TRUE, scale = FALSE))

  # get variance components
    se2 <- as.numeric( object$sigma$units )
    su<-  unlist(object$sigma[-which(names(object$sigma)=="units")],
           recursive = F)

  # get estimates
    b.hat <- as.numeric( object$Beta$Estimate )[-1]
    S.b.hat <- object$VarBeta[-1, -1]
    u.tilde <-  lapply(object$U, function(x) x$y)

  # variance-covariance matrix for  BLUPs
    var.u <-  lapply(object$VarU, function(x) x$y)

  h1 <- compute_h1(Xc = X, Z = Z, su, se2)

  # decomposition works with centered matrices!
  deco <- decomp(X = X, Z = Z,
                 se2 = se2, su = su,
                 b.hat = b.hat, S.b.hat = S.b.hat,
                 u.tilde = u.tilde,
                 var.u =var.u, h1 = h1
  )
  deco= structure(c(var.y = var(model.response(model.frame(object$call$fixed,
                                                           data = object$dataOriginal)
                                               )), 
                    deco), 
                  class = "varExp")
  return(deco)
}
