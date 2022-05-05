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


#' @importFrom lme4 lFormula ranef fixef VarCorr getME
#' @importFrom Matrix t
#' @importFrom stats sigma var vcov
#' @export
#' @rdname varianceExplained
varianceExplained.lmerMod <- varianceExplained.lmerModLmerTest <- function(object, ...){

  # get matrices
    X <-getME(object,"X")[,-1, drop = FALSE]
    X <- scale(X, center = TRUE, scale = FALSE)
    Z <-lapply(getME(object,"Ztlist"), Matrix::t)
    Z <- lapply(Z, function(x) scale( x, center = TRUE, scale = FALSE))
    
    # check whether correltated random intercept and random slope for same grouping variable
      correlatedREs <- sapply(object@cnms, function(x) length(x)>1)
      if (any (correlatedREs)) stop("Correlated random effects currently not supported. 
                                    Use two vertical bars (||) to specify multiple uncorrelated random effects for the same grouping variable, 
                                    e.g. (x||group)
                                    or equivantly (1|group) + (0+x|group).")
    # folgender code war um correlated random effects mit einem | zuszulassen, muss aber angepasst werden  
    #  Z0 <-lapply(lmod$reTrms$Ztlist,Matrix::t)
    # # check whether correltated random intercept and random slope for same grouping variable
    #   multipleREs <- lapply(object@cnms, function(x) length(x))
    # 
    # # lme4 combines all random effects of one grouping variable into a single Z, following code separates this
    #   Z=list()
    #   j=1
    #   for (i in 1:length(Z0)){
    #     if (multipleREs[[i]]==1) {
    #       Z[[j]] <- Z0[[i]]
    #       names(Z)[j] <- paste(names(object@cnms)[i],object@cnms[[i]][1], sep=".")
    #       j=j+1
    #     }  else {
    #       for (k in 1: multipleREs[[i]]){
    #         Z[[j]] <- Z0[[i]][, seq(k,
    #                                           by=multipleREs[[i]],
    #                                           length.out=ncol(Z0[[i]])/multipleREs[[i]]) ]
    #         names(Z)[j] <- paste(names(object@cnms)[i],object@cnms[[i]][k], sep=".")
    #         j = j +1
    # 
    #       }
    #     }
    #   }
    #   Z <- lapply(Z, function(x) scale( x, center = TRUE, scale = FALSE))
    #   rm(Z0)

      # alternative funktionen die hilfreich sein könnnen  
        # reTrms <- mkReTrms(findbars(lme4:::RHSForm(mc$formula)), object@frame)
        # 
        # data("Pixel", package="nlme")
        # mform <- pixel ~ day + I(day^2) + (day || Dog) 
        # (bar.f <- findbars(mform)) # list with 3 terms
        # mf <- model.frame(subbars(mform),data=Pixel)
        # rt <- mkReTrms(bar.f,mf)
        # names(rt)
        # expandDoubleVerts(mc$formula)
        # 
        # mform <- pixel ~ day + I(day^2) + (day | Dog) 
        # (bar.f <- findbars(mform)) 
        # lapply(bar.f,lme4:::mkBlist,fr=mf)
        # 
        
      
      
      
  # get variance components
    se2 <-sigma(object)^2
    su<-unlist(lapply(VarCorr(object), function(x) attr(x, "stddev")^2),
                   recursive = F)

  # get estimates
    b.hat <- fixef(object)[-1]
    S.b.hat <- vcov(object, full = FALSE)[-1, -1]
    u.tilde <- unlist(ranef(object,condVar = FALSE), recursive = F)

  # assign names to Z matrices
    # CHECK WHETHER NECESSARY!!!
    names(su)  <- names(u.tilde)
    

  h1 <- compute_h1(Xc = X, Z = Z, su, se2)

  # variance-covariance matrix for  BLUPs
    var.u <- sapply(names(su),
                    function(id) su[id]^2 * t(Z[[id]]) %*% h1 %*% Z[[id]] / se2,
                    simplify = FALSE
                    )

  # decomposition works with centered matrices!
    deco <- decomp(X = X, Z = Z,
                   se2 = se2, su = su,
                   b.hat = b.hat, S.b.hat = S.b.hat,  u.tilde = u.tilde,
                   var.u =var.u, h1 = h1
                   )
  var.y <- var(getME(object, "y"))
  deco <- structure(c(var.y = var.y,
                      deco,
                      error = var.y- deco$se2 - 
                        deco$Rx - 
                        sum(deco$Rz.1 + deco$Rz.2) -  
                        sum(deco$Rz.pairs, na.rm = TRUE) - 
                        sum(deco$Rxz)
  ), 
                    class = "VarExp")
  return(deco)
}






#' @param X Design matrix for fixed effects
#' @param Z *named* List of desgin matrices for random effects (one list element for each effect).
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
  var.y <- var(model.response(model.frame(object$call$fixed,
                                          data = object$dataOriginal)
                              ))
  deco= structure(c(var.y = var.y, 
                    deco,
                    error = var.y - deco$se2 - 
                      deco$Rx - 
                      sum(deco$Rz.1 + deco$Rz.2) -  
                      sum(deco$Rz.pairs, na.rm = TRUE) - 
                      sum(deco$Rxz)
  ), 
  class = "VarExp")
  return(deco)
}

#' @export
varianceExplained.matrix <- function(X, 
                                  Z, #list
                                  fixef, #fixed effects excluding intercept
                                  ranef,
                                  vcov, #excluding intercept
                                  var.u,
                                  var.eps,
                                  var.y,...
){   
  # center matrices
  X <- scale(X, center = TRUE, scale = FALSE)
  Z <- lapply(Z, function(x) scale( x, center = TRUE, scale = FALSE))
  
  # get variance components
  se2 <- var.eps #as.numeric( object$sigma$units )
  su<-  var.u #unlist(object$sigma[-which(names(object$sigma)=="units")],  recursive = F)
  
  # get estimates
  b.hat <- fixef #fixef(object)[-1]
  S.b.hat <- vcov #vcov(object, full = FALSE)[-1, -1]
  u.tilde <- ranef # unlist(ranef(object,condVar = FALSE), recursive = F)
  names(su)  <- names(u.tilde)
  
  h1 <- compute_h1(Xc = X, Z = Z, su, se2)
  # variance-covariance matrix for  BLUPs
  var.u <- sapply(names(su),
                  function(id) su[id]^2 * t(Z[[id]]) %*% h1 %*% Z[[id]] / se2,
                  simplify = FALSE
  )
  
  
  # decomposition works with centered matrices!
  deco <- decomp(X = X, Z = Z,
                 se2 = se2, su = su,
                 b.hat = b.hat, S.b.hat = S.b.hat,
                 u.tilde = u.tilde,
                 var.u =var.u, h1 = h1
  )
  deco= structure(c(var.y = var.y, 
                    deco,
                    error = var.y - deco$se2 - 
                      deco$Rx - 
                      sum(deco$Rz.1 + deco$Rz.2) -  
                      sum(deco$Rz.pairs, na.rm = TRUE) - 
                      sum(deco$Rxz)
  ), 
  class = "VarExp")
  return(deco)
}

