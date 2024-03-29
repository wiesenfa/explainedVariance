% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrap.R
\name{bootVarianceExplained}
\alias{bootVarianceExplained}
\alias{bootVarianceExplained.default}
\alias{bootVarianceExplained.lmerMod}
\alias{bootVarianceExplained.mmer}
\title{Derive bootstrap confidence intervals for variance decomposition (based on percentile method)}
\usage{
bootVarianceExplained(object, ...)

\method{bootVarianceExplained}{default}(object, ...)

\method{bootVarianceExplained}{lmerMod}(object, ...)

\method{bootVarianceExplained}{mmer}(
  object,
  X,
  Z,
  nsim = 1000,
  parallel = FALSE,
  progress = "time",
  ...
)
}
\arguments{
\item{object}{a \code{lmerMod} or \code{lmerModLmerTest} object created by \code{\link[lme4:lmer]{lme4::lmer()}} or \code{\link[lmerTest:lmer]{lmerTest::lmer()}}, respectively, or a \code{mmer} object created by  \code{\link[sommer:mmer]{sommer::mmer()}}.}

\item{...}{arguments passed to  \code{\link[lme4:bootMer]{lme4::bootMer()}}, in particular \code{nsim} for the number of simulations, the type of bootstrap and arguments for parallel computing.  In case of mmer() objects arguments  X, Z need to be provided}

\item{parallel}{only mmer objects. TRUE for parallelization. Initiallize parallization e.g. with \code{library(doParallel); registerDoParallel(cores = 8)} before}

\item{progress}{only mmer objects. passed to plyr::llply}
}
\description{
Derive bootstrap confidence intervals for variance decomposition (based on percentile method)
}
