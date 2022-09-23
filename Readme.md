Decomposition of the Explained Variation in the Linear Mixed Model
================

-   <a href="#introduction" id="toc-introduction">Introduction</a>
-   <a href="#installation" id="toc-installation">Installation</a>
-   <a href="#usage" id="toc-usage">Usage</a>
-   <a href="#reference" id="toc-reference">Reference</a>

# Introduction

The concept of variation explained is widely used to assess the
relevance of factors in the analysis of variance. In the linear model,
it is the main contribution to the coefficient of determination which is
widely used to assess the proportion of variation explained, to
determine model goodness-of-fit and to compare models with different
covariables. There has not been a consensus on a similar concept of
explained variation for the class of linear mixed models yet. Based on
the restricted maximum likelihood equations, we prove a full
decomposition of the sum of squares of the dependent variable in the
context of the variance components form of the linear mixed model. This
decomposition is dimensionless relative to the variation of the
dependent variable, has an intuitive and simple definition in terms of
variance explained, is additive for several random effects and reduces
to the decomposition in the linear model. Our result leads us to propose
a natural extension of the well-known adjusted coefficient of
determination to the linear mixed model. To this end, we introduce novel
measures for the explained variation which we allocate to specific
contributions of covariates associated with fixed and random effects.
These partial explained variations constitute easily interpretable
quantities, quantifying relevance of covariates associated with both
fixed and random effects on a common scale, and thus allowing to rank
their importance.

# Installation

To get the latest released version (master branch) of the R package from
GitHub:

``` r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("wiesenfa/explainedVariance", dependencies = TRUE)
```

# Usage

``` r
library(lme4)
library(explainedVariance)

fm2 <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy)

ve <- varianceExplained(fm2)
summary(ve)

# including bootstrap confidence intervals
bb <- bootVarianceExplained(fm2, nsim=500)
summary(bb)
# equivalently 
bb <- boot(ve, nsim=500)
summary(bb)

# normalized to the variance of Reaction
summary(proportionOf(varianceExplained(fm2), 
                     type = "dataset-specific"))
summary(proportionOf(bb))

# add simplified output to common model summary including coefficients, standard errors and p-values
expandResults(summary(proportionOf(bb)))

# Similarly applied to mmer objects from package sommer.
```

# Reference

Schreck, N., Wiesenfarth, M. (2022). Decomposition of the Explained
Variation in the Linear Mixed Model.
