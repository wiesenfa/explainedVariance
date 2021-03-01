print.varExp <-function(x,...)  cat("Variance decomposition object. Sum of variance components deviates from var(y) by, ", proportionOf(x)$error,".")

print.varExpProp <-function(x,...)   cat("Proportion of variance explained object. Sum of variance components deviates from var(y) by, ", proportionOf(x)$error/x$var.y*100,"%.")
