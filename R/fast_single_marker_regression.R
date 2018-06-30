# libraries
library(Rcpp)

# functions
sourceCpp("mylm.cpp")

# extract pvalues
pvalue <- function(object, ...) {
  tval <- object$coefficients/object$stderr
  p.value = 2*pt(-abs(tval), df=object$df.residual)
  return(p.value[2,1])
}

# extract coefficients
# = mean(B) - mean(A)
# here Berkut is A, so negative means Berkut increase the values
coeff <- function(object, ...) {
  object$coefficients[2,1]
}

###### examples
marker[marker=="H"] = NA # only two alleles: A and B
marker2 = data.frame(t(marker))
marker3 = sapply(marker2, as.numeric)
rownames(marker3) = rownames(marker2)
marker3[1:10,1:3]
marker3 = marker3 - 1 # seems not necessary for pvalue, but may be for coefficient

system.time(aa <- lapply(dd[trait.col], function(x) mylm(x,marker3)))
pp = sapply(aa, function(x) {
  p1 = sapply(x, pvalue)
  names(p1) = colnames(marker3)
  return(p1)
})

head(pp[,1:6])

cc = sapply(aa, function(x) {
  p1 = sapply(x, coeff)
  names(p1) = colnames(marker3)
  return(p1)
})

head(cc[,1:6])


############ more complex models
## but slower
## test fastLm in
#library(RcppEigen)
# default contrast to type 3, same as SAS
# myanova2 = function(formula, data){# y ~ x + z
#   dd = na.omit(data[all.vars(formula)])
#   mm = RcppEigen::fastLmPure(model.matrix(formula[-2], dd), dd[,1])
#   pval = 2*pt(abs(mm$coefficients/mm$se), mm$df.residual, lower.tail=FALSE)
#   return(pval)
# }
#library(RcppEigen)
myanova2  <-  function(formula, data){# y ~ x + z
  mm  <-  RcppEigen::fastLm(formula, data)
  pval  <-  2*pt(abs(mm$coefficients/mm$se), mm$df.residual, lower.tail=FALSE)
  return(pval)
}


