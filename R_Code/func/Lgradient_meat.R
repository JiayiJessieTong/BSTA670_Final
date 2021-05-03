######### Function 5: meat of the sandwich variance estimator ########
## Input: X, n*p matrix (n: total number of patients. p: dimension of covariates) 
## Output: meat of the sandwich variance estimator 
Lgradient_meat = function(beta,X,Y){
  design = cbind(1,X)
  Z = expit(design%*%beta)
  t(c(Y-Z)*design)%*%(c(Y-Z)*design)/nrow(X)
}
##### ------------------------------------ #####