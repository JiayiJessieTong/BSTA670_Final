######### Function 4: 2nd order gradient of log-likelihood ########
## Input: X, n*p matrix (n: total number of patients. p: dimension of covariates) 
##        beta, p*1 vector (coefficients of the covariates)
## Output: 2nd order gradient of log-likelihood
Lgradient2 =function(beta,X){
  design = cbind(1,X)
  Z=expit(design%*%beta)
  t(c(-1*Z*(1-Z))*design)%*%design/nrow(X)
}
##### ------------------------------------ #####