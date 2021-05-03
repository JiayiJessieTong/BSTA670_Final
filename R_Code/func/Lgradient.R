######### Function 3: 1st order gradient of log-likelihood ########
## Input: X, n*p matrix (n: total number of patients. p: dimension of covariates) 
##        Y, n*1 binary vector (1 indicates case and 0 indicates control)
##        beta, p*1 vector (coefficients of the covariates)
## Output: 1st order gradient of log-likelihood
Lgradient = function(beta,X,Y){
  design = cbind(1,X)
  t(Y-expit(design%*%t(t(beta))))%*%design/length(Y)
}
##### ------------------------------------ #####