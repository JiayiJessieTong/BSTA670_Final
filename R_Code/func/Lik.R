######### Function 2: likelihood function ########
## Input: X, n*p matrix (n: total number of patients, p: dimension of covariates) 
##        Y, n*1 binary vector (1 indicates case and 0 indicates control)
##        beta, p*1 vector (coefficients of the covariates)
## Output: log-likelihood value
Lik = function(beta,X,Y){
  design = cbind(1,X)
  sum(Y*(design%*%t(t(beta)))-log(1+exp(design%*%t(t(beta)))))/length(Y)
}
##### ------------------------------------ #####
