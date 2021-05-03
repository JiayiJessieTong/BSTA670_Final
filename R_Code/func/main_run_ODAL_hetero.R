
######### Function 11: main_ODAL hetero function ########
## add the conditional regression model as the gold standard
## Input: Data, N*(ID + Y + p covariates)
## Output: MSE of pooled, ODAL, and local estimators
main_run_ODAL_hetero <-function(Data, beta_true){
  
  N = dim(Data)[1]
  p = dim(Data)[2]- 2
  K = length(unique(Data$ID))
  
  ########### METHOD 0: pooled analysis (naive) ##############
  # overall data
  X = as.matrix(Data[,-c(5,6)])
  Y = Data$Y
  
  # pooled beta and var
  fit_pooled = summary(glm(Y~X, family = "binomial"(link = "logit")))
  est_pooled = fit_pooled$coefficients[,1]
  sd_pooled = fit_pooled$coefficients[,2]
  var_pooled =sd_pooled^2
  ##### ------------------------------------ #####
  
  ########### METHOD 1: pooled analysis (clogit) ##############
  # pooled beta and var
  fit_clogit = summary(clogit(Y~X + strata(Data$ID)), method = "exact")
  est_clogit = fit_clogit$coefficients[,1]
  sd_clogit = fit_clogit$coefficients[,2]
  var_clogit =sd_clogit^2
  ##### ------------------------------------ #####
  
  ########### METHOD 2: local analysis ##############
  # extract the local data. Treat ID = 1 as local site
  local_ind = which(Data$ID == 1)
  local_data = Data[local_ind,]
  Xlocal = as.matrix(local_data[,-c(5,6)])
  Ylocal = local_data$Y
  
  # Local beta and var
  fit_local = summary(glm(Ylocal~Xlocal, family = "binomial"(link = "logit")))
  est_local = fit_local$coefficients[,1]
  sd_local = fit_local$coefficients[,2]
  var_local = sd_local^2
  ##### ------------------------------------ #####
  
  ########### METHOD 3: ODAL method ##############
  
  ######### Function 6: surrogate likelihood function ########
  ## Input: beta, p*1 vector (coefficients of the covariates)
  ##       Xlocal: local covariates
  ##       Ylocal: local outcome
  ## Output: surrogate likelihood value
  SL = function(beta){
    -Lik(beta,Xlocal,Ylocal) - L%*%beta
  }
  ##### ------------------------------------ #####
  
  # use local point estimate as initial value
  L = Lgradient(est_local,X,Y)
  
  # ODAL point estimate
  est_ODAL = optim(est_local,SL,control = list(maxit = 10000,reltol = 1e-10))$par
  
  # var of the estimator
  cov_ODAL = Sandwich(est_ODAL,Xlocal,Ylocal,N)
  var_ODAL = diag(cov_ODAL)
  ##### ------------------------------------ #####
  
  bias_pooled = (est_pooled - beta_true)^2
  MSE_pooled = bias_pooled + var_pooled
  
  bias_clogit = (est_pooled_clogit - beta_true[-1])^2
  MSE_clogit = bias_clogit + var_clogit
  
  bias_local = (est_local - beta_true)^2
  MSE_local = bias_local + var_local
  
  bias_ODAL = (est_ODAL - beta_true)^2
  MSE_ODAL = bias_ODAL + var_ODAL
  
  
  return(list(MSE_pooled = bias_pooled,
              MSE_clogit = bias_clogit,
              MSE_local = bias_local,
              MSE_ODAL = bias_ODAL))
}
##### ------------------------------------ #####
