##################################################
###### Rcode for BSTA670 final project ###########
##################################################
setwd("/Users/jiayito/Dropbox/000_Penn_Phd/1_2021_5_Spring/BSTA670/Final/BSTA670_final/R_Code")

## Step I: setup the enviroment
## libraries
library(mvtnorm)
library(ggplot2)
library(matlib)

## Step II: load the functions
source("BSTA670_final_functions.R")



##################################################
###Distributed algorithm for logistic regression##
###----------------- ODAL ------------------ #####
##################################################

set.seed(4321)
### ------------ setup values ----------- ###
beta_true = c(-1,1,-1,1,-1)


### ------ Replication of Setting A ------ ###
### fix the number of sites ###
### sites have the same number of patients
### increase the number of patients in each site ###
K_1 = 10 # total number of sites
n_1_list = seq(100,1000,by=100)  # number in one site
N_1_list = K_1 * n_1_list # overall patients across K_1 sites


### run the ODAL function ###
Nsim = 50
MSE_result_pooled = MSE_result_local = MSE_result_ODAL = rep(0, length(N_1_list))
for (i in 1:length(N_1_list)){
  print(paste("The total sample size across",K_1,"sites is", N_1_list[i]))
  MSE_result_pooled_tmp = 
    MSE_result_local_tmp = 
    MSE_result_ODAL_tmp = matrix(NA, ncol = 5, nrow = Nsim)
  for (iter in 1:Nsim){
    # generate data
    data = data_generator(N_1_list[i], beta_true, K_1, n_1_list[i]) 
    
    # run ODAL
    out = main_run_ODAL(data, beta_true)
    MSE_result_pooled_tmp[iter,] =  out$MSE_pooled
    MSE_result_local_tmp[iter,] =  out$MSE_local
    MSE_result_ODAL_tmp[iter,] =  out$MSE_ODAL
  }
  MSE_result_pooled[i] = sum(apply(MSE_result_pooled_tmp, 2, mean))/length(beta_true)
  MSE_result_local[i] = sum(apply(MSE_result_local_tmp, 2, mean))/length(beta_true)
  MSE_result_ODAL[i] = sum(apply(MSE_result_ODAL_tmp, 2, mean))/length(beta_true)
}

plot(N_1_list, MSE_result_local,
     type="b", col=rgb(26/255,133/255,172/255,0.5), lwd=2, pch=15, lty = 1,
     xlab="Total number of patients across 10 sites", ylab="MSE", ylim = c(0, 1.1))
lines(N_1_list, MSE_result_pooled, lty = 2,
      type="b", col=rgb(26/255,133/255,0/255,0.5), lwd=2, pch=17)
lines(N_1_list, MSE_result_ODAL, lty = 3,
      type="b", col=rgb(255/255,101/255,80/255,0.5), lwd=2, pch=19)
legend(6000, 1, legend = c("Local", "Pooled", "ODAL"), 
       lwd=2, lty = c(1,2,3), col=c(rgb(26/255,133/255,172/255,0.5), 
                                    rgb(26/255,133/255,0/255,0.5), 
                                    rgb(255/255,101/255,80/255,0.5)), 
       pch=c(15, 17, 19), bty = "n")










