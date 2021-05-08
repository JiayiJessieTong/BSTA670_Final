##################################################
###### Rcode for BSTA670 final project ###########
##################################################

## !! This path is needed to be changed 
path = "/Users/jiayito/Dropbox/000_Penn_Phd/1_2021_5_Spring/BSTA670/Final/BSTA670_final/R_Code"
setwd(path)

## Step I: setup the enviroment
## libraries
library(ggplot2)
library(matlib)
library(survival)
library(parallel)
library(miceadds)

## Step II: load the functions
## !!!!! This path is needed to be changed 
path_func = "/Users/jiayito/Dropbox/000_Penn_Phd/1_2021_5_Spring/BSTA670/Final/BSTA670_final/R_Code/func"
source.all(path_func, grepstring="\\.R",  print.source=TRUE, file_sep="__")


## Step III: run ODAL
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
     xlab="Total number of patients across 10 sites", ylab="MSE", ylim = c(0, 1.1), 
     xaxt='n', main = "Replication A: fixed site number, increasing site size")
axis(side = 1,at =N_1_list,
     labels=N_1_list,lwd.ticks = TRUE)
lines(N_1_list, MSE_result_pooled, lty = 2,
      type="b", col=rgb(26/255,133/255,0/255,0.5), lwd=2, pch=17)
lines(N_1_list, MSE_result_ODAL, lty = 3,
      type="b", col=rgb(255/255,101/255,80/255,0.5), lwd=2, pch=19)
legend(6000, 1, legend = c("Local", "Pooled", "ODAL"), 
       lwd=2, lty = c(1,2,3), col=c(rgb(26/255,133/255,172/255,0.5), 
                                    rgb(26/255,133/255,0/255,0.5), 
                                    rgb(255/255,101/255,80/255,0.5)), 
       pch=c(15, 17, 19), bty = "n")
##### ------------------------------------ #####





### ------ Replication of Setting B ------ ###
### fix patients in each site ###
### sites have the same number of patients
### increase the number of patients in each site ###
K_1_list = c(2, 5, 10, 15, 30, 40, 50, 75, 90, 100) # total number of sites
n_1 = 1000  # number in one site
N_1_list = K_1_list * n_1 # overall patients across K_1 sites


### run the ODAL function ###
Nsim = 50
MSE_result_pooled = MSE_result_local = MSE_result_ODAL = rep(0, length(N_1_list))
for (i in 1:length(N_1_list)){
  print(paste("The total sample size across",K_1_list[i],"sites is", N_1_list[i]))
  MSE_result_pooled_tmp = 
    MSE_result_local_tmp = 
    MSE_result_ODAL_tmp = matrix(NA, ncol = 5, nrow = Nsim)
  for (iter in 1:Nsim){
    # generate data
    data = data_generator(N_1_list[i], beta_true, K_1_list[i], n_1) 
    
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

plot(K_1_list, MSE_result_local,
     type="b", col=rgb(26/255,133/255,172/255,0.5), lwd=2, pch=15, lty = 1,
     xlab="Total number of sites", ylab="MSE", ylim = c(0, 0.2), 
     xaxt='n', main = "Replication B: fixed site size, increasing site number")
axis(side = 1,at =K_1_list,
     labels=K_1_list,lwd.ticks = TRUE)
lines(K_1_list, MSE_result_pooled, lty = 2,
      type="b", col=rgb(26/255,133/255,0/255,0.5), lwd=2, pch=17)
lines(K_1_list, MSE_result_ODAL, lty = 3,
      type="b", col=rgb(255/255,101/255,80/255,0.5), lwd=2, pch=19)
legend(60, 0.15, legend = c("Local", "Pooled", "ODAL"), 
       lwd=2, lty = c(1,2,3), col=c(rgb(26/255,133/255,172/255,0.5), 
                                    rgb(26/255,133/255,0/255,0.5), 
                                    rgb(255/255,101/255,80/255,0.5)), 
       pch=c(15, 17, 19), bty = "n")
##### ------------------------------------ #####




### -------------- Extension 1 -------------- ###
### different size sites
K_1 = 10 # total number of sites
n_1_list = matrix(NA, ncol = K_1, nrow = 10)
for (i in 1:K_1){
  tmp = sample(seq(200,1000,1),10)
  n_1_list[i,] = tmp # K_1 different size sites
}
N_1_list = apply(n_1_list, 1, sum) # overall patients across K_1 sites (10 settings)


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
    data = data_generator(N_1_list[i], beta_true, K_1, n_1_list[i,]) 
    
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

order_index = order(N_1_list)
plot(N_1_list[order_index], MSE_result_local[order_index],
     type="b", col=rgb(26/255,133/255,172/255,0.5), lwd=2, pch=15, lty = 1,
     xlab="Total number of patients across 10 sites", ylab="MSE", ylim = c(0, 0.35),
     xaxt='n', main = "Extension 1: 10 sites with different sizes")
axis(side = 1,at =N_1_list[order_index],
     labels=N_1_list[order_index],lwd.ticks = TRUE)
lines(N_1_list[order_index], MSE_result_pooled[order_index], lty = 2,
      type="b", col=rgb(26/255,133/255,0/255,0.5), lwd=2, pch=17)
lines(N_1_list[order_index], MSE_result_ODAL[order_index], lty = 3,
      type="b", col=rgb(255/255,101/255,80/255,0.5), lwd=2, pch=19)
legend(6000, 0.30, legend = c("Local", "Pooled", "ODAL"), 
       lwd=2, lty = c(1,2,3), col=c(rgb(26/255,133/255,172/255,0.5), 
                                    rgb(26/255,133/255,0/255,0.5), 
                                    rgb(255/255,101/255,80/255,0.5)), 
       pch=c(15, 17, 19), bty = "n")
##### ------------------------------------ #####



### -------------- Extension 2 -------------- ###
### heterogeneous covariates
K_1 = 10 # total number of sites
n_1_list = seq(100,1000,by=100)  # number in one site
N_1_list = K_1 * n_1_list # overall patients across K_1 sites

### run the ODAL function ###
Nsim = 10
MSE_result_pooled = MSE_result_clogit = MSE_result_local = MSE_result_ODAL = rep(0, length(N_1_list))
for (i in 1:length(N_1_list)){
  print(paste("The total sample size across",K_1,"sites is", N_1_list[i]))
  MSE_result_pooled_tmp = 
    MSE_result_local_tmp = 
    MSE_result_ODAL_tmp = matrix(NA, ncol = 5, nrow = Nsim)
  MSE_result_clogit_tmp =  matrix(NA, ncol = 4, nrow = Nsim)
  for (iter in 1:Nsim){
    # generate data
    data = data_generator_hetero(N_1_list[i], beta_true, K_1, n_1_list[i]) 
    
    # run ODAL
    out = main_run_ODAL_hetero(data, beta_true)
    MSE_result_pooled_tmp[iter,] =  out$MSE_pooled
    MSE_result_clogit_tmp[iter,] = out$MSE_clogit
    MSE_result_local_tmp[iter,] =  out$MSE_local
    MSE_result_ODAL_tmp[iter,] =  out$MSE_ODAL
  }
  MSE_result_pooled[i] = sum(apply(MSE_result_pooled_tmp, 2, mean))/length(beta_true)
  MSE_result_clogit[i] = sum(apply(MSE_result_clogit_tmp, 2, mean))/(length(beta_true)-1)
  MSE_result_local[i] = sum(apply(MSE_result_local_tmp, 2, mean))/length(beta_true)
  MSE_result_ODAL[i] = sum(apply(MSE_result_ODAL_tmp, 2, mean))/length(beta_true)
}

plot(N_1_list, MSE_result_local,
     type="b", col=rgb(26/255,133/255,172/255,0.5), lwd=2, pch=15, lty = 1,
     xlab="Total number of patients across 10 sites", ylab="Bias", ylim = c(0, 1.1), 
     xaxt='n', main = "Extension B: heterogenous prevalence across fixed K sites")
axis(side = 1,at =N_1_list,
     labels=N_1_list,lwd.ticks = TRUE)
lines(N_1_list, MSE_result_pooled, lty = 2,
      type="b", col=rgb(26/255,133/255,0/255,0.5), lwd=2, pch=17)
lines(N_1_list, MSE_result_ODAL, lty = 3,
      type="b", col=rgb(255/255,101/255,80/255,0.5), lwd=2, pch=19)
lines(N_1_list, MSE_result_clogit, lty = 4,
      type="b", col=rgb(172/255,85/255,255/255,0.5), lwd=2, pch=18)
legend(7000, 1.1, legend = c("Local", "Pooled","ODAL","clogit"), 
       lwd=2, lty = c(1,2,3,4),
       col=c(rgb(26/255,133/255,172/255,0.5), 
             rgb(26/255,133/255,0/255,0.5), 
             rgb(255/255,101/255,80/255,0.5),
             rgb(172/255,85/255,255/255,0.5)), 
       pch=c(15, 17, 19, 18), bty = "n")
##### ------------------------------------ #####





