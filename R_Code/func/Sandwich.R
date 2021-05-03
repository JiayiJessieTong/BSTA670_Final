######### Function 8: sandwich var est ########
## Input:  beta, X, Y, N
## Output: sandwich variance estimate of est_ODAL
Sandwich <- function(beta,X,Y,N){
  mat_L1 = Lgradient_meat(beta,X,Y)
  mat_L2 = Lgradient2(beta,X)
  inv_L2 = solve.default(mat_L2)
  out = inv_L2%*%mat_L1%*%inv_L2/N
  return(out)
}
##### ------------------------------------ #####