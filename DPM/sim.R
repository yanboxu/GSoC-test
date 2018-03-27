require(mvnfast)
sim_dpm <- function(dimI, datafile){
  K = 3
  dimP = 2
  
  dataX = array(0, dim = c(dimI, dimP))
  Mu = array(0, dim=c(dimP,K))
  Sigma = array(0, dim=c(dimP, dimP,K))
  
  Mu[,1] = c(1,1)
  Mu[,2] = c(10,10)
  Mu[,3] = c(10,1)
  
  Sigma[,,1] = diag(1, dimP)
  Sigma[,,2] = diag(4, dimP)
  Sigma[,,3] = diag(1, dimP)
  
  Z = sample(1:K, dimI, replace = TRUE, prob = c(0.5, 0.5, 0.25))
  for(i in 1:dimI)
      dataX[i,] = rmvn(1, t(Mu[,Z[i]]), Sigma[,,Z[i]])
  
  output = list("dimI" = dimI, "dimP" = dimP, "dataX" = dataX, "K" = K, "Mu" = Mu, "Sigma" = Sigma, "Z" = Z)
  return(output)
}