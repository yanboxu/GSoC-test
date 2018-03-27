# setwd("~/Dropbox/Yanbo/LDTR/experiments/C_sampling")
require(MASS)
require(reshape2)
require(statmod)
require(mnormt)
require(matrixcalc)
require(fGarch)
require(stats)
require(mvnfast)
require(sn)
require(MCMCpack)
require(Rcpp)
#require(inline)
require(RcppArmadillo)
#require(rbenchmark)

logit <- function(p){
  return(log(p/(1-p)));
}

invlogit <- function(x){
  return(exp(x)/(1+exp(x)));
}

sourceCpp("sampling_C.cpp")

dpm_fit <- function(tol_iter, thin, K, sim_data){
  tol_iter = tol_iter
  thin = thin
  K = K
  
  dimI = sim_data$dimI
  dimP = sim_data$dimP
  
  dataX = sim_data$dataX
  
  paramZ = sim_data$Z
  paramMu = sim_data$Mu
  paramS = sim_data$Sigma
  
  # Initialization
  # NIW params for Mu_b, Sigma_b
  k0 = 1;
  v0 = dimP+1;
  beta0 = rep(0,dimP);
  S0 = diag(1,dimP);
  
  # DPM parameter M for clustering
  M_c = 2;
  M_d = 2;
  
  # cluster indicator
  Z = sample(0:(K-1),dimI,replace = TRUE);
  
  # mu[k], Sigma[k], hyperparameters for kth cluster for Beta
  Mu = array(0,dim=c(dimP,K));
  Sigma = array(0,dim=c(dimP,dimP,K));
  for(k in 1:K){
    # Mu[,k] = paramMu[,k];
    Sigma[,,k] = diag(1,dimP); #xxxxxxxxxxxxxxxxxxx 0.2
    # Sigma[,,k] = paramS[,,k]; #xxxxxxxxxxxxxxxxxxx 0.2
  }
  
  # M
  M = rgamma(1, M_c, M_d);
  
  Sigma_3D = array(0, dim=c(dimP, dimP * K));
  for(k in 1:K)
    for(p in 1:dimP)
      for(pp in 1:dimP)
        Sigma_3D[p, (k-1)*dimP + pp] = Sigma[p,pp,k];
  
  modelFit_DPM = sampling_dpm_C(tol_iter, thin, K,
                                dimI, dimP, dataX,
                                Z, Mu, Sigma_3D, M,
                                k0, v0, beta0, S0, M_c, M_d,
                                mvnfast::rmvn, mvnfast::dmvn, rinvgamma, riwish)
  
  return(modelFit_DPM);
}