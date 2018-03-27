#include <sampling_utlities.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::export]]
// ********* Declaration ***********//

void Update_Z(int, int, int, double,
              NumericMatrix, NumericMatrix, NumericMatrix, 
              IntegerVector *, IntegerVector *,
              NumericMatrix *, NumericMatrix *,
              Function);

void Update_mu_u(int, int,
                 double, double, arma::vec, arma::mat,
                 IntegerVector,
                 NumericMatrix, NumericMatrix,
                 NumericMatrix *, NumericMatrix *,
                 Function,
                 Function); // Gibbs to update mu

void Update_M(int,
               double, double,
               int,
               IntegerVector,
               double *); // Gibbs to update M1


//*********** Defination ****************//

// Gibbs to update Z
void Update_Z(int i, int dimP,
                int K, double M,
                NumericMatrix dataX, NumericMatrix Mu, NumericMatrix Sigma_3D, 
                IntegerVector * Z, IntegerVector * N,
                NumericMatrix * sum_beta, NumericMatrix * sum_S3D,
                Function dmvn){
  (*N)((*Z)(i))--;
  NumericVector hatP(K);
  
  // stick-breaking
  NumericVector v(K);
  NumericVector pi(K);
  
  for(int k = 0; k < K; k++){
    if(k == K - 1){
      v(k) = 1;
    }else{
      double tmp_sum = 0;
      for(int kk = k+1; kk < K; kk++)
        tmp_sum += (*N)(kk);
      v(k) = rbeta(1, 1 + (*N)(k), M + tmp_sum)(0);
    }
    
    if(k == 0){
      pi(k) = v(k);
    }else{
      double tmp_prod = 1;
      for(int kk = 0; kk < k; kk++)
        tmp_prod *= (1-v(kk));
      pi(k) = v(k) * tmp_prod;
    }
    
    arma::mat tmpS = arma::zeros(dimP,dimP);
    arma::vec tmpB = arma::zeros(dimP);
    arma::vec tmpMu = arma::zeros(dimP);
    for(int p = 0; p < dimP; p++){
      tmpB(p) = dataX(i,p);
      tmpMu(p) = Mu(p,k);
      tmpS(p,p) = Sigma_3D(p, k*dimP + p);
      for(int pp = p + 1; pp < dimP; pp++){
        tmpS(p,pp) = Sigma_3D(p, k * dimP + pp);
        tmpS(pp,p) = tmpS(p,pp);
      }
    }
    
    hatP(k) = log(pi(k));
    hatP(k) += as<double>(dmvn(tmpB.t(), tmpMu.t(), tmpS,true));
  }
  
  // Update Z
  hatP = exp(hatP-max(hatP));
  if(sum(hatP) != 0){
    hatP = hatP/(sum(hatP));
    double tmp_rand = runif(1,0,1)(0);
    double tmp_sum = 0;
    int s = 0;
    for(s = 0; s < K; s++){
      tmp_sum += hatP(s);
      if(tmp_sum > tmp_rand)
        break;
    }
    (*Z)(i) = s;
  }
  
  (*N)((*Z)(i))++;
  
  for(int p = 0; p < dimP; p++){
    (*sum_beta)(p,(*Z)(i)) += dataX(i,p);
    (*sum_S3D)(p, (*Z)(i) * dimP + p) += dataX(i,p) * dataX(i,p);
    for(int pp = p+1; pp < dimP; pp++){
      (*sum_S3D)(p, (*Z)(i) * dimP + pp) += dataX(i,p) * dataX(i,pp);
      (*sum_S3D)(pp, (*Z)(i) * dimP + p) += dataX(i,pp) * dataX(i,p);
    }
  }
}

// Gibbs to update mu_logit_rho_u, mu_log_sd2_u
void Update_mu_u(int dimP, int K,
                 double k0, double v0, arma::vec beta0, arma::mat S0,
                 IntegerVector N,
                 NumericMatrix sum_beta, NumericMatrix sum_S3D,
                 NumericMatrix * Mu, NumericMatrix * Sigma_3D,
                 Function rmvn,
                 Function riwish){
  for(int k = 0; k < K; k++){
    
    // Update Mu_b[k], Sigma_b[k]
    double hat_k = k0 + N[k];
    double hat_v = v0 + N[k];
    
    arma::vec sum_beta_k = arma::zeros(dimP);
    arma::mat sum_S_k = arma::zeros(dimP, dimP);
    for(int p = 0; p < dimP; p++){
      sum_beta_k(p) = sum_beta(p,k);
      sum_S_k(p,p) = sum_S3D(p, k*dimP + p);
      for(int pp = p+1; pp < dimP; pp++){
        sum_S_k(p,pp) = sum_S3D(p,k*dimP + pp);
        sum_S_k(pp,p) = sum_S_k(p,pp);
      }
    }
    
    arma::vec hat_beta = (k0 * beta0 + sum_beta_k)/hat_k;
    arma::mat hat_S = S0 + sum_S_k + k0 * beta0 * beta0.t() - hat_k * hat_beta * hat_beta.t();
    
    NumericMatrix new_Sigma = as<NumericMatrix>(riwish(hat_v,hat_S));
    //arma::mat new_Sigma = arma::eye(dimP, dimP) * 0.09;
    NumericVector new_Mu_b = as<NumericVector>(rmvn(1,hat_beta,new_Sigma/hat_k));
    for(int p = 0; p < dimP; p++){
      (*Mu)(p,k) = new_Mu_b(p);
      (*Sigma_3D)(p, k*dimP+p) = new_Sigma(p,p);
      for(int pp = p+1; pp < dimP; pp++){
        (*Sigma_3D)(p, k*dimP+pp) = new_Sigma(p,pp);
        (*Sigma_3D)(pp, k*dimP+p) = new_Sigma(pp,p);
      }
    }
  }
}

// Gibbs to update M1
void Update_M(int dimI,
               double M_c, double M_d,
               int K,
               IntegerVector Z,
               double * M){
  double eta = R::rbeta(*M + 1, dimI);
  IntegerVector k_ind(K);
  for(int i = 0; i < dimI; i++)
    k_ind(Z(i)) = 1;
  int k = sum(k_ind);
  
  arma::vec hatP = arma::zeros(2);
  hatP(0) = (M_c + k -1)/(M_c + k - 1 + dimI*(M_d - log(eta)));
  hatP(1) = dimI*(M_d - log(eta))/(M_c + k - 1 + dimI*(M_d - log(eta)));
  
  arma::vec tmp_seq = arma::zeros(2);
  tmp_seq(1) = 1;
  int ind = RcppArmadillo::sample(tmp_seq,1,false,hatP)(0);
  if(ind == 0){
    *M = rgamma(1, M_c+k, M_d-log(eta))(0);
  }else{
    *M = rgamma(1, M_c+k-1, M_d-log(eta))(0);
  }
}


