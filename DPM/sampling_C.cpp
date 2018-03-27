#include <sampling_update.h>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::export]]


List sampling_dpm_C(const int tol_iter, const int thin, const int K,
                    const int dimI, const int dimP, const NumericMatrix dataX,
                    IntegerVector Z, NumericMatrix Mu, NumericMatrix Sigma_3D, double M,
                    const int k0, const int v0, const NumericVector beta0, const NumericMatrix S0, 
                    const double M_c, const double M_d,
                    Function rmvn, Function dmvn, Function rinvgamma, Function riwish){
  
  int _Thined_tol_iter = (int) (tol_iter/thin) + 1;
  
  NumericMatrix __param_3DBeta(dimP, dimI * _Thined_tol_iter);
  NumericMatrix __param_3DMu(dimP, K * _Thined_tol_iter);
  NumericMatrix __param_4DSigma(dimP, dimP * K * _Thined_tol_iter);
  
  NumericMatrix __param_Z(dimI, _Thined_tol_iter);
  NumericVector __param_M(_Thined_tol_iter);
  
  NumericVector __trainll(_Thined_tol_iter);
  
  IntegerVector N(K);
  for(int i = 0; i < dimI; i++){
    N(Z(i))++;
  }
  
  int __thined_iter = 0;
  for(int iter = 0; iter < tol_iter; iter++){
    NumericMatrix sum_beta(dimP, K);
    NumericMatrix sum_S3D(dimP, dimP * K);
    
    __trainll(__thined_iter) = 0;
    //Rcout << std::endl << "Iter " << iter; 
    for(int i = 0; i < dimI; i++){
      Update_Z(i, dimP,
                 K, M,
                 dataX, Mu, Sigma_3D,
                 & Z, & N, & sum_beta, & sum_S3D,
                 dmvn);
      
      
      if(iter % thin == 0){
        // calculate log-likelihood
        NumericMatrix tmpS(dimP,dimP);
        NumericVector tmpMu(dimP);
        for(int p=0; p < dimP; p++){
          tmpMu(p) = Mu(p,Z(i));
          for(int pp=0; pp < dimP; pp++)
            tmpS(p,pp) = Sigma_3D(p, Z(i)*dimP + pp);
        }
        __trainll(__thined_iter) += as<double>(dmvn(dataX(i,_),tmpMu, tmpS, true));
        if(i % 20 == 0){
          Rcout << ".";
        }
        __param_Z(i, __thined_iter) = Z(i);
      }
    }
    
    arma::vec a_beta0 = as<arma::vec>(beta0);
    arma::mat a_S0 = as<arma::mat>(S0);
    Update_mu_u(dimP, K,
                k0, v0, a_beta0, a_S0, N,
                sum_beta, sum_S3D,
                & Mu, & Sigma_3D,
                rmvn, riwish);

    Update_M(dimI,
              M_c, M_d, K, Z,
              & M);
    
    if(iter % thin == 0){
      Rcout << std::endl << "thined iter: " << __thined_iter;
      for(int k = 0; k < K; k++){
        for(int p = 0; p < dimP; p++){
          __param_3DMu(p, __thined_iter * K + k) = Mu(p,k);
          __param_4DSigma(_, __thined_iter * (K * dimP) + k * dimP + p) = Sigma_3D(_, k * dimP + p);
        }
      }
      __param_M(__thined_iter) = M;
      __thined_iter++;
    }
  }
  
  List output;
  __param_3DBeta.attr("dim") = IntegerVector::create(dimP,dimI,_Thined_tol_iter);
  __param_3DMu.attr("dim") = IntegerVector::create(dimP,K,_Thined_tol_iter);
  __param_4DSigma.attr("dim") = IntegerVector::create(dimP,dimP,K,_Thined_tol_iter);
  
  output["Mu"] = __param_3DMu;
  output["Sigma"] = __param_4DSigma;
  
  output["Z"] = __param_Z; 
  output["M"] = __param_M;
  output["ll"] = __trainll;
  return(output);
}