#include <RcppArmadillo.h>
#include "mrf2d.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double log_pl_mrf(const IntegerMatrix Z,
                  const IntegerMatrix R,
                  const arma::fcube theta){

  int N = Z.nrow(); int M = Z.ncol();
  int n_R = R.nrow();
  int C = theta.n_rows - 1;
  double log_pl = 0.0;
  double this_cond_prob = 0.0;
  IntegerVector position(2);
  int zij = -1;

  for(int i = 0; i < N; i++){
    for(int j = 0; j < M; j++){
      zij = Z(i,j);
      position[0] = i+1; position[1] = j+1;
      this_cond_prob = conditional_probabilities_mrf(Z, position, R, theta, N, M, n_R, C)(zij);
      log_pl += log(this_cond_prob);
    }
  }
  return(log_pl);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double log_pl_mrf_sub(const IntegerMatrix Z,
                      const LogicalMatrix sub_mat,
                      const IntegerMatrix R,
                      const arma::fcube theta){

  int N = Z.nrow(); int M = Z.ncol();
  int n_R = R.nrow();
  int C = theta.n_rows - 1;
  double log_pl = 0.0;
  double this_cond_prob = 0.0;
  IntegerVector position(2);
  int zij = -1;

  for(int i = 0; i < N; i++){
    for(int j = 0; j < M; j++){
      if(sub_mat(i,j)){
        zij = Z(i,j);
        position[0] = i+1; position[1] = j+1;
        this_cond_prob = conditional_probabilities_mrf_sub(Z, sub_mat, position, R, theta, N, M, n_R, C)(zij);
        log_pl += log(this_cond_prob);
      }
    }
  }
  return(log_pl);
}
