#include <RcppArmadillo.h>
#include "mrf2d.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double log_pl_mrf(IntegerMatrix Z,  IntegerMatrix R, const arma::fcube theta){
  int N = Z.nrow(); const int M = Z.ncol();
  int n_R = R.nrow();
  int C = theta.n_rows - 1;
  double log_pl = 0.0;
  IntegerVector position;
  int zij;

  for(int i = 0; i < N; i++){
    for(int j = 0; j < M; j++){
      zij = Z(i,j);
      position[0] = i+1; position[1] = j+1;
      log_pl += log(conditional_probabilities_mrf(Z, position, R, theta, N, M, n_R, C)[zij]);
    }
  }
  return(log_pl);
}
