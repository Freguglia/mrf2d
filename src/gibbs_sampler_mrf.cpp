#include <RcppArmadillo.h>
#include "mrf2d.h"

using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
IntegerMatrix gibbs_sampler_mrf2d(IntegerMatrix init_Z, IntegerMatrix R, const arma::fcube theta, int n_steps){
  int N = init_Z.nrow(); int M = init_Z.ncol();
  int C = theta.n_rows - 1;
  IntegerMatrix Z(N,M); Z = Rcpp::clone(init_Z);
  IntegerVector values = seq_len(C+1) - 1;
  IntegerVector order = seq_len(N*M);
  int x, y;
  int n_R = R.nrow();
  int num = N*M;
  IntegerVector position(2);
  NumericVector cprobs(C+1);

  for(int step = 0; step < n_steps; step++){
    order = sample(order, num, false);
    for(int i = 0; i < num; i++){
      x = (order[i] / M) + 1; y = (order[i] % M) + 1;
      position[0] = x; position[1] = y;
      cprobs = conditional_probabilities_mrf(Z, position, R, theta, N, M, n_R, C);
      Z(x-1, y-1) = sample(values, 1, false, cprobs)[0];
    }
  }
  return(Z);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
IntegerMatrix gibbs_sampler_mrf2d_sub(IntegerMatrix init_Z,
                                      LogicalMatrix sub_mat,
                                      LogicalMatrix fix_mat, IntegerMatrix R,
                                      const arma::fcube theta, int n_steps){
  int N = init_Z.nrow(); int M = init_Z.ncol();
  int C = theta.n_rows - 1;
  IntegerMatrix Z(N,M); Z = Rcpp::clone(init_Z);
  IntegerVector values = seq_len(C+1) - 1;
  IntegerVector order = seq_len(N*M);
  int x, y;
  int n_R = R.nrow();
  int num = N*M;
  IntegerVector position(2);
  NumericVector cprobs(C+1);

  for(int step = 0; step < n_steps; step++){
    order = sample(order, num, false);
    for(int i = 0; i < num; i++){
      x = (order[i] / M) + 1; y = (order[i] % M) + 1;
      if(sub_mat(x-1,y-1) & !fix_mat(x-1,y-1)){
        position[0] = x; position[1] = y;
        cprobs = conditional_probabilities_mrf_sub(Z, sub_mat, position, R, theta, N, M, n_R, C);
        Z(x-1, y-1) = sample(values, 1, false, cprobs)[0];
      }
    }
  }
  return(Z);
}
