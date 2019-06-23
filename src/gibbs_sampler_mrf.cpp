#include <RcppArmadillo.h>
#include "mrf2d.h"

using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
IntegerMatrix gibbs_sampler_mrf2d(IntegerMatrix init_Z, IntegerMatrix R, const arma::fcube theta, int n_steps){
  int N = init_Z.nrow(); const int M = init_Z.ncol();
  int C = theta.n_rows - 1;
  IntegerMatrix Z(N,M); Z = Rcpp::clone(init_Z);
  IntegerVector values = seq_len(C+1) - 1;
  IntegerVector order = seq_len(N*M);
  int x, y;
  int n_R = R.nrow();
  const int num = N*M;
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
