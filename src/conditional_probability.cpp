#include <RcppArmadillo.h>
#include "mrf2d.h"

using namespace Rcpp;


// Computes the conditional probabilities in a specific pixel given an image
// and MRF parameters.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector conditional_probabilities_mrf(const IntegerMatrix &Z,
                                            IntegerVector position,
                                            IntegerMatrix R,
                                            const arma::fcube &theta,
                                            int N, int M,
                                            int n_R, int C){

  IntegerVector this_pos(2);
  NumericVector probs(C+1);
  float this_prob;
  int dx, dy;
  int x = position[0] -1; int y = position[1] -1;

  for(int value = 0; value <= C; value++){
    this_prob = 0;
    for(int i = 0; i < n_R; i++){
      dx = R(i,0); dy = R(i,1);
      if(0 <= x+dx && x+dx < N && 0 <= y+dy && y+dy < M){
        this_prob = this_prob + theta(value, Z(x+dx, y+dy), i);}
      if(0 <= x-dx && x-dx < N && 0 <= y-dy && y-dy < M){
        this_prob = this_prob + theta(Z(x-dx, y-dy), value, i);}
    }
    probs[value] = exp(this_prob);
  }
  return(probs/sum(probs));
}
