#include <RcppArmadillo.h>
#include "mrf2d.h"

using namespace Rcpp;


// Computes the conditional probabilities in a specific pixel given an image
// and MRF parameters.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector conditional_probabilities_mrf(const IntegerMatrix &Z,
                                            const IntegerVector position,
                                            const IntegerMatrix R,
                                            const arma::fcube &theta,
                                            const int N, const int M,
                                            const int n_R, const int C){

  NumericVector probs(C+1);
  double this_prob;
  int dx, dy;
  int x = position[0] -1; int y = position[1] -1;

  for(int value = 0; value <= C; value++){
    this_prob = 0.0;
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

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector conditional_probabilities_mrf_sub(const IntegerMatrix &Z,
                                                const LogicalMatrix &sub_mat,
                                                const IntegerVector position,
                                                const IntegerMatrix R,
                                                const arma::fcube &theta,
                                                const int N, const int M,
                                                const int n_R, const int C){

  NumericVector probs(C+1);
  double this_prob;
  int dx, dy;
  int x = position[0] -1; int y = position[1] -1;

  for(int value = 0; value <= C; value++){
    this_prob = 0.0;
    for(int i = 0; i < n_R; i++){
      dx = R(i,0); dy = R(i,1);
      if(0 <= x+dx && x+dx < N && 0 <= y+dy && y+dy < M){
        if(sub_mat(x+dx, y+dy)){
          this_prob = this_prob + theta(value, Z(x+dx, y+dy), i);}}
      if(0 <= x-dx && x-dx < N && 0 <= y-dy && y-dy < M){
        if(sub_mat(x-dx, y-dy)){
          this_prob = this_prob + theta(Z(x-dx, y-dy), value, i);}}
    }
    probs[value] = exp(this_prob);
  }
  return(probs/sum(probs));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
IntegerMatrix icm_restoration_cpp(const IntegerMatrix &init_Z,
                                  const IntegerMatrix R,
                                  const arma::fcube &theta,
                                  const double corr_prob,
                                  const int cycles){
  int N = init_Z.nrow(); const int M = init_Z.ncol();
  int C = theta.n_rows - 1;
  IntegerMatrix Z(N,M); Z = Rcpp::clone(init_Z);
  IntegerVector values = seq_len(C+1) - 1;
  IntegerVector order = seq_len(N*M);
  int x, y;
  int n_R = R.nrow();
  long num = N*M;
  IntegerVector position(2);
  NumericVector cprobs(C+1);
  NumericVector postprobs(C+1);

  for(int step = 0; step < cycles; step++){
    order = sample(order, num, false);
    for(int i = 0; i < num; i++){
      x = (order[i] / M) + 1; y = (order[i] % M) + 1;
      position[0] = x; position[1] = y;
      cprobs = conditional_probabilities_mrf(Z, position, R, theta, N, M, n_R, C);
      postprobs = cprobs * corr_prob/C;
      postprobs[Z(x-1,y-1)] = cprobs[Z(x-1,y-1)] * (1 - corr_prob);
      Z(x-1, y-1) = which_max(postprobs);
    }
  }
  return(Z);
}

// Computes conditional probabilities in one specific position.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector cprob_ghm_one(const IntegerMatrix &Z,
                            const IntegerVector position,
                            const IntegerMatrix R, const arma::fcube &theta,
                            const int N, const int M,
                            const int n_R, const int C,
                            const NumericVector mus, const NumericVector sigmas,
                            const NumericMatrix &Y){

  NumericVector probs(C+1);
  NumericVector y(1);
  y[0] = Y(position[0] - 1, position[1] - 1);
  probs = conditional_probabilities_mrf(Z, position, R, theta, N, M, n_R, C);
  for(int i = 0; i <= C; i++){
    probs[i] = probs[i]*dnorm(y, mus[i], sigmas[i])[0];
  }
  return(probs/sum(probs));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector cprob_ghm_one_sub(const IntegerMatrix &Z,
                                const LogicalMatrix &sub_mat,
                                const IntegerVector position,
                                const IntegerMatrix R, const arma::fcube &theta,
                                const int N, const int M,
                                const int n_R, const int C,
                                const NumericVector mus,
                                const NumericVector sigmas,
                                const NumericMatrix Y){

  NumericVector probs(C+1);
  NumericVector y(1);
  y[0] = Y(position[0] - 1, position[1] - 1);
  probs = conditional_probabilities_mrf_sub(Z, sub_mat, position, R, theta,
                                            N, M, n_R, C);
  for(int i = 0; i <= C; i++){
    probs[i] = probs[i]*dnorm(y, mus[i], sigmas[i])[0];
  }
  return(probs/sum(probs));
}

// Computes conditional probabilities for all positions. Returns an array with
// probabilities for (position i,position j, class).
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::fcube cprob_ghm_all(const IntegerMatrix &Z,
                          const IntegerMatrix R,
                          const arma::fcube &theta,
                          const NumericVector mus,
                          const NumericVector sigmas,
                          const NumericMatrix &Y){

  int N = Z.nrow();
  int M = Z.ncol();
  int C = theta.n_rows - 1;
  int n_R = theta.n_slices;
  arma::fcube result(N, M, C+1);

  IntegerVector position(2);
  NumericVector cprobs(C+1);

  for(int i = 0; i < N; i++){
    for(int j = 0; j < M; j++){
      position[0] = i+1; position[1] = j+1;
      cprobs = cprob_ghm_one(Z, position, R, theta, N, M, n_R, C, mus, sigmas, Y);
      for(int k = 0; k<= C; k++){
        result(i,j,k) = cprobs[k];
      }
    }
  }
  return(result);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::fcube cprob_ghm_all_sub(const IntegerMatrix &Z,
                              const LogicalMatrix &sub_mat,
                              const IntegerMatrix R,
                              const arma::fcube &theta,
                              const NumericVector mus,
                              const NumericVector sigmas,
                              const NumericMatrix &Y){

  int N = Z.nrow();
  int M = Z.ncol();
  int C = theta.n_rows - 1;
  int n_R = theta.n_slices;
  arma::fcube result(N, M, C+1);

  IntegerVector position(2);
  NumericVector cprobs(C+1);

  for(int i = 0; i < N; i++){
    for(int j = 0; j < M; j++){
      if(sub_mat(i,j)){
        position[0] = i+1; position[1] = j+1;
        cprobs = cprob_ghm_one_sub(Z, sub_mat, position, R, theta,
                                   N, M, n_R, C, mus, sigmas, Y);
        for(int k = 0; k<= C; k++){
          result(i,j,k) = cprobs[k];
        }
      }
      else {
        for(int k = 0; k<= C; k++){
          result(i,j,k) = NA_REAL;
        }
      }
    }
  }
  return(result);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
IntegerMatrix icm_gaussian_cpp(const NumericMatrix &Y,
                               const IntegerMatrix R,
                               const IntegerMatrix &init_Z,
                               const arma::fcube &theta,
                               const NumericVector mus,
                               const NumericVector sigmas,
                               const int cycles){
  int N = Y.nrow(); int M = Y.ncol();
  int C = theta.n_rows - 1;
  IntegerMatrix Z(N,M); Z = Rcpp::clone(init_Z);
  IntegerVector values = seq_len(C+1) - 1;
  IntegerVector order = seq_len(N*M);
  int x, y;
  int n_R = R.nrow();
  long num = N*M;
  IntegerVector position(2);
  NumericVector cprobs(C+1);
  NumericVector pix(1);

  for(int step = 0; step < cycles; step++){
    order = sample(order, num, false);
    for(int i = 0; i < num; i++){
      x = (order[i] / M) + 1; y = (order[i] % M) + 1;
      position[0] = x; position[1] = y;
      pix[0] = Y(x - 1, y - 1);
      cprobs = cprob_ghm_one(Z, position, R, theta, N, M, n_R, C, mus, sigmas, Y);
      Z(x-1, y-1) = which_max(cprobs);
    }
  }
  return(Z);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
IntegerMatrix icm_gaussian_cpp_sub(const NumericMatrix &Y,
                                   const LogicalMatrix &sub_mat,
                                   const IntegerMatrix R,
                                   const IntegerMatrix &init_Z,
                                   const arma::fcube &theta,
                                   const NumericVector mus,
                                   const NumericVector sigmas,
                                   const int cycles){
  int N = Y.nrow(); int M = Y.ncol();
  int C = theta.n_rows - 1;
  IntegerMatrix Z(N,M); Z = Rcpp::clone(init_Z);
  IntegerVector values = seq_len(C+1) - 1;
  IntegerVector order = seq_len(N*M);
  int x, y;
  int n_R = R.nrow();
  long num = N*M;
  IntegerVector position(2);
  NumericVector cprobs(C+1);
  NumericVector pix(1);

  for(int step = 0; step < cycles; step++){
    order = sample(order, num, false);
    for(int i = 0; i < num; i++){
      x = (order[i] / M) + 1; y = (order[i] % M) + 1;
      if(sub_mat(x-1, y-1)){
        position[0] = x; position[1] = y;
        pix[0] = Y(x - 1, y - 1);
        cprobs = cprob_ghm_one_sub(Z, sub_mat, position, R, theta,
                                   N, M, n_R, C, mus, sigmas, Y);
        Z(x-1, y-1) = which_max(cprobs);
      }
    }
  }
  return(Z);
}


