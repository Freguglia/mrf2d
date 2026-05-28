#include <RcppArmadillo.h>
#include "mrf2d.h"

using namespace Rcpp;

// Computes log P(Z_ij = zij | Z_rest) directly in log-space:
//   log P = s_zij - log_sum_exp(s_0, ..., s_C)
// where s_v = sum of theta potentials for value v with all neighbors.
// Avoids: per-pixel NumericVector allocation, redundant exp/log roundtrip,
// and numerical overflow via the log-sum-exp trick.
static inline double log_cond_pixel(const IntegerMatrix &Z,
                                    const IntegerMatrix &R,
                                    const arma::fcube &theta,
                                    const int N, const int M,
                                    const int n_R, const int C,
                                    const int i, const int j,
                                    std::vector<double> &log_probs){
  int dx, dy;
  for(int v = 0; v <= C; v++){
    double lp = 0.0;
    for(int r = 0; r < n_R; r++){
      dx = R(r, 0); dy = R(r, 1);
      if(i+dx >= 0 && i+dx < N && j+dy >= 0 && j+dy < M)
        lp += theta(v, Z(i+dx, j+dy), r);
      if(i-dx >= 0 && i-dx < N && j-dy >= 0 && j-dy < M)
        lp += theta(Z(i-dx, j-dy), v, r);
    }
    log_probs[v] = lp;
  }

  // log-sum-exp trick: subtract max before exponentiation to prevent overflow
  double max_lp = *std::max_element(log_probs.begin(), log_probs.end());
  double sum_exp = 0.0;
  for(int v = 0; v <= C; v++)
    sum_exp += std::exp(log_probs[v] - max_lp);

  return log_probs[Z(i, j)] - max_lp - std::log(sum_exp);
}

static inline double log_cond_pixel_sub(const IntegerMatrix &Z,
                                        const LogicalMatrix &sub_mat,
                                        const IntegerMatrix &R,
                                        const arma::fcube &theta,
                                        const int N, const int M,
                                        const int n_R, const int C,
                                        const int i, const int j,
                                        std::vector<double> &log_probs){
  int dx, dy;
  for(int v = 0; v <= C; v++){
    double lp = 0.0;
    for(int r = 0; r < n_R; r++){
      dx = R(r, 0); dy = R(r, 1);
      if(i+dx >= 0 && i+dx < N && j+dy >= 0 && j+dy < M && sub_mat(i+dx, j+dy))
        lp += theta(v, Z(i+dx, j+dy), r);
      if(i-dx >= 0 && i-dx < N && j-dy >= 0 && j-dy < M && sub_mat(i-dx, j-dy))
        lp += theta(Z(i-dx, j-dy), v, r);
    }
    log_probs[v] = lp;
  }

  double max_lp = *std::max_element(log_probs.begin(), log_probs.end());
  double sum_exp = 0.0;
  for(int v = 0; v <= C; v++)
    sum_exp += std::exp(log_probs[v] - max_lp);

  return log_probs[Z(i, j)] - max_lp - std::log(sum_exp);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double log_pl_mrf(const IntegerMatrix Z,
                  const IntegerMatrix R,
                  const arma::fcube theta){

  int N = Z.nrow(); int M = Z.ncol();
  int n_R = R.nrow();
  int C = theta.n_rows - 1;
  double log_pl = 0.0;
  std::vector<double> log_probs(C + 1);

  for(int i = 0; i < N; i++){
    for(int j = 0; j < M; j++){
      log_pl += log_cond_pixel(Z, R, theta, N, M, n_R, C, i, j, log_probs);
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
  std::vector<double> log_probs(C + 1);

  for(int i = 0; i < N; i++){
    for(int j = 0; j < M; j++){
      if(sub_mat(i,j)){
        log_pl += log_cond_pixel_sub(Z, sub_mat, R, theta, N, M, n_R, C, i, j, log_probs);
      }
    }
  }
  return(log_pl);
}
