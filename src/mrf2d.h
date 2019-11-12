#ifndef MRF2D_H
#define MRF2D_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
Rcpp::NumericVector conditional_probabilities_mrf(const Rcpp::IntegerMatrix &Z, const Rcpp::IntegerVector position, const Rcpp::IntegerMatrix R, const arma::fcube &theta, const int N, const int M, const int n_R, const int C);
Rcpp::NumericVector conditional_probabilities_mrf_sub(const Rcpp::IntegerMatrix &Z, const Rcpp::LogicalMatrix &sub_mat, const Rcpp::IntegerVector position, const Rcpp::IntegerMatrix R, const arma::fcube &theta, const int N, const int M, const int n_R, const int C);
arma::dmat table_relative(const Rcpp::IntegerMatrix &Z, const Rcpp::IntegerVector r, const int C, const bool prop);



#endif
