#ifndef MRF2D_H
#define MRF2D_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
Rcpp::NumericVector conditional_probabilities_mrf(const Rcpp::IntegerMatrix &Z, Rcpp::IntegerVector position, Rcpp::IntegerMatrix R, const arma::fcube &theta, int N, int M, int n_R, int C);
Rcpp::NumericVector conditional_probabilities_mrf_sub(const Rcpp::IntegerMatrix &Z, const Rcpp::LogicalMatrix &sub_mat, Rcpp::IntegerVector position, Rcpp::IntegerMatrix R, const arma::fcube &theta, int N, int M, int n_R, int C);
arma::dmat table_relative(const Rcpp::IntegerMatrix &Z, Rcpp::IntegerVector r, int C, bool prop);



#endif
