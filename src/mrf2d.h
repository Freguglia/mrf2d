#ifndef MRF2D_H
#define MRF2D_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
Rcpp::NumericVector conditional_probabilities_multinomial(const Rcpp::IntegerMatrix &Z, Rcpp::IntegerVector position, Rcpp::IntegerMatrix R, const arma::fcube &theta, int N, int M, int n_R, int C);
arma::dmat table_relative(const Rcpp::IntegerMatrix &Z, Rcpp::IntegerVector r, int C, bool prop);
arma::dmat table_relative(const Rcpp::IntegerMatrix &Z, Rcpp::IntegerVector r, int C, bool prop);

#endif
