#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


// Creates a matrix where each element (i,j) represents the counts of ocurrences
// i,j in relative position r.
// i: the 'central' pixel value; j: the pixel value in relative position r.
// [[Rcpp::export]]
arma::dmat table_relative(const IntegerMatrix &Z,
                          IntegerVector r,
                          int C){
  const int N = Z.nrow();
  const int M = Z.ncol();
  const int dx = r[0];
  const int dy = r[1];
  int a,b;
  arma::mat res = arma::zeros<arma::mat>(C+1,C+1);

  for(int i = std::max(0, -dx); i < std::min(N, N - dx); i++){
    for(int j = std::max(0, -dy); j < std::min(M, M - dy); j++){
      a = Z(i,j);
      b = Z(i+dx,j+dy); // NA's are stored as the lowest integer.
      if(a >= 0 && b >= 0){
        res(a,b) = res(a,b) + 1;
      }
    }
  }

  return(res);
}


// Creates an array where each slice is a matrix with the counts in a relative
// position identified by the row of R.
// [[Rcpp::export]]
arma::dcube table_relative_3d(const IntegerMatrix &Z, IntegerMatrix R, int C){
  int n_R = R.nrow();
  IntegerMatrix r(2);
  arma::dcube tab3d = arma::zeros<arma::dcube>(C+1, C+1, n_R);

  for(int i=0; i < n_R; i++){
    r[0] = R(i,0); r[1] = R(i,1);
    tab3d.slice(i) = table_relative(Z, r, C);
  }

  return(tab3d);
}
