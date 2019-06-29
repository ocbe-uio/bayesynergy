
#include "RcppArmadillo.h"
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat computeDelta(List Delta_list, arma::vec x1, arma::vec x2){
  
  double gamma0 = Delta_list["gamma0"];
  double gamma1 = Delta_list["gamma1"];
  double gamma2 = Delta_list["gamma2"];
  arma::mat B = Delta_list["B"];
  int n1 = B.n_rows;
  int n2 = B.n_cols;
  arma::mat Delta(n1,n2);
  
    for(int i = 0; i < n1; i++){
      for(int j = 0; j < n2; j++){
        Delta(i,j) = gamma0 + gamma1 * x1(i) + gamma2 * x2(j) + B(i,j);
      }
    }
  return Delta;
}