
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double logaccept(arma::mat y_mat, arma::vec p_ij_vec, arma::vec p_ij_vec_new, double s2_eps, double is_accept){
  double log_accept = 0;
  int n_rep = y_mat.n_cols;
  arma::vec z, z_new;
  
  for(int r = 0;  r < n_rep; r++){
    z = (y_mat.col(r) - p_ij_vec);
    z_new = (y_mat.col(r) - p_ij_vec_new);
    log_accept = log_accept - 0.5 * (accu(z_new % z_new) - is_accept * accu(z % z))/s2_eps;
     // (sum( (pow(y_mat.col(r) - p_ij_vec_new,2) - is_accept * pow(y_mat.col(r) - p_ij_vec,2))/s2_eps ));
  }
  
  return log_accept;
}