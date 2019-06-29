#ifndef LOGACCEPT_H
#define LOGACCEPT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

double logaccept(arma::mat y_mat, arma::vec p_ij_vec, arma::vec p_ij_vec_new, double s2_eps, double is_accept);

#endif