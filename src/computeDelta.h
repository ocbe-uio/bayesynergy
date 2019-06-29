#ifndef COMPUTE_DELTA_H
#define COMPUTE_DELTA_H

#include <RcppArmadillo.h>
using namespace Rcpp;

arma::mat computeDelta(List Delta_list, arma::vec x1, arma::vec x2);

#endif