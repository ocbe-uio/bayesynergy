#ifndef GP_KERNEL_H
#define GP_KERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

arma::mat GPKernel(arma::mat x_vec, List kern_spec);

#endif
