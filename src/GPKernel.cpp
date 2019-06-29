// Function to compute the Gaussian process kernel

// [[Rcpp::depends(BH)]]
#include "RcppArmadillo.h"
#define BOOST_DISABLE_ASSERTS
#include <boost/math/special_functions/bessel.hpp>
using namespace Rcpp;
using namespace boost;

// [[Rcpp::export]]
arma::mat GPKernel(arma::mat x_vec, List model_spec){
  int n3 = x_vec.n_rows;
  double ell = model_spec["ell"], nu = 0, alpha = 0, r_arg = 0;
  int type  = model_spec["type"];
  
  //Empty kernel matrix
  arma::mat Kxx(n3,n3,arma::fill::eye), eye_n3(n3,n3,arma::fill::eye);
  
  if(type == 2){ // double exponential
    for(int i = 0; i < n3; i++){
      for(int j = i; j < n3; j++){
        Kxx(i,j) = exp(-0.5 * accu((x_vec.row(i) - x_vec.row(j)) % (x_vec.row(i) - x_vec.row(j)))/pow(ell,2));
        Kxx(j,i) = Kxx(i,j);
      }
    }
  }
  if(type == 3){ // Matern(nu)
    nu = model_spec["nu"];
    for(int i = 0; i < n3; i++){
      for(int j = i+1; j < n3; j++){
        r_arg = sqrt(2*nu) * sqrt(accu((x_vec.row(i) - x_vec.row(j)) % (x_vec.row(i) - x_vec.row(j))))/ell;
        if(r_arg > 0){
          Kxx(i,j) = pow(2, 1 - nu) / tgamma(nu) * pow(r_arg, nu) * boost::math::cyl_bessel_k(nu, r_arg);
        }else{
          Kxx(i,j) = 1;
        }
        Kxx(j,i) = Kxx(i,j);
      }
    }
  }
  if(type == 4){ // rational quadratic
    alpha = model_spec["alpha"];
    for(int i = 0; i < n3; i++){
      for(int j = i; j < n3; j++){
        Kxx(i,j) = pow(1 + (0.5 * accu((x_vec.row(i) - x_vec.row(j)) % (x_vec.row(i) - x_vec.row(j)))/pow(ell,2)/alpha), -alpha);
        Kxx(j,i) = Kxx(i,j);
      }
    }
  }
  
  //For bad conditioning
  Kxx = Kxx + pow(10,-10) * eye_n3;
  return Kxx;
}