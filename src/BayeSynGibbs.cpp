// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "RcppArmadillo.h"
#include "GPKernel.h"
#include "computeDelta.h"
#include "logaccept.h"
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;

// [[Rcpp::export]]

List BayeSynGibbs(arma::mat y_mat, arma::mat x_mat, List model_spec, List Alg_param, List Hyper_param) {
  
  
 
  // Find missing values and create a working y_mat
  arma::mat missing_y(y_mat.n_rows,y_mat.n_cols,arma::fill::zeros);
  int n_missing = 0;
  for (int i = 0; i < y_mat.n_rows; i++){
    for (int j = 0; j < y_mat.n_cols; j++){
      if (Rcpp::NumericVector::is_na(y_mat(i,j))){
        missing_y(i,j) = 1;
        n_missing++;
      }
    }
  }
  
  arma::mat wy_mat(y_mat.n_rows,y_mat.n_cols,arma::fill::zeros);
  
  
  //INITIALIZE COMMON PARAMETERS AND OUTPUT ARRAYS
  //The difference is in the parameters of the splines/GP
  
  arma::vec x1 = unique(x_mat.col(0));
  arma::vec x2 = unique(x_mat.col(1));
  
  int n1 = x1.n_elem;
  int n2 = x2.n_elem;
  //Needed for the log_accept function
  int n_rep = y_mat.n_cols;
  unsigned int n3 = y_mat.n_rows;
  if(n3 != x_mat.n_rows){
    Rcout << "Error: y and x must have the same length" << "\n";
    return 0;
  }
  
  int n_var_par = Hyper_param.length(), var_prior;
  if(n_var_par == 20){//IG priors
    var_prior = 1;
  }
  if(n_var_par == 14){//HC priors
    var_prior = 2;
  }
  
  
  
  //Gamma_0
  double gamma_0 = 0;
  double gamma_0_accept = 0;
  double gamma_0_count = 0;
  double s_gamma_0 = 0.01;
  
  //Gamma_1
  double gamma_1 = 0;
  double gamma_1_accept = 0;
  double gamma_1_count = 0;
  double s_gamma_1 = 0.01;
  
  //Gamma_2
  double gamma_2 = 0;
  double gamma_2_accept = 0;
  double gamma_2_count = 0;
  double s_gamma_2 = 0.01;
  
  //Log-Logistic curve f_1
  //Slope
  double Slope_1 = 1;
  double Slope_1_accept = 0;
  double Slope_1_count = 0;
  double s_Slope_1 = 0.01;
  double a_Slope_1 = Hyper_param["a_Slope_1"];
  double b_Slope_1 = Hyper_param["b_Slope_1"];
  //EC50
  double Ec50_1 = 0;
  double Ec50_1_accept = 0;
  double Ec50_1_count = 0;
  double s_Ec50_1 = 0.01;
  //LowerAsymptote
  double La_1 = 0.5;
  double La_1_accept = 0;
  double La_1_count = 0;
  double s_La_1 = 0.01;
  double a_La_1 = Hyper_param["a_La_1"];
  double b_La_1 = Hyper_param["b_La_1"];
  

  
  arma::colvec f_1(n1,arma::fill::ones);

  
  //Log-Logistic curve f_2
  //Slope
  double Slope_2 = 1;
  double Slope_2_accept = 0;
  double Slope_2_count = 0;
  double s_Slope_2 = 0.01;
  double a_Slope_2 = Hyper_param["a_Slope_2"];
  double b_Slope_2 = Hyper_param["b_Slope_2"];
  //EC50
  double Ec50_2 = 0;
  double Ec50_2_accept = 0;
  double Ec50_2_count = 0;
  double s_Ec50_2 = 0.01;
  //LowerAsymptote
  double La_2 = 0.5;
  double La_2_accept = 0;
  double La_2_count = 0;
  double s_La_2 = 0.01;
  double a_La_2 = Hyper_param["a_La_2"];
  double b_La_2 = Hyper_param["b_La_2"];
  
  arma::colvec f_2(n2,arma::fill::ones);
  for(int j = 1; j < n2; j++){
    // f_2(j) = pow(1 + pow(10,Slope_2 * (x2(j) - Ec50_2)),-1);
    f_2(j) = La_2 + (1-La_2)*pow(1 + pow(10,Slope_2 * (x2(j) - Ec50_2)),-1);
  }
  
  //Compute baseline level
  arma::mat p0 = f_1 * f_2.t();
  
  //For variance updates
  //IG
  double a_s2_gamma_0 = 0, b_s2_gamma_0 = 0, a_s2_gamma_1 = 0, b_s2_gamma_1 = 0, a_s2_gamma_2 = 0, b_s2_gamma_2 = 0, a_s2_eps = 0, b_s2_eps = 0, a_s2_Ec50_1 = 0, b_s2_Ec50_1 = 0, a_s2_Ec50_2 = 0, b_s2_Ec50_2 = 0;
  //HC
  double s2_gamma_0_new, s2_gamma_1_new, s2_gamma_2_new, s2_eps_new, s2_Ec50_1_new, s2_Ec50_2_new;
  double h_s2_gamma_0 = 0, h_s2_gamma_1 = 0, h_s2_gamma_2 = 0, h_s2_eps = 0, h_s2_Ec50_1 = 0, h_s2_Ec50_2 = 0;
  double s_s2_gamma_0 = 0, s_s2_gamma_1 = 0, s_s2_gamma_2 = 0, s_s2_eps = 0, s_s2_Ec50_1 = 0, s_s2_Ec50_2 = 0;
  double s2_gamma_0_accept = 0, s2_gamma_1_accept = 0, s2_gamma_2_accept = 0, s2_eps_accept = 0, s2_Ec50_1_accept = 0, s2_Ec50_2_accept = 0;
  double s2_gamma_0_count = 0, s2_gamma_1_count = 0, s2_gamma_2_count = 0, s2_eps_count = 0, s2_Ec50_1_count = 0, s2_Ec50_2_count = 0;
  if(var_prior == 1){//Inverse-Gamma
    a_s2_gamma_0 = Hyper_param["a_s2_gamma_0"];
    b_s2_gamma_0 = Hyper_param["b_s2_gamma_0"];
    a_s2_gamma_1 = Hyper_param["a_s2_gamma_1"];
    b_s2_gamma_1 = Hyper_param["b_s2_gamma_1"];
    a_s2_gamma_2 = Hyper_param["a_s2_gamma_2"];
    b_s2_gamma_2 = Hyper_param["b_s2_gamma_2"];
    a_s2_eps = Hyper_param["a_s2_eps"];
    b_s2_eps = Hyper_param["b_s2_eps"];
    a_s2_Ec50_1 = Hyper_param["a_s2_Ec50_1"];
    b_s2_Ec50_1 = Hyper_param["b_s2_Ec50_1"];
    a_s2_Ec50_2 = Hyper_param["a_s2_Ec50_2"];
    b_s2_Ec50_2 = Hyper_param["b_s2_Ec50_2"];
  }else{//Half-Cauchy
    h_s2_gamma_0 = Hyper_param["h_s2_gamma_0"];
    h_s2_gamma_1 = Hyper_param["h_s2_gamma_1"];
    h_s2_gamma_2 = Hyper_param["h_s2_gamma_2"];
    h_s2_eps = Hyper_param["h_s2_eps"];
    h_s2_Ec50_1 = Hyper_param["h_s2_Ec50_1"];
    h_s2_Ec50_2 = Hyper_param["h_s2_Ec50_2"];
  }
  double s2_gamma_0 = 1;
  double s2_gamma_1 = 1;
  double s2_gamma_2 = 1;
  double s2_eps = 1;
  double s2_Ec50_1 = 1;
  double s2_Ec50_2 = 1;
  
  
  //Some objects needed outside if loop
  
  //Splines
  int K1 = 0, K2 = 0, K3 = 0;
  double s_d_C = 0, s_d_B = 0, C_accept = 0, B_accept = 0, aux_s2eps = 0, C_count = 0, B_count = 0;
  arma::vec vec_C, vec_B, vec_C_new, vec_B_new, sum_C, sum_B, b_new, mu_B(n3,arma::fill::zeros), mu_C;
  arma::mat C, B(n1,n2,arma::fill::zeros), C_new, s_C, D_C, U_C, sigma_aux_C, prod_C, eye_K1, eye_K2, eye_K3, B_K1(n1,K1,arma::fill::zeros), B_K2(n2,K2,arma::fill::zeros);
  arma::field<arma::mat> B_K3(n1,n2); //field object of matrices along dimensions n1 and n2 (drug doses)
  
  //GP
  bool update_ell = false, update_nu = false, update_alpha = false, update_sigma2_f = false;
  
  //x_distances for GP kernel calculations
  arma::mat x_dist(n3,n3,arma::fill::eye);
  for(int i = 0; i < n3; i++){
    for(int j = i; j < n3; j++){
      x_dist(i,j) = accu((x_mat.row(i) - x_mat.row(j)) % (x_mat.row(i) - x_mat.row(j)));
      x_dist(j,i) = x_dist(i,j);
    }
  }
  
  double a_ell = 0, b_ell = 0, ell = 0, ell_accept = 0, s_ell = 0, a_nu = 0, b_nu = 0, nu = 0, nu_accept = 0, s_nu = 0, a_alpha = 0, b_alpha = 0, alpha = 0, alpha_accept = 0, s_alpha = 0, logdet_Kxx = 0, sign = 1, a_sigma2_f = 0, b_sigma2_f = 0, h_sigma2_f = 0, sigma2_f = 0, s_sigma2_f = 0, sigma2_f_accept = 0, ell_count = 0, nu_count = 0, alpha_count = 0, sigma2_f_count = 0;
  int sigma2_f_prior = 0;
  arma::vec p0_vec(n3,arma::fill::zeros), z(n3), z_new(n3);
  arma::mat eye_n3(n3,n3,arma::fill::eye), eye_2(2,2,arma::fill::eye), Kxx(n3,n3), Kxx_new(n3,n3), s_B(n3,n3,arma::fill::zeros), prod_B(n3,n3,arma::fill::zeros), sigma_aux_B(n3,n3,arma::fill::zeros), sigma_aux_B_new(n3,n3);
  
  List model_spec_new;
  int type = model_spec["type"];
  bool update_LA = model_spec["lower_asymptote"];
  if(!update_LA){
    La_1 = 0;
    La_2 = 0;
  }
  if(type == 1){//Spline model
    
    //Spline of degree 3
    double deg = 3; 
    double ndx = 5;
    
    //Dim 1
    double xl1 = min(x1) - 0.1;
    double xr1 = max(x1) + 0.1;
    double dx1 = (xr1 - xl1) / ndx;
    // Construct a B-spline basis of degree 'deg'
    double start1 = xl1 - dx1 * deg;
    double end1 = xr1 + dx1 * deg;
    arma::vec knots_1 = arma::regspace(start1, dx1, end1);
    K1 = knots_1.n_elem;
    arma::mat P1(n1,K1,arma::fill::zeros), D1(K1,K1,arma::fill::zeros);
    // Truncated p-th power function
    for(int i1 = 0; i1 < n1; i1 ++){
      for(int j1 = 0; j1 < K1; j1 ++){
        if(x1(i1) > knots_1(j1)){
          P1(i1,j1) = pow(x1(i1) - knots_1(j1), deg);
        }
      }
    }
    arma::mat eye_K1_pre(K1,K1,arma::fill::eye);
    D1 = arma::diff(eye_K1_pre, deg + 1, 0) / (tgamma(deg + 1) * pow(dx1, deg));
    B_K1 = pow(-1, deg + 1) * P1 * D1.t();
    K1 = B_K1.n_cols;
    arma::mat eye_K1(K1,K1,arma::fill::eye);
    
    //Dim 2
    double xl2 = min(x2) - 0.1;
    double xr2 = max(x2) + 0.1;
    double dx2 = (xr2 - xl2) / ndx;
    // Construct a B-spline basis of degree 'deg'
    double start2 = xl2 - dx2 * deg;
    double end2 = xr2 + dx2 * deg;
    arma::vec knots_2 = arma::regspace(start2, dx2, end2);
    K2 = knots_2.n_elem;
    arma::mat P2(n2,K2,arma::fill::zeros), D2(K2,K2,arma::fill::zeros);
    // Truncated p-th power function
    for(int i2 = 0; i2 < n2; i2 ++){
      for(int j2 = 0; j2 < K2; j2 ++){
        if(x2(i2) > knots_2(j2)){
          P2(i2,j2) = pow(x2(i2) - knots_2(j2), deg);
        }
      }
    }
    arma::mat eye_K2_pre(K2,K2,arma::fill::eye);
    D2 = arma::diff(eye_K2_pre, deg + 1, 0) / (tgamma(deg + 1) * pow(dx2, deg));
    B_K2 = pow(-1, deg + 1) * P2 * D2.t();
    K2 = B_K2.n_cols;
    arma::mat eye_K2(K2,K2,arma::fill::eye);
    
    //Spline for interaction
    K3 = K1*K2;
    //This splines are product of the previous ones
    for(int j1 = 0; j1 < n1; j1 ++){
      for(int j2 = 0; j2 < n2; j2 ++){
        arma::mat aux_mat(K1,K2,arma::fill::zeros);
        for(int i1 = 0; i1 < K1; i1 ++){
          for(int i2 = 0; i2 < K2; i2 ++){
            aux_mat(i1,i2) = B_K1(j1,i1) * B_K2(j2,i2);
          }
        }
        B_K3(j1,j2) = aux_mat;
      }
    }
    //kronecker product matrices
    D_C = 2 * eye_K1;
    U_C = 2 * eye_K2;
    arma::mat eye_K3(K3,K3,arma::fill::eye);
    for(int i1 = 1; i1 < K1; i1 ++){
      D_C(i1, i1 - 1) = -1;
    }
    for(int i1 = 0; i1 < K1-1; i1 ++){
      D_C(i1, i1 + 1) = -1;
    }
    for(int i2 = 1; i2 < K2; i2 ++){
      U_C(i2, i2 - 1) = -1;
    }
    for(int i2 = 0; i2 < K2-1; i2 ++){
      U_C(i2, i2 + 1) = -1;
    }
    
    //Matrix of coefficients for the spline term
    //Differentiate between B and the matrix C of coefficients!
    arma::mat C(K1,K2,arma::fill::zeros);

    C_accept = 0;
    C_count = 0;
    s_C = 0.01*eye_K3;
    vec_C = vectorise(C,0);
    mu_C = vec_C;
    //Adaptive
    s_d_C = pow(2.4,2)/K3;
    sum_C = vec_C;
    prod_C = vec_C * vec_C.t();
    
    //We need to pass the spline already multiplied to avoid problems with array dimensions etc...
    for(int i = 0; i < n1; i ++){
      for(int j = 0; j < n2; j ++){
        B(i,j) = arma::accu(C % B_K3(i,j));
      }
    }
    
    //Needed for density calculations
    sigma_aux_C = arma::trans(arma::chol(arma::inv(arma::kron(U_C, D_C))));
    z = sigma_aux_C * (vec_C - mu_C);  
    
  }else{//GP model
    
    p0_vec = vectorise(p0,0);
    //Matrix used for the GP
    B_accept = 0;
    B_count = 0;
    vec_B = vectorise(B);
    s_B = 0.01*eye_n3;
    //Adaptive
    s_d_B = pow(2.4,2)/n3;
    sum_B = vec_B;
    prod_B = vec_B * vec_B.t();
    
    NumericVector ell_param = model_spec["ell"];
    if(ell_param.length() > 1){// prior on ell
      a_ell = ell_param(0);
      b_ell = ell_param(1);
      ell = 1.0;
      ell_accept = 0.0;
      ell_count = 0;
      s_ell = 0.01;
      update_ell = true;
      model_spec["ell"] = ell;
      model_spec_new = clone(model_spec);
    }else{
      ell = model_spec["ell"];
    }
    
    //GP's with additional parameters
    if(type == 3){//Matern kernel
      NumericVector nu_param = model_spec["nu"];
      if(nu_param.length() > 1){// random nu
        a_nu = nu_param(0);
        b_nu = nu_param(1);
        nu = 1.0;
        nu_accept = 0.0;
        nu_count = 0;
        s_nu = 0.01;
        update_nu = true;
        model_spec["nu"] = nu;
        model_spec_new = clone(model_spec);
      }else{
        nu = model_spec["nu"];
      }
    }
    
    if(type == 4){//Rational quadratic kernel
      NumericVector alpha_param = model_spec["alpha"];
      if(alpha_param.length() > 1){// random alpha
        a_alpha = alpha_param(0);
        b_alpha = alpha_param(1);
        alpha = 1.0;
        alpha_accept = 0.0;
        alpha_count = 0;
        s_alpha = 0.01;
        update_alpha = true;
        model_spec["alpha"] = alpha;
        model_spec_new = clone(model_spec);
      }else{
        alpha = model_spec["alpha"];
      }
      
    }
    
    Kxx = GPKernel(x_dist, model_spec);
    //sigma2_f
    sigma2_f_prior = model_spec["sigma2_f_prior"];
    if(sigma2_f_prior == 0){//sigma2_f is fixed
      sigma2_f = model_spec["sigma2_f"];
      update_sigma2_f = false;
    }
    if(sigma2_f_prior == 1){//sigma2_f is IG
      a_sigma2_f = model_spec["a_sigma2_f"];;
      b_sigma2_f = model_spec["b_sigma2_f"];;
      sigma2_f = 1.0;
      update_sigma2_f = true;
      model_spec["sigma2_f"] = sigma2_f;
    }
    if(sigma2_f_prior == 2){//sigma2_f is HC
      h_sigma2_f = model_spec["h_sigma2_f"];
      sigma2_f = 1.0;
      sigma2_f_accept = 0.0;
      sigma2_f_count = 0;
      s_sigma2_f = 0.01;
      update_sigma2_f = true;
      model_spec["sigma2_f"] = sigma2_f;
    }
    model_spec_new = clone(model_spec);
    
    
    //Needed for density calculations
    sigma_aux_B = arma::trans(arma::chol(arma::inv(Kxx)));
    arma::log_det(logdet_Kxx, sign, Kxx);
    z = sigma_aux_B * (vec_B - mu_B);    
  }
  
  //Interaction Delta before transformation
  List Delta_list = List::create(Named("gamma0") = gamma_0, Named("gamma1") = gamma_1, Named("gamma2") = gamma_2, Named("B") = B);
  arma::mat Delta = computeDelta(Delta_list, x1, x2);
  
  //Joint update of b
  arma::colvec b(2,arma::fill::ones);
  double a_b1 = Hyper_param["a_b1"];
  double b_b1 = Hyper_param["b_b1"];
  double a_b2 = Hyper_param["a_b2"];
  double b_b2 = Hyper_param["b_b2"];
  double b_accept = 0;
  int b_count = 0;
  eye_2.eye();
  arma::mat s_b = 0.01 * eye_2;
  //Adaptive
  double s_d_b = pow(2.4,2)/2;
  arma::vec sum_b = log(b);
  arma::mat prod_b = sum_b * sum_b.t();
  
  //Transformation for interaction matrix Delta
  arma::mat Delta_trans = - p0 % pow(1 + exp(b(0)*Delta),-1) + (1 - p0) % pow(1 + exp(-b(1)*Delta),-1);
  //Interaction only for doses > 0
  arma::mat id(n1,n2,arma::fill::ones);
  id.row(0) = arma::rowvec(n2,arma::fill::zeros);
  id.col(0) = arma::colvec(n1,arma::fill::zeros);
  Delta_trans = Delta_trans % id;
  
  arma::mat p_ij = p0 + Delta_trans;
  arma::vec p_ij_vec = vectorise(p_ij,0);
  
  //MCMC parameters
  double g0 = Alg_param["g0"], wg = Alg_param["wg"], opt_rate = Alg_param["opt_rate"], eps = Alg_param["eps"];
  double n_burn = Alg_param["n_burn"], thin = Alg_param["thin"];
  int n_save = Alg_param["n_save"], G = n_burn + thin * n_save, iter;
  
  if(g0 > n_burn){
    Rcout << "The adaptive burn-in (g0) should be smaller than the MCMC burn-in (n_burn)." << "\n";
    return 0;
  }
  
  //Output matrices
  arma::vec EC50_1(n_save,arma::fill::zeros), SLOPE_1(n_save,arma::fill::zeros), LA_1(n_save,arma::fill::zeros), EC50_2(n_save,arma::fill::zeros), SLOPE_2(n_save,arma::fill::zeros),LA_2(n_save,arma::fill::zeros);
  arma::vec GAMMA_0(n_save,arma::fill::zeros), GAMMA_1(n_save,arma::fill::zeros), GAMMA_2(n_save,arma::fill::zeros);
  
  if(type != 1){//GP does not need the matrix C_OUT, but we intialise it or we get an error
    K3 = 0;
  }
  
  arma::mat C_OUT(n_save,K3,arma::fill::zeros), B_OUT(n_save,n3,arma::fill::zeros);
  arma::vec ELL(n_save), NU(n_save), ALPHA(n_save), SIGMA2_F(n_save);
  
  arma::vec B1(n_save,arma::fill::zeros), B2(n_save,arma::fill::zeros);
  arma::vec S2_GAMMA_0(n_save,arma::fill::zeros), S2_GAMMA_1(n_save,arma::fill::zeros), S2_GAMMA_2(n_save,arma::fill::zeros), S2_EC50_1(n_save,arma::fill::zeros), S2_EC50_2(n_save,arma::fill::zeros), S2_EPS(n_save,arma::fill::zeros);
  
  //Performance evaluation (LPML)
  arma::cube CPO(n_rep,n1,n2,arma::fill::zeros);
  
  //Additional variables
  double log_accept, accept, Ec50_1_new, Ec50_2_new, Slope_1_new, Slope_2_new, La_1_new, La_2_new, gamma_0_new, gamma_1_new, gamma_2_new, ell_new, nu_new, alpha_new, logdet_Kxx_new, sigma2_f_new;
  arma::mat Delta_new(n1,n2), Delta_trans_new(n1,n2), p0_new(n1,n2), p_ij_new(n1,n2), B_new(n1,n2);
  arma::vec f_1_new(n1,arma::fill::ones), f_2_new(n2,arma::fill::ones), p_ij_vec_new(n3);
  List Delta_list_new;
  arma::rowvec y_rowtmp; 
  
  ///////////////////
  //  START GIBBS  //
  ///////////////////
  Progress progr(G, true);
  
  for(int g = 0; g < G; g++){
    
    // Start each iteration by filling in missing values
    // This is done by sampling iid N(0,\sigma^2) and adding in row-wise means p_ij
    // p_ij
    for (int i = 0; i < y_mat.n_rows; i++){
      y_rowtmp = y_mat.row(i);
      y_rowtmp.replace(arma::datum::nan,0);
      y_rowtmp = y_rowtmp+missing_y.row(i)%(p_ij_vec(i) + sqrt(s2_eps)*arma::randn(1,y_mat.n_cols));
      wy_mat.row(i) = y_rowtmp;
    }
    
    //Slope1
    //Propose new value
    Slope_1_new = Slope_1 * exp(R::rnorm(0,1) * sqrt(s_Slope_1));
    
    for(int i = 1; i < n1; i ++){
      f_1_new(i) = La_1 + (1-La_1)*pow(1 + pow(10,Slope_1_new * (x1(i) - Ec50_1)),-1);
    }
    p0_new = f_1_new * f_2.t();
    
    Delta_trans_new = - p0_new % pow(1 + exp(b(0)*Delta),-1) + (1 - p0_new) % pow(1 + exp(-b(1)*Delta),-1);
    Delta_trans_new = Delta_trans_new % id;
    
    p_ij_new = p0_new + Delta_trans_new;
    p_ij_vec_new = vectorise(p_ij_new,0);
    
    //Computing MH ratio:
    log_accept = a_Slope_1 *(log(Slope_1_new) - log(Slope_1)) - b_Slope_1 * (Slope_1_new - Slope_1);
    log_accept = log_accept + logaccept(wy_mat, p_ij_vec, p_ij_vec_new, s2_eps, 1);
    accept = 1;
    if( arma::is_finite(log_accept) ){
      if(log_accept < 0){
        accept = exp(log_accept);
      }
    }else{
      accept = 0;
    }
    
    Slope_1_accept = Slope_1_accept + accept;
    Slope_1_count = Slope_1_count + 1;
    
    if( R::runif(0,1) < accept ){
      Slope_1 = Slope_1_new;
      f_1 = f_1_new;
      p0 = p0_new;
      Delta_trans = Delta_trans_new;
      p_ij = p_ij_new;
      p_ij_vec = p_ij_vec_new;
    }
    
    s_Slope_1 = s_Slope_1 + pow(g+1,-wg)*(accept - opt_rate);
    if(s_Slope_1 > exp(50)){
      s_Slope_1 = exp(50);
    }else{
      if(s_Slope_1 < exp(-50)){
        s_Slope_1 = exp(-50);
      }
    }
    
    
    //Slope2
    //Propose new value
    Slope_2_new = Slope_2 * exp(R::rnorm(0,1) * sqrt(s_Slope_2));
    
    for(int j = 1; j < n2; j ++){
      f_2_new(j) = La_2 + (1-La_2)*pow(1 + pow(10,Slope_2_new * (x2(j) - Ec50_2)),-1);
    }
    p0_new = f_1 * f_2_new.t();
    
    Delta_trans_new = - p0_new % pow(1 + exp(b(0)*Delta),-1) + (1 - p0_new) % pow(1 + exp(-b(1)*Delta),-1);
    Delta_trans_new = Delta_trans_new % id;
    
    p_ij_new = p0_new + Delta_trans_new;
    p_ij_vec_new = vectorise(p_ij_new,0);
    
    //Computing MH ratio:
    log_accept = a_Slope_2 *(log(Slope_2_new) - log(Slope_2)) - b_Slope_2 * (Slope_2_new - Slope_2);
    log_accept = log_accept + logaccept(wy_mat, p_ij_vec, p_ij_vec_new, s2_eps, 1);
    
    accept = 1;
    if( arma::is_finite(log_accept) ){
      if(log_accept < 0){
        accept = exp(log_accept);
      }
    }else{
      accept = 0;
    }
    
    Slope_2_accept = Slope_2_accept + accept;
    Slope_2_count = Slope_2_count + 1;
    
    if( R::runif(0,1) < accept ){
      Slope_2 = Slope_2_new;
      f_2 = f_2_new;
      p0 = p0_new;
      Delta_trans = Delta_trans_new;
      p_ij = p_ij_new;
      p_ij_vec = p_ij_vec_new;
    }
    
    s_Slope_2 = s_Slope_2 + pow(g+1,-wg)*(accept - opt_rate);
    if(s_Slope_2 > exp(50)){
      s_Slope_2 = exp(50);
    }else{
      if(s_Slope_2 < exp(-50)){
        s_Slope_2 = exp(-50);
      }
    }
    
    
    //Ec501
    //Propose new value
    Ec50_1_new = Ec50_1 + R::rnorm(0,1) * sqrt(s_Ec50_1);
    
    for(int i = 1; i < n1; i ++){
      f_1_new(i) = La_1 + (1-La_1)*pow(1 + pow(10,Slope_1 * (x1(i) - Ec50_1_new)),-1);
    }
    p0_new = f_1_new * f_2.t();
    
    Delta_trans_new = - p0_new % pow(1 + exp(b(0)*Delta),-1) + (1 - p0_new) % pow(1 + exp(-b(1)*Delta),-1);
    Delta_trans_new = Delta_trans_new % id;
    
    p_ij_new = p0_new + Delta_trans_new;
    p_ij_vec_new = vectorise(p_ij_new,0);
    
    //Computing MH ratio:
    log_accept = - 0.5 * (pow(Ec50_1_new,2) - pow(Ec50_1,2))/s2_Ec50_1;
    log_accept = log_accept + logaccept(wy_mat, p_ij_vec, p_ij_vec_new, s2_eps, 1);
    
    accept = 1;
    if( arma::is_finite(log_accept) ){
      if(log_accept < 0){
        accept = exp(log_accept);
      }
    }else{
      accept = 0;
    }
    
    Ec50_1_accept = Ec50_1_accept + accept;
    Ec50_1_count = Ec50_1_count + 1;
    
    if( R::runif(0,1) < accept ){
      Ec50_1 = Ec50_1_new;
      f_1 = f_1_new;
      p0 = p0_new;
      Delta_trans = Delta_trans_new;
      p_ij = p_ij_new;
      p_ij_vec = p_ij_vec_new;
    }
    
    s_Ec50_1 = s_Ec50_1 + pow(g+1,-wg)*(accept - opt_rate);
    if(s_Ec50_1 > exp(50)){
      s_Ec50_1 = exp(50);
    }else{
      if(s_Ec50_1 < exp(-50)){
        s_Ec50_1 = exp(-50);
      }
    }
    
    
    //Ec502
    //Propose new value
    Ec50_2_new = Ec50_2 + R::rnorm(0,1) * sqrt(s_Ec50_2);
    
    for(int j = 1; j < n2; j ++){
      f_2_new(j) = La_2 + (1-La_2)*pow(1 + pow(10,Slope_2 * (x2(j) - Ec50_2_new)),-1);
    }
    p0_new = f_1 * f_2_new.t();
    
    Delta_trans_new = - p0_new % pow(1 + exp(b(0)*Delta),-1) + (1 - p0_new) % pow(1 + exp(-b(1)*Delta),-1);
    Delta_trans_new = Delta_trans_new % id;
    
    p_ij_new = p0_new + Delta_trans_new;
    p_ij_vec_new = vectorise(p_ij_new,0);
    
    //Computing MH ratio:
    log_accept = - 0.5 * (pow(Ec50_2_new,2) - pow(Ec50_2,2))/s2_Ec50_2;
    log_accept = log_accept + logaccept(wy_mat, p_ij_vec, p_ij_vec_new, s2_eps, 1);
    
    accept = 1;
    if( arma::is_finite(log_accept) ){
      if(log_accept < 0){
        accept = exp(log_accept);
      }
    }else{
      accept = 0;
    }
    
    Ec50_2_accept = Ec50_2_accept + accept;
    Ec50_2_count = Ec50_2_count + 1;
    
    if( R::runif(0,1) < accept ){
      Ec50_2 = Ec50_2_new;
      f_2 = f_2_new;
      p0 = p0_new;
      Delta_trans = Delta_trans_new;
      p_ij = p_ij_new;
      p_ij_vec = p_ij_vec_new;
    }
    
    s_Ec50_2 = s_Ec50_2 + pow(g+1,-wg)*(accept - opt_rate);
    if(s_Ec50_2 > exp(50)){
      s_Ec50_2 = exp(50);
    }else{
      if(s_Ec50_2 < exp(-50)){
        s_Ec50_2 = exp(-50);
      }
    }
    
    if (update_LA){
      //La_1
      //Propose new value
      La_1_new = pow(1+exp(-(La_1+R::rnorm(0,1)*sqrt(s_La_1))),-1);
      
      for(int i = 1; i < n1; i++){
        f_1_new(i) = La_1_new + (1-La_1_new)*pow(1 + pow(10,Slope_1 * (x1(i) - Ec50_1)),-1);
      }
      p0_new = f_1_new * f_2.t();
      
      Delta_trans_new = - p0_new % pow(1 + exp(b(0)*Delta),-1) + (1 - p0_new) % pow(1 + exp(-b(1)*Delta),-1);
      Delta_trans_new = Delta_trans_new % id;
      
      p_ij_new = p0_new + Delta_trans_new;
      p_ij_vec_new = vectorise(p_ij_new,0);
      
      //Computing MH ratio:
      log_accept = a_La_1*(log(La_1_new)-log(La_1))+b_La_1*(log(1-La_1_new)-log(1-La_1));
      log_accept = log_accept + logaccept(wy_mat, p_ij_vec, p_ij_vec_new, s2_eps, 1);
      
      accept = 1;
      if( arma::is_finite(log_accept) ){
        if(log_accept < 0){
          accept = exp(log_accept);
        }
      }else{
        accept = 0;
      }
      
      La_1_accept = La_1_accept + accept;
      La_1_count = La_1_count+1;
      if( R::runif(0,1) < accept ){
        La_1 = La_1_new;
        f_1 = f_1_new;
        p0 = p0_new;
        Delta_trans = Delta_trans_new;
        p_ij = p_ij_new;
        p_ij_vec = p_ij_vec_new;
      }

      s_La_1 = s_La_1 + pow(g+1,-wg)*(accept - opt_rate);
      if(s_La_1 > exp(50)){
        s_La_1 = exp(50);
      }else{
        if(s_La_1 < exp(-50)){
          s_La_1 = exp(-50);
        }
      }
      
      //La_2
      //Propose new value
      La_2_new = pow(1+exp(-(La_2+R::rnorm(0,1)*sqrt(s_La_2))),-1);
      for(int j = 1; j < n2; j++){
        // f_2(j) = pow(1 + pow(10,Slope_2 * (x2(j) - Ec50_2)),-1);
        f_2_new(j) = La_2_new + (1-La_2_new)*pow(1 + pow(10,Slope_2 * (x2(j) - Ec50_2)),-1);
      }
      p0_new = f_1 * f_2_new.t();
      
      Delta_trans_new = - p0_new % pow(1 + exp(b(0)*Delta),-1) + (1 - p0_new) % pow(1 + exp(-b(1)*Delta),-1);
      Delta_trans_new = Delta_trans_new % id;
      
      p_ij_new = p0_new + Delta_trans_new;
      p_ij_vec_new = vectorise(p_ij_new,0);
      
      //Computing MH ratio:
      // Note! The prior becomes something else, due to transform
      log_accept = a_La_2*(log(La_2_new)-log(La_2))+b_La_2*(log(1-La_2_new)-log(1-La_2));
      log_accept = log_accept + logaccept(wy_mat, p_ij_vec, p_ij_vec_new, s2_eps, 1);
      
      accept = 1;
      if( arma::is_finite(log_accept) ){
        if(log_accept < 0){
          accept = exp(log_accept);
        }
      }else{
        accept = 0;
      }
      
      La_2_accept = La_2_accept + accept;
      La_2_count = La_2_count+1;
      if( R::runif(0,1) < accept ){
        La_2 = La_2_new;
        f_2 = f_2_new;
        p0 = p0_new;
        Delta_trans = Delta_trans_new;
        p_ij = p_ij_new;
        p_ij_vec = p_ij_vec_new;
      }
      
      s_La_2 = s_La_2 + pow(g+1,-wg)*(accept - opt_rate);
      if(s_La_2 > exp(50)){
        s_La_2 = exp(50);
      }else{
        if(s_La_2 < exp(-50)){
          s_La_2 = exp(-50);
        }
      }
    }
    
    
    
    //gamma0
    //Propose new value
    gamma_0_new = gamma_0 + R::rnorm(0,1) * sqrt(s_gamma_0);
    
    Delta_list_new = clone(Delta_list);
    Delta_list_new["gamma0"] = gamma_0_new;
    
    Delta_new = computeDelta(Delta_list_new, x1, x2);
    
    Delta_trans_new = - p0 % pow(1 + exp(b(0)*Delta_new),-1) + (1 - p0) % pow(1 + exp(-b(1)*Delta_new),-1);
    Delta_trans_new = Delta_trans_new % id;
    
    p_ij_new = p0 + Delta_trans_new;
    p_ij_vec_new = vectorise(p_ij_new,0);
    
    //Computing MH ratio:
    log_accept = - 0.5 * (pow(gamma_0_new,2) - pow(gamma_0,2))/s2_gamma_0;
    log_accept = log_accept + logaccept(wy_mat, p_ij_vec, p_ij_vec_new, s2_eps, 1);
    
    accept = 1;
    if( arma::is_finite(log_accept) ){
      if(log_accept < 0){
        accept = exp(log_accept);
      }
    }else{
      accept = 0;
    }
    
    gamma_0_accept = gamma_0_accept + accept;
    gamma_0_count = gamma_0_count + 1;
    
    if( R::runif(0,1) < accept ){
      gamma_0 = gamma_0_new;
      Delta = Delta_new;
      Delta_trans = Delta_trans_new;
      Delta_list = clone(Delta_list_new);
      p_ij = p_ij_new;
      p_ij_vec = p_ij_vec_new;
    }
    
    s_gamma_0 = s_gamma_0 + pow(g+1,-wg)*(accept - opt_rate);
    if(s_gamma_0 > exp(50)){
      s_gamma_0 = exp(50);
    }else{
      if(s_gamma_0 < exp(-50)){
        s_gamma_0 = exp(-50);
      }
    }
    
    
    //gamma1
    //Propose new value
    gamma_1_new = gamma_1 + R::rnorm(0,1) * sqrt(s_gamma_1);
    
    Delta_list_new = clone(Delta_list);
    Delta_list_new["gamma1"] = gamma_1_new;
    
    Delta_new = computeDelta(Delta_list_new, x1, x2);
    
    Delta_trans_new = - p0 % pow(1 + exp(b(0)*Delta_new),-1) + (1 - p0) % pow(1 + exp(-b(1)*Delta_new),-1);
    Delta_trans_new = Delta_trans_new % id;
    
    p_ij_new = p0 + Delta_trans_new;
    p_ij_vec_new = vectorise(p_ij_new,0);
    
    //Computing MH ratio:
    log_accept = - 0.5 * (pow(gamma_1_new,2) - pow(gamma_1,2))/s2_gamma_1;
    log_accept = log_accept + logaccept(wy_mat, p_ij_vec, p_ij_vec_new, s2_eps, 1);
    
    accept = 1;
    if( arma::is_finite(log_accept) ){
      if(log_accept < 0){
        accept = exp(log_accept);
      }
    }else{
      accept = 0;
    }
    
    gamma_1_accept = gamma_1_accept + accept;
    gamma_1_count = gamma_1_count + 1;
    
    if( R::runif(0,1) < accept ){
      gamma_1 = gamma_1_new;
      Delta = Delta_new;
      Delta_trans = Delta_trans_new;
      Delta_list = clone(Delta_list_new);
      p_ij = p_ij_new;
      p_ij_vec = p_ij_vec_new;
    }
    
    s_gamma_1 = s_gamma_1 + pow(g+1,-wg)*(accept - opt_rate);
    if(s_gamma_1 > exp(50)){
      s_gamma_1 = exp(50);
    }else{
      if(s_gamma_1 < exp(-50)){
        s_gamma_1 = exp(-50);
      }
    }
    
    
    //gamma2
    //Propose new value
    gamma_2_new = gamma_2 + R::rnorm(0,1) * sqrt(s_gamma_2);
    
    Delta_list_new = clone(Delta_list);
    Delta_list_new["gamma2"] = gamma_2_new;
    
    Delta_new = computeDelta(Delta_list_new, x1, x2);
    
    Delta_trans_new = - p0 % pow(1 + exp(b(0)*Delta_new),-1) + (1 - p0) % pow(1 + exp(-b(1)*Delta_new),-1);
    Delta_trans_new = Delta_trans_new % id;
    
    p_ij_new = p0 + Delta_trans_new;
    p_ij_vec_new = vectorise(p_ij_new,0);
    
    //Computing MH ratio:
    log_accept = - 0.5 * (pow(gamma_2_new,2) - pow(gamma_2,2))/s2_gamma_2;
    log_accept = log_accept + logaccept(wy_mat, p_ij_vec, p_ij_vec_new, s2_eps, 1);
    
    accept = 1;
    if( arma::is_finite(log_accept) ){
      if(log_accept < 0){
        accept = exp(log_accept);
      }
    }else{
      accept = 0;
    }
    
    gamma_2_accept = gamma_2_accept + accept;
    gamma_2_count = gamma_2_count + 1;
    
    if( R::runif(0,1) < accept ){
      gamma_2 = gamma_2_new;
      Delta = Delta_new;
      Delta_trans = Delta_trans_new;
      Delta_list = clone(Delta_list_new);
      p_ij = p_ij_new;
      p_ij_vec = p_ij_vec_new;
    }
    
    s_gamma_2 = s_gamma_2 + pow(g+1,-wg)*(accept - opt_rate);
    if(s_gamma_2 > exp(50)){
      s_gamma_2 = exp(50);
    }else{
      if(s_gamma_2 < exp(-50)){
        s_gamma_2 = exp(-50);
      }
    }
    
    //Update matrix of coefficients (C or B, depending on the model)
    if(type == 1){//Splines
      
      //Propose new value (matrix-normal)
      vec_C_new = vec_C + arma::trans(arma::chol(s_C)) * arma::randn(K3,1);
      C_new = reshape(vec_C_new,K1,K2);
      
      for(int i = 0; i < n1; i ++){
        for(int j = 0; j < n2; j ++){
          B_new(i,j) = arma::accu(C_new % B_K3(i,j));
        }
      }
      
      Delta_list_new = clone(Delta_list);
      Delta_list_new["B"] = B_new;
      Delta_new = computeDelta(Delta_list_new, x1, x2);
      
      Delta_trans_new = - p0 % pow(1 + exp(b(0)*Delta_new),-1) + (1 - p0) % pow(1 + exp(-b(1)*Delta_new),-1);
      Delta_trans_new = Delta_trans_new % id;
      
      p_ij_new = p0 + Delta_trans_new;
      p_ij_vec_new = vectorise(p_ij_new,0);
      
      //Computing MH ratio:
      z_new = sigma_aux_C * (vec_C_new - mu_C);    
      log_accept = - 0.5 * (accu(z_new % z_new) - accu(z % z));     
      log_accept = log_accept + logaccept(wy_mat, p_ij_vec, p_ij_vec_new, s2_eps, 1);
      
      accept = 1;
      if( arma::is_finite(log_accept) ){
        if(log_accept < 0){
          accept = exp(log_accept);
        }
      }else{
        accept = 0;
      }
      
      C_accept = C_accept + accept;
      C_count = C_count + 1;
      
      if( R::runif(0,1) < accept ){
        C = C_new;
        vec_C = vec_C_new;
        B = B_new;
        z = z_new;
        Delta = Delta_new;
        Delta_trans = Delta_trans_new;
        Delta_list = clone(Delta_list_new);
        p_ij = p_ij_new;
        p_ij_vec = p_ij_vec_new;
      }
      
      sum_C = sum_C + vec_C;
      prod_C = prod_C + vec_C * vec_C.t();
      
      s_d_C = s_d_C + pow(g+1,-wg)*(accept - opt_rate);
      if(s_d_C > exp(50)){
        s_d_C = exp(50);
      }else{
        if(s_d_C < exp(-50)){
          s_d_C = exp(-50);
        }
      }
      if(g > g0){
        arma::mat eye_K3_adapt(K3,K3,arma::fill::eye);
        s_C = s_d_C/(g-1) * (prod_C - sum_C * sum_C.t()/g) + s_d_C * eps * eye_K3_adapt;
      }
    }else{//GP
      //Propose new value (mv-normal)
      vec_B_new = vec_B + arma::trans(arma::chol(s_B)) * arma::randn(n3,1);
      B_new = reshape(vec_B_new,n1,n2);
      
      Delta_list_new = clone(Delta_list);
      Delta_list_new["B"] = B_new;
      Delta_new = computeDelta(Delta_list_new, x1, x2);
      
      Delta_trans_new = - p0 % pow(1 + exp(b(0)*Delta_new),-1) + (1 - p0) % pow(1 + exp(-b(1)*Delta_new),-1);
      Delta_trans_new = Delta_trans_new % id;
      
      p_ij_new = p0 + Delta_trans_new;
      p_ij_vec_new = vectorise(p_ij_new,0);
      
      //Computing MH ratio:
      z_new = sigma_aux_B * (vec_B_new - mu_B);    
      log_accept = - 0.5 * (accu(z_new % z_new) - accu(z % z))/sigma2_f;     
      log_accept = log_accept + logaccept(wy_mat, p_ij_vec, p_ij_vec_new, s2_eps, 1);
      
      accept = 1;
      if( arma::is_finite(log_accept) ){
        if(log_accept < 0){
          accept = exp(log_accept);
        }
      }else{
        accept = 0;
      }
      
      B_accept = B_accept + accept;
      B_count = B_count + 1;
      
      if( R::runif(0,1) < accept ){
        vec_B = vec_B_new;
        B = B_new;
        z = z_new;
        Delta = Delta_new;
        Delta_trans = Delta_trans_new;
        Delta_list = clone(Delta_list_new);
        p_ij = p_ij_new;
        p_ij_vec = p_ij_vec_new;
      }
      
      sum_B = sum_B + vec_B;
      prod_B = prod_B + vec_B * vec_B.t();
      
      s_d_B = s_d_B + pow(g+1,-wg)*(accept - opt_rate);
      if(s_d_B > exp(50)){
        s_d_B = exp(50);
      }else{
        if(s_d_B < exp(-50)){
          s_d_B = exp(-50);
        }
      }
      if(g > g0){
        s_B = s_d_B/(g-1) * (prod_B - sum_B * sum_B.t()/g) + s_d_B * eps * eye_n3;
      }
      
      
      if(update_ell){
        //Propose new value
        ell_new = ell * exp(R::rnorm(0,1) * sqrt(s_ell));
        
        model_spec_new = clone(model_spec);
        model_spec_new["ell"] = ell_new;
        Kxx_new = GPKernel(x_dist, model_spec_new);
        arma::log_det( logdet_Kxx_new, sign, Kxx_new );
        sigma_aux_B_new = arma::trans(arma::chol(arma::inv(Kxx_new)));
        
        //Computing MH ratio:
        log_accept = a_ell *(log(ell_new) - log(ell)) - b_ell * (ell_new - ell);
        z_new = sigma_aux_B_new * (vec_B - mu_B);    
        log_accept = log_accept - 0.5 * ( (accu(z_new % z_new) - accu(z % z))/sigma2_f + logdet_Kxx_new - logdet_Kxx);     
        
        accept = 1;
        if( arma::is_finite(log_accept) ){
          if(log_accept < 0){
            accept = exp(log_accept);
          }
        }else{
          accept = 0;
        }
        
        ell_accept = ell_accept + accept;
        ell_count = ell_count + 1;
        
        if( R::runif(0,1) < accept ){
          ell = ell_new;
          Kxx = Kxx_new;
          logdet_Kxx = logdet_Kxx_new;
          sigma_aux_B = sigma_aux_B_new;
          z = z_new;
          model_spec = clone(model_spec_new);
        }
        
        s_ell = s_ell + pow(g+1,-wg)*(accept - opt_rate);
        if(s_ell > exp(50)){
          s_ell = exp(50);
        }else{
          if(s_ell < exp(-50)){
            s_ell = exp(-50);
          }
        }
      }
      
      
      if(update_nu){
        //Propose new value
        nu_new = nu * exp(R::rnorm(0,1) * sqrt(s_nu));
        
        model_spec_new = clone(model_spec);
        model_spec_new["nu"] = nu_new;
        Kxx_new = GPKernel(x_dist, model_spec_new);
        arma::log_det( logdet_Kxx_new, sign, Kxx_new );
        sigma_aux_B_new = arma::trans(arma::chol(arma::inv(Kxx_new)));
        
        //Computing MH ratio:
        log_accept = a_nu *(log(nu_new) - log(nu)) - b_nu * (nu_new - nu);
        z_new = sigma_aux_B_new * (vec_B - mu_B);    
        log_accept = log_accept - 0.5 * ( (accu(z_new % z_new) - accu(z % z))/sigma2_f + logdet_Kxx_new - logdet_Kxx);     
        
        accept = 1;
        if( arma::is_finite(log_accept) ){
          if(log_accept < 0){
            accept = exp(log_accept);
          }
        }else{
          accept = 0;
        }
        
        nu_accept = nu_accept + accept;
        nu_count = nu_count + 1;
        
        if( R::runif(0,1) < accept ){
          nu = nu_new;
          Kxx = Kxx_new;
          logdet_Kxx = logdet_Kxx_new;
          sigma_aux_B = sigma_aux_B_new;
          z = z_new;
          model_spec = clone(model_spec_new);
        }
        
        s_nu = s_nu + pow(g+1,-wg)*(accept - opt_rate);
        if(s_nu > exp(50)){
          s_nu = exp(50);
        }else{
          if(s_nu < exp(-50)){
            s_nu = exp(-50);
          }
        }
      }
      
      if(update_alpha){
        //Propose new value
        alpha_new = alpha * exp(R::rnorm(0,1) * sqrt(s_alpha));
        
        model_spec_new = clone(model_spec);
        model_spec_new["alpha"] = alpha_new;
        Kxx_new = GPKernel(x_dist, model_spec_new);
        arma::log_det( logdet_Kxx_new, sign, Kxx_new );
        sigma_aux_B_new = arma::trans(arma::chol(arma::inv(Kxx_new)));
        
        //Computing MH ratio:
        log_accept = a_alpha *(log(alpha_new) - log(alpha)) - b_alpha * (alpha_new - alpha);
        z_new = sigma_aux_B_new * (vec_B - mu_B);    
        log_accept = log_accept - 0.5 * ( (accu(z_new % z_new) - accu(z % z))/sigma2_f + logdet_Kxx_new - logdet_Kxx);     
        
        accept = 1;
        if( arma::is_finite(log_accept) ){
          if(log_accept < 0){
            accept = exp(log_accept);
          }
        }else{
          accept = 0;
        }
        
        alpha_accept = alpha_accept + accept;
        alpha_count = alpha_count + 1;
        
        if( R::runif(0,1) < accept ){
          alpha = alpha_new;
          Kxx = Kxx_new;
          logdet_Kxx = logdet_Kxx_new;
          sigma_aux_B = sigma_aux_B_new;
          z = z_new;
          model_spec = clone(model_spec_new);
        }
        
        s_alpha = s_alpha + pow(g+1,-wg)*(accept - opt_rate);
        if(s_alpha > exp(50)){
          s_alpha = exp(50);
        }else{
          if(s_alpha < exp(-50)){
            s_alpha = exp(-50);
          }
        }
      }
      
      if(update_sigma2_f){
        if(sigma2_f_prior == 1){//IG(3,2)
          sigma2_f = 1/R::rgamma(a_sigma2_f + n3/2, 1/(b_sigma2_f + 0.5 * accu(z % z) ));
        }else{//HC
          //Propose new value
          sigma2_f_new = sqrt(sigma2_f) * exp(R::rnorm(0,1) * sqrt(s_sigma2_f));
          sigma2_f_new = pow(sigma2_f_new,2);
          
          //Computing MH ratio:
          log_accept = log(sigma2_f + pow(h_sigma2_f,2)) - log(sigma2_f_new + pow(h_sigma2_f,2));
          log_accept = log_accept + (1 - n3)/2 * (log(sigma2_f_new) - log(sigma2_f)) - 0.5 * accu(z % z) * (1/sigma2_f_new - 1/sigma2_f);   
          
          accept = 1;
          if( arma::is_finite(log_accept) ){
            if(log_accept < 0){
              accept = exp(log_accept);
            }
          }else{
            accept = 0;
          }
          
          sigma2_f_accept = sigma2_f_accept + accept;
          sigma2_f_count = sigma2_f_count + 1;
          
          if( R::runif(0,1) < accept ){
            sigma2_f = sigma2_f_new;
          }
          
          s_sigma2_f = s_sigma2_f + pow(g+1,-wg)*(accept - opt_rate);
          if(s_sigma2_f > exp(50)){
            s_sigma2_f = exp(50);
          }else{
            if(s_sigma2_f < exp(-50)){
              s_sigma2_f = exp(-50);
            }
          }
        }
        model_spec["sigma2_f"] = sigma2_f;
      }
    }
    
    
    //Update b
    //Propose new values (log-normal)
    b_new = exp(log(b) + arma::chol(s_b,"lower") * arma::randn(2,1));
    // Rcout << "b_new: " << b_new.t() << "\n";
  
    Delta_trans_new = - p0 % pow(1 + exp(b_new(0)*Delta),-1) + (1 - p0) % pow(1 + exp(-b_new(1)*Delta),-1);
    Delta_trans_new = Delta_trans_new % id;
    
    p_ij_new = p0 + Delta_trans_new;
    p_ij_vec_new = vectorise(p_ij_new,0);
    
    //Computing MH ratio:
    log_accept = a_b1 * (log(b_new(0)) - log(b(0))) - b_b1 * (b_new(0) - b(0)) + a_b2 * (log(b_new(1)) - log(b(1))) - b_b2 * (b_new(1) - b(1));
    log_accept = log_accept + logaccept(wy_mat, p_ij_vec, p_ij_vec_new, s2_eps, 1);
    
    accept = 1;
    if( arma::is_finite(log_accept) ){
      if(log_accept < 0){
        accept = exp(log_accept);
      }
    }else{
      accept = 0;
    }
    b_accept = b_accept + accept;
    b_count = b_count + 1;
    
    if( R::runif(0,1) < accept ){
      // Rcout << "accept!" <<"\n";
      b = b_new;
      Delta_trans = Delta_trans_new;
      p_ij = p_ij_new;
      p_ij_vec = p_ij_vec_new;
    }
    // LR 10.04.20 added logs here, since covariance estimation should happen with transformed parameters
    sum_b = sum_b + log(b);
    prod_b = prod_b + log(b) * log(b).t();

    
    s_d_b = s_d_b + pow(g+1,-wg)*(accept - opt_rate);
    if(s_d_b > exp(50)){
      s_d_b = exp(50);
    }else{
      if(s_d_b < exp(-50)){
        s_d_b = exp(-50);
      }
    }
    if(g > g0){
      s_b = s_d_b/(g-1) * (prod_b - sum_b * sum_b.t()/g) + s_d_b * eps * eye_2;
      // Rcout << "s_b:" << s_b << "\n";
    } 

    //Update variances
    aux_s2eps =  -2 * logaccept(wy_mat, p_ij_vec, p_ij_vec, 1, 0);
    
    if(var_prior == 1){//IG update
      
      s2_eps = 1/R::rgamma(a_s2_eps + (n_rep * n1 * n2)/2, 1/(b_s2_eps + 0.5 *aux_s2eps ));
      
      s2_gamma_0 = 1/R::rgamma(a_s2_gamma_0 + 0.5 , 1/(b_s2_gamma_0 + 0.5 * pow(gamma_0,2) ));
      
      s2_gamma_1 = 1/R::rgamma(a_s2_gamma_1 + 0.5 , 1/(b_s2_gamma_1 + 0.5 * pow(gamma_1,2) ));
      
      s2_gamma_2 = 1/R::rgamma(a_s2_gamma_2 + 0.5 , 1/(b_s2_gamma_2 + 0.5 * pow(gamma_2,2) ));
      
      s2_Ec50_1 = 1/R::rgamma(a_s2_Ec50_1 + 0.5 , 1/(b_s2_Ec50_1 + 0.5 * pow(Ec50_1,2) ));
      
      s2_Ec50_2 = 1/R::rgamma(a_s2_Ec50_2 + 0.5 , 1/(b_s2_Ec50_2 + 0.5 * pow(Ec50_2,2) ));
      
    }else{//HC update
      //s2_eps
      //Propose new value
      s2_eps_new = sqrt(s2_eps) * exp(R::rnorm(0,1) * sqrt(s_s2_eps));
      s2_eps_new = pow(s2_eps_new,2);
      
      //Computing MH ratio:
      log_accept = log(s2_eps + pow(h_s2_eps,2)) - log(s2_eps_new + pow(h_s2_eps,2));
      log_accept = log_accept + (1 - (n_rep * n1 * n2))/2 * (log(s2_eps_new) - log(s2_eps)) - 0.5 * aux_s2eps * (1/s2_eps_new - 1/s2_eps);    
      
      accept = 1;
      if( arma::is_finite(log_accept) ){
        if(log_accept < 0){
          accept = exp(log_accept);
        }
      }else{
        accept = 0;
      }
      
      s2_eps_accept = s2_eps_accept + accept;
      s2_eps_count = s2_eps_count + 1;
      
      if( R::runif(0,1) < accept ){
        s2_eps = s2_eps_new;
      }
      
      s_s2_eps = s_s2_eps + pow(g+1,-wg)*(accept - opt_rate);
      if(s_s2_eps > exp(50)){
        s_s2_eps = exp(50);
      }else{
        if(s_s2_eps < exp(-50)){
          s_s2_eps = exp(-50);
        }
      }
      
      
      
      
      //s2_gamma_0
      //Propose new value
      s2_gamma_0_new = sqrt(s2_gamma_0) * exp(R::rnorm(0,1) * sqrt(s_s2_gamma_0));
      s2_gamma_0_new = pow(s2_gamma_0_new,2);
      
      //Computing MH ratio:
      log_accept = log(s2_gamma_0 + pow(h_s2_gamma_0,2)) - log(s2_gamma_0_new + pow(h_s2_gamma_0,2));
      log_accept = log_accept - 0.5 * pow(gamma_0,2) * (1/s2_gamma_0_new - 1/s2_gamma_0);    
      
      accept = 1;
      if( arma::is_finite(log_accept) ){
        if(log_accept < 0){
          accept = exp(log_accept);
        }
      }else{
        accept = 0;
      }
      
      s2_gamma_0_accept = s2_gamma_0_accept + accept;
      s2_gamma_0_count = s2_gamma_0_count + 1;
      
      if( R::runif(0,1) < accept ){
        s2_gamma_0 = s2_gamma_0_new;
      }
      
      s_s2_gamma_0 = s_s2_gamma_0 + pow(g+1,-wg)*(accept - opt_rate);
      if(s_s2_gamma_0 > exp(50)){
        s_s2_gamma_0 = exp(50);
      }else{
        if(s_s2_gamma_0 < exp(-50)){
          s_s2_gamma_0 = exp(-50);
        }
      }
      
      
      
      
      //s2_gamma_1
      //Propose new value
      s2_gamma_1_new = sqrt(s2_gamma_1) * exp(R::rnorm(0,1) * sqrt(s_s2_gamma_1));
      s2_gamma_1_new = pow(s2_gamma_1_new,2);
      
      //Computing MH ratio:
      log_accept = log(s2_gamma_1 + pow(h_s2_gamma_1,2)) - log(s2_gamma_1_new + pow(h_s2_gamma_1,2));
      log_accept = log_accept - 0.5 * pow(gamma_1,2) * (1/s2_gamma_1_new - 1/s2_gamma_1);    
      
      accept = 1;
      if( arma::is_finite(log_accept) ){
        if(log_accept < 0){
          accept = exp(log_accept);
        }
      }else{
        accept = 0;
      }
      
      s2_gamma_1_accept = s2_gamma_1_accept + accept;
      s2_gamma_1_count = s2_gamma_1_count + 1;
      
      if( R::runif(0,1) < accept ){
        s2_gamma_1 = s2_gamma_1_new;
      }
      
      s_s2_gamma_1 = s_s2_gamma_1 + pow(g+1,-wg)*(accept - opt_rate);
      if(s_s2_gamma_1 > exp(50)){
        s_s2_gamma_1 = exp(50);
      }else{
        if(s_s2_gamma_1 < exp(-50)){
          s_s2_gamma_1 = exp(-50);
        }
      }
      
      
      
      
      //s2_gamma_2
      //Propose new value
      s2_gamma_2_new = sqrt(s2_gamma_2) * exp(R::rnorm(0,1) * sqrt(s_s2_gamma_2));
      s2_gamma_2_new = pow(s2_gamma_2_new,2);
      
      //Computing MH ratio:
      log_accept = log(s2_gamma_2 + pow(h_s2_gamma_2,2)) - log(s2_gamma_2_new + pow(h_s2_gamma_2,2));
      log_accept = log_accept - 0.5 * pow(gamma_2,2) * (1/s2_gamma_2_new - 1/s2_gamma_2);      
      
      accept = 1;
      if( arma::is_finite(log_accept) ){
        if(log_accept < 0){
          accept = exp(log_accept);
        }
      }else{
        accept = 0;
      }
      
      s2_gamma_2_accept = s2_gamma_2_accept + accept;
      s2_gamma_2_count = s2_gamma_2_count + 1;
      
      if( R::runif(0,1) < accept ){
        s2_gamma_2 = s2_gamma_2_new;
      }
      
      s_s2_gamma_2 = s_s2_gamma_2 + pow(g+1,-wg)*(accept - opt_rate);
      if(s_s2_gamma_2 > exp(50)){
        s_s2_gamma_2 = exp(50);
      }else{
        if(s_s2_gamma_2 < exp(-50)){
          s_s2_gamma_2 = exp(-50);
        }
      }
      
      
      
      
      //s2_Ec50_1
      //Propose new value
      s2_Ec50_1_new = sqrt(s2_Ec50_1) * exp(R::rnorm(0,1) * sqrt(s_s2_Ec50_1));
      s2_Ec50_1_new = pow(s2_Ec50_1_new,2);
      
      //Computing MH ratio:
      log_accept = log(s2_Ec50_1 + pow(h_s2_Ec50_1,2)) - log(s2_Ec50_1_new + pow(h_s2_Ec50_1,2));
      log_accept = log_accept - 0.5 * pow(Ec50_1,2) * (1/s2_Ec50_1_new - 1/s2_Ec50_1);      
      
      accept = 1;
      if( arma::is_finite(log_accept) ){
        if(log_accept < 0){
          accept = exp(log_accept);
        }
      }else{
        accept = 0;
      }
      
      s2_Ec50_1_accept = s2_Ec50_1_accept + accept;
      s2_Ec50_1_count = s2_Ec50_1_count + 1;
      
      if( R::runif(0,1) < accept ){
        s2_Ec50_1 = s2_Ec50_1_new;
      }
      
      s_s2_Ec50_1 = s_s2_Ec50_1 + pow(g+1,-wg)*(accept - opt_rate);
      if(s_s2_Ec50_1 > exp(50)){
        s_s2_Ec50_1 = exp(50);
      }else{
        if(s_s2_Ec50_1 < exp(-50)){
          s_s2_Ec50_1 = exp(-50);
        }
      }
      
      
      
      //s2_Ec50_2
      //Propose new value
      s2_Ec50_2_new = sqrt(s2_Ec50_2) * exp(R::rnorm(0,1) * sqrt(s_s2_Ec50_2));
      s2_Ec50_2_new = pow(s2_Ec50_2_new,2);
      
      //Computing MH ratio:
      log_accept = log(s2_Ec50_2 + pow(h_s2_Ec50_2,2)) - log(s2_Ec50_2_new + pow(h_s2_Ec50_2,2));
      log_accept = log_accept - 0.5 * pow(Ec50_2,2) * (1/s2_Ec50_2_new - 1/s2_Ec50_2);   
      
      accept = 1;
      if( arma::is_finite(log_accept) ){
        if(log_accept < 0){
          accept = exp(log_accept);
        }
      }else{
        accept = 0;
      }
      
      s2_Ec50_2_accept = s2_Ec50_2_accept + accept;
      s2_Ec50_2_count = s2_Ec50_2_count + 1;
      
      if( R::runif(0,1) < accept ){
        s2_Ec50_2 = s2_Ec50_2_new;
      }
      
      s_s2_Ec50_2 = s_s2_Ec50_2 + pow(g+1,-wg)*(accept - opt_rate);
      if(s_s2_Ec50_2 > exp(50)){
        s_s2_Ec50_2 = exp(50);
      }else{
        if(s_s2_Ec50_2 < exp(-50)){
          s_s2_Ec50_2 = exp(-50);
        }
      }
    }
    
    
    //Save output
    if( (g + 1 > n_burn) & (((g + 1 - n_burn) / thin - floor((g + 1 - n_burn) / thin)) == 0 )){
      
      iter = (g + 1 - n_burn)/thin - 1;
      
      SLOPE_1(iter) = Slope_1;
      EC50_1(iter) = Ec50_1;
      LA_1(iter) = La_1;
      SLOPE_2(iter) = Slope_2;
      EC50_2(iter) = Ec50_2;
      LA_2(iter) = La_2;
      
      GAMMA_0(iter) = gamma_0;
      GAMMA_1(iter) = gamma_1;
      GAMMA_2(iter) = gamma_2;
      
      if(type == 1){
        C_OUT.row(iter) = vec_C.t();
      }else{
        B_OUT.row(iter) = vec_B.t();
        ELL(iter) = ell;
        if(type == 3){
          NU(iter) = nu;
        }
        if(type == 4){
          ALPHA(iter) = alpha;
        }
        SIGMA2_F(iter) = sigma2_f;
      }
      
      B1(iter) = b(0);
      B2(iter) = b(1);
      
      S2_EPS(iter) = s2_eps;
      S2_GAMMA_0(iter) = s2_gamma_0;
      S2_GAMMA_1(iter) = s2_gamma_1;
      S2_GAMMA_2(iter) = s2_gamma_2;
      S2_EC50_1(iter) = s2_Ec50_1;
      S2_EC50_2(iter) = s2_Ec50_2;
      
      //Performance evaluation (LPML)
      for(int r = 0; r < n_rep; r++){
        for(int i = 0; i < n1; i++){
          for(int j = 0; j < n2; j++){
            CPO(r,i,j) = CPO(r,i,j) + 1 / R::dnorm(wy_mat(j*n1 + i,r), p_ij(i,j), sqrt(s2_eps), 0);
          }
        }
      }
      
    }
    
    progr.increment(); 
  }
  
  
  double LPML = -accu(log( CPO / n_save ));
  
  List MCMC_Output_List, Monotherapies_List, Output_List;
  
  if(type == 1){//Splines
    Monotherapies_List = List::create(Named("SLOPE_1") = SLOPE_1, Named("EC50_1") = EC50_1, Named("SLOPE_2") = SLOPE_2, Named("EC50_2") = EC50_2,
                                            Named("LA_1") = LA_1, Named("LA_2") = LA_2);
    MCMC_Output_List = List::create(Named("MONO") = Monotherapies_List,
                                          Named("GAMMA_0") = GAMMA_0, Named ("GAMMA_1") = GAMMA_1, Named("GAMMA_2") = GAMMA_2,
                                                Named("C") = C_OUT, Named("B_K1") = B_K1, Named("B_K2") = B_K2, Named("B_K3") = B_K3, Named("B1") = B1, Named("B2") = B2,
                                                      Named("S2_EC50_1") = S2_EC50_1, Named("S2_EC50_2") = S2_EC50_2, Named("S2_GAMMA_0") = S2_GAMMA_0, Named("S2_GAMMA_1") = S2_GAMMA_1, Named("S2_GAMMA_2") = S2_GAMMA_2,
                                                            Named("S2_EPS") = S2_EPS, Named("LPML") = LPML);
  }
  
  if(type == 2){//Double exponential
    Monotherapies_List = List::create(Named("SLOPE_1") = SLOPE_1, Named("EC50_1") = EC50_1, Named("SLOPE_2") = SLOPE_2, Named("EC50_2") = EC50_2,
                                            Named("LA_1") = LA_1, Named("LA_2") = LA_2);
    MCMC_Output_List = List::create(Named("MONO") = Monotherapies_List,
                                            Named("GAMMA_0") = GAMMA_0, Named ("GAMMA_1") = GAMMA_1, Named("GAMMA_2") = GAMMA_2,
                                                     Named("B") = B_OUT, Named("ELL") = ELL, Named("SIGMA2_F") = SIGMA2_F, Named("B1") = B1, Named("B2") = B2,
                                                      Named("S2_EC50_1") = S2_EC50_1, Named("S2_EC50_2") = S2_EC50_2, Named("S2_GAMMA_0") = S2_GAMMA_0, Named("S2_GAMMA_1") = S2_GAMMA_1, Named("S2_GAMMA_2") = S2_GAMMA_2,
                                                            Named("S2_EPS") = S2_EPS, Named("LPML") = LPML);
  }
  
  if(type == 3){//Matern
    Monotherapies_List = List::create(Named("SLOPE_1") = SLOPE_1, Named("EC50_1") = EC50_1, Named("SLOPE_2") = SLOPE_2, Named("EC50_2") = EC50_2,
                                            Named("LA_1") = LA_1, Named("LA_2") = LA_2);
    MCMC_Output_List = List::create(Named("MONO") = Monotherapies_List,
                                          Named("GAMMA_0") = GAMMA_0, Named ("GAMMA_1") = GAMMA_1, Named("GAMMA_2") = GAMMA_2,
                                                Named("B") = B_OUT, Named("ELL") = ELL, Named("NU") = NU, Named("SIGMA2_F") = SIGMA2_F, Named("B1") = B1, Named("B2") = B2,
                                                Named("S2_EC50_1") = S2_EC50_1, Named("S2_EC50_2") = S2_EC50_2, Named("S2_GAMMA_0") = S2_GAMMA_0, Named("S2_GAMMA_1") = S2_GAMMA_1, Named("S2_GAMMA_2") = S2_GAMMA_2,
                                                      Named("S2_EPS") = S2_EPS, Named("LPML") = LPML);
  }
  
  if(type == 4){//Rational quadratic
    Monotherapies_List = List::create(Named("SLOPE_1") = SLOPE_1, Named("EC50_1") = EC50_1, Named("SLOPE_2") = SLOPE_2, Named("EC50_2") = EC50_2,
                                            Named("LA_1") = LA_1, Named("LA_2") = LA_2);
    MCMC_Output_List = List::create(Named("MONO") = Monotherapies_List,
                                          Named("GAMMA_0") = GAMMA_0, Named ("GAMMA_1") = GAMMA_1, Named("GAMMA_2") = GAMMA_2,
                                                Named("B") = B_OUT, Named("ELL") = ELL, Named("ALPHA") = ALPHA, Named("SIGMA2_F") = SIGMA2_F, Named("B1") = B1, Named("B2") = B2,
                                                      Named("S2_EC50_1") = S2_EC50_1, Named("S2_EC50_2") = S2_EC50_2, Named("S2_GAMMA_0") = S2_GAMMA_0, Named("S2_GAMMA_1") = S2_GAMMA_1, Named("S2_GAMMA_2") = S2_GAMMA_2,
                                                            Named("S2_EPS") = S2_EPS, Named("LPML") = LPML);
  }
  
  
  //Print acceptance rates
  Rcout << "Ec50_1 acc. rate = " << Ec50_1_accept/Ec50_1_count << "\n";
  Rcout << "Ec50_2 acc. rate = " << Ec50_2_accept/Ec50_2_count << "\n";
  Rcout << "Slope_1 acc. rate = " << Slope_1_accept/Slope_1_count << "\n";
  Rcout << "Slope_2 acc. rate = " << Slope_2_accept/Slope_2_count << "\n";
  if(update_LA){
    Rcout << "La_1 acc. rate = " << La_1_accept/La_1_count << "\n";
    Rcout << "La_2 acc. rate = " << La_2_accept/La_2_count << "\n";
  }
  Rcout << "gamma_0 acc. rate = " << gamma_0_accept/gamma_0_count << "\n";
  Rcout << "gamma_1 acc. rate = " << gamma_1_accept/gamma_1_count << "\n";
  Rcout << "gamma_2 acc. rate = " << gamma_2_accept/gamma_2_count << "\n";
  Rcout << "b acc. rate = " << b_accept/b_count << "\n";
  
  if(var_prior == 2){
    Rcout << "s2 Ec50_1 acc. rate = " << s2_Ec50_1_accept/s2_Ec50_1_count << "\n";
    Rcout << "s2 Ec50_2 acc. rate = " << s2_Ec50_2_accept/s2_Ec50_2_count << "\n";
    Rcout << "s2 gamma_0 acc. rate = " << s2_gamma_0_accept/s2_gamma_0_count << "\n";
    Rcout << "s2 gamma_1 acc. rate = " << s2_gamma_1_accept/s2_gamma_1_count << "\n";
    Rcout << "s2 gamma_2 acc. rate = " << s2_gamma_2_accept/s2_gamma_2_count << "\n"; 
    Rcout << "s2 eps acc. rate = " << s2_eps_accept/s2_eps_count << "\n";
  }
  
  if(type == 1){
    Rcout << "C acc. rate = " << C_accept/C_count << "\n";
  }else{
    Rcout << "B acc. rate = " << B_accept/B_count << "\n";
    if((update_sigma2_f) & (sigma2_f_prior == 2)){
      Rcout << "sigma2_f acc. rate = " << sigma2_f_accept/sigma2_f_count << "\n";
    }
    if(update_ell){
      Rcout << "ell acc. rate = " << ell_accept/ell_count << "\n";
    }
    if(update_nu){
      Rcout << "nu acc. rate = " << nu_accept/nu_count << "\n";
    }
    if(update_alpha){
      Rcout << "alpha acc. rate = " << alpha_accept/alpha_count << "\n";
    }
  }
  
  Output_List = List::create(Named("MCMC_Output") = MCMC_Output_List);
  
  return Output_List;
}

