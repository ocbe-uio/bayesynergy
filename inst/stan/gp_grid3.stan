functions{
<<<<<<< HEAD
  #include "/include/kron_mvprod.stan"
  #include "/include/VUS.stan"
  #include "/include/VUS3.stan"
  #include "/include/dss.stan"
  #include "/include/matrix_min.stan"
  #include "/include/matrix_max.stan"
  #include "/include/matrix_array_min.stan"
  #include "/include/matrix_array_max.stan"
  #include "/include/one_minus_matrix_array.stan"
  #include "/include/sum_matrix_array.stan"
=======
#include "/include/kron_mvprod.stan"
#include "/include/VUS.stan"
#include "/include/VUS3.stan"
#include "/include/dss.stan"
#include "/include/matrix_min.stan"
#include "/include/matrix_max.stan"
#include "/include/matrix_array_min.stan"
#include "/include/matrix_array_max.stan"
#include "/include/one_minus_matrix_array.stan"
#include "/include/sum_matrix_array.stan"
>>>>>>> 47c44f2 (finished basic functionality for three-drug combinations)
}
data {
  int<lower=1> n1;                                        // No. of concentrations for drug 1
  int<lower=1> n2;                                        // No. of concentrations for drug 2
  int<lower=1> n3;                                        // No. of concentrations for drug 3
  int<lower=0> nmissing;                                  // No. of missing grid locations
  int<lower=1> nrep;                                      // No. of replicates
  vector[(1 + n1 + n2 + n3 + n1*n2 + n1*n3 + n2*n3 + n1*n2*n3)*nrep-nmissing] y; // Observed drug-response
  int ii_obs[(1 + n1 + n2 + n3 + n1*n2 + n1*n3 + n2*n3 + n1*n2*n3)*nrep-nmissing];              // Indices of observed variables
  real x1[n1];                                            // Concentrations of drug 1
  real x2[n2];                                            // Concentrations of drug 2
  real x3[n3];                                            // Concentrations of drug 2
  int est_la;                                             // Boolean, estimating lower-asymptotes?
  real lambda;                                            // Real, standard deviation ratio between positive and negative controls
}
transformed data{
  int N = (1 + n1 + n2 + n3 + n1*n2 + n1*n3 + n2*n3 + n1*n2*n3)*nrep-nmissing; // Total number of observations
  real c11 = min(x1);                       // Integration limits dss_1
  real c12 = max(x1);                       // Integration limits dss_1
  real c21 = min(x2);                       // Integration limits dss_2
  real c22 = max(x2);                       // Integration limits dss_2
  real c31 = min(x2);                       // Integration limits dss_3
  real c32 = max(x2);                       // Integration limits dss_3
}
parameters {
  // For Non-Interaction
  real<lower=0,upper=1> la_1[est_la ? 1 : 0];
  real<lower=0,upper=1> la_2[est_la ? 1 : 0];
  real<lower=0,upper=1> la_3[est_la ? 1 : 0];
  real log10_ec50_1;
  real log10_ec50_2;
  real log10_ec50_3;
  real theta_1;
  real theta_2;
  real theta_3;
  real<lower=0> slope_1;
  real<lower=0> slope_2;
  real<lower=0> slope_3;
  
  // For Interaction transformation
  real<lower=0> b1;
  real<lower=0> b2;
  
  // For the GP
  real<lower=0> ell_1;
  real<lower=0> ell_2;
  real<lower=0> ell_3;
  real<lower=0> sigma_f;
  matrix[n2,n1] z_12;
  matrix[n3,n1] z_13;
  matrix[n3,n2] z_23;
  matrix[n2,n1] z_123[n3];
  
  
  // Variances
  real<lower=0> s;
  real<lower=0> s2_log10_ec50_1;
  real<lower=0> s2_log10_ec50_2;
  real<lower=0> s2_log10_ec50_3;
}
transformed parameters{
  // Monotherapies
  vector<lower=0, upper=1>[n1] p01;          // Monotherapy drug 1
  vector<lower=0, upper=1>[n2] p02;          // Monotherapy drug 2
  vector<lower=0, upper=1>[n3] p03;          // Monotherapy drug 2
  // Noninteractions
  matrix<lower=0, upper=1> [n2,n1] p0_12; // Non-interaction between drugs 1 and 2
  matrix<lower=0, upper=1> [n3,n1] p0_13; // Non-interaction between drugs 1 and 3
  matrix<lower=0, upper=1> [n3,n2] p0_23; // Non-interaction between drugs 2 and 3
  
  matrix<lower=0, upper=1>[n2,n1] p0_123[n3];        // Non-interaction between drugs 1,2 and 3
  
  // Interactions
  matrix<lower=-1,upper=1>[n2,n1] Delta_12; // Interaction between drugs 1 and 2
  matrix<lower=-1,upper=1>[n3,n1] Delta_13; // Interaction between drugs 1 and 3
  matrix<lower=-1,upper=1>[n3,n2] Delta_23; // Interaction between drugs 2 and 3
  matrix<lower=-1,upper=1>[n2,n1] Delta_123[n3];  // Interaction between drugs 1 and 2 and 3
  
  
  {
    // The GP itself
    matrix[n2,n1] GP_12;
    matrix[n3,n1] GP_13;
    matrix[n3,n2] GP_23;
    matrix[n2,n1] GP_123[n3];
    real la_1_param;
    real la_2_param;
    real la_3_param;
    
    // The kernel construction
    matrix[n1,n1] cov1;
    matrix[n2,n2] cov2;
    matrix[n3,n3] cov3;
    matrix[n1,n1] L_cov1;
    matrix[n2,n2] L_cov2;
    matrix[n3,n3] L_cov3;
    
    cov1 = gp_matern32_cov(x1,sqrt(sigma_f),ell_1);
    cov2 = gp_matern32_cov(x2,1,ell_2);
    cov3 = gp_matern32_cov(x3,1,ell_3);
    
    L_cov1 = cholesky_decompose(cov1);
    L_cov2 = cholesky_decompose(cov2);
    L_cov3 = cholesky_decompose(cov3);
    
    // Sampling the GP
    GP_12 = kron_mvprod(L_cov1,L_cov2,z_12);
    GP_13 = kron_mvprod(L_cov1,L_cov3,z_13);
    GP_23 = kron_mvprod(L_cov2,L_cov3,z_23);
    
    for (k in 1:n3){
      GP_123[k] = kron_mvprod(L_cov1,L_cov2, z_123[k]);
    }
    for (j in 1:n1){
      for (i in 1:n2){
        GP_123[,i,j] = to_array_1d(L_cov3 * to_vector(z_123[,i,j]));
      }
    }
    
    if (est_la){
      la_1_param = la_1[1];
      la_2_param = la_2[1];
      la_3_param = la_3[1];
    } else {
      la_1_param = 0;
      la_2_param = 0;
      la_3_param = 0;
    }
    
    for (j in 1:n1){
      p01[j] = la_1_param+(1-la_1_param)/(1+10^(slope_1*(x1[j]-log10_ec50_1)));
      for (i in 1:n2){
        p02[i] = la_2_param+(1-la_2_param)/(1+10^(slope_2*(x2[i]-log10_ec50_2)));
        p0_12[i,j] = p01[j]*p02[i];
        
        Delta_12[i,j] = -p0_12[i,j]/(1+exp(b1*GP_12[i,j]+log(p0_12[i,j]/(1-p0_12[i,j]))))+(1-p0_12[i,j])/(1+exp(-b2*GP_12[i,j]-log(p0_12[i,j]/(1-p0_12[i,j]))));
        for (k in 1:n3){
          p03[k] = la_3_param+(1-la_3_param)/(1+10^(slope_3*(x3[k]-log10_ec50_3)));
          p0_13[k,j] = p01[j]*p03[k];
          p0_23[k,i] = p02[i]*p03[k];
          
          p0_123[k,i,j] = p01[j]*p02[i]*p03[k];
          
          Delta_13[k,j] = -p0_13[k,j]/(1+exp(b1*GP_13[k,j]+log(p0_13[k,j]/(1-p0_13[k,j]))))+(1-p0_13[k,j])/(1+exp(-b2*GP_13[k,j]-log(p0_13[k,j]/(1-p0_13[k,j]))));
          Delta_23[k,i] = -p0_23[k,i]/(1+exp(b1*GP_23[k,i]+log(p0_23[k,i]/(1-p0_23[k,i]))))+(1-p0_23[k,i])/(1+exp(-b2*GP_23[k,i]-log(p0_23[k,i]/(1-p0_23[k,i]))));
          
          Delta_123[k,i,j] = -p0_123[k,i,j]/(1+exp(b1*GP_123[k,i,j]+log(p0_123[k,i,j]/(1-p0_123[k,i,j]))))+(1-p0_123[k,i,j])/(1+exp(-b2*GP_123[k,i,j]-log(p0_123[k,i,j]/(1-p0_123[k,i,j]))));
        }
        
        
      }
    }
  }
}
model {
  // Constructing the dose-response over the full (n2+1)x(n1+1) grid
  matrix[n2+1,n1+1] f[(n3+1)];              // The dose-response function
  vector[N] fobs;                           // The observed dose response
  vector[N] noise;                          // The observation noise
  vector[(n2+1)*(n1+1)*(n3+1)] fobs_full;   // Vector version of f
  for (k in 1:(n3+1)){
    if (k == 1){
      f[k,1,1] = 1;                               // At (0,0,0) dose response is one
      f[k,1,2:(n1+1)] = to_row_vector(p01);         // At (0,0,:) dose response is mono1
      f[k,2:(n2+1),1] = p02;                      // At (0,:,0) dose response is mono2
      f[k,2:(n2+1),2:(n1+1)] = p0_12 + Delta_12;      // At the interior, dose response is non-interaction + interaction between drug 1 and 2
    } else {
      f[k,1,1] = p03[(k-1)]; // At (k,0,0) dose response is equal to corresponding value for mono3
      f[k,1,2:(n1+1)] = p0_13[(k-1),] + Delta_13[(k-1),];
      f[k,2:(n2+1),1] = to_vector(p0_23[(k-1),] + Delta_23[(k-1),]);
      f[k,2:(n2+1),2:(n1+1)] = p0_123[(k-1),,] + Delta_123[(k-1),,];
    }
    fobs_full[(k-1)*(n2+1)*(n1+1)+1:(k)*(n2+1)*(n1+1)] = to_vector(f[k]);
  }
  
  
  
  
  
  
  // Variances
  target += cauchy_lpdf(s | 0, 1) - cauchy_lccdf(0 | 0, 1);
  // target += inv_gamma_lpdf(s | 5, 1);
  target += inv_gamma_lpdf(s2_log10_ec50_1 | 3,2);
  target += inv_gamma_lpdf(s2_log10_ec50_2 | 3,2);
  target += inv_gamma_lpdf(s2_log10_ec50_3 | 3,2);
  
  
  // Monotherapies
  target += beta_lpdf(la_1 | 1,1.25);
  target += beta_lpdf(la_2 | 1,1.25);
  target += beta_lpdf(la_3 | 1,1.25);
  target += gamma_lpdf(slope_1 | 1,1);
  target += gamma_lpdf(slope_2 | 1,1);
  target += gamma_lpdf(slope_3 | 1,1);
  target += std_normal_lpdf(theta_1);
  target += std_normal_lpdf(theta_2);
  target += std_normal_lpdf(theta_3);
  target += normal_lpdf(log10_ec50_1 | theta_1,sqrt(s2_log10_ec50_1));
  target += normal_lpdf(log10_ec50_2 | theta_2,sqrt(s2_log10_ec50_2));
  target += normal_lpdf(log10_ec50_3 | theta_3,sqrt(s2_log10_ec50_3));
  
  // Interaction transformation
  target += normal_lpdf(b1 | 1,0.1);
  target += normal_lpdf(b2 | 1,0.1);
  
  // Interaction
  
  target += inv_gamma_lpdf(ell_1 | 5,5);
  target += inv_gamma_lpdf(ell_2 | 5,5);
  target += inv_gamma_lpdf(ell_3 | 5,5);
  target += lognormal_lpdf(sigma_f | 1,1);
  
  target += std_normal_lpdf(to_vector(z_12));
  target += std_normal_lpdf(to_vector(z_13));
  target += std_normal_lpdf(to_vector(z_23));
  for (k in 1:n3){
    target += std_normal_lpdf(to_vector(z_123[k]));
  }
  
  // Response
  fobs = fobs_full[ii_obs];
  noise = s*sqrt((fobs+lambda));
  
  target += normal_lpdf(y | fobs,noise);
  
}
generated quantities {
  real ec50_1;
  real ec50_2;
  real ec50_3;
  real dss_1 = 0;                             // DSS drug 1
  real dss_2 = 0;                             // DSS drug 2
  real dss_3 = 0;                             // DSS drug 2
  real rVUS_f_12 = 0;                            // Overall efficacy
  real rVUS_f_13 = 0;                            // Overall efficacy
  real rVUS_f_23 = 0;                            // Overall efficacy
  real rVUS_f_123 = 0;                            // Overall efficacy
  real rVUS_p0_12 = 0;                           // Noninteraction Efficacy
  real rVUS_p0_13 = 0;                           // Noninteraction Efficacy
  real rVUS_p0_23 = 0;                           // Noninteraction Efficacy
  real rVUS_p0_123 = 0;                           // Noninteraction Efficacy
  real VUS_Delta_12 = 0;                        // Overall interaction
  real VUS_Delta_13 = 0;                        // Overall interaction
  real VUS_Delta_23 = 0;                        // Overall interaction
  real VUS_Delta_123 = 0;                        // Overall interaction
  real VUS_syn_12 = 0;                          // Synergy
  real VUS_syn_13 = 0;                          // Synergy
  real VUS_syn_23 = 0;                          // Synergy
  real VUS_syn_123 = 0;                          // Synergy
  real VUS_ant_12 = 0;                          // Antagonism
  real VUS_ant_13 = 0;                          // Antagonism
  real VUS_ant_23 = 0;                          // Antagonism
  real VUS_ant_123 = 0;                          // Antagonism
  {
    // Constructing the dose-response over the full (n2+1)x(n1+1) grid
    matrix[n2+1,n1+1] f[(n3+1)];              // The dose-response function
    real la_1_param;                          // lower_asymptotes
    real la_2_param;                          // lower_asymptotes
    real la_3_param;                          // lower_asymptotes
    for (k in 1:(n3+1)){
      if (k == 1){
        f[k,1,1] = 1;                               // At (0,0,0) dose response is one
        f[k,1,2:(n1+1)] = to_row_vector(p01);         // At (0,0,:) dose response is mono1
        f[k,2:(n2+1),1] = p02;                      // At (0,:,0) dose response is mono2
        f[k,2:(n2+1),2:(n1+1)] = p0_12 + Delta_12;      // At the interior, dose response is non-interaction + interaction between drug 1 and 2
      } else {
        f[k,1,1] = p03[(k-1)]; // At (k,0,0) dose response is equal to corresponding value for mono3
        f[k,1,2:(n1+1)] = p0_13[(k-1),] + Delta_13[(k-1),];
        f[k,2:(n2+1),1] = to_vector(p0_23[(k-1),] + Delta_23[(k-1),]);
        f[k,2:(n2+1),2:(n1+1)] = p0_123[(k-1),,] + Delta_123[(k-1),,];
      }
    }
    
    // Getting EC50 on original scale
    ec50_1 = 10^log10_ec50_1;
    ec50_2 = 10^log10_ec50_2;
    ec50_3 = 10^log10_ec50_3;
    
    // Calculating DSS
    if (est_la){
      la_1_param = la_1[1];
      la_2_param = la_2[1];
      la_3_param = la_3[1];
    } else {
      la_1_param = 0;
      la_2_param = 0;
      la_3_param = 0;
    }
    dss_1 = dss(c11,c12,la_1_param,slope_1,log10_ec50_1);
    dss_2 = dss(c21,c22,la_2_param,slope_2,log10_ec50_2);
    dss_3 = dss(c31,c32,la_3_param,slope_3,log10_ec50_3);
    
    // Calculating drug combination scores
    // Drug 1 vs. Drug 2
    rVUS_f_12 = VUS(1-(p0_12+Delta_12),x1,x2);
    rVUS_p0_12 = VUS(1-p0_12,x1,x2);
    VUS_Delta_12 = VUS(Delta_12,x1,x2);
    VUS_syn_12 = VUS(matrix_min(Delta_12,0),x1,x2);
    VUS_ant_12 = VUS(matrix_max(Delta_12,0),x1,x2);
    
    // Drug 1 vs. Drug 3
    rVUS_f_13 = VUS(1-(p0_13+Delta_13),x1,x3);
    rVUS_p0_13 = VUS(1-p0_13,x1,x3);
    VUS_Delta_13 = VUS(Delta_13,x1,x3);
    VUS_syn_13 = VUS(matrix_min(Delta_13,0),x1,x3);
    VUS_ant_13 = VUS(matrix_max(Delta_13,0),x1,x3);
    
    // Drug 2 vs. Drug 3
    rVUS_f_23 = VUS(1-(p0_23+Delta_23),x2,x3);
    rVUS_p0_23 = VUS(1-p0_23,x2,x3);
    VUS_Delta_23 = VUS(Delta_23,x2,x3);
    VUS_syn_23 = VUS(matrix_min(Delta_23,0),x2,x3);
    VUS_ant_23 = VUS(matrix_max(Delta_23,0),x2,x3);
    
    // Drug 1 vs. Drug 2 vs. Drug 3
    rVUS_f_123 = VUS3(one_minus_matrix_array(sum_matrix_array(p0_123,Delta_123)),x1,x2,x3);
    rVUS_p0_123 = VUS3(one_minus_matrix_array(p0_123),x1,x2,x3);
    VUS_Delta_123 = VUS3(Delta_123,x1,x2,x3);
    VUS_syn_123 = VUS3(matrix_array_min(Delta_123,0),x1,x2,x3);
    VUS_ant_123 = VUS3(matrix_array_max(Delta_123,0),x1,x2,x3);
  }
}
