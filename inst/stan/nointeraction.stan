// A version of gp_grid.stan that only estimates f = p0. Used for computing the Bayes factor
functions{
#include /include/VUS.stan
#include /include/lptn.stan
#include /include/dss.stan
}
data {
  int<lower=1> n1;                                        // No. of concentrations for drug 1
  int<lower=1> n2;                                        // No. of concentrations for drug 2
  int<lower=0> nmissing;                                  // No. of missing grid locations
  int<lower=1> nrep;                                      // No. of replicates
  vector[(n1+n2+n1*n2+1)*nrep-nmissing] y;                // Observed drug-response
  int ii_obs[(n1+n2+n1*n2+1)*nrep-nmissing];              // Indices of observed variables
  real x1[n1];                                            // Concentrations of drug 1
  real x2[n2];                                            // Concentrations of drug 2
  int est_la;                                             // Boolean, estimating lower-asymptotes?
  int heteroscedastic;                                    // Boolean, are we assuming heteroscedastic measurement error?
  int robust;                                             // Boolean, using mixture likelihood?
  int pcprior;                                            // Boolean, using PC prior for Matern covariance?
  vector[4] pcprior_hypers;                               // Hyperparameters for the PC prior
  real rho;                                               // Hyperparameter for LPTN
  real lambda;                                            // Real, standard deviation ratio between positive and negative controls
}
transformed data{
  int N = (n1+n2+n1*n2+1)*nrep-nmissing;                  // Total number of observations
  real c11 = min(x1);                       // Integration limits dss_1
  real c12 = max(x1);                       // Integration limits dss_1
  real c21 = min(x2);                       // Integration limits dss_2
  real c22 = max(x2);                       // Integration limits dss_2
  real tauLPTN = inv_Phi((1+rho)/2);
  real lambdaLPTN = 2*inv(1-rho)*exp(normal_lpdf(tauLPTN | 0,1))*tauLPTN*log(tauLPTN);
}
parameters {
  // For Non-Interaction
  real<lower=0,upper=1> la_1[est_la ? 1 : 0];
  real<lower=0,upper=1> la_2[est_la ? 1 : 0];
  real log10_ec50_1;
  real log10_ec50_2;
  real theta_1;
  real theta_2;
  real<lower=0> slope_1;
  real<lower=0> slope_2;

  // Variances
  real<lower=0> s;
  real<lower=0> s2_log10_ec50_1;
  real<lower=0> s2_log10_ec50_2;

}
transformed parameters{
  matrix<lower=0, upper=1>[n2,n1] p0;        // Non-interaction
  row_vector<lower=0, upper=1>[n1] p01;      // Monotherapy drug 1
  vector<lower=0, upper=1>[n2] p02;          // Monotherapy drug 2

  {
    real la_1_param;
    real la_2_param;

    if (est_la){
      la_1_param = la_1[1];
      la_2_param = la_2[1];
    } else {
      la_1_param = 0;
      la_2_param = 0;
    }

    for (j in 1:n1){
      p01[j] = la_1_param+(1-la_1_param)/(1+10^(slope_1*(x1[j]-log10_ec50_1)));
      for (i in 1:n2){
        p02[i] = la_2_param+(1-la_2_param)/(1+10^(slope_2*(x2[i]-log10_ec50_2)));
        p0[i,j] = p01[j]*p02[i];
      }
    }
  }
}
model {
  // Constructing the dose-response over the full (n2+1)x(n1+1) grid
  matrix[n2+1,n1+1] f;                      // The dose-response function
  vector[N] fobs;                           // The observed dose response
  vector[N] noise;                          // The observation noise
  f[1,1] = 1;                               // At (-Inf,-Inf) dose response is one
  f[1,2:(n1+1)] = p01;                      // At (.,0) dose response is mono1
  f[2:(n2+1),1] = p02;                      // At (0,.) dose response is mono2
  f[2:(n2+1),2:(n1+1)] = p0;              // At the interior, dose response is non-interaction

  // Variances
  target += cauchy_lpdf(s | 0, 1)  - cauchy_lccdf(0 | 0, 1);
  target += inv_gamma_lpdf(s2_log10_ec50_1 | 3,2);
  target += inv_gamma_lpdf(s2_log10_ec50_2 | 3,2);


  // Monotherapies
  target += beta_lpdf(la_1 | 1,1.25);
  target += beta_lpdf(la_2 | 1,1.25);
  target += gamma_lpdf(slope_1 | 1,1);
  target += gamma_lpdf(slope_2 | 1,1);
  target += std_normal_lpdf(theta_1);
  target += std_normal_lpdf(theta_2);
  target += normal_lpdf(log10_ec50_1 | theta_1,sqrt(s2_log10_ec50_1));
  target += normal_lpdf(log10_ec50_2 | theta_2,sqrt(s2_log10_ec50_2));

  // Response
  fobs = to_vector(f)[ii_obs];
  noise = s*sqrt((fobs+lambda));
  if (heteroscedastic){
    if (robust){
      for (n in 1:N) {
        target += lptn(rho,(y[n]-fobs[n])/noise[n],lambdaLPTN,tauLPTN) - log(noise[n]);
      }
    } else{
      target += normal_lpdf(y | fobs,noise);
    }

  } else{
    if (robust){
      for (n in 1:N) {
        target += lptn(rho,(y[n]-fobs[n])/s,lambdaLPTN,tauLPTN) - log(s);
      }
    } else{
      target += normal_lpdf(y | fobs,s);
    }

  }

}
generated quantities {
  real ec50_1;
  real ec50_2;
  vector[N] CPO = rep_vector(0,N);            // Model diagnostics
  real dss_1 = 0;                             // DSS drug 1
  real dss_2 = 0;                             // DSS drug 2
  real rVUS_f = 0;                            // Overall efficacy
  real rVUS_p0 = 0;                           // Noninteraction Efficacy

  {
    matrix[n2+1,n1+1] f;                      // The dose-response function
    matrix[n2,n1] fc_interior;                // The complement of dose-response function interior
    vector[N] fobs;                           // The observed dose response
    vector[N] noise;                          // The observation noise
    real la_1_param;                          // lower_asymptotes
    real la_2_param;                          // lower_asymptotes

    // Setting up drug response function
    f[1,1] = 1;                               // At (-Inf,-Inf) dose response is one
    f[1,2:(n1+1)] = p01;                      // At (.,0) dose response is mono1
    f[2:(n2+1),1] = p02;                      // At (0,.) dose response is mono2
    f[2:(n2+1),2:(n1+1)] = p0;              // At the interior, dose response is non-interaction
    fc_interior[1:n2,1:n1] = 1 - (p0);        // At the interior, dose response is non-interaction


    // Calculating CPO
    fobs = to_vector(f)[ii_obs];
    noise = s*sqrt((fobs+lambda));
    for (i in 1:N){
      if (heteroscedastic){
        CPO[i] = exp(-normal_lpdf(y[i] | fobs[i], noise[i]));
      } else {
        CPO[i] = exp(-normal_lpdf(y[i] | fobs[i], s));
      }
    }
    // Getting EC50 on original scale
    ec50_1 = 10^log10_ec50_1;
    ec50_2 = 10^log10_ec50_2;
    // Calculating DSS
    if (est_la){
      la_1_param = la_1[1];
      la_2_param = la_2[1];
    } else {
      la_1_param = 0;
      la_2_param = 0;
    }

    dss_1 = dss(c11,c12,la_1_param,slope_1,log10_ec50_1);
    dss_2 = dss(c21,c22,la_2_param,slope_2,log10_ec50_2);

    // Calculating drug combination scores
    rVUS_f = VUS(fc_interior,x1,x2);
    rVUS_p0 = VUS(1-p0,x1,x2);

    // Fixing zero valued integrals for so that sampler doesn't complain too much
    if (rVUS_f == 0){rVUS_f = uniform_rng(1e-6,1e-4);}
    if (rVUS_p0 == 0){rVUS_p0 = uniform_rng(1e-6,1e-4);}
  }
}


