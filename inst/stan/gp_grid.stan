functions{
  #include /include/kron_mvprod.stan
  #include /include/VUS.stan
  #include /include/dss.stan
  #include /include/matrix_min.stan
  #include /include/matrix_max.stan
  #include /include/pc_prior.stan
  #include /include/lptn.stan
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
  int robust;                                             // Boolean, using mixture likelihood?
  int pcprior;                                            // Boolean, using PC prior for Matern covariance?
  int heteroscedastic;                                    // Boolean, are we assuming heteroscedastic measurement error?
  vector[4] pcprior_hypers;                               // Hyperparameters for the PC prior
  real rho;                                               // Hyperparameter for LPTN
  real lambda;                                            // Real, standard deviation ratio between positive and negative controls
  int kernel;                                             // Indicator, which kernel are we using?
  int nu_matern;                                          // Specification of nu parameter for matern kernel
  int est_alpha;                                          // Are we estimating alpha for RQ kernel?
}
transformed data{
  int N = (n1+n2+n1*n2+1)*nrep-nmissing;                  // Total number of observations
  real c11 = min(x1);                       // Integration limits dss_1
  real c12 = max(x1);                       // Integration limits dss_1
  real c21 = min(x2);                       // Integration limits dss_2
  real c22 = max(x2);                       // Integration limits dss_2
  real d_half = 2.0/2.0;
  real tauLPTN = inv_Phi((1+rho)/2);
  real lambdaLPTN = 2*inv(1-rho)*exp(normal_lpdf(tauLPTN | 0,1))*tauLPTN*log(tauLPTN);
  matrix[n1,n1] x1dist;
  matrix[n2,n2] x2dist;
  matrix[n1,n1] x1dist_squared;
  matrix[n2,n2] x2dist_squared;
  real lambda1 = -log(pcprior_hypers[2])*pow(pcprior_hypers[1],d_half);
  real lambda2 = -log(pcprior_hypers[4])*pow(pcprior_hypers[3],-1);

  // Pairwise distances for kernel creation
  for (i in 1:n1){
    for (j in i:n1){
      x1dist[i,j] = sqrt((x1[i]-x1[j])^2);
      x1dist[j,i] = x1dist[i,j];
      x1dist_squared[i,j] = x1dist[i,j]^2;
      x1dist_squared[j,i] = x1dist_squared[i,j];
    }
  }
  for (i in 1:n2){
    for (j in i:n2){
      x2dist[i,j] = sqrt((x2[i]-x2[j])^2);
      x2dist[j,i] = x2dist[i,j];
      x2dist_squared[i,j] = x2dist[i,j]^2;
      x2dist_squared[j,i] = x2dist_squared[i,j];
    }
  }
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

  // For Interaction transformation
  real<lower=0> b1;
  real<lower=0> b2;

  // For the GP
  real<lower=0> ell;
  real<lower=0> sigma_f;
  real<lower=0> alpha[est_alpha ? 1 : 0];
  matrix[n2,n1] z;


  // Variances
  real<lower=0> s;
  real<lower=0> s2_log10_ec50_1;
  real<lower=0> s2_log10_ec50_2;
}
transformed parameters{
  matrix<lower=0, upper=1>[n2,n1] p0;        // Non-interaction
  row_vector<lower=0, upper=1>[n1] p01;      // Monotherapy drug 1
  vector<lower=0, upper=1>[n2] p02;          // Monotherapy drug 2
  matrix<lower=-1,upper=1>[n2,n1] Delta;     // Interaction
  matrix[n2,n1] GP; // The GP itself

  {
    real la_1_param;
    real la_2_param;




    matrix[n2,n1] B; // Shorthand for what goes into g()

    // The kernel construction
    matrix[n1,n1] cov1;
    matrix[n2,n2] cov2;
    matrix[n1,n1] L_cov1;
    matrix[n2,n2] L_cov2;
    if (kernel==1){ // RBF
    cov1 =  sigma_f*exp(- x1dist_squared ./ (2*ell^2))  + diag_matrix(rep_vector(1e-10, n1));
    cov2 = exp(- x2dist_squared ./ (2*ell^2))  + diag_matrix(rep_vector(1e-10, n2));
    }
    else if (kernel==2){ // Matern class
    matrix[n1,n1] poly1;
    matrix[n2,n2] poly2;
    if (nu_matern==1){ // nu = 1/2
    poly1 = rep_matrix(1,n1,n1);
    poly2 = rep_matrix(1,n2,n2);
    cov1 = sigma_f*(poly1 .* exp(-x1dist ./ ell)) + diag_matrix(rep_vector(1e-10, n1));
    cov2 = (poly2 .* exp(-x2dist ./ ell)) + diag_matrix(rep_vector(1e-10, n2));
    }
    else if (nu_matern==2){ // nu = 3/2
    poly1 = (1+sqrt(3)*(x1dist ./ ell));
    poly2 = (1+sqrt(3)*(x2dist ./ ell));
    cov1 = sigma_f*(poly1 .* exp(-sqrt(3)*x1dist ./ ell)) + diag_matrix(rep_vector(1e-10, n1));
    cov2 = poly2 .* exp(-sqrt(3)*x2dist ./ ell) + diag_matrix(rep_vector(1e-10, n2));
    }
    else if (nu_matern==3){ // nu = 5/2
    poly1 = (1+sqrt(5)*(x1dist ./ ell)+(5./3.)*(x1dist_squared ./ (ell^2)));
    poly2 = (1+sqrt(5)*(x2dist ./ ell)+(5./3.)*(x2dist_squared ./ (ell^2)));
    cov1 = sigma_f*(poly1 .* exp(-sqrt(5)*x1dist ./ ell)) + diag_matrix(rep_vector(1e-10, n1));
    cov2 = poly2 .* exp(-sqrt(5)*x2dist ./ ell) + diag_matrix(rep_vector(1e-10, n2));
    }
    }
    else if (kernel==3){ // Rational quadratic
    cov1 = sigma_f*(exp(-alpha[1]*log(1 + (x1dist_squared ./ (2*alpha[1]*ell^2))))) + diag_matrix(rep_vector(1e-10, n1));
    cov2 = exp(-alpha[1]*log(1 + (x2dist_squared ./ (2*alpha[1]*ell^2)))) + diag_matrix(rep_vector(1e-10, n2));
    }
    L_cov1 = cholesky_decompose(cov1);
    L_cov2 = cholesky_decompose(cov2);
    GP = kron_mvprod(L_cov1,L_cov2,z);

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
        B[i,j] = GP[i,j];
        Delta[i,j] = -p0[i,j]/(1+exp(b1*B[i,j]+log(p0[i,j]/(1-p0[i,j]))))+(1-p0[i,j])/(1+exp(-b2*B[i,j]-log(p0[i,j]/(1-p0[i,j]))));
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
  f[2:(n2+1),2:(n1+1)] = p0 + Delta;      // At the interior, dose response is non-interaction + interaction


  // Variances
  target += cauchy_lpdf(s | 0, 1) - cauchy_lccdf(0 | 0, 1);
  // target += inv_gamma_lpdf(s | 5, 1);
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

  // Interaction transformation
  target += normal_lpdf(b1 | 1,0.1);
  target += normal_lpdf(b2 | 1,0.1);

  // Interaction
  if (pcprior){
    target += pc_prior(ell,sigma_f,lambda1,lambda2);

  } else{
    target += inv_gamma_lpdf(ell | 5,5);
    target += lognormal_lpdf(sigma_f | 1,1);
  }

  if (est_alpha){
    target += gamma_lpdf(alpha | 1,1);
  }
  target += std_normal_lpdf(to_vector(z));

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
  real VUS_Delta = 0;                        // Overall interaction
  real VUS_syn = 0;                          // Synergy
  real VUS_ant = 0;                          // Antagonism
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
    f[2:(n2+1),2:(n1+1)] = p0 + Delta;      // At the interior, dose response is non-interaction + interaction
    fc_interior[1:n2,1:n1] = 1 - (p0 + Delta);// At the interior, dose response is non-interaction + interaction


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
    VUS_Delta = VUS(Delta,x1,x2);
    VUS_syn = VUS(matrix_min(Delta,0),x1,x2);
    VUS_ant = VUS(matrix_max(Delta,0),x1,x2);

    // Fixing zero valued integrals so that sampler doesn't complain too much
    if (rVUS_f == 0){rVUS_f = uniform_rng(1e-6,1e-4);}
    if (rVUS_p0 == 0){rVUS_p0 = uniform_rng(1e-6,1e-4);}
    if (VUS_syn == 0){VUS_syn = uniform_rng(1e-6,1e-4);}
    if (VUS_ant == 0){VUS_ant = uniform_rng(1e-6,1e-4);}
  }
}


