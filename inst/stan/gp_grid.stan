
functions{
  matrix kron_mvprod(matrix A, matrix B, matrix V) {
    // return (A kron_prod B) v where:
      // A is n1 x n1, B = n2 x n2, V = n2 x n1 = reshape(v,n2,n1)
      return transpose(A * transpose(B * V));
  }
}
data {
  int<lower=1> n1;                                        // No. of concentrations for drug 1
  int<lower=1> n2;                                        // No. of concentrations for drug 2
  int<lower=0> nmissing;                                  // No. of missing grid locations
  int<lower=1> nrep;                                      // No. of replicates
  vector[(n1+n2+n1*n2+1)*nrep-nmissing] pij;              // Observed drug-response
  int ii_obs[(n1+n2+n1*n2+1)*nrep-nmissing];              // Indices of observed variables
  real x1[n1];                                            // Concentrations of drug 1
  real x2[n2];                                            // Concentrations of drug 2
  int est_la;                                             // Boolean, estimating lower-asymptotes?
  int kernel;                                             // Indicator, which kernel are we using?
  int nu_matern;                                          // Specification of nu parameter for matern kernel
  int est_alpha;                                          // Are we estimating alpha for RQ kernel?
}
transformed data{
  int N = (n1+n2+n1*n2+1)*nrep-nmissing;                  // Total number of observations
  matrix[n1,n1] x1dist;
  matrix[n2,n2] x2dist;
  matrix[n1,n1] x1dist_squared;
  matrix[n2,n2] x2dist_squared;
 
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
    real ec50_1;
    real ec50_2;
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
    real<lower=0> s2;
    real<lower=0> s2_ec50_1;
    real<lower=0> s2_ec50_2; 
  }
  transformed parameters{
    matrix<lower=0, upper=1>[n2,n1] pij_0;        // Non-interaction
    row_vector<lower=0, upper=1>[n1] pij_01;      // Monotherapy drug 1
    vector<lower=0, upper=1>[n2] pij_02;          // Monotherapy drug 2
    matrix<lower=-1,upper=1>[n2,n1] Delta_ij;     // Interaction
    
    
      {
        real la_1_param;
        real la_2_param;
        
        
        matrix[n2,n1] GPij; // The GP itself
        matrix[n2,n1] Bij; // Shorthand for what goes into g()
        
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
        GPij = kron_mvprod(L_cov1,L_cov2,z);
        
        if (est_la){
          la_1_param = la_1[1];
          la_2_param = la_2[1];
        } else {
          la_1_param = 0;
          la_2_param = 0;
        }
        for (j in 1:n1){
          pij_01[j] = la_1_param+(1-la_1_param)/(1+10^(slope_1*(x1[j]-ec50_1)));
          for (i in 1:n2){
            pij_02[i] = la_2_param+(1-la_2_param)/(1+10^(slope_2*(x2[i]-ec50_2)));
            pij_0[i,j] = pij_01[j]*pij_02[i];
            Bij[i,j] = GPij[i,j];
            Delta_ij[i,j] = -pij_0[i,j]/(1+exp(b1*Bij[i,j]+log(pij_0[i,j]/(1-pij_0[i,j]))))+(1-pij_0[i,j])/(1+exp(-b2*Bij[i,j]-log(pij_0[i,j]/(1-pij_0[i,j]))));
          }
        }
      }
  }
  model {
    // Constructing the dose-response over the full (n2+1)x(n1+1) grid
    matrix[n2+1,n1+1] f;                      // The dose-response function
    f[1,1] = 1;                               // At (-Inf,-Inf) dose response is one
    f[1,2:(n1+1)] = pij_01;                   // At (.,0) dose response is mono1
    f[2:(n2+1),1] = pij_02;                   // At (0,.) dose response is mono2
    f[2:(n2+1),2:(n1+1)] = pij_0 + Delta_ij;  // At the interior, dose response is non-interaction + interaction
    
    // Variances
    s2 ~ inv_gamma(3,0.5);
    s2_ec50_1 ~ inv_gamma(3,2);
    s2_ec50_2 ~ inv_gamma(3,2);
    
    // Monotherapies
    la_1 ~ beta(.5,.5);
    la_2 ~ beta(.5,.5);
    slope_1 ~ gamma(1,1);
    slope_2 ~ gamma(1,1);
    ec50_1 ~ normal(0,sqrt(s2_ec50_1));
    ec50_2 ~ normal(0,sqrt(s2_ec50_2));
    
    // Interaction transformation
    b1 ~ gamma(1,1);
    b2 ~ gamma(1,1);
    
    // Interaction
    ell ~ inv_gamma(5,5);
    sigma_f ~ lognormal(0,1);
    if (est_alpha){
      alpha ~ gamma(1,1);
    }
    to_vector(z) ~ std_normal();
    
    // Response
    pij ~ normal(to_vector(f)[ii_obs],sqrt(s2));
    
  }
  generated quantities {
    vector[N] CPO = rep_vector(0,N);            // Model diagnostics
    real dss_1 = 0;                                 // DSS drug 1
    real dss_2 = 0;                                 // DSS drug 2
    real rVUS_p = 0;                            // Overall efficacy
    real rVUS_Delta = 0;                        // Overall interaction
    real rVUS_syn = 0;                          // Synergy
    real rVUS_ant = 0;                          // Antagonism
    {
      matrix[n2+1,n1+1] f;                      // The dose-response function
      matrix[n2,n1] f_interior;                 // The dose-response function interior
      real la_1_param;                          // lower_asymptotes
      real la_2_param;                          // lower_asymptotes
      real eps = 0.05;                          // for integration limits
      real c11;                                 // Integration limits dss_1
      real c12 = max(x1);                       // Integration limits dss_1
      real c21;                                 // Integration limits dss_1
      real c22 = max(x2);                       // Integration limits dss_1
      vector[n2] B_rVUS;                        // Placeholder for trapezoidal rule
      vector[n2] B_Delta;                       // Placeholder for trapezoidal rule
      vector[n2] B_syn;                         // Placeholder for trapezoidal rule
      vector[n2] B_ant;                         // Placeholder for trapezoidal rule
      // Setting up drug response function
      f[1,1] = 1;                               // At (-Inf,-Inf) dose response is one
      f[1,2:(n1+1)] = pij_01;                   // At (.,0) dose response is mono1
      f[2:(n2+1),1] = pij_02;                   // At (0,.) dose response is mono2
      f[2:(n2+1),2:(n1+1)] = pij_0 + Delta_ij;  // At the interior, dose response is non-interaction + interaction
      f_interior[1:n2,1:n1] = 1 - (pij_0 + Delta_ij);  // At the interior, dose response is non-interaction + interaction
      for (i in 1:N){
        CPO[i] = exp(-normal_lpdf(pij[i] | to_vector(f)[ii_obs[i]],sqrt(s2)));
      }
      // Calculating DSS
      if (est_la){
          la_1_param = la_1[1];
          la_2_param = la_2[1];
        } else {
          la_1_param = 0;
          la_2_param = 0;
        }
      
     if ((1 - eps/2) > la_1_param){
        c11 = (1/slope_1)*log10((1-la_1_param)/((1 - eps/2)-la_1_param)-1)+ec50_1;
        if (c11 < max(x1)){
          dss_1 = (c12-c11)+(la_1_param-1)/slope_1*(log10(1+10^(slope_1*(c12-ec50_1))) - log10(1+10^(slope_1*(c11-ec50_1))));
          dss_1 = 100 * (1-dss_1/((1-eps/2)*(c12-c11))); // Normalizing
        } else {dss_1 = uniform_rng(1e-16,1e-15);}
      } else {dss_1 = uniform_rng(1e-16,1e-15);}
      if ((1 - eps/2) > la_2_param){
        c21 = (1/slope_2)*log10((1-la_2_param)/((1 - eps/2)-la_2_param)-1)+ec50_2;
        if (c21 < max(x2)){
          dss_2 = (c22-c21)+(la_2_param-1)/slope_2*(log10(1+10^(slope_2*(c22-ec50_2))) - log10(1+10^(slope_2*(c21-ec50_2))));
          dss_2 = 100 * (1-dss_2/((1-eps/2)*(c22-c21))); // Normalizing
        } else {dss_2 = uniform_rng(1e-16,1e-15);}
      } else {dss_2 = uniform_rng(1e-16,1e-15);}
      
      // Calculating drug combination scores
       for (i in 1:n2){
         real b_rVUS = 0;
         real b_Delta = 0;
         real b_syn = 0;
         real b_ant = 0;
         for (j in 2:n1){
           b_rVUS += (x1[j]-x1[(j-1)])*(f_interior[i,j]+f_interior[i,(j-1)])/2;
           b_Delta += (x1[j]-x1[(j-1)])*(fabs(Delta_ij[i,j])+fabs(Delta_ij[i,(j-1)]))/2;
           b_syn += (x1[j]-x1[(j-1)])*(fabs(fmin(Delta_ij[i,j],0))+fabs(fmin(Delta_ij[i,(j-1)],0)))/2;
           b_ant += (x1[j]-x1[(j-1)])*(fabs(fmax(Delta_ij[i,j],0))+fabs(fmax(Delta_ij[i,(j-1)],0)))/2;
         }
         B_rVUS[i] = b_rVUS;
         B_Delta[i] = b_Delta;
         B_syn[i] = b_syn;
         B_ant[i] = b_ant;
         if (i > 1){
           rVUS_p += (x2[i]-x2[(i-1)])*(B_rVUS[i]+B_rVUS[(i-1)]) / 2;
           rVUS_Delta += (x2[i]-x2[(i-1)])*(B_Delta[i]+B_Delta[(i-1)]) / 2;
           rVUS_syn += (x2[i]-x2[(i-1)])*(B_syn[i]+B_syn[(i-1)]) / 2;
           rVUS_ant += (x2[i]-x2[(i-1)])*(B_ant[i]+B_ant[(i-1)]) / 2;
         }
       }
       // Normalizing
       rVUS_p = 100 * rVUS_p / ((max(x1)-min(x1))*(max(x2)-min(x2)));
       rVUS_Delta = 100 * rVUS_Delta / (((max(x1)-min(x1))*(max(x2)-min(x2)))*fmax(max(pij_0),max(1-pij_0)));
       rVUS_syn = 100 * rVUS_syn / (((max(x1)-min(x1))*(max(x2)-min(x2)))*max(pij_0));
       rVUS_ant = 100 * rVUS_ant / (((max(x1)-min(x1))*(max(x2)-min(x2)))*max(1-pij_0));
    }
  }
  
      
