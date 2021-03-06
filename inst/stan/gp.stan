
functions{
}
data {
  int<lower=1> n1;                                        // No. of concentrations for drug 1
  int<lower=1> n2;                                        // No. of concentrations for drug 2
  int<lower=0> nmissing;                                  // No. of missing grid locations
  int<lower=1> nrep;                                      // No. of replicates
  vector[(n1+n2+n1*n2+1)*nrep-nmissing] y;                   // Observed drug-response
  int ii_obs[(n1+n2+n1*n2+1)*nrep-nmissing];                   // Indices of observed variables
  real x1[n1];                                            // Concentrations of drug 1
  real x2[n2];                                            // Concentrations of drug 2
  int est_la;                                             // Boolean, estimating lower-asymptotes?
  int kernel;                                             // Indicator, which kernel are we using?
  int nu_matern;                                          // Specification of nu parameter for matern kernel
  int est_alpha;                                          // Are we estimating alpha for RQ kernel?
}
transformed data{
  int N = (n1+n2+n1*n2+1)*nrep-nmissing;                       // Total number of observations
  matrix[n1*n2,n1*n2] xdist;
  matrix[n1*n2,n1*n2] xdist_squared;
  // First we expand.grid to get full matrix X of all points
  matrix[n1*n2,2] X;
  int k = 1;
  for (j in 1:n1){
    for (i in 1:n2){
      X[k,1] = x1[j];
      X[k,2] = x2[i];
      k += 1;
    }
  }
  // Pairwise distances for kernel creation
  for(j in 1:(n1*n2)){
    for (i in j:(n1*n2)){
      xdist[i,j] = distance(X[i,],X[j,]);
      xdist[j,i] = xdist[i,j];
      xdist_squared[i,j] = squared_distance(X[i,],X[j,]);
      xdist_squared[j,i] = xdist_squared[i,j];
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
    vector[n1*n2] z;
    
    
    // Variances
    real<lower=0> s2;
    real<lower=0> s2_ec50_1;
    real<lower=0> s2_ec50_2; 
  }
  transformed parameters{
    matrix<lower=0, upper=1>[n2,n1] p0;        // Non-interaction
    row_vector<lower=0, upper=1>[n1] p01;      // Monotherapy drug 1
    vector<lower=0, upper=1>[n2] p02;          // Monotherapy drug 2
    matrix<lower=-1,upper=1>[n2,n1] Delta;     // Interaction
    
    
      {
        real la_1_param;
        real la_2_param;
        
        matrix[n2,n1] B; // Shorthand for what goes into g()
        matrix[n2,n1] GP; // The GP itself
        
        // The kernel construction
        matrix[(n1*n2),(n1*n2)] cov;
        matrix[(n1*n2),(n1*n2)] L_cov;
        if (kernel==1){ // RBF
          cov =  sigma_f*exp(- xdist_squared ./ (2*ell^2))  + diag_matrix(rep_vector(1e-10, (n1*n2)));
        }
        else if (kernel==2){ // Matern class
          matrix[(n1*n2),(n1*n2)] poly;
          if (nu_matern==1){ // nu = 1/2
            poly = rep_matrix(1,(n1*n2),(n1*n2));
            cov = sigma_f*(poly .* exp(-xdist ./ ell)) + diag_matrix(rep_vector(1e-10, (n1*n2)));
          }
          else if (nu_matern==2){ // nu = 3/2
            poly = (1+sqrt(3)*(xdist ./ ell));
            cov = sigma_f*(poly .* exp(-sqrt(3)*xdist ./ ell)) + diag_matrix(rep_vector(1e-10, (n1*n2)));
          }
          else if (nu_matern==3){ // nu = 5/2
            poly = (1+sqrt(5)*(xdist ./ ell)+(5./3.)*(xdist_squared ./ (ell^2)));
            cov = sigma_f*(poly .* exp(-sqrt(5)*xdist ./ ell)) + diag_matrix(rep_vector(1e-10, (n1*n2)));
          }
        }
        else if (kernel==3){ // Rational quadratic
          cov = sigma_f*(exp(-alpha[1]*log(1 + (xdist_squared ./ (2*alpha[1]*ell^2))))) + diag_matrix(rep_vector(1e-10, (n1*n2)));
        }
        L_cov = cholesky_decompose(cov);
        GP = to_matrix(L_cov*z,n2,n1);
        
        if (est_la){
          la_1_param = la_1[1];
          la_2_param = la_2[1];
        } else {
          la_1_param = 0;
          la_2_param = 0;
        }
        
        
        for (j in 1:n1){
          p01[j] = (la_1_param+(1-la_1_param)/(1+10^(slope_1*(x1[j]-ec50_1))));
          for (i in 1:n2){
            p02[i] = (la_2_param+(1-la_2_param)/(1+10^(slope_2*(x2[i]-ec50_2))));
            p0[i,j] = p01[j]*p02[i];
            B[i,j] = GP[i,j];
            Delta[i,j] = -p0[i,j]/(1+exp(b1*B[i,j]+log(p0[i,j]/(1-p0[i,j]))))+(1-p0[i,j])/(1+exp(-b2*B[i,j]-log(p0[i,j]/(1-p0[i,j]))));    //  Modified Andreas version
          }
        }
      }
  }
  model {
    // Constructing the dose-response over the full (n2+1)x(n1+1) grid
    matrix[n2+1,n1+1] f;                      // The dose-response function
    f[1,1] = 1;                               // At (-Inf,-Inf) dose response is one
    f[1,2:(n1+1)] = p01;                   // At (.,0) dose response is mono1
    f[2:(n2+1),1] = p02;                   // At (0,.) dose response is mono2
    f[2:(n2+1),2:(n1+1)] = p0 + Delta;  // At the interior, dose response is non-interaction + interaction
    
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
    ell ~ gamma(1,1);
    sigma_f ~ lognormal(1,1);
    if (est_alpha){
      alpha ~ gamma(1,1);
    }
    z ~ std_normal();
    
    // Response
    y ~ normal(to_vector(f)[ii_obs],sqrt(s2));
    
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
      real c11;                       // Integration limits dss_1
      real c12 = max(x1);                       // Integration limits dss_1
      real c21;                       // Integration limits dss_1
      real c22 = max(x2);                       // Integration limits dss_1
      vector[n2] B_rVUS;                        // Placeholder for trapezoidal rule
      vector[n2] B_Delta;                       // Placeholder for trapezoidal rule
      vector[n2] B_syn;                         // Placeholder for trapezoidal rule
      vector[n2] B_ant;                         // Placeholder for trapezoidal rule
      // Setting up drug response function
      f[1,1] = 1;                               // At (-Inf,-Inf) dose response is one
      f[1,2:(n1+1)] = p01;                   // At (.,0) dose response is mono1
      f[2:(n2+1),1] = p02;                   // At (0,.) dose response is mono2
      f[2:(n2+1),2:(n1+1)] = p0 + Delta;  // At the interior, dose response is non-interaction + interaction
      f_interior[1:n2,1:n1] = 1 - (p0 + Delta);  // At the interior, dose response is non-interaction + interaction
      for (i in 1:N){
        CPO[i] = exp(-normal_lpdf(y[i] | to_vector(f)[ii_obs[i]],sqrt(s2)));
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
           b_Delta += (x1[j]-x1[(j-1)])*(fabs(Delta[i,j])+fabs(Delta[i,(j-1)]))/2;
           b_syn += (x1[j]-x1[(j-1)])*(fabs(fmin(Delta[i,j],0))+fabs(fmin(Delta[i,(j-1)],0)))/2;
           b_ant += (x1[j]-x1[(j-1)])*(fabs(fmax(Delta[i,j],0))+fabs(fmax(Delta[i,(j-1)],0)))/2;
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
       rVUS_Delta = 100 * rVUS_Delta / (((max(x1)-min(x1))*(max(x2)-min(x2)))*fmax(max(p0),max(1-p0)));
       rVUS_syn = 100 * rVUS_syn / (((max(x1)-min(x1))*(max(x2)-min(x2)))*max(p0));
       rVUS_ant = 100 * rVUS_ant / (((max(x1)-min(x1))*(max(x2)-min(x2)))*max(1-p0));
    }
  }
  
      
