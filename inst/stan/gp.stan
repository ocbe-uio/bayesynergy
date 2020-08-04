
functions{
}
data {
  int<lower=1> n1;                                        // No. of concentrations for drug 1
  int<lower=1> n2;                                        // No. of concentrations for drug 2
  int<lower=0> nmissing;                                  // No. of missing grid locations
  vector[(n1+n2+n1*n2+1-nmissing)] pij;                   // Observed drug-response
  int ii_obs[(n1+n2+n1*n2+1-nmissing)];                   // Indices of observed variables
  real x1[n1];                                            // Concentrations of drug 1
  real x2[n2];                                            // Concentrations of drug 2
  int est_la;                                             // Boolean, estimating lower-asymptotes?
  int kernel;                                             // Indicator, which kernel are we using?
  int nu_matern;                                          // Specification of nu parameter for matern kernel
  int est_alpha;                                          // Are we estimating alpha for RQ kernel?
}
transformed data{
  int N = (n1+n2+n1*n2+1-nmissing);                       // Total number of observations
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
    real<lower=0> s2_gamma0;
    real<lower=0> s2_gamma1;
    real<lower=0> s2_gamma2;
  }
  transformed parameters{
    matrix<lower=0, upper=1>[n2,n1] pij_0;        // Non-interaction
    row_vector<lower=0, upper=1>[n1] pij_01;      // Monotherapy drug 1
    vector<lower=0, upper=1>[n2] pij_02;          // Monotherapy drug 2
    matrix<lower=-1,upper=1>[n2,n1] Delta_ij;     // Interaction
    
    
      {
        real la_1_param;
        real la_2_param;
        
        matrix[n2,n1] Bij; // Shorthand for what goes into g()
        matrix[n2,n1] GPij; // The GP itself
        
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
        GPij = to_matrix(L_cov*z,n2,n1);
        
        if (est_la){
          la_1_param = la_1[1];
          la_2_param = la_2[1];
        } else {
          la_1_param = 0;
          la_2_param = 0;
        }
        
        
        for (j in 1:n1){
          pij_01[j] = (la_1_param+(1-la_1_param)/(1+10^(slope_1*(x1[j]-ec50_1))));
          for (i in 1:n2){
            pij_02[i] = (la_2_param+(1-la_2_param)/(1+10^(slope_2*(x2[i]-ec50_2))));
            pij_0[i,j] = pij_01[j]*pij_02[i];
            Bij[i,j] = GPij[i,j];
            Delta_ij[i,j] = -pij_0[i,j]/(1+exp(b1*Bij[i,j]+log(pij_0[i,j]/(1-pij_0[i,j]))))+(1-pij_0[i,j])/(1+exp(-b2*Bij[i,j]-log(pij_0[i,j]/(1-pij_0[i,j]))));    //  Modified Andreas version
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
    s2 ~ inv_gamma(3,2);
    s2_ec50_1 ~ inv_gamma(3,2);
    s2_ec50_2 ~ inv_gamma(3,2);
    s2_gamma0 ~ inv_gamma(3,2);
    s2_gamma1 ~ inv_gamma(3,2);
    s2_gamma2 ~ inv_gamma(3,2);
    
    // Monotherapies
    if (est_la){
      la_1 ~ beta(.5,.5);
      la_2 ~ beta(.5,.5);
    }
    slope_1 ~ gamma(1,1);
    slope_2 ~ gamma(1,1);
    ec50_1 ~ normal(0,sqrt(s2_ec50_1));
    ec50_2 ~ normal(0,sqrt(s2_ec50_2));
    
    // Interaction transformation
    b1 ~ gamma(1,1);
    b2 ~ gamma(1,1);
    
    // Interaction
    // ell ~ cauchy(0,2.5);
    ell ~ gamma(1,1);
    // sigma_f ~ inv_gamma(3,2);
    // sigma_f ~ normal(0,5);
    sigma_f ~ lognormal(0,1);
    if (est_alpha){
      alpha ~ gamma(1,1);
    }
    z ~ std_normal();
    
    // Response
    pij ~ normal(to_vector(f)[ii_obs],sqrt(s2));
    
  }
  generated quantities {
    vector[N] CPO = rep_vector(0,N);
    {
      matrix[n2+1,n1+1] f;                      // The dose-response function
      f[1,1] = 1;                               // At (-Inf,-Inf) dose response is one
      f[1,2:(n1+1)] = pij_01;                   // At (.,0) dose response is mono1
      f[2:(n2+1),1] = pij_02;                   // At (0,.) dose response is mono2
      f[2:(n2+1),2:(n1+1)] = pij_0 + Delta_ij;  // At the interior, dose response is non-interaction + interaction
    
      for (i in 1:N){
        CPO[i] += exp(-normal_lpdf(pij[i] | to_vector(f)[ii_obs[i]],sqrt(s2)));
      }
    }
  }
  
      
