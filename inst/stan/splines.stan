
functions{
  matrix kron_mvprod(matrix A, matrix B, matrix V) {
    // return (A kron_prod B) v where:
      // A is n1 x n1, B = n2 x n2, V = n2 x n1 = reshape(v,n2,n1)
      return transpose(A * transpose(B * V));
  }
  // Building b spline basis
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order);
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    // INPUTS:
      //    t:          the points at which the b_spline is calculated
    //    ext_knots:  the set of extended knots
    //    ind:        the index of the b_spline
    //    order:      the order of the b-spline
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order==1)
      for (i in 1:size(t)) // B-splines of order 1 are piece-wise constant
    b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]);
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) /
          (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) /
          (ext_knots[ind+order] - ext_knots[ind+1]);
      // Calculating the B-spline recursively as linear interpolation of two lower-order splines
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) +
        w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  }
}
data {
  int<lower=1> n1;                                        // No. of concentrations for drug 1
  int<lower=1> n2;                                        // No. of concentrations for drug 2
  int<lower=1> n_knots1;                                  // No. of knots for drug 1
  int<lower=1> n_knots2;                                  // No. of knots for drug 2
  int degree;                                             // Degree of splines
  int<lower=0> nmissing;                                  // No. of missing grid locations
  vector[(n1+n2+n1*n2+1-nmissing)] pij;                   // Observed drug-response
  int ii_obs[(n1+n2+n1*n2+1-nmissing)];                   // Indices of observed variables
  real x1[n1];                                            // Concentrations of drug 1
  real x2[n2];                                            // Concentrations of drug 2
  vector[n_knots1] t1;                                    // Knots for drug 1
  vector[n_knots2] t2;                                    // Knots for drug 2
  int est_la;                                             // Boolean, estimating lower-asymptotes
}
transformed data{
  int N = (n1+n2+n1*n2+1-nmissing);
  
  // Construction of B spline basis
  int num_basis1 = n_knots1 + degree - 1;              // total number of B-splines drug 1
  matrix[num_basis1, n1] B1;                           // matrix of B-splines drug 1
  int num_basis2 = n_knots2 + degree - 1;              // total number of B-splines drug 2
  matrix[num_basis2, n2] B2;                           // matrix of B-splines drug 2
  
  matrix[num_basis2, num_basis1] B[n2,n1];             // Combined Spline array
  
  matrix[num_basis1,num_basis1] U = 2*diag_matrix(rep_vector(1,num_basis1)); // Second order difference
  matrix[num_basis2,num_basis2] V = 2*diag_matrix(rep_vector(1,num_basis2)); // Second order difference
  
  matrix[num_basis1,num_basis1] L_cov1;
  matrix[num_basis2,num_basis2] L_cov2;
  
  
  
  {
    vector[degree + n_knots1] ext_knots_temp1;
    vector[2*degree + n_knots1] ext_knots1;               // set of extended knots
    vector[degree + n_knots2] ext_knots_temp2;
    vector[2*degree + n_knots2] ext_knots2;               // set of extended knots
    ext_knots_temp1 = append_row(rep_vector(t1[1], degree), t1);
    ext_knots1 = append_row(ext_knots_temp1, rep_vector(t1[n_knots1], degree));
    ext_knots_temp2 = append_row(rep_vector(t2[1], degree), t2);
    ext_knots2 = append_row(ext_knots_temp2, rep_vector(t2[n_knots2], degree));
    for (i in 1:num_basis1){
      B1[i,:] = to_row_vector(build_b_spline(x1, to_array_1d(ext_knots1), i, degree + 1));
    }
    for (i in 1:num_basis2){
      B2[i,:] = to_row_vector(build_b_spline(x2, to_array_1d(ext_knots2), i, degree + 1));
    }
    B1[n_knots1 + degree - 1, n1] = 1;
    B2[n_knots2 + degree - 1, n2] = 1;
    
    for (i in 2:num_basis1){
      U[i,i-1] = -1;
    }
    for (i in 1:num_basis1-1){
      U[i,i+1] = -1;
    }
    for (i in 2:num_basis2){
      V[i,i-1] = -1;
    }
    for (i in 1:num_basis2-1){
      V[i,i+1] = -1;
    }
    
    for (i in 1:n2){
      for (j in 1:n1){
        B[i,j,:,:] = B2[:,i]*B1[:,j]';
          }
        }
        
      }

    L_cov1 = cholesky_decompose(U);
    L_cov2 = cholesky_decompose(V);
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
    real gamma0;
    real gamma1;
    real gamma2;
    real<lower=0> b1;
    real<lower=0> b2;
    
    // For the interaction splines
    matrix[num_basis2, num_basis1] z;
    
    
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
        matrix[num_basis2,num_basis1] C; // Spline coefficients
         matrix[n2,n1] Bij; // Shorthand for what goes into g() i.e. the splines here
        
        C = kron_mvprod(L_cov1, L_cov2,z);
        
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
            Bij[i,j] = gamma0 + gamma1*x1[j] + gamma2*x2[i] + sum(C .* B[i,j,:,:]); // SPLINES
            Delta_ij[i,j] = -pij_0[i,j]/(1+exp(b1*Bij[i,j]))+(1-pij_0[i,j])/(1+exp(-b2*Bij[i,j]));    //  Modified Andreas version
          }
        }
      }
  }
  model {
    // Constructing the dose-response over the fulle (n2+1)x(n1+1) grid
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
    gamma0 ~ normal(0,sqrt(s2_gamma0));
    gamma1 ~ normal(0,sqrt(s2_gamma1));
    gamma2 ~ normal(0,sqrt(s2_gamma2));
    b1 ~ gamma(1,1);
    b2 ~ gamma(1,1);
    
    // Interaction (defined via latent z, because covariance is fixed)
    to_vector(z) ~ std_normal();
    
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
  
      
