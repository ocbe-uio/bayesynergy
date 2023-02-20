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
  vector[(n1+n2+n1*n2+1)*nrep-nmissing] y;                // Observed drug-response
  int ii_obs[(n1+n2+n1*n2+1)*nrep-nmissing];              // Indices of observed variables
  real x1[n1];                                            // Concentrations of drug 1
  real x2[n2];                                            // Concentrations of drug 2
  int est_la;                                             // Boolean, estimating lower-asymptotes?
  real lambda;                                            // Real, standard deviation ratio between positive and negative controls
}
transformed data{
  int N = (n1+n2+n1*n2+1)*nrep-nmissing;                  // Total number of observations
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

    cov1 = gp_matern32_cov(x1,sqrt(sigma_f),ell);
    cov2 = gp_matern32_cov(x1,1,ell);

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

  target += inv_gamma_lpdf(ell | 5,5);
  target += lognormal_lpdf(sigma_f | 1,1);

  target += std_normal_lpdf(to_vector(z));

  // Response
  fobs = to_vector(f)[ii_obs];
  noise = s*sqrt((fobs+lambda));

  target += normal_lpdf(y | fobs,noise);

}
generated quantities {

}


