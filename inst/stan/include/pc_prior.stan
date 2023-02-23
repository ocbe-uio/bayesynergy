real pc_prior(real ell, real sigmaf, real lambda1, real lambda2){
    // Custom lpdf for the PC prior of Matern covariance function
    real d_half = 2.0/2.0;
    real lprob = log(d_half) + log(lambda1) + log(lambda2) + (-d_half-1)*log(ell)-lambda1*pow(ell,-d_half)  - lambda2*sqrt(sigmaf) -log(sqrt(sigmaf));         ;
    return lprob;
  }
