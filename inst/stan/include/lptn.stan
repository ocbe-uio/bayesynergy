real lptn(real rho, real z, real lambda, real tau){
  // Returns the log-pdf of the log-pareto tailed normal distribution for robust regression
  if (fabs(z) <= tau){
    return std_normal_lpdf(z);
  }
  else {
    real logt = log(tau);
    real logl = log(lambda);
    // real logabsz = log(fabs(z));
    real logabsz = log(fmax(fabs(z),1.0001)); // Small hack to avoid numeric instability in next line
    // real lpdf = (-tau^2/2) + logt - logabsz + (lambda+1)*(log(logt)-log(logabsz));
    real lpdf = std_normal_lpdf(tau) + logt - logabsz + (lambda+1)*(log(logt)-log(logabsz));
    return lpdf;
  }
}
