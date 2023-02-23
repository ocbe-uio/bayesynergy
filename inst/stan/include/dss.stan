real dss(real c1,real c2,real la,real slope,real log10_ec50){
  real val = (c2-c1)+(la-1)/slope*(log10(1+10^(slope*(c2-log10_ec50))) - log10(1+10^(slope*(c1-log10_ec50))));
  val = 100 * (1-val/(c2-c1)); // Normalizing
  return val;
}
