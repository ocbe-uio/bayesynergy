real VUS3(matrix[] Z, real[] x1, real[] x2, real[] x3){
  // Function calculating the "volume under the surface" metric for 3D dose-response functions
  int d[3] = dims(Z);
  int n1 = d[3];
  int n2 = d[2];
  int n3 = d[1];
  real rVUS = 0;
  vector[n3] slicesum;
  for (k in 1:n3){
    real matrixsum = 0;
    vector[n2] colsums = rep_vector(0,n2);
    for (j in 2:n2){
      // For each new row, need a variable summing up the trap rule
      real colsum = 0;
      for (i in 2:n3){
        colsum += (x1[i]-x1[(i-1)])*(Z[k,j,i]+Z[k,j,(i-1)])/2;
      }
      // And store this value for the jth row
      colsums[j] = colsum;
      // Next we integrate over all the colsums if j > 2
      // To get a single number for all k slices
      if (j > 2){
        matrixsum += (x2[j]-x2[(j-1)])*(colsums[j] + colsums[(j-1)])/2;
      }
    }
    slicesum[k] = matrixsum;
    // Then once slicesum starts getting built, we integrate over this as well
    if (k > 2){
      rVUS += (x3[k]-x3[(k-1)])*(slicesum[k] + slicesum[(k-1)])/2;
    }
  }
  // And this is the integral that's normalized and returned
  rVUS = 100 * rVUS / ((max(x1)-min(x1))*(max(x2)-min(x2))*(max(x3)-min(x3)));
  return rVUS;
}
