real VUS(matrix Z, real[] x1, real[] x2){
  // Function calculating the "volume under the surface" metric for 2D dose-response functions
    int d[2] = dims(Z);
    int n1 = d[2];
    int n2 = d[1];
    real rVUS = 0;
    vector[n2] B_rVUS;
    for (i in 1:n2){
      real b_rVUS = 0;
      for (j in 2:n1){
        b_rVUS += (x1[j]-x1[(j-1)])*(Z[i,j]+Z[i,(j-1)])/2;
      }
      B_rVUS[i] = b_rVUS;
      if (i > 1){
        rVUS += (x2[i]-x2[(i-1)])*(B_rVUS[i]+B_rVUS[(i-1)]) / 2;
      }
    }
    rVUS = 100 * rVUS / ((max(x1)-min(x1))*(max(x2)-min(x2)));
    return rVUS;
  }
