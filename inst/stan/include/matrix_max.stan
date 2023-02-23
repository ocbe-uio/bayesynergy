matrix matrix_max(matrix A, real a){
    // Elementwise max of matrix -- not yet in Stan
    int d[2] = dims(A);
    matrix[d[1],d[2]] Z;
    for (i in 1:d[1]){
      for (j in 1:d[2]){
        Z[i,j] = fmax(A[i,j],a);
      }
    }
    return Z;
  }
