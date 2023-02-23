matrix[] one_minus_matrix_array(matrix[] A){
  // Elementwise max of matrix -- not yet in Stan
  int d[3] = dims(A);
  matrix[d[2],d[3]] Z[d[1]];
  for (k in 1:d[1]){
    for (j in 1:d[2]){
      for (i in 1:d[3]){
        Z[k,j,i] = 1-A[k,j,i];
      }
    }
  }
  return Z;
}
