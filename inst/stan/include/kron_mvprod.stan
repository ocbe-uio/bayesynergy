matrix kron_mvprod(matrix A, matrix B, matrix V) {
    // return (A kron_prod B) v where:
    // A is n1 x n1, B = n2 x n2, V = n2 x n1 = reshape(v,n2,n1)
    return transpose(A * transpose(B * V));
  }
