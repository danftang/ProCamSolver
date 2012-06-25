#include <iostream>
#include "numerical.h"

int main() {
  MatDoub 	A(3,3);
  MatDoub	myNullsp;

  A[0][0] = 0.0; A[0][1] = 0.0; A[0][2] = 0.0;
  A[1][0] = 0.0; A[1][1] = 1.0; A[1][2] = 0.0;
  A[2][0] = 0.0; A[2][1] = 0.0; A[2][2] = 1.0;

  SVD	mySVD(A);

  A = mySVD.nullspace();

  int i;
  int j;

  for(i = 0; i<A.nrows(); ++i) {
    for(j=0; j<A.ncols(); ++j) {
      std::cout << A[i][j] << " ";
    }
    std::cout << std::endl;
  }

  return(0);
}
