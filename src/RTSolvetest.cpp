#include "RTSolver.cpp"
#include "transformMatrix.cpp"

int main() {
  transformMatrix F;
  transformMatrix M(intrinsic, 1.0, 1.0, 0.0, 0.0);
  transformMatrix MP(intrinsic, 1.0, 1.2, 0.01, 0.02);
  transformMatrix R(yrotate, atan(0.5));
  transformMatrix NPT(inverse_intrinsic, 1.0,1.2,0.01,0.02);
  transformMatrix N(inverse_intrinsic, 1.0,1.0,0.0,0.0);
  transformMatrix T(cross_prod,1.0,0.0,0.0);
  coord		  t(3);
  transformMatrix I(identity);
  transformMatrix P;

  NPT.transpose();
  t[0] = -1.0; t[1] = 0.0; t[2] = 0.0;

  F[0][0] = 0.0; F[0][1] = 0.44721; F[0][2] = 0.0;
  F[1][0] = 0.0; F[1][1] = 0.0;	    F[1][2] = -0.83333;
  F[2][0] = 0.0; F[2][1] = 0.88996; F[2][2] = 0.016667;

  P = RTSolver(F,M,MP);

  std::cout << MP*R*(I|t);
  std::cout << "---" << std::endl;
  std::cout << P;
  std::cout << "---" << std::endl;
  std::cout << P - MP*R*(I|t);
  return(0);
}
