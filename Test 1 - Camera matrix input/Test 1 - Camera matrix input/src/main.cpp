///////////////////////////////////////////////////////////////////////////////
//
// Test of Fundamental matrix solver and camera matrix extractor
//
///////////////////////////////////////////////////////////////////////////////

#include "FSolver.h"
#include "RTSolver.h"
#include "transformMatrix.h"

int main() {
  transformMatrix 	F;
  transformMatrix 	M(intrinsic, 1.0, 1.0, 0.0, 0.0);
  transformMatrix 	MP(intrinsic, 1.0, 1.2, 0.01, 0.02);
  transformMatrix 	P;
  FSolver		myFSolver;
  coord			zero(3);

  // --- Load correspondences and solve for the fundamental matrix
  // --- putting the result in F
  // ---------------------------
  myFSolver.load("data.dat");
  myFSolver.solve(F);

  std::cout << "Fundamental matrix is:" << std::endl;
  std::cout << F;

  // --- now extract the camera matrix for given fundamental martix (F),
  // --- and given intrinsic matrices for the reference camera (M)
  // --- and the other camera (MP)
  // ------------------------------------------------------
  P = RTSolver(F,M,MP);

  std::cout << std::endl;
  std::cout << "Camera matrix is:" << std::endl;
  std::cout << P;
  std::cout << std::endl;
  std::cout << "Reference camera matrix is:" << std::endl;
  std::cout << (M|zero);

  return(0);
}
