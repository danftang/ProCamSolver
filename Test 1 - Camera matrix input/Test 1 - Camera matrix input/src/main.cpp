///////////////////////////////////////////////////////////////////////////////
//
// Test of Fundamental matrix solver and camera matrix extractor
//
///////////////////////////////////////////////////////////////////////////////

#include "FSolver.h"
#include "RTSolver.h"
#include "transformMatrix.h"

int main() {
#if defined(_MSC_VER)
	//do vc++ related
	transformMatrix	I(matrixType::identity);
#else
	//do gcc related
	transformMatrix	I(identity);
#endif
	coord		t(3);			// translation vector
	t[0] = 1.0f; t[1] = 0.0f; t[2] = 0.0f;
	transformMatrix referenceRotate(yrotate, atan(0.5));
	Matrix<double> referenceTranslate = I|t;

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
  std::cout << "ViewProjetion matrix estimation:" << std::endl;
  std::cout << P;
  std::cout << std::endl;

  std::cout << "Reference Projection matrix is:" << std::endl;
  std::cout << (MP|zero);
  std::cout << std::endl;

  std::cout << "Reference Rotate matrix is:" << std::endl;
  std::cout << referenceRotate;
  std::cout << std::endl;

  std::cout << "Reference Translate matrix is:" << std::endl;
  std::cout << referenceTranslate;
  std::cout << std::endl;

  std::cout << "Reference ViewProjection matrix is:" << std::endl;
  std::cout << (MP * referenceRotate * referenceTranslate);
  std::cout << std::endl;


  return(0);
}
