#include <iostream>
#include <cstdlib>

#include "MeasurementMatrix.h"
#include "ShapeMatrix.h"
#include "MotionMatrix.h"

using namespace Eigen;

int main() {
  const int 			N = 10;			// number of cameras
  const int 			F = 0;			// number of fixed cams
  MotionMatrix<N>		reconstruction;		// camera matrices
  ShapeMatrix			shape;			// 3D-points
  ImageTransform<N>		Ks,Rs;			// camera matrices
  MeasurementMatrix<N>		correspondences;	// correspondences
  ImageTransform<F>		fixedIntrinsics;
  Affine2d			K;
  Matrix3d			R;

  // --- synthesise a vector of pixel correspondences
  // ------------------------------------------------
  Eigen::Matrix4d H;
  Eigen::VectorXd err;
  int i;

  try {

    correspondences.load("/home/daniel/data/data/correspondence-test");

    reconstruction.svd_solve(correspondences);

    std::cout << "SVD reconstructed camera matrix =" << std::endl;
    std::cout << reconstruction << std::endl;

    H = reconstruction.diac_euclidean_lift(3);
    reconstruction *= H;
    std::cout << "Projective transform =" << std::endl;
    std::cout << H << std::endl;

    reconstruction.reprojection_err(correspondences, err);
    std::cout << "SVD reprojection err = " << std::endl;
    std::cout << err.norm()/(correspondences.cols()*N*2.0*(1.0-0.08))<<std::endl;

    std::cout << "Reconstructed camera matrix =" << std::endl;
    std::cout << reconstruction << std::endl;
    reconstruction.KR_decompose(Ks,Rs);
    std::cout << "Intrinsics are" << std::endl;
    std::cout << Ks << std::endl;
    std::cout << "Rotations are" << std::endl;
    std::cout << Rs << std::endl;
  }
  catch(const char *err) {
    std::cout << "Caught error: " << err << std::endl;
  }
  /******
  reconstruction = cameras;
  reconstruction += MotionMatrix<N>::Random()*0.001;
  reconstruction.view(0) = cameras.view(0);

  std::cout << "Starting Bundle Adjustment" << std::endl;

  reconstruction.bundle_adjust(1, correspondences);

  shape.solve(correspondences, reconstruction);

  reconstruction.reprojection_err(correspondences, err);
  std::cout << "Adjusted reprojection err = " << std::endl;
  std::cout << err.norm()/(correspondences.cols()*N*2.0*(1.0-0.08))<<std::endl;

  std::cout << "Adjusted camera matrix =" << std::endl;
  std::cout << reconstruction << std::endl;
  skew = reconstruction.KR_decompose(Ks,Rs);
  std::cout << "Intrinsics are" << std::endl;
  std::cout << Ks << std::endl;
  std::cout << "Rotations are" << std::endl;
  std::cout << Rs << std::endl;
  *****/

  return(0);
}
