#include <iostream>
#include <cstdlib>

#include "MeasurementMatrix.h"
#include "ShapeMatrix.h"
#include "MotionMatrix.h"

using namespace Eigen;

int main() {
  const int 			N = 5;			// number of cameras
  const int 			F = 0;			// number of fixed cams
  MotionMatrix<N>		cameras;		// camera matrices
  MotionMatrix<N>		reconstruction;		// camera matrices
  ShapeMatrix			shape;			// 3D-points
  ImageTransform<N>		Ks,Rs;			// camera matrices
  MeasurementMatrix<N>		correspondences;	// correspondences
  double			skew;
  ImageTransform<F>		fixedIntrinsics;
  Affine2d			K;
  Matrix3d			R;

  // --- Set up a scene with N cameras
  // ---------------------------------
  K = Scaling(1.234);
  cameras.M(0) = 		K.matrix();
  cameras.T(0) = 		Vector3d(0.0, 0.0, 0.0);
  //fixedIntrinsics.view(0) = 	K.matrix();

  K = Scaling(1.2);
  R = AngleAxisd(0.1,Vector3d::UnitY());
  cameras.M(1) = 		K * R;
  cameras.T(1) = 		Vector3d(1.0, 0.0, 0.0);
  //fixedIntrinsics.view(1) = 	K.matrix();

  K = Scaling(1.0345);
  R = AngleAxisd(0.2,Vector3d::UnitY());
  cameras.M(2) = 		K * R;
  cameras.T(2) = 		Vector3d(2.0, 1.0, -0.1);

  K = Translation2d(0,0.25) * Scaling(1.0123);
  R = AngleAxisd(0.1,Vector3d::UnitX()) *
      AngleAxisd(0.1,Vector3d::UnitY()) * 
      AngleAxisd(0.1,Vector3d::UnitZ());
  cameras.M(3) = 		K * R;
  cameras.T(3) = 		Vector3d(3.0, -1.0, 0.1);

  K = Translation2d(0,0.25) * Scaling(1.0123);
  R = AngleAxisd(-0.1,Vector3d::UnitX()) *
      AngleAxisd(0.1,Vector3d::UnitY()) * 
      AngleAxisd(-0.1,Vector3d::UnitZ());
  cameras.M(4) = 		K * R;
  cameras.T(4) = 		Vector3d(4.0, -1.0, -0.1);

  std::cout << "Original MotionMatrix = " << std::endl;
  std::cout << cameras << std::endl;

  // --- synthesise a vector of pixel correspondences
  // ------------------------------------------------
  Eigen::Matrix4d H;
  Eigen::VectorXd err;
  int i;

  correspondences.synthesise_measurements(cameras, 2000, 0.001);
  correspondences.synthesise_occlusions(0.25);

  std::cout << "Original Correspondence norm = " << std::endl;
  // std::cout << correspondences << std::endl;
  std::cout << correspondences.rowwise().norm() << std::endl;

  reconstruction.svd_solve(correspondences);

  std::cout << "SVD reconstructed camera matrix =" << std::endl;
  std::cout << reconstruction << std::endl;

  //Eigen::Vector3d p;
  //p = reconstruction.plane_at_infinity(3);
  //H = reconstruction.euclidean_lift(p);
  H = reconstruction.diac_euclidean_lift(3);
  reconstruction *= H;
  std::cout << "Projective transform =" << std::endl;
  std::cout << H << std::endl;

  reconstruction.reprojection_err(correspondences, err);
  std::cout << "SVD reprojection err = " << std::endl;
  std::cout << err.norm()/sqrt(err.rows())<<std::endl;

  std::cout << "Reconstructed camera matrix =" << std::endl;
  std::cout << reconstruction << std::endl;
  reconstruction.KR_decompose(Ks,Rs);
  std::cout << "Intrinsics are" << std::endl;
  std::cout << Ks << std::endl;
  std::cout << "Rotations are" << std::endl;
  std::cout << Rs << std::endl;

  //reconstruction = cameras;
  //reconstruction += MotionMatrix<N>::Random()*0.001;
  //reconstruction.view(0) = cameras.view(0);

  std::cout << "Starting Bundle Adjustment" << std::endl;

  reconstruction.bundle_adjust(0, correspondences);

  reconstruction.reprojection_err(correspondences, err);
  std::cout << "Adjusted reprojection err = " << std::endl;
  std::cout << err.norm()/sqrt(err.rows())<<std::endl;

  std::cout << "Adjusted camera matrix =" << std::endl;
  std::cout << reconstruction << std::endl;
  skew = reconstruction.KR_decompose(Ks,Rs);
  std::cout << "Intrinsics are" << std::endl;
  std::cout << Ks << std::endl;
  std::cout << "Rotations are" << std::endl;
  std::cout << Rs << std::endl;

  return(0);
}
