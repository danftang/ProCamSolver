#include <iostream>
#include <cstdlib>

#include "MeasurementMatrix.h"
#include "ShapeMatrix.h"
#include "MotionMatrix.h"

int main() {
  const int 			M = 3;			// number of cameras
  const int 			N = 10;			// number of 3D points
  MotionMatrix<M>		cameras;		// camera matrices
  MeasurementMatrix<M>		correspondences;	// correspondences
  ViewProjectionMatrix		view;			// single view

  // --- Set up a scene with N cameras
  // ---------------------------------
  view.fx = 1.0;
  view.fy = 1.0;
  view.cx = 0.0;
  view.cy = 0.0;
  view.rot(0) = 0.0;
  view.rot(1) = 0.0;
  view.rot(2) = 0.0;
  view.t(0) = 0.0;
  view.t(1) = 0.0;
  view.t(2) = 0.0;
  cameras.view(0) = view.P();

  view.fx = 1.0;
  view.fy = 1.2;
  view.cx = 0.01;
  view.cy = 0.02;
  view.rot(0) = 0.0;
  view.rot(1) = 0.1;
  view.rot(2) = 0.0;
  view.t(0) = 1.0;
  view.t(1) = 0.0;
  view.t(2) = 0.0;
  cameras.view(1) = view.P();

  view.fx = 1.0123;
  view.fy = 1.0345;
  view.cx = 0.0678;
  view.cy = 0.0789;
  view.rot(0) = 0.0;
  view.rot(1) = 0.2;
  view.rot(2) = 0.0;
  view.t(0) = 2.0;
  view.t(1) = 1.0;
  view.t(2) = 0.0;
  cameras.view(2) = view.P();

  // --- synthesise a vector of pixel correspondences
  // ------------------------------------------------
  correspondences.synthesise_measurements(cameras, 100, 0.0);

  // --- create Fundamental matrices
  // -------------------------------------
  int i,j,p;
  FundamentalMatrix F;
  double err;

  for(i = 0; i<M; ++i) {
    for(j = 0; j<M; ++j) {
      if(i == j) continue;
      correspondences.eight_point_algorithm(F,i,j);
      std::cout << "F" << i << j << " = " << std::endl;
      std::cout << F << std::endl;
      std::cout << "e_ij =" << F.e_ij().transpose() << std::endl;
      std::cout << "e_ji =" << F.e_ji().transpose() << std::endl;
      // --- calculate error
      err = 0.0;
      for(p=0; p<N; ++p) {
	err += correspondences.pixel(i,p).transpose()*F*
	  correspondences.pixel(j,p);
      }
      std::cout << "err = " << err << std::endl;
      std::cout << "e_ij^TF = " << F.e_ij().transpose()*F << std::endl;
      std::cout << "Fe_ji = " << (F*F.e_ji()).transpose() << std::endl;
      std::cout << std::endl;
    }
  }

  return(0);
}
