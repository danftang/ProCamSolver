///////////////////////////////////////////////////////////////////////////////
//
// Class for storing matrices. Specialises Matrix<> by providing constructors
// for various matrices used for epi-polar geometry.
//
// By default 3x3 identity matrix, but may also have other geometry
// see .cpp for constructors
//
///////////////////////////////////////////////////////////////////////////////
#ifndef TRANSFORM_MATRIX_H
#define TRANSFORM_MATRIX_H

#include "stdincludes.h"

typedef Eigen::Vector3d coord;

enum matrixType {
  intrinsic,			// Camera intrinsics matrix
  inverse_intrinsic,		// inverse of camera intrinsics
  transpose_inverse_intrinsic,	// transpose of inverse of camera intrinsics
  xrotate,			// clockwise rotation about x-axis (left/right)
  yrotate,			// rotation about y-axis
  zrotate,			// rotation about z-axis
  identity,
  cross_prod,			// cross product matrix Mv = m x v
  camera,			// camera matrix
  radial_distortion,		// 3x4 radial distortion for lifted coords
  transpose_radial_distortion,	// 3x4 radial distortion for lifted coords
  zero
};


class TransformMatrix : public Eigen::Matrix3d {
public:

  TransformMatrix(matrixType=identity, 
		  double=0.0, double=0.0, double=0.0, double=0.0);


  TransformMatrix & 	operator *=(TransformMatrix &);
  TransformMatrix & 	operator =(const Eigen::Matrix3d &);
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline TransformMatrix & TransformMatrix::operator *=(TransformMatrix &M) {
  Eigen::Matrix3d::operator *=(M);
  return(*this);
}

inline TransformMatrix & TransformMatrix::operator =(const Eigen::Matrix3d &M){
  Eigen::Matrix3d::operator =(M);
  return(*this);
}

#endif
