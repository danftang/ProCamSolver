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

#include <cmath>
#include <valarray>
#include <ostream>
#include "matrix.h"

typedef std::valarray<double> coord;

enum matrixType {
  intrinsic,			// Camera intrinsics matrix
  inverse_intrinsic,		// inverse of camera intrinsics
  xrotate,			// clockwise rotation about x-axis (left/right)
  yrotate,			// rotation about y-axis
  zrotate,			// rotation about z-axis
  identity,			
  cross_prod,			// cross product matrix Mv = m x v
  camera			// camera matrix
};


class transformMatrix : public Matrix<double> {
public:

  transformMatrix(matrixType=matrixType::identity, 
		  double=0.0, double=0.0, double=0.0, double=0.0);


  transformMatrix & 	operator *=(transformMatrix &);
  transformMatrix & 	operator =(const Matrix<double> &);
  transformMatrix & 	operator =(double);
};

std::ostream &operator <<(std::ostream &, const std::valarray<double> &);

inline transformMatrix & transformMatrix::operator *=(transformMatrix &M) {
  Matrix<double>::operator *=(M);
  return(*this);
}

inline transformMatrix & transformMatrix::operator =(const Matrix<double> &M){
  Matrix<double>::operator =(M);
  return(*this);
}

inline transformMatrix & transformMatrix::operator =(double c) {
  Matrix<double>::operator =(c);
  return(*this);
}


#endif
