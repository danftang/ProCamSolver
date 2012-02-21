///////////////////////////////////////////////////////////////////////////////
//
// Extract the camera matrix from a fundamental matrix for given intrinsic
// matrices of both cameras, on the assumption that the camera matrix of
// the reference camera is of the form P=M[I|0] and the other camera matrix
// is of the form P'=M'R[I|-t]
// where M and M' are intrinsic matrices of the form
//     (f  0  cx )
//     (0  af cy )
//     (0  0  1  )
// R = rotation matrix
// I = identity
// t = translation of centre of projection
//
// The constructor takes a triple (F,M,M') where
// F  = the fundamental matrix
// M  = the intrinsic matrix of the reference camera
// M' = the intrinsic matrix of the second camera
//
// on construction, the solution is calculated and 'this' is set to the
// camera matrix of the second camera.
//
///////////////////////////////////////////////////////////////////////////////
#ifndef RTSOLVER_H
#define RTSOLVER_H

#include <valarray>
#include <vector>
#include "nrlib/nr3.h"
#include "transformMatrix.h"

class RTSolver : public transformMatrix {
public:
  RTSolver(Matrix<double> &, Matrix<double> &, Matrix<double> &);

  double	operator()(VecDoub_I &);
  void		df(VecDoub_I &, VecDoub_O &);

protected:
  double elemental_dot_prod(const Matrix<double> &, const Matrix<double> &);

protected:
  transformMatrix	F;	// Fundamental matrix
  transformMatrix	N;	// Inverse reference intrinsics M^{-1}
  transformMatrix	NPT;	// Transpose inverse intrinsics primed M'^{-1T}
  transformMatrix	R;	// temporary rotation matrix

  static const double		PI;
};

#endif
