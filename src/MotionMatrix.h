///////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2012 Daniel Tang.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing,
//   software distributed under the License is distributed on an "AS
//   IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
//   express or implied.  See the License for the specific language
//   governing permissions and limitations under the License.
//
///////////////////////////////////////////////////////////////////////////////
#ifndef MOTIONMATRIX_H
#define MOTIONMATRIX_H

#include "stdincludes.h"
#include "ImageTransform.h"

template<int M> class MeasurementMatrix;

///////////////////////////////////////////////////////////////////////////////
/// Class to represent a set of projective camera matrices. The
/// matrices are placed one on top of the other into a 'projective
/// motion matrix' as described in Sturm, P. and Triggs, B. 1996: A
/// factorization based algorithm for multi-image projective structure
/// and motion. European conference on computer vision, 1996: 709-720.
///
/// The template parameter gives the number of cameras/projectors
/// represented in this matrix.
//////////////////////////////////////////////////////////////////////////////
template<int VIEWS>
class MotionMatrix : public Eigen::Matrix<double, 3*VIEWS, 4> {
public:
  typedef typename Eigen::Matrix<double,3*VIEWS,4>	Base;
  typedef typename Eigen::Block<Base,3,4>		Block34;
  typedef typename Eigen::Block<Base,3,3>		Block33;
  typedef typename Eigen::Block<Base,3,1>		Block31;
  typedef const typename Eigen::Block<const Base,3,4>	constBlock34;
  typedef const typename Eigen::Block<const Base,3,3>	constBlock33;
  typedef const typename Eigen::Block<const Base,3,1>	constBlock31;

  Block34		view(int);
  constBlock34		view(int) const;
  Block33		M(int);
  constBlock33		M(int) const;
  Block31		T(int);
  constBlock31		T(int) const;

  void			svd_solve(const MeasurementMatrix<VIEWS> &);
  double 		bundle_adjust(int, const MeasurementMatrix<VIEWS> &);

  Eigen::Matrix4d	diac_euclidean_lift(int);
  template<int F>
  Eigen::Matrix4d	diac_euclidean_lift(int, const ImageTransform<F> &);
  Eigen::Matrix4d	euclidean_lift(const Eigen::Vector3d &);
  Eigen::Vector3d	plane_at_infinity(int);
  double		KR_decompose(ImageTransform<VIEWS> &,
				     ImageTransform<VIEWS> &) const;
  double		KR_decompose(int, Eigen::Affine2d &, 
				          Eigen::AngleAxisd &) const;

  void			reprojection_err(const MeasurementMatrix<VIEWS> &,
					 Eigen::VectorXd &);

  template<class D>
  MotionMatrix &	operator=(const Eigen::MatrixBase<D> &);

protected:
  void			quadratic_prod(const Eigen::Vector4d &,
				       const Eigen::Vector4d &,
		       typename Eigen::Matrix<double,Eigen::Dynamic,10>::RowXpr);
  void			quadratic_prod(const Eigen::Vector3d &,
				       const Eigen::Vector3d &,
		       typename Eigen::Matrix<double,2*VIEWS+1,6>::RowXpr);
};

#include "MotionMatrix.hpp"

#endif
