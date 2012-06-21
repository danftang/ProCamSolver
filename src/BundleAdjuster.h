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
#ifndef BUNDLEADJUSTER_H
#define BUNDLEADJUSTER_H

#include "stdincludes.h"
#include "LevMarSolver.h"

template<int M> class MotionMatrix;
template<int M> class MeasurementMatrix;

///////////////////////////////////////////////////////////////////////////////
///
/// This class finds the MotionMatrix that gives the least-squares
/// reprojection error compared to a given MeasurementMatrix. It uses
/// the Levenberg-Marquardt method to perform bundle adjustment.
/// If any camera intrinsic matrices are already known, these can be
/// held fixed.
///
/// This uses Eigen's Levenberg-Marquardt module. Since this doesn't
/// account for sparse Jacobians it would be impractical to include
/// the shape matrix as an independent variable. However, the shape
/// matrix that minimises reprojection can be calculated very quickly
/// given the motion matrix, so here we take the Euclidean motion
/// matrix as the only ndependent variable and calculate the shape on the
/// fly.
///
/// This isn't as fast as a full sparse Lev-Mar implementation, but is
/// nevertheless a reasonably fast way to do full bundle adjustment
/// without the need for a sparse Lev-Mar solver.
///
///////////////////////////////////////////////////////////////////////////////
template<int M>
class BundleAdjuster : public LevMarSolver<BundleAdjuster<M> > {
public:
  typedef LevMarSolver<BundleAdjuster<M> >		Base;
  typedef Eigen::Matrix3d     	 			Matrix3d;
  typedef Eigen::Block<Eigen::VectorXd,3,1>		Block3;
  typedef Eigen::Block<Eigen::VectorXd,2,1>		Block2;
  typedef Eigen::Block<const Eigen::VectorXd,3,1>	constBlock3;
  typedef Eigen::Block<const Eigen::VectorXd,2,1>	constBlock2;

  using Base::levmar_solve;

  BundleAdjuster(const MeasurementMatrix<M> &, MotionMatrix<M> &);
  int				levmar_solve(int);
  constBlock3			translation(int,  const Eigen::VectorXd &);
  constBlock2			principal_pt(int, const Eigen::VectorXd &);
  const double &		focal_len(int,    const Eigen::VectorXd &);
  constBlock3			rotation(int,     const Eigen::VectorXd &);
  Block3			translation(int,  Eigen::VectorXd &);
  Block2			principal_pt(int, Eigen::VectorXd &);
  double &			focal_len(int,    Eigen::VectorXd &);
  Block3			rotation(int,     Eigen::VectorXd &);

  Eigen::AngleAxisd		rotation_to_angleaxis(constBlock3);
  Eigen::Vector3d		angleaxis_to_rotation(const Eigen::AngleAxisd &);

  int 				operator ()(const Eigen::VectorXd &,
					    Eigen::VectorXd &);

protected:
  void				params_to_matrices(const Eigen::VectorXd &);
  void				matrices_to_params();

protected:
  Eigen::VectorXd		variableParams;
  std::vector<Eigen::Affine2d>	fixedIntrinsics;
  const MeasurementMatrix<M> &	measurement;
  MotionMatrix<M> &		motion;
};

#include "BundleAdjuster.hpp"

#endif
