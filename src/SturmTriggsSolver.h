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
#ifndef STURMTRIGGSSOLVER_H
#define STURMTRIGGSSOLVER_H

#include "stdincludes.h"

class ShapeMatrix;
template<int M> class MeasurementMatrix;
template<int M> class MotionMatrix;

///////////////////////////////////////////////////////////////////////////////
/// Implementation of the algorithm described in Sturm and Triggs
/// (1996) which approximates the shape and motion matrices for any number
/// of cameras from a measurement matrix containing a set of pixel
/// correspondences.
///
/// This is a very fast solver, but the result does not minimise the
/// reprojection error, instead it minimises the Frobenius norm of the error
/// in the weighted measurement matrix. Also, the resulting motion matrix is
/// not guaranteed to be liftable to a Euclidean form.
///////////////////////////////////////////////////////////////////////////////
template<int M>
class SturmTriggsSolver {
public:
  SturmTriggsSolver(MotionMatrix<M> &, ShapeMatrix &);

  void	solve(const MeasurementMatrix<M> &);

  void	approx_scale();
  void  factorise_with_occlusions(MotionMatrix<M> &, ShapeMatrix &);
  void	create_connectivity_graph(Eigen::Matrix<double,M,M> &);
  void  form_subspace(std::vector<int>::iterator, MeasurementMatrix<M> &);

  bool operator ()(int, int);

protected:
  void	random_permute(std::vector<int> &);
  void	set_measurement_matrix(const MeasurementMatrix<M> &);

protected:
  MeasurementMatrix<M> 		measurement;
  MotionMatrix<M> & 		motion;
  ShapeMatrix &			shape;
  Eigen::Matrix<int,M,1>	spanning_tree;      // tree of views
  Eigen::Matrix<int,M,1>	tree_traversal;      // tree of views
};

#include "SturmTriggsSolver.hpp"

#endif
