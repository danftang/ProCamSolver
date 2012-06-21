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
#ifndef SHAPEMATRIX_H
#define SHAPEMATRIX_H

#include "stdincludes.h"

template<int M> class MeasurementMatrix;
template<int M> class MotionMatrix;

///////////////////////////////////////////////////////////////////////////////
/// Class to represent the homogeneous coordinates of a set of 3D points
///////////////////////////////////////////////////////////////////////////////
class ShapeMatrix : public Eigen::Matrix4Xd {
public:

  template<int V>
  void		solve(const MeasurementMatrix<V> &, const MotionMatrix<V> &);

  template<class D>
  ShapeMatrix &operator =(const Eigen::MatrixBase<D> &);
};

#include "ShapeMatrix.hpp"

#endif
