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
#include "FundamentalMatrix.h"


///////////////////////////////////////////////////////////////////////////////
/// Calculates the epipoles of this matrix and enforces the epipolar
/// constraint. After execution the epipoles are stored in 'eij' and 'eji'.
///////////////////////////////////////////////////////////////////////////////
void FundamentalMatrix::calc_epipoles_and_constrain()
{
  Eigen::JacobiSVD<Eigen::Matrix3d>  svd;

  svd.compute(*this, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::Vector3d D;
  D = svd.singularValues();
  D(2) = 0.0;
  (*this) = svd.matrixU()* D.asDiagonal() * svd.matrixV().transpose();
  eij = svd.matrixU().col(2);
  eji = svd.matrixV().col(2);
}


