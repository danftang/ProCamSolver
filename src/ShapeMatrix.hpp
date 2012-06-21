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
#include "MotionMatrix.h"
#include "MeasurementMatrix.h"


///////////////////////////////////////////////////////////////////////////////
/// Sets *this to the shape that minimises the reprojection error for a
/// given motion matrix and (un-scaled, occluded) measurement matrix.
/// Works in two stages. For each column in the measurement matrix:
/// 1) Approximate the projective depths by least squares fitting 
///    (P_n0:1 - pi_n0:1P_n2)X = 0 then setting lambda_n = P_n2X
/// 2) Do least squares fit of lambda_n^-1(P_n0:1 - pi_n0:1P_n2)X = 0, for
///    0 <= n < M
///
/// The result is not the exact minimum reprojection since the initial
/// guess at lambda does not minimise error. The result is that the
/// errors of different views will have slightly differing weights. However,
/// the difference will be very small indeed.
///
///////////////////////////////////////////////////////////////////////////////
template<int V>
void ShapeMatrix::solve(const MeasurementMatrix<V> &M, const MotionMatrix<V> &P) {
  double			lambda;	// projective depth
  double			px, py;
  int 				i,j;
  int				visiblePixels;
  Eigen::Matrix<double,2*V,3>	coeffs;	// matrix for least squares solution
  Eigen::Matrix<double,2*V,1>	y;	// values for least squares solution

  resize(4,M.cols());
  row(3).setOnes();
  
  for(j = 0; j<cols(); ++j) {
    // --- Do linear least squares.
    // --- For each non-occluded pixel, pi, constraint is
    // --- l^-1(M.row(0:1) - pi.row(0:1)M.row(2))\chi
    // --- = -l^-1(T.row(0:1) - pi.row(0:1)*T.row(2))
    // --- Where P = [M|T] and X = [chi|1]^T
    // ---------------------------------------
    visiblePixels = 0;
    for(i=0; i<V; ++i) {
      if(!M.pixel_is_occluded(i,j)) {
	px = M.pixel(i,j)(0) / M.pixel(i,j)(2);
	py = M.pixel(i,j)(1) / M.pixel(i,j)(2);
	coeffs.row(2*i)   = P.M(i).row(0) - px * P.M(i).row(2); 
	coeffs.row(2*i+1) = P.M(i).row(1) - py * P.M(i).row(2);
	y(2*i)   = px * P.T(i)(2) - P.T(i)(0);
	y(2*i+1) = py * P.T(i)(2) - P.T(i)(1);
	++visiblePixels;
      } else {
	coeffs.middleRows<2>(2*i).setZero();
	y.middleRows<2>(2*i).setZero();
      }
    }
    

    if(visiblePixels < 2) {
      // --- not enough pixels to constrain X
      // --- so just choose a null vector
      block<3,1>(0,j) =
	coeffs.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(y);
      continue;
    }

    // --- first estimate lambda
    // --- and use to set weights

    block<3,1>(0,j) =
      (coeffs.transpose() * coeffs).ldlt().solve(coeffs.transpose()*y);

    //std::cout << "Approx X = " << col(j) << std::endl;
    //std::cout << "Error = " << std::endl << coeffs * block<3,1>(0,j) - y << std::endl;

    for(i=0; i<V; ++i) {
      lambda = P.view(i).row(2)*col(j);
      coeffs.middleRows<2>(2*i) /= lambda;
      y.middleRows<2>(2*i)   	/= lambda;
    }

    // --- now re-solve with weights
    /****
    block<3,1>(0,j) = 
      coeffs.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(y);
    ****/
    block<3,1>(0,j) =
      (coeffs.transpose() * coeffs).ldlt().solve(coeffs.transpose()*y);

    //std::cout << "Min X = " << col(j) << std::endl;
    //std::cout << "Shape/Motion error = " << std::endl << coeffs * block<3,1>(0,j) - y << std::endl;
  }
  //std::cout << "Reconstruction = " << std::endl;
  //std::cout << P*(*this) << std::endl;
  //std::cout << "Measurement = " << std::endl;
  //std::cout << M << std::endl;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
template<class D>
ShapeMatrix &ShapeMatrix::operator =(const Eigen::MatrixBase<D> &o) {
  Eigen::Matrix4Xd::operator =(o);
  return(*this);
}

