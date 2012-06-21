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
#ifndef FUNDAMENTALMATRIX_H
#define FUNDAMENTALMATRIX_H

#include "stdincludes.h"

//////////////////////////////////////////////////////////////////////////////
/// Class to represent the fundamental matrix between two views
//////////////////////////////////////////////////////////////////////////////
class FundamentalMatrix : public Eigen::Matrix3d {
public:
  FundamentalMatrix();

  template<class DERIVED1, class DERIVED2>
  void  eight_point_algorithm(const Eigen::MatrixBase<DERIVED1> &,
			      const Eigen::MatrixBase<DERIVED2> &);
  void	calc_epipoles_and_constrain();
  const Eigen::Vector3d & e_ij() 			{return(eij);}
  const Eigen::Vector3d & e_ji()			{return(eji);}

  template<class D>
  FundamentalMatrix &operator =(const Eigen::MatrixBase<D> &o);

protected:
  Eigen::Vector3d	eij;
  Eigen::Vector3d	eji;
};


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline FundamentalMatrix::FundamentalMatrix()
{
}

///////////////////////////////////////////////////////////////////////////////
/// Sets '*this' fundamental matrix from a set of corresponding pixels
/// from two separate views using the 8-point algorithm. The inputs
/// 'a' and 'b' are 3xN matrices whose columns contain pixel
/// coordinates in the form (s.x,s.y,s)^T where s is an arbitrary
/// scale. Corresponding columns in 'a' and 'b' represent
/// corresponding pixels in the two views. The number of
/// correspondences must be at least 8 but may be more. Pixels with a scale
/// of zero are assumed to be occluded and so are ignored.
///
/// In order to ensure stability, the pixel coordinates should have
/// previously been normalised to have zero average and unit standard
/// deviation. For more details, see Hartley, R.I. 1997: In defense of
/// the eight-point algorithm. IEEE Trans.
/// Patt. Anal. Mach. Int. 19:580-593.
///
/// \param a the pixels from the first camera in the form of a 3xN matrix
/// \param b the pixels from the second camera in the form of a 3xN matrix
///
///////////////////////////////////////////////////////////////////////////////
template<class DERIVED1, class DERIVED2>
void FundamentalMatrix::
eight_point_algorithm(const Eigen::MatrixBase<DERIVED1> &a,
		      const Eigen::MatrixBase<DERIVED2> &b){
  // --- A is the lifted matrix of correspondences
  Eigen::Matrix<double,Eigen::Dynamic,9>			A;
  Eigen::Matrix<double,9,9>					ATA;
  Eigen::JacobiSVD<Eigen::Matrix<double,Eigen::Dynamic,9> >  	Adecomposition;
  int 	point;
  int 	validCorrespondences;
  int 	totalCorrespondences;	
  int 	n;

  // --- Sanity checks
  if(a.rows() != 3 || b.rows() != 3) 
    throw("Inputs to eight_point_algorithm must have 3 rows");

  // --- form normalised correspondence matrix
  // -----------------------------------------
  validCorrespondences = 0;
  totalCorrespondences = (a.cols()<b.cols()?a.cols():b.cols());
  for(point = 0; point<totalCorrespondences; ++point) {
    if(a(2,point) != 0.0 && b(2,point) != 0.0) {
      ++validCorrespondences;
    }
  }
  A.resize(validCorrespondences,9);
  n = 0;
  for(point = 0; point<totalCorrespondences; ++point) {
    if(a(2,point) != 0.0 && b(2,point) != 0.0) {
      A(n,0) = b(0,point)*a(0,point);
      A(n,1) = b(0,point)*a(1,point);
      A(n,2) = b(0,point);
      A(n,3) = b(1,point)*a(0,point);
      A(n,4) = b(1,point)*a(1,point);
      A(n,5) = b(1,point);
      A(n,6) = a(0,point);
      A(n,7) = a(1,point);
      A(n,8) = 1.0;
      ++n;
    }
  }

  // --- Approximately solve AF = 0 by SVD
  // -------------------------------------
  ATA = A.transpose()*A;
  Adecomposition.compute(ATA, Eigen::ComputeFullV);
  col(0) = Adecomposition.matrixV().block(0,8,3,1);
  col(1) = Adecomposition.matrixV().block(3,8,3,1);
  col(2) = Adecomposition.matrixV().block(6,8,3,1);

  // --- Enforce zero determinant of F
  // ---------------------------------
  calc_epipoles_and_constrain();

  std::cout << "Residual in eight_point_algorithm = " << std::endl;
  std::cout << (a.transpose() * (*this) * b).diagonal() << std::endl;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
template<class D>
FundamentalMatrix &FundamentalMatrix::operator =(const Eigen::MatrixBase<D> &o) {
  Eigen::Matrix3d::operator =(o);
  return(*this);
}


#endif
