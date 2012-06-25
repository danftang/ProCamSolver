///////////////////////////////////////////////////////////////////////////////
#ifndef FUNDAMENTALMATRIX_H
#define FUNDAMENTALMATRIX_H

#include "stdincludes.h"
#include "ViewProjectionMatrix.h"

//////////////////////////////////////////////////////////////////////////////
/// Class to represent a projective camera matrix
//////////////////////////////////////////////////////////////////////////////
class FundamentalMatrix : public Eigen::Matrix3d {
public:
  FundamentalMatrix();
//  FundamentalMatrix(const ViewProjectionMatrix &,const ViewProjectionMatrix &);

//  void	set(const ViewProjectionMatrix &, const ViewProjectionMatrix &);  
  const Eigen::Vector3d & e_ij() 			{return(eij);}
  const Eigen::Vector3d & e_ji()			{return(eji);}
  void	calc_epipoles_and_constrain();

  template<class D>
  FundamentalMatrix &operator =(const Eigen::MatrixBase<D> &o);

protected:
  Eigen::Vector3d	eij;
  Eigen::Vector3d	eji;
};


//////////////////////////////////////////////////////////////////////////////
/// Construct fundamental matrix between two cameras with known (Euclidean)
/// projection matrices
//////////////////////////////////////////////////////////////////////////////
/***
inline FundamentalMatrix::FundamentalMatrix(const ViewProjectionMatrix &pi,
					    const ViewProjectionMatrix &pj)
{
  set(pi,pj);
}
****/

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline FundamentalMatrix::FundamentalMatrix()
{
}


///////////////////////////////////////////////////////////////////////////////
///
/// Calculates the epipoles of this matrix and enforces the epipolar
/// constraint. After execution the epipoles are stored in 'eij' and 'eji'.
///
///////////////////////////////////////////////////////////////////////////////
inline void FundamentalMatrix::calc_epipoles_and_constrain()
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


/****
//////////////////////////////////////////////////////////////////////////////
/// Set *this to the fundamental matrix between two cameras with known
/// (Euclidean) projection matrices
//////////////////////////////////////////////////////////////////////////////
inline void FundamentalMatrix::set(const ViewProjectionMatrix &pi,
			    const ViewProjectionMatrix &pj) {
  (*this) = (pj.M1T() * pj.R() * (pj.T() - pi.T()) * pi.R1() * pi.M1());
}
*****/

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
template<class D>
FundamentalMatrix &FundamentalMatrix::operator =(const Eigen::MatrixBase<D> &o) {
  Eigen::Matrix3d::operator =(o);
  return(*this);
}



#endif
