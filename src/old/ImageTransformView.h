#ifndef IMAGETRANSFORMVIEW_H
#define IMAGETRANSFORMVIEW_H

#include "stdincludes.h"

///////////////////////////////////////////////////////////////////////////////
/// Class to apply various normalisations to Measurement Matrices and
/// Motion Matrices by left multiplication
//////////////////////////////////////////////////////////////////////////////

enum transformType {
  intrinsic,                    // Camera intrinsics matrix
  inverse_intrinsic,            // inverse of camera intrinsics
  transpose_inverse_intrinsic,  // transpose of inverse of camera intrinsics
  xrotate,                      // clockwise rotation about x-axis
  yrotate,                      // rotation about y-axis
  zrotate,                      // rotation about z-axis
  xyzrotate,
  cross_prod	                // cross product matrix Mv = m x v
};


template<class CONTAINER>
class ImageTransformView : public Eigen::Block<CONTAINER,3,3> {
public:
  typedef Eigen::Block<CONTAINER,3,3>	Base;

  ImageTransformView(const Base &);

  void	set(transformType, const double &, const double &, const double &);

  template<class D>
  ImageTransformView &operator =(const Eigen::DenseBase<D> &);
};


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
template<class CONTAINER>
ImageTransformView<CONTAINER>::ImageTransformView(const Base &o) : Base(o) {  
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
template<class CONTAINER>
template<class DERIVED>
ImageTransformView<CONTAINER> ImageTransformView<CONTAINER>::operator =(const Eigen::DenseBase<DERIVED> &o) {
  Base::operator=(o);
  return(*this);
}


//////////////////////////////////////////////////////////////////////////////
///
//////////////////////////////////////////////////////////////////////////////
template<int M>
template<int J>
inline void ImageTransform<M>::set(transformType t,
				   const double &a,
				   const double &b,
				   const double &c) {
  fill(0.0);
 
  switch(t) {
  case intrinsic:
    (*this)(0,0) = a;
    (*this)(1,1) = a;
    (*this)(0,2) = b;
    (*this)(1,2) = c;
    (*this)(2,2) = 1.0;
    break;
  case inverse_intrinsic:
    (*this)(0,0) = 1.0/a;
    (*this)(1,1) = 1.0/a;
    (*this)(0,2) = -b/a;
    (*this)(1,2) = -c/a;
    (*this)(2,2) = 1.0;
    break;
  case transpose_inverse_intrinsic:
    (*this)(0,0) = 1.0/a;
    (*this)(1,1) = 1.0/a;
    (*this)(2,0) = -b/a;
    (*this)(2,1) = -c/a;
    (*this)(2,2) = 1.0;
    break;
  case xrotate:
    (*this)(2,2) = (*this)(1,1) = std::cos(a);
    (*this)(2,1) = std::sin(a);
    (*this)(1,2) = -(*this)(2,1);
    (*this)(0,0) = 1.0;
    break;
  case yrotate:
    (*this)(2,2) = (*this)(0,0) = std::cos(a);
    (*this)(0,2) = std::sin(a);
    (*this)(2,0) = -(*this)(0,2);
    (*this)(1,1) = 1.0;
  break;
  case zrotate:
    (*this)(1,1) = (*this)(0,0) = std::cos(a);
    (*this)(1,0) = std::sin(a);
    (*this)(0,1) = -(*this)(1,0);
    (*this)(2,2) = 1.0;
    break;
  case xyzrotate:
    (*this) =     
      Eigen::AngleAxisd(a,Vector3d::UnitX())*
      Eigen::AngleAxisd(b,Vector3d::UnitY())*
      Eigen::AngleAxisd(c,Vector3d::UnitZ());
    break;
  case cross_prod:
    (*this)(0,0) = 0.0; (*this)(0,1) = -c;  (*this)(0,2) = b;
    (*this)(1,0) = c;   (*this)(1,1) = 0.0; (*this)(1,2) = -a;
    (*this)(2,0) = -b;  (*this)(2,1) = a;   (*this)(2,2) = 0.0;
    break;
  case zero:
    (*this)(0,0) = 0.0;
    (*this)(1,1) = 0.0;
    (*this)(2,2) = 0.0;
    break;
  }
}


#endif

