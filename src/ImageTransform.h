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
#ifndef IMAGETRANSFORM_H
#define IMAGETRANSFORM_H

#include "stdincludes.h"

///////////////////////////////////////////////////////////////////////////////
/// Class to apply various transformations to Measurement Matrices and
/// Motion Matrices. The transformations are of the form X' = TX where
/// T is a sparse, RxR matrix and R is the number of rows of X.
///
/// T is zero everywhere except for in the 3x3 blocks whose top left corner
/// are at (3n,3n). This has the effect of applying each 3x3 transformation
/// to the respective views of X.
///////////////////////////////////////////////////////////////////////////////
template<int M>
class ImageTransform {
public:
  typedef Eigen::Matrix3d	TransformType;

  ImageTransform();
  ImageTransform(const ImageTransform &);

  TransformType &		view(int v);
  const TransformType &		view(int v) const;
  template<class D>  void	apply_to(Eigen::MatrixBase<D> &);
  template<class D>  void	apply_inverse_to(Eigen::MatrixBase<D> &);
  
protected:
  TransformType			views[M];
};


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
template<int M>
ImageTransform<M>::ImageTransform() {
  int i;
  for(i = 0; i<M; ++i) {
    view(i) = Eigen::Matrix3d::Identity();
  }
}

//////////////////////////////////////////////////////////////////////////////
/// returns a reference to the 'v''th 3x3 block.
//////////////////////////////////////////////////////////////////////////////
template<int M>
typename ImageTransform<M>::TransformType &
ImageTransform<M>::view(int v) {
  return(views[v]);
}

template<int M>
const typename ImageTransform<M>::TransformType &
ImageTransform<M>::view(int v) const {
  return(views[v]);
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
template<int M>
ImageTransform<M>::ImageTransform(const ImageTransform<M> &o) {
  int i;
  for(i = 0; i<M; ++i) {
    view(i) = o.view(i);
  }
}


//////////////////////////////////////////////////////////////////////////////
/// transforms X so that X -> TX. Multiplication syntax will be added
/// soon(ish)!
//////////////////////////////////////////////////////////////////////////////
template<int M>
template<class D>
void ImageTransform<M>::apply_to(Eigen::MatrixBase<D> &X) {
  int i;

  for(i = 0; i< M; ++i) {
    X.middleRows(3*i,3) = view(i) * X.middleRows(3*i,3);
  }
}


//////////////////////////////////////////////////////////////////////////////
/// transforms X so that TX -> X. i.e. X -> T^-1X
//////////////////////////////////////////////////////////////////////////////
template<int M>
template<class D>
void ImageTransform<M>::apply_inverse_to(Eigen::MatrixBase<D> &X) {
  int i;

  for(i = 0; i< M; ++i) {
    X.middleRows(3*i,3) = view(i).inverse() * X.middleRows(3*i,3);
  }
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
template<int M>
std::ostream &operator <<(std::ostream &out, const ImageTransform<M> &trans) {
  int i;
  for(i = 0; i<M; ++i) {
    out << trans.view(i) << std::endl;
  }
  return(out);
}


#endif

