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
#ifndef MEASUREMENTMATRIX_H
#define MEASUREMENTMATRIX_H

#include "stdincludes.h"
#include "ShapeMatrix.h"
#include "MotionMatrix.h"
#include "Correspondence.h"
#include "ImageTransform.h"
#include "CorrespondenceSet.h"

//////////////////////////////////////////////////////////////////////////////
/// Class to represent the measurements of N 3D points from M views.
/// As defined in Sturm, P. and Triggs, B. 1996: A factorization based
/// algorithm for multi-image projective structure and
/// motion. European conference on computer vision, 1996: 709-720.
//////////////////////////////////////////////////////////////////////////////
template<int M>
class MeasurementMatrix : public Eigen::Matrix<double, 3*M, Eigen::Dynamic> {
public:
  typedef Eigen::Matrix<double, 3*M,Eigen::Dynamic> 		Base;
  typedef typename Eigen::Block<Base,3,1>			PixelRef;
  typedef const typename Eigen::Block<const Base,3,1>		constPixelRef;
  typedef typename Base::ColXpr					ColXpr;
  typedef typename Base::ColsBlockXpr				ColsBlockXpr;
  typedef typename Base::RowsBlockXpr				ViewRef;
  using Base::col;
  using Base::row;
  using Base::cols;
  using Base::rows;
  using Base::middleRows;
  using Base::middleCols;

  MeasurementMatrix()				 		{}
  MeasurementMatrix(const CorrespondenceSet &);
  int	add_correspondence(const Correspondence &);
  bool	load(const char *);
  void	normalise(ImageTransform<M> &, bool fixedAspect=false);
  void	unscale();
  void	scale_and_fill(const MotionMatrix<M> &, const ShapeMatrix &);
  void	synthesise_measurements(const MotionMatrix<M> &, int, double);
  void	synthesise_occlusions(double);
  void	order_cols_by_occlusions();
  
  PixelRef	pixel(int,int);
  constPixelRef pixel(int,int) const;
  bool		pixel_is_occluded(int,int) const;
  ViewRef	view(int);
  
  void 	insert_pixel(int, const GPixel &);
  void 	merge_columns(int, int);
  int  	new_column();
  int  	new_columns(int);
  int  	insert_columns(ColsBlockXpr);
  int  	insert_column_and_expand(ColXpr);
  void	delete_column(int);
  void	delete_columns(int,int);
  int	column_containing(const GPixel &);
  int	rank(double = 0.0);
  bool	is_fully_expanded(int);

  template<class D>
  MeasurementMatrix &operator=(const Eigen::MatrixBase<D> &);
  
};

#include "MeasurementMatrix.hpp"

#endif
