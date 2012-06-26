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


///////////////////////////////////////////////////////////////////////////////
/// Loads a set of pixel correspondences from a sorted file in the format
///  <id1> <id2> <x1> <y1> <x2> <y2>
///
/// where
///
/// id1/2      = id of the projector / camera
/// x1/2, y1/2 = pixel coordinate of the projector / camera
///
///////////////////////////////////////////////////////////////////////////////
template<int M>
bool MeasurementMatrix<M>::load(const char *filename) {
  std::ifstream 	myFile(filename);
  Correspondence 	c;
  GPixel		currentProjPixel;
  int			points;
  int			currentCol;
  const double		p = 0.002; // probability of loading point
  const int		cameras = 5;

  if(!myFile) return(false);
  std::cout << "Loading correspondences" << std::endl;

  // --- count number of 3D points
  currentProjPixel.id = -1;
  points = 0;
  while(myFile >> c) {
    if(c.i != currentProjPixel) {
      currentProjPixel = c.i;
      ++points;
    }
  }
  std::cout <<"Found "<< points <<" points" << std::endl;
  Base::resize(3*M, points);

  // ### scan back to beginning of file ###
  myFile.clear();
  myFile.seekg(0, std::ios::beg);

  // --- form columns one by one
  // --- assuming file is ordered
  // --- by projector pixels.
  currentProjPixel.id = -1;
  currentCol = 0;
  myFile >> c;
  while(myFile) {
    if((rand()*1.0/RAND_MAX) < p) {
      currentProjPixel = c.i;

      c.i.id += cameras; // test

      pixel(c.i.id, currentCol)(0) = c.i.x;
      pixel(c.i.id, currentCol)(1) = c.i.y;
      pixel(c.i.id, currentCol)(2) = 1.0;

      //      std::cout << currentProjPixel << std::endl;

      do {
	c.j.y *= 0.75; 	// test
	c.j.id -= 1;    // test
	
	pixel(c.j.id, currentCol)(0) = c.j.x;
	pixel(c.j.id, currentCol)(1) = c.j.y;
	pixel(c.j.id, currentCol)(2) = 1.0;

	myFile >> c;
      } while(myFile && c.i == currentProjPixel);
      ++currentCol;
    } else {
      // skip 3D point
      while(myFile >> c && c.i == currentProjPixel) {;}
    }
  }
  Base::conservativeResize(3*M, currentCol);
  std::cout <<"Loaded "<< currentCol <<" points" << std::endl;
  myFile.close();
  return(true);
}


///////////////////////////////////////////////////////////////////////////////
/// Normalises *this by setting the scales of all un-occluded pixels
/// to 1.0, offsetting the pixels so that the mean of pixels in
/// each view is (0,0) and scaling the (x,y) values so that the
/// standard deviation in each view is 1. If 'fixedAspect' is false
/// (the default) then the standard deviation in both x and y
/// directions is set to 1, otherwise the sum of x and y variances is 1
///
/// \param denormalisation This is set to the inverse transform 
/// for subsequent de-normalisation.
///
/// \param fixedAspect If true, aspect ratio is held fixed (default is false).
///////////////////////////////////////////////////////////////////////////////
template<int M>
void MeasurementMatrix<M>::normalise(ImageTransform<M> &denormalisation,
				     bool fixedAspect) {
  Eigen::Vector3d mean;
  Eigen::Vector2d sd;
  Eigen::Matrix3d K;
  int v;
  int point;

  unscale();
  K.setZero();
  for(v = 0; v < M; ++v) {
    mean = view(v).rowwise().sum();
    sd   = view(v).topRows(2).rowwise().squaredNorm();
    if(mean(2) == 0.0) mean(2) = 1.0; // no data for view
    sd   /= mean(2);
    mean /= mean(2);
    sd = (sd - mean.cwiseProduct(mean).topRows(2));
    if(fixedAspect) {
      sd(0) = sd(1) = (sd(0) + sd(1))*0.5;
    }
    sd = sd.cwiseSqrt();

    if(sd(0) == 0.0) sd(0) = 1.0;
    if(sd(1) == 0.0) sd(1) = 1.0;

    K(0,0) = 1.0/sd(0);
    K(1,1) = 1.0/sd(1);
    K(0,2) = -mean(0)/sd(0);
    K(1,2) = -mean(1)/sd(1);
    K(2,2) = 1.0;    
    view(v) = K*view(v);

    // --- form inverse
    K(0,0) = sd(0);
    K(1,1) = sd(1);
    K(0,2) = mean(0);
    K(1,2) = mean(1);
    denormalisation.view(v) = K;
  }
}


///////////////////////////////////////////////////////////////////////////////
/// Scales the pixels of *this and fills in occluded pixels
/// based on the given motion matrix and shape matrix. Scale is calculated
/// as lambda = P.row(2)X
///////////////////////////////////////////////////////////////////////////////
template<int M>
void MeasurementMatrix<M>::scale_and_fill(const MotionMatrix<M> &P, 
					  const ShapeMatrix &X) {
  int i,j;

  for(j = 0; j<cols(); ++j) {
    for(i = 0; i<M; ++i) {
      if(pixel_is_occluded(i,j)) {
	pixel(i,j) = P.view(i) * X.col(j);
      } else {
	pixel(i,j) *= (P.view(i).row(2) * X.col(j))/pixel(i,j)(2);
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
/// Adds a pixel correspondence to '*this'
///
/// \return the index of the column that the correspondence was added to 
///////////////////////////////////////////////////////////////////////////////
template<int M>
int MeasurementMatrix<M>::add_correspondence(const Correspondence &corr) {
  int icolumn = column_containing(corr.i);
  int jcolumn = column_containing(corr.j);

  if(icolumn < 0) {
    if(jcolumn < 0) {                    // neither exist yet
      icolumn = new_column();
      insert_pixel(icolumn, corr.i);
      insert_pixel(icolumn, corr.j);      
    } else {                            // j exists, i doesn't
      insert_pixel(jcolumn, corr.i);
      return(jcolumn);
    }
  } else {
    if(jcolumn < 0) {                    // i exists, j doesn't
      insert_pixel(icolumn, corr.j);
    } else {                            // both exist
      merge_columns(icolumn, jcolumn);
    }
  }
  return(icolumn);

}


///////////////////////////////////////////////////////////////////////////////
/// Returns the index of the column contining the given global pixel, G, or
/// -1 if no such column exists
///
/// To avoid problems with rounding error, pixel coordinates that are closer
/// than 10^-5 are treated as equal
///////////////////////////////////////////////////////////////////////////////
template<int M>
int MeasurementMatrix<M>::column_containing(const GPixel &G) {
  int j;
  for(j = 0; j<cols(); ++j) {
    if(fabs(pixel(G.id, j)(0) - G.x) < 1e-5 &&
       fabs(pixel(G.id, j)(1) - G.y) < 1e-5) {
      return(j);
    }
  }
  return(-1);
}


///////////////////////////////////////////////////////////////////////////////
/// Returns all pixel measurements from a given view, in the form of
/// a 3 x N matrix, where N is the number of 3D points and each column
/// contains a pixel. Some columns may be zero, meaning there are no
/// measurements for this point from the given view.
///////////////////////////////////////////////////////////////////////////////
template<int M>
typename MeasurementMatrix<M>::ViewRef MeasurementMatrix<M>::view(int i) {
  return(middleRows(3*i,3));
}


///////////////////////////////////////////////////////////////////////////////
/// Synthesises a set of pixel correspondences based on a random set of
/// 'pts' 3D-points and a given motion matrix. A random noise of amplitude
/// given by 'noise' is added to the pixels.
///
/// \param P The motion matrix defining the position and orientation
/// of the cameras
///
/// \param pts The number of 3D points to synthesize data for
///
/// \param noise The amplitude of noise to add to the pixel measurements
///
///////////////////////////////////////////////////////////////////////////////
template<int M>
void MeasurementMatrix<M>::synthesise_measurements(const MotionMatrix<M> &P, 
						  int pts, double noise) {
  ShapeMatrix	 	X;
  Eigen::Matrix4d	shift;
  int 			i,j;

  // --- fill X with homogeneous 3D points centred around origin
  X.resize(4,pts);
  for(i = 0; i<3; ++i) {
    for(j = 0; j<X.cols(); ++j) {
      X(i,j) = rand()*2.0/RAND_MAX - 1.0;
    }
  }
  X.row(3).fill(1.0);
  
  // --- shift centre to (0,0,1.5)
  shift = Eigen::Matrix4d::Identity();
  shift(2,3) = 1.5;
  X = shift*X;

  // --- form matrix and add noise
  (*this) = P*X;
  unscale();
  for(i = 0; i<M; ++i) {
    for(j = 0; j<cols(); ++j) {
      (*this)(3*i,j) += 2.0*rand()*noise/RAND_MAX - noise;
      (*this)(3*i+1,j) += 2.0*rand()*noise/RAND_MAX - noise;
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
/// Synthesises occlusions in this measurement matrix by randomly
/// replacing some of the pixels with zero. Each pixel has a
/// 'proportion' probability of being occluded.
///////////////////////////////////////////////////////////////////////////////
template<int M>
void MeasurementMatrix<M>::synthesise_occlusions(double proportion) {
  int i,j;
  for(i = 0; i<M; ++i) {
    for(j=0; j<cols(); ++j) {
      if(rand()*1.0/RAND_MAX <= proportion) {
	pixel(i,j).fill(0.0);
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
/// Removes scaling from *this by multiplying each pixel, p, by 1/p(2)
/// i.e. by setting the homogeneos coord to z=1.0. Pixels with scale = 0.0
/// are assumed to be occluded and so are unchanged
///////////////////////////////////////////////////////////////////////////////
template<int M>
void MeasurementMatrix<M>::unscale() {
  int i,j;
  
  for(i = 0; i<M; ++i) {
    for(j = 0; j<cols(); ++j) {
      if(!pixel_is_occluded(i,j)) {
	pixel(i,j)(0) /= pixel(i,j)(2);
	pixel(i,j)(1) /= pixel(i,j)(2);
	pixel(i,j)(2) = 1.0;
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
/// Provides access to individual pixels within this matrix.
///
/// \param view The view of the pixel (i.e. row index/3)
/// \param point The 3D point of the pixel (i.e. column index)
///
/// \return reference to a 3x1 block containing the pixel
///////////////////////////////////////////////////////////////////////////////
template<int M>
typename MeasurementMatrix<M>::PixelRef
MeasurementMatrix<M>::pixel(int view, int point) {
  return(Base::template block<3,1>(view*3,point));
}
template<int M>
typename MeasurementMatrix<M>::constPixelRef
MeasurementMatrix<M>::pixel(int view, int point) const {
  return(Base::template block<3,1>(view*3,point));
}


///////////////////////////////////////////////////////////////////////////////
/// true if the given pixel of the given view is occluded (i.e. the
/// scale is zero).
///////////////////////////////////////////////////////////////////////////////
template<int M>
bool MeasurementMatrix<M>::pixel_is_occluded(int view, int point) const {
  double scale = pixel(view, point)(2);
  return(scale == 0.0);// && std::signbit(scale) == 0);
}


///////////////////////////////////////////////////////////////////////////////
/// Adds a new column to '*this', corresponding to measurements of a
/// new 3D point
///
/// \return the index of the new column
///////////////////////////////////////////////////////////////////////////////
template<int M>
int MeasurementMatrix<M>::new_column() {
  Base::conservativeResize(Eigen::NoChange, cols()+1);
  return(cols()-1);
}


///////////////////////////////////////////////////////////////////////////////
/// inserts 'n' new columns into *this.
///
/// \return index of the left-most column inserted
///////////////////////////////////////////////////////////////////////////////
template<int M>
int MeasurementMatrix<M>::new_columns(int n) {
  Base::conservativeResize(Eigen::NoChange, cols()+n);
  return(cols()-n);
}


///////////////////////////////////////////////////////////////////////////////
/// Inserts the column 'c' into *this and applies a number of
/// expansions to form a subspace matrix, as described in Martinec,
/// D. and Pajdla, T. 2002: Structire from many perspective images
/// with occlusions. ECCV 2002, LNCS 2351:355-369.
///
/// There are two types of expansion, each of which can be applied to
/// a view (i.e. 3 rows, corresponding to the 3 rows of a pixel)
///
/// - Pixel expansion: a pixel in 'c' is moved into a new column
/// that has all other pixels set to zero, the original column is set
/// to zero in that pixel position
///
/// - Full expansion: Three new columns are added that have zeroes in
/// all pixels apart from one, which have the values (1,0,0)^T,
/// (0,1,0)^T and (0,0,1)^T respectively.
///   
/// Full expansion is applied to each view in 'c' that has a missing
/// pixel (flagged by a 0.0 in the third row of the pixel), unless
/// that view has already been fully expanded. Pixel expansion is
/// applied to each view in 'c' that has an unweighted pixel
/// (flagged by a 1.0 in the third row of the pixel), unless that view
/// has already been fully expanded.
///
/// \return the index of the leftmost column that was inserted. The inserted
/// columns are always the rightmost columns.
///////////////////////////////////////////////////////////////////////////////
template<int M>
int MeasurementMatrix<M>::insert_column_and_expand(ColXpr c) {
  int 	i,j,k;
  int 	extracols;
  bool	fullyExpand[M];
  bool	pixelExpand[M];

  // --- First calculate how many cols to insert
  extracols = 1;
  for(i = 0; i<M; ++i) {
    if(c(3*i+2) == 0.0 && !is_fully_expanded(i)) {
      if(std::signbit(c(3*i+2)) == 0) {
	// --- full expansion
	fullyExpand[i] = true;
	pixelExpand[i] = false;
	extracols += 3;
      } else {
	// --- Pixel expansion
	fullyExpand[i] = false;
	pixelExpand[i] = true;
	++extracols;
      }
    } else {
      // --- no expansion
      fullyExpand[i] = false;
      pixelExpand[i] = false;
    }
  }

  // --- now insert cols and fill
  k = j = new_columns(extracols);
  col(j) = c;
  ++j;
  for(i = 0; i<M; ++i) {
    if(fullyExpand[i]) {
      middleCols(j,3).fill(0.0);
      (*this)(3*i   , j  ) = 1.0;
      (*this)(3*i+1 , j+1) = 1.0;
      (*this)(3*i+2 , j+2) = 1.0;
      j += 3;
    }
    else if(pixelExpand[i]) {
      col(j).fill(0.0);
      col(j)(3*i+2) = 1.0;
      col(j).middleRows(3*i,2) = c.middleRows(3*i,2);
      col(k).middleRows(3*i,3).fill(0.0);
      ++j;
    }
  }
  return(k);
}


///////////////////////////////////////////////////////////////////////////////
// reorder the columns of *this so that columns with the same occluded
// pixels are adjacent.
///////////////////////////////////////////////////////////////////////////////
template<int M>
void MeasurementMatrix<M>::order_cols_by_occlusions() {
}


///////////////////////////////////////////////////////////////////////////////
/// Inserts a number of columns into *this
///
/// \return index of the first column added
///////////////////////////////////////////////////////////////////////////////
template<int M>
int MeasurementMatrix<M>::insert_columns(ColsBlockXpr c) {
  int initialCols = cols();

  conservativeResize(Eigen::NoChange, initialCols + c.cols());  
  middleCols(initialCols, c.cols()) = c;
  return(initialCols);
}


///////////////////////////////////////////////////////////////////////////////
/// Removes column 'i' from this matrix
///////////////////////////////////////////////////////////////////////////////
template<int M>
void MeasurementMatrix<M>::delete_column(int i) {
  while(i < cols()-1) {
    col(i) = col(i+1);
  }
  conservativeResize(Eigen::NoChange, cols()-1);
}


///////////////////////////////////////////////////////////////////////////////
/// Removes columns 'i' to i+n from this matrix
///////////////////////////////////////////////////////////////////////////////
template<int M>
void MeasurementMatrix<M>::delete_columns(int i, int n) {
  while(i < cols()-n) {
    col(i) = col(i+n);
  }
  conservativeResize(Eigen::NoChange, cols()-n);
}


///////////////////////////////////////////////////////////////////////////////
/// Inserts a Gpixel in the correct row of the 'col'th column
/// \param col the column (3d point) to insert the pixel in
/// \param p the global pixel to insert
///////////////////////////////////////////////////////////////////////////////
template<int M>
void MeasurementMatrix<M>::insert_pixel(int col, const GPixel &p) {
  pixel(p.id, col)(0) = p.x;
  pixel(p.id, col)(1) = p.y;
  pixel(p.id, col)(2) = 1.0;
}

///////////////////////////////////////////////////////////////////////////////
/// Merges measurements from two different 3D points into one 3D point
/// by copying pixels in column j2 into the undefined portions of j1,
/// then deleting j2
///////////////////////////////////////////////////////////////////////////////
template<int M>
void MeasurementMatrix<M>::merge_columns(int j1, int j2) {
  // --- copy pixels from j2 to j1
  for(int i = 0; i<M; ++i) {
    if(pixel(i,j1)(2) == 0.0 && pixel(i,j2)(2) != 0.0) {
      pixel(i,j1) = pixel(i,j2);
    }
  }
  delete_column(j2);
}


///////////////////////////////////////////////////////////////////////////////
/// Returns the number of singular values of *this whose absolute
/// value is greater than s1*'threshold' where s1 is the absolute value
/// of the first singular value.
///////////////////////////////////////////////////////////////////////////////
template<int M>
int MeasurementMatrix<M>::rank(double threshold) {
  int i;
  Eigen::JacobiSVD<Base> svd(*this);

  if(svd.singularValues().size() == 0) return(0);

  threshold *= fabs(svd.singularValues()(0));
  i = 1;
  while(i < svd.singularValues().size() && 
	fabs(svd.singularValues()(i)) > threshold)
    ++i;
  return(i);
}


///////////////////////////////////////////////////////////////////////////////
/// True if the given view in *this has been fully expanded during
/// formation of a subspace matrix when accounting for
/// occlusions. Works by looking for a (1,x,0)^T pixel in 'view',
/// which shouldn't exist in a subspace matrix unless expanded.
///////////////////////////////////////////////////////////////////////////////
template<int M>
bool MeasurementMatrix<M>::is_fully_expanded(int view) {
  int j;

  for(j = 0; j<cols(); ++j) {
    if(pixel(view, j)(2) == 0.0 && pixel(view, j)(0) == 1.0) return(true);
  }
  return(false);
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
template<int M>
template<class D>
MeasurementMatrix<M> &MeasurementMatrix<M>::operator=(const Eigen::MatrixBase<D> &o) {
  Base::operator =(o);
  return(*this);
}
