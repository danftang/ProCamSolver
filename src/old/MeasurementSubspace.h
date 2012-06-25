#ifndef MEASUREMENTSUBSPACE_H
#define MEASUREMENTSUBSPACE_H

#include "MeasurementMatrix.h"

//////////////////////////////////////////////////////////////////////////////
/// Class to represent a subspace of the space of measurements of N 3D
/// points from M views. As defined in Martinec and Pajdla (2002)
//////////////////////////////////////////////////////////////////////////////
template<int M>
class MeasurementSubspace : public MeasurementMatrix<M> {
public:
  int  	insert_column_and_expand(ColXpr);

  
protected:
  int	rank(double = 0.0);
  bool	is_fully_expanded(int);
};


///////////////////////////////////////////////////////////////////////////////
/// Loads a set of correspondences from a file in the format
///  <id1> <id2> <x1> <y1> <x2> <y2>
///
///////////////////////////////////////////////////////////////////////////////
template<int M>
bool MeasurementMatrix<M>::load(const char *filename) {
  std::ifstream myFile(filename);
  int i;
  Correspondence c;

  if(!myFile) return(false);
  i = 10000;
  while(i >0 && (myFile >> c)) {
    add_correspondence(c);
    //    std::cout << c << std::endl;
    --i;
  }
  myFile.close();
  return(true);
}


///////////////////////////////////////////////////////////////////////////////
/// Sets the mean of all pixels in each view to zero and the standard deviation
/// to 1 in both x and y directions.
///
/// Saves the offset and re-scaling to the matrix 'denormalisation'
/// for subsequent de-normalisation using denormalise().
///////////////////////////////////////////////////////////////////////////////
template<int M>
void MeasurementMatrix<M>::normalise() {
  Eigen::Vector3d mean;
  Eigen::Vector2d sd;
  Eigen::Matrix3d K;
  int v;
  int point;

  if(normalised) return;
  K.fill(0.0);
  for(v = 0; v < M; ++v) {
    mean = view(v).rowwise().sum();
    sd   = view(v).topRows(2).rowwise().squaredNorm();
    sd   /= mean(2);
    mean /= mean(2);
    sd = (sd - mean.cwiseProduct(mean).topRows(2)).cwiseSqrt();

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
  normalised = true;
}

template<int M>
void MeasurementMatrix<M>::denormalise() {
  if(!normalised) return;
  denormalisation.transform(*this);
  normalised = false;
}


///////////////////////////////////////////////////////////////////////////////
/// Implementation of the normalised eight-point algorithm to
/// calculate the fundamental matrix from a set of correspondences.
///
/// See Hartley (1997) for details of the 8-point algorithm.
///
/// \param F variable to hold the resulting fundamental matrix
/// \param i id of the first camera
/// \param j id of the second camera
///////////////////////////////////////////////////////////////////////////////
template<int M>
void MeasurementMatrix<M>::eight_point_algorithm(FundamentalMatrix &F, 
						 int i, int j) {
  // A is the lifted matrix of correspondences
  Eigen::Matrix<double,Eigen::Dynamic,9>			A;
  Eigen::Matrix<double,9,9>					ATA;
  Eigen::JacobiSVD<Eigen::Matrix<double,Eigen::Dynamic,9> >  	Adecomposition;
  Eigen::JacobiSVD<Eigen::Matrix3d>                   		Fdecomposition;
  int point;
  int n;

  normalise();

  // --- form normalised correspondence matrix
  // -----------------------------------------
  A.resize(points(),9);
  n = 0;
  for(point = 0; point<points(); ++point) {
    if(pixel(i,point)(2) != 0.0 && pixel(j,point)(2) != 0.0) {
      A(n,0) = pixel(j,point)(0)*pixel(i,point)(0);
      A(n,1) = pixel(j,point)(0)*pixel(i,point)(1);
      A(n,2) = pixel(j,point)(0);
      A(n,3) = pixel(j,point)(1)*pixel(i,point)(0);
      A(n,4) = pixel(j,point)(1)*pixel(i,point)(1);
      A(n,5) = pixel(j,point)(1);
      A(n,6) = pixel(i,point)(0);
      A(n,7) = pixel(i,point)(1);
      A(n,8) = 1.0;
      ++n;
    }
  }
  ATA = A.transpose()*A;

  // --- Approximately solve AF = 0 by SVD
  // -------------------------------------

  Adecomposition.compute(ATA, Eigen::ComputeFullV);
  F.col(0) = Adecomposition.matrixV().block(0,8,3,1);
  F.col(1) = Adecomposition.matrixV().block(3,8,3,1);
  F.col(2) = Adecomposition.matrixV().block(6,8,3,1);

  // --- Enforce zero determinant of F
  // ---------------------------------
  F.calc_epipoles_and_constrain();

  /****
  Fdecomposition.compute(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::Vector3d D;
  D = Fdecomposition.singularValues();
  D(2) = 0.0;
  F = Fdecomposition.matrixU()* D.asDiagonal() *Fdecomposition.matrixV().transpose();
  *****/
}


///////////////////////////////////////////////////////////////////////////////
/// Adds a pair of corresponding pixels to '*this'
///
/// \return the column that the correspondence was added to 
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
      icolumn = merge_columns(icolumn, jcolumn);
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
  for(j = 0; j<points(); ++j) {
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
  return(Base::middleRows(3*i,3));
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
  Eigen::RowVectorXd	noiseVec;
  int 			i;

  // --- fill X with homogeneous 3D points centred around origin
  X.resize(4,pts);
  X = ShapeMatrix::Random(4,pts);
  X.row(3).fill(1.0);
  
  // --- shift centre to (0,0,1.5)
  shift = Eigen::Matrix4d::Identity();
  shift(2,3) = 1.5;
  X = shift*X;

  // --- form matrix and add noise
  (*this) = P*X;
  unscale();
  noiseVec.resize(points());
  for(i = 0; i<M; ++i) {
    noiseVec = Eigen::RowVectorXd::Random(points())*noise;
    Base::row(3*i) += noiseVec;
    noiseVec = Eigen::RowVectorXd::Random(points())*noise;
    Base::row(3*i+1) += noiseVec;
  }
}


///////////////////////////////////////////////////////////////////////////////
/// Removes scaling from *this by multiplying each pixel, p, by 1/p(2)
/// i.e. by setting the homogeneos coord to z=1.0
///////////////////////////////////////////////////////////////////////////////
template<int M>
void MeasurementMatrix<M>::unscale() {
  int i;
  
  for(i = 0; i<3*M; i += 3) {
    Base::row(i) 	= Base::row(i).cwiseQuotient(Base::row(i+2));
    Base::row(i+1) 	= Base::row(i+1).cwiseQuotient(Base::row(i+2));
    Base::row(i+2).fill(1.0);
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
  return(Base::middleRows(view*3,3).col(point));
}


///////////////////////////////////////////////////////////////////////////////
/// Adds a new column to '*this', corresponding to measurements of a
/// new 3D point
///
/// \return the index of the new column
///////////////////////////////////////////////////////////////////////////////
template<int M>
int MeasurementMatrix<M>::new_column() {
  conservativeResize(Eigen::NoChange, points()+1);
  return(points()-1);
}


///////////////////////////////////////////////////////////////////////////////
/// inserts 'n' new columns into *this.
///
/// \return index of the left-most column inserted
///////////////////////////////////////////////////////////////////////////////
template<int M>
int MeasurementMatrix<M>::new_columns(int n) {
  conservativeResize(Eigen::NoChange, points()+n);
  return(points()-n);
}


///////////////////////////////////////////////////////////////////////////////
/// Inserts 'c' into *this and applies a number of expansions.
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
    if(c(2+3*i) == 0.0 && !is_fully_expanded(i)) {
      // --- full expansion
      fullyExpand[i] = true;
      pixelExpand[i] = false;
      extracols += 3;
    }
    else if(c(2+3*i) == 1.0 && !is_fully_expanded(i)) {
      // --- Pixel expansion
      fullyExpand[i] = false;
      pixelExpand[i] = true;
      ++extracols;
    } else {
      fullyExpand[i] = false;
      pixelExpand[i] = false;
    }
  }

  // --- now insert cols and fill
  k = j = new_columns(extracols);
  Base::col(j) = c;
  ++j;
  for(i = 0; i<M; ++i) {
    if(fullyExpand[i]) {
      Base::middleCols(j,3).fill(0.0);
      (*this)(j   ,3*i) = 1.0;
      (*this)(j+1 ,3*i+1) = 1.0;
      (*this)(j+2 ,3*i+2) = 1.0;
      j += 3;
    }
    else if(pixelExpand[i]) {
      Base::col(j).fill(0.0);
      Base::col(j).middleRows(3*i,3) = c.middleRows(3*i,3);
      Base::col(k).middleRows(3*i,3).fill(0.0);
      ++j;
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
/// Inserts a number of columns into *this
///
/// \return index of the first column added
///////////////////////////////////////////////////////////////////////////////
template<int M>
int MeasurementMatrix<M>::insert_columns(ColsBlockXpr cols) {
  int initialCols = cols();

  conservativeResize(Eigen::NoChange, initialCols + cols.cols());  
  middleCols(initialCols, cols.cols()) = cols;
  return(initialCols);
}


///////////////////////////////////////////////////////////////////////////////
/// Removes column 'i' from this matrix
///////////////////////////////////////////////////////////////////////////////
template<int M>
void MeasurementMatrix<M>::delete_column(int i) {
  while(i < points()-1) {
    Base::col(i) = Base::col(i+1);
  }
  conservativeResize(Eigen::NoChange, points()-1);
}


///////////////////////////////////////////////////////////////////////////////
/// Removes columns 'i' to i+n from this matrix
///////////////////////////////////////////////////////////////////////////////
template<int M>
void MeasurementMatrix<M>::delete_columns(int i, int n) {
  while(i < points()-n) {
    Base::col(i) = Base::col(i+n);
  }
  conservativeResize(Eigen::NoChange, points()-n);
}


///////////////////////////////////////////////////////////////////////////////
/// Inserts a Gpixel in the correct row of the 'col'th column
/// \param col the column (3d point) to insert the pixel in
/// \param p the global pixel to insert
///////////////////////////////////////////////////////////////////////////////
template<int M>
void MeasurementMatrix<M>::insert_pixel(int col, GPixel &p) {
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
  for(int i = 0; i<views(); ++i) {
    if(pixel(i,j1)(2) == 0.0 && pixel(i,j2)(2) != 0.0) {
      pixel(i,j1) = pixel(i,j2);
    }
  }
  delete_column(j2);
}


///////////////////////////////////////////////////////////////////////////////
/// Sets G to the graph whose vertices represent views and whose edge
/// weights are a measure of the expected error in the fundamental
/// matrix between views. G is in the form of an adjacency matrix.
///
/// 
///////////////////////////////////////////////////////////////////////////////
template<int M>
void
MeasurementMatrix<M>::create_connectivity_graph(Eigen::Matrix<double,M,M> &G) {
  int 		v1,v2,j;
  int		connections;
  const double 	k = -0.2;

  G.fill(1e6);
  for(v1 = 0; v1<M; ++v1) {
    for(v2 = v1+1; v2<M; ++v2) {
      j = 0;
      connections = 0;
      while(j < points() && connections < 50) {
	if(pixel(v1,j)(2) > 0 && pixel(v2,j)(2) > 0) {
	  ++connections;
	}
	++j;
      }
      if(connections > 7) {
	G(v1,v2) = G(v2,v1) = exp(k*connections);
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
template<int M>
template<class D>
MeasurementMatrix<M> &MeasurementMatrix<M>::operator=(const Eigen::MatrixBase<D> &o) {
  Base::operator =(o);
  return(*this);
}


///////////////////////////////////////////////////////////////////////////////
/// Approximates the scaling of '*this' using the method described in Sturm
/// and Triggs (1996). The views are ordered by the minimum depth spanning
/// tree of the connectivity graph
///////////////////////////////////////////////////////////////////////////////
template<int M>
void MeasurementMatrix<M>::approx_scale() {
  AdjacencyMatrix<M> 		view_connections;
  int				tree_root;
  Eigen::Matrix<int,M,1>	spanning_tree;
  Eigen::Matrix<bool,M,1>	scaled;
  Eigen::Vector3d		cross_prod, pix;
  int				v,pt,parent;
  bool				updated;
  double			scale;
  FundamentalMatrix		F[M];

  create_connectivity_graph(view_connections);
  view_connections.shortest_spanning_tree(tree_root, spanning_tree);
  scaled.fill(false);
  scaled(tree_root) = true;


  std::cout << "Spanning tree = " << std::endl << spanning_tree << std::endl; 

  // --- create fundamental matrices for
  // --- each edge in the spanning tree 
  // ----------------------------------
  for(v = 0; v<M; ++v) {
    if(v != tree_root) {
      eight_point_algorithm(F[v],v,spanning_tree(v));      
    }
  }

  updated = true;
  while(updated) {
    updated = false;
    for(v = 0; v<M; ++v) {
      parent = spanning_tree(v);
      if(!scaled(v) && scaled(parent)) {
	// --- scale between v and parent of v
	for(pt = 0; pt<points(); ++pt) {
	  pix = pixel(v,pt);
	  cross_prod = F[v].e_ij().cross(pix);
	  scale = cross_prod.dot(F[v]*pixel(parent,pt))/
	          cross_prod.dot(cross_prod);
	  //	  std::cout << "Scaling pixel " << pixel(v,pt).transpose() << std::endl;
	  pixel(v,pt) *= scale;
	  //std::cout << "Scaled pixel " << pixel(v,pt).transpose() << std::endl;
	}
	scaled(v) = true;
	updated = true;
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
/// Decompose *this into MotionMatrix and ShapeMatrix
///////////////////////////////////////////////////////////////////////////////
template<int M>
void MeasurementMatrix<M>::factorise(Eigen::Matrix<double,3*M,4> &P) {
  Eigen::JacobiSVD<Base> svd;

  approx_scale();
  svd.compute(*this, Eigen::ComputeThinV | Eigen::ComputeThinU);
  P = svd.matrixU().leftCols(4);
  P *= svd.singularValues().topRows(4).asDiagonal();

  std::cout << "Singular values are" << std::endl;
  std::cout << svd.singularValues() << std::endl;
  std::cout << "U = " << std::endl;
  std::cout << svd.matrixU() << std::endl;
  std::cout << "V = " << std::endl;
  std::cout << svd.matrixV() << std::endl;

}


template<int M>
void MeasurementMatrix<M>::factorise(Eigen::Matrix<double,3*M,4> &P, Eigen::Matrix4Xd &X) {
  Eigen::JacobiSVD<Base> svd;

  approx_scale();
  svd.compute(*this, Eigen::ComputeThinV | Eigen::ComputeThinU);
  P = svd.matrixU().leftCols(4);
  P *= svd.singularValues().topRows(4).asDiagonal();
  X = svd.matrixV().leftCols(4).transpose();

  std::cout << "Singular values are" << std::endl;
  std::cout << svd.singularValues() << std::endl;
  std::cout << "U = " << std::endl;
  std::cout << svd.matrixU() << std::endl;
  std::cout << "V = " << std::endl;
  std::cout << svd.matrixV() << std::endl;

}


///////////////////////////////////////////////////////////////////////////////
/// Returns the number of singular values of *this whose absolute
/// value is greater than 'threshold'
///////////////////////////////////////////////////////////////////////////////
template<int M>
int MeasurementMatrix<M>::rank(double threshold) {
  Eigen::JacobiSVD<Base> svd(*this);

  int i = 0;
  while(i < svd.singularValues().size() && 
	fabs(svd.singularValues()(i)) > threshold)
    ++i;
  return(i);
}


///////////////////////////////////////////////////////////////////////////////
/// True if this view has been fully expanded. Works by looking for a pixel
/// (1,x,0)^T, which should never occur other than in a full expansion
///////////////////////////////////////////////////////////////////////////////
template<int M>
bool MeasurementMatrix<M>::is_fully_expanded(int view) {
  int j;

  for(j = 0; j<points(); ++j) {
    if(pixel(view, j)(2) == 0.0 && pixel(view, j)(0) == 1.0) return(true);
  }
  return(false);
}


#endif
