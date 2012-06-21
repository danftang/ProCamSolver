
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
/// Normalises *this by setting the scales of all un-occluded pixels
/// to 1.0, the mean of all pixels in each view to zero and the
/// standard deviation in each view to 1 in both x and y directions.
///
/// \param denormalisation This is set to the inverse transform 
/// for subsequent de-normalisation.
///////////////////////////////////////////////////////////////////////////////
template<int M>
void MeasurementMatrix<M>::normalise(ImageTransform<M> &denormalisation) {
  Eigen::Vector3d mean;
  Eigen::Vector2d sd;
  Eigen::Matrix3d K;
  int v;
  int point;

  unscale();
  K.fill(0.0);
  for(v = 0; v < M; ++v) {
    mean = view(v).rowwise().sum();
    sd   = view(v).topRows(2).rowwise().squaredNorm();
    if(mean(2) == 0.0) mean(2) = 1.0; // no data for view
    sd   /= mean(2);
    mean /= mean(2);
    sd = (sd - mean.cwiseProduct(mean).topRows(2)).cwiseSqrt();
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
/// Implementation of the eight-point algorithm to calculate the
/// fundamental matrix from a set of correspondences. This matrix
/// should be normalised before calling this method to ensure
/// stability.
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
  int validCorrespondences;
  int n;

  // --- form normalised correspondence matrix
  // -----------------------------------------
  validCorrespondences = 0;
  for(point = 0; point<cols(); ++point) {
    if(pixel(i,point)(2) != 0.0 && pixel(j,point)(2) != 0.0) {
      ++validCorrespondences;
    }
  }
  A.resize(validCorrespondences,9);
  n = 0;
  for(point = 0; point<cols(); ++point) {
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
/// Fills in missing pixels and weights using the method described in Martinec
/// and Pajdla (2002)
///////////////////////////////////////////////////////////////////////////////
template<int M>
void MeasurementMatrix<M>::factorise_with_occlusions(MotionMatrix<M> &motion, 
						     ShapeMatrix &shape) {
  MeasurementMatrix<M>	notMotion; 	// nullspace of this is motion matrix
  MeasurementMatrix<M> 	subspace;	// nullspace of this is a constraint
                                        // on motionMatrix
  int 		i, j, j0;
  int 		c;
  int		srank;
  int		nmRank;
  const double	singular_cutoff=0.001;	// cutoff value for rank calculation

  // --- calculate motion matrix
  // --- from subspace constraints
  // -----------------------------
  c = 0;
  nmRank = 0;
  while(c<cols()) {
    // --- create constraint
    subspace.delete_columns(0,subspace.cols());
    for(j = 0; c<cols() && j<4; ++j) {
      // --- add a column to subspace
      j0 = subspace.insert_column_and_expand(col(c));
      ++c;
      srank = subspace.rank(singular_cutoff);
      while(c<cols() && srank < subspace.cols() && srank < 3*M) {
	subspace.delete_columns(j0, subspace.cols()-j0);
	j0 = subspace.insert_column_and_expand(col(c));
	++c;
      }
    }

    std::cout << "Adding constraint:" << std::endl << subspace << std::endl;

    // --- add nullspace of constraint to notMotion
    if(j == 4 && subspace.cols() < 3*M) {
      notMotion.new_columns(3*M-subspace.cols());
      notMotion.rightCols(3*M-subspace.cols()) = subspace.jacobiSvd(Eigen::ComputeFullU).matrixU().rightCols(3*M-subspace.cols());
      nmRank = notMotion.rank(singular_cutoff);
    }
    std::cout << "notmotion =" << std::endl << notMotion << std::endl;
    std::cout << "notMotion singular values = " << std::endl;
    std::cout << notMotion.jacobiSvd().singularValues() << std::endl;
  }
  
  if(nmRank < 3*M-4) {
    throw("Too many occlusions to fill in the Measurement matrix.");
  }


  motion = notMotion.jacobiSvd(Eigen::ComputeFullU).matrixU().rightCols(4);
  std::cout << "motion =" << std::endl << motion << std::endl;

  shape.solve(*this, motion);

  /************
  Eigen::Vector4d	point;
  MotionMatrix<M>	maskedMotion; // zeroes inserted in unknown pixels
  Eigen::Matrix<double,3*M,1>	maskedPixels;

  for(j = 0; j<cols(); ++j) {
    // --- calculate 3D point by least squares
    maskedMotion = motion;
    maskedPixels = col(j);
    for(i = 0; i<M; ++i) {
      if(pixel(i,j)(2) == 0.0) {
	if(std::signbit(pixel(i,j)(2)) == 0) {
	  // --- missing point
	  maskedMotion.view(i).fill(0.0);
	} else {
	  // --- missing scale, use constraints:
	  // --- (P0/x - P2)X = 0
	  // --- (P1/y - P2)X = 0
	  maskedMotion.view(i).row(0) = 
	    maskedMotion.view(i).row(0) / maskedPixels(3*i) - 
	    maskedMotion.view(i).row(2);
	  maskedMotion.view(i).row(1) = 
	    maskedMotion.view(i).row(1) / maskedPixels(3*i+1) - 
	    maskedMotion.view(i).row(2);
	  maskedMotion.view(i).row(2).fill(0.0); 
	}
	maskedPixels.middleRows(3*i,3).fill(0.0);
      }
    }

    if(maskedPixels.squaredNorm() == 0.0) {
      // all pixels are unscaled, so set scale
      maskedPixels(3*M-1) = 1.0;
      maskedMotion.row(3*M-1) = motion.row(3*M-1);
    }

    std::cout << "MaskedPixels = " << std::endl << maskedPixels << std::endl;
    std::cout << "MaskedMotion = " << std::endl << maskedMotion << std::endl;

    point = maskedMotion.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(maskedPixels);

    std::cout << "Least squares error = " << std::endl << maskedMotion*point - maskedPixels << std::endl;
    
    // --- fill in pixels from 3D-point
    double scale;
    for(i = 0; i<M; ++i) {
      if(pixel(i,j)(2) == 0.0) {
	if(std::signbit(pixel(i,j)(2)) == 0) {
	  // --- missing pixel
	  pixel(i,j) = motion.view(i)*point;
	} else {
	  // --- unweighted pixel
	  // --- find least squares scale
	  pixel(i,j)(2) = 1.0;
	  scale = 
	    pixel(i,j).dot(motion.view(i)*point) / 
	    pixel(i,j).squaredNorm();
	  pixel(i,j) *= scale;
	}
      }
    }
  }
  *************/
}


///////////////////////////////////////////////////////////////////////////////
/// scales the un-scaled pixels of *this and fills in occluded pixels
/// based on the given motion matrix and shape matrix. Scale is calculated
/// as lambda = P_2X.
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
  noiseVec.resize(cols());
  for(i = 0; i<M; ++i) {
    noiseVec = Eigen::RowVectorXd::Random(cols())*noise;
    row(3*i) += noiseVec;
    noiseVec = Eigen::RowVectorXd::Random(cols())*noise;
    row(3*i+1) += noiseVec;
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
/// scale is positive-zero).
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
  conservativeResize(Eigen::NoChange, cols()+1);
  return(cols()-1);
}


///////////////////////////////////////////////////////////////////////////////
/// inserts 'n' new columns into *this.
///
/// \return index of the left-most column inserted
///////////////////////////////////////////////////////////////////////////////
template<int M>
int MeasurementMatrix<M>::new_columns(int n) {
  conservativeResize(Eigen::NoChange, cols()+n);
  return(cols()-n);
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
      while(j < cols() && connections < 50) {
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
/// Approximates the scaling of '*this' by forming the fundamental
/// matrices, as described in Sturm and Triggs (1996). The views are
/// ordered by the minimum depth spanning tree of the connectivity
/// graph. Unscaled pixels have their scale set to -0.0
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
  double			parent_scale, child_scale;
  FundamentalMatrix		F[M];
  ImageTransform<M>		denormalisation;

  // --- create spanning tree of views
  // ---------------------------------
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

  // --- Do scaling over tree
  // ------------------------
  updated = true;
  while(updated) {
    updated = false;
    for(v = 0; v<M; ++v) {
      parent = spanning_tree(v);
      if(!scaled(v) && scaled(parent)) {
	// --- scale between v and parent of v
	for(pt = 0; pt<cols(); ++pt) {
	  if(pixel(v,pt)(2) != 0.0) {
	    if(pixel(parent,pt)(2) != 0.0) {
	      pix = pixel(v,pt);
	      cross_prod = F[v].e_ij().cross(pix);
	      scale = cross_prod.dot(F[v]*pixel(parent,pt))/
		cross_prod.dot(cross_prod);
	      pixel(v,pt) *= scale;
	    } else {
	      // --- parent missing so all descendants are unscaled
	      pixel(v,pt)(2) = -0.0;
	    }
	  }
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
void MeasurementMatrix<M>::svd_factorise(Eigen::Matrix<double,3*M,4> &P) {
  Eigen::JacobiSVD<Base> svd;

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
void MeasurementMatrix<M>::svd_factorise(Eigen::Matrix<double,3*M,4> &P, Eigen::Matrix4Xd &X) {
  Eigen::JacobiSVD<Base> svd;

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
  if(cols() == 0) return(0);
  Eigen::JacobiSVD<Base> svd(*this);

  int i = 0;
  while(i < svd.singularValues().size() && 
	fabs(svd.singularValues()(i)) > threshold)
    ++i;
  return(i);
}


///////////////////////////////////////////////////////////////////////////////
/// True if this view has been fully expanded. Works by looking for a
/// (1,x,0)^T pixel in 'view', which shouldn't exist IN A SUBSPACE
/// MATRIX if not expanded.
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


///////////////////////////////////////////////////////////////////////////////
// precompiled
///////////////////////////////////////////////////////////////////////////////
//template class MeasurementMatrix<4>;
