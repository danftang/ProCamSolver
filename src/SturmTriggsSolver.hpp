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

#include "AdjacencyMatrix.h"
#include "FundamentalMatrix.h"
#include "MotionMatrix.h"
#include "MeasurementMatrix.h"
#include "ShapeMatrix.h"

///////////////////////////////////////////////////////////////////////////////
template<int M>
SturmTriggsSolver<M>::SturmTriggsSolver(MotionMatrix<M> &P,
					ShapeMatrix &S) :
  motion(P),
  shape(S)
{
}


///////////////////////////////////////////////////////////////////////////////
/// Approximately factorises the supplied measurement matrix into
/// shape and motion matrices. The results are put into the motion and
/// shape matrices supplied at construction.
///////////////////////////////////////////////////////////////////////////////
template<int M>
void SturmTriggsSolver<M>::solve(const MeasurementMatrix<M> &E) {
  Eigen::JacobiSVD<typename MeasurementMatrix<M>::Base> 	svd;
  ImageTransform<M>				denormalisation;

  set_measurement_matrix(E);
  measurement.normalise(denormalisation);

  std::cout << "Normalised measurement norm =" << std::endl;
  std::cout << measurement.rowwise().norm()/sqrt(measurement.cols()) << std::endl;

  approx_scale();
  factorise_with_occlusions(motion, shape);

  std::cout << "First guess at motion matrix = " << std::endl;
  std::cout << motion << std::endl;
  std::cout << "shapenorm = " << std::endl;
  std::cout << shape.rowwise().norm()/sqrt(shape.cols()) << std::endl;

  Eigen::VectorXd err;
  motion.reprojection_err(measurement, err);
  std::cout << "Normalised reprojection err = " << std::endl;
  std::cout << err.norm()/sqrt(err.rows())<<std::endl;

  measurement.scale_and_fill(motion, shape);
  svd.compute(measurement, Eigen::ComputeThinV | Eigen::ComputeThinU);
  motion = svd.matrixU().leftCols(4);
  // motion *= svd.singularValues().topRows(4).asDiagonal();

  std::cout << "Sturm-Triggs singular values are:" << std::endl;
  std::cout << svd.singularValues() << std::endl;

  denormalisation.apply_to(motion);
}


///////////////////////////////////////////////////////////////////////////////
/// Finds an approximate factorisation of the occluded measurement
/// matrix by forming a subspace of possible motion matrices, as
/// described in Martinec and Pajdla (2002). Once the motion matrix is
/// constrained, the shape that minimises reprojection error is calculated.
///////////////////////////////////////////////////////////////////////////////
template<int M>
void SturmTriggsSolver<M>::factorise_with_occlusions(MotionMatrix<M> &motion, 
						     ShapeMatrix &shape) {
  MeasurementMatrix<M>	notMotion; 	// nullspace of this is motion matrix
  MeasurementMatrix<M> 	subspace;	// nullspace of this is a constraint
                                        // on motionMatrix
  Eigen::JacobiSVD<typename MeasurementMatrix<M>::Base> subspaceSvd;
  std::vector<int>	permutation(measurement.cols());
  int 			i, j, j0;
  int 			c;
  int			srank;
  int			nmRank;
  const double	singular_cutoff=0.01;	// cutoff value for rank calculation

  // --- calculate the best permutation of
  // --- columns of the measurement matrix
  // -------------------------------------
  for(i=0; i<permutation.size(); ++i) {
    permutation[i] = i;
  }
  std::cout << "Permuting..." << std::endl;
  random_permute(permutation);
  std::cout << "Sorting..." << std::endl;
  std::sort(permutation.begin(), permutation.end(), colCompare(measurement));
  std::cout << "Forming subspaces..." << std::endl;

  // --- calculate motion matrix
  // --- from subspace constraints
  // -----------------------------
  c = 0;	// current position in permutation
  nmRank = 0;
  while(c < permutation.size()-4) {

    form_subspace(permutation.begin() + c, subspace);

    if(subspace.cols() < 3*M-1) {
      subspaceSvd.compute(subspace, Eigen::ComputeFullU);
      if(subspaceSvd.singularValues()(subspace.cols()-1)/subspaceSvd.singularValues()(0) > singular_cutoff) {
	// -- found full rank subspace
	notMotion.new_columns(3*M-subspace.cols());
	notMotion.rightCols(3*M-subspace.cols()) = 
	  subspaceSvd.matrixU().rightCols(3*M-subspace.cols());

	//	std::cout << "Subspace singular values = " << std::endl 
	//  << subspaceSvd.singularValues() << std::endl;
	//std::cout << "NotMotion singular values = " << std::endl 
	//<< notMotion.jacobiSvd().singularValues() << std::endl;

	//	c += 4;
	c += 1;
      } else {
	//	c += 2;
	c += 1;
      }
    } else {
      // c += 2;
      c += 1;
    }
  }

  if(notMotion.cols() < 3*M-4) {
    throw("Couldn't gather enough information to fill in the Measurement matrix.");
  }

  subspaceSvd.compute(notMotion, Eigen::ComputeFullU);
  if(subspaceSvd.singularValues()(3*M-5) < singular_cutoff) {
    throw("Too many occlusions to fill in the Measurement matrix.");
  }

  std::cout << "NotMotion singular values = " << std::endl 
	    << subspaceSvd.singularValues() << std::endl;


  // --- motion is approximate nullspace of the complement
  motion = subspaceSvd.matrixU().rightCols(4);

  shape.solve(measurement, motion);
}


///////////////////////////////////////////////////////////////////////////////
/// Sets G to the graph whose vertices represent views and whose edge
/// weights are a measure of the expected error in the fundamental
/// matrix between views. G is in the form of an adjacency matrix. 
///////////////////////////////////////////////////////////////////////////////
template<int M>
void SturmTriggsSolver<M>::
create_connectivity_graph(Eigen::Matrix<double,M,M> &G) {
  int 		v1,v2,j;
  int		connections;
  const double 	k = -0.2;

  G.fill(1e6);
  for(v1 = 0; v1<M; ++v1) {
    for(v2 = v1+1; v2<M; ++v2) {
      j = 0;
      connections = 0;
      while(j < measurement.cols() && connections < 50) {
	if(measurement.pixel(v1,j)(2) > 0 && measurement.pixel(v2,j)(2) > 0) {
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
void SturmTriggsSolver<M>::approx_scale() {
  Eigen::Vector3d		cross_prod, pix;    
  int				v,pt,parent,i;
  double			scale;
  FundamentalMatrix		F[M];

  // --- create fundamental matrices for
  // --- each edge in the spanning tree 
  // ----------------------------------
  for(v = 0; v<M; ++v) {
    if(v != tree_traversal(0)) {
      F[v].eight_point_algorithm(measurement.view(v),
				 measurement.view(spanning_tree(v)));
      F[v].remove_outliers(measurement.view(v),
			   measurement.view(spanning_tree(v)), 3.0);
    }
  }

  // --- Do scaling over tree
  // ------------------------
  for(v = 1; v<M; ++v) {
    i = tree_traversal(v);
    parent = spanning_tree(i);
    // --- scale between i and parent of i
    for(pt = 0; pt<measurement.cols(); ++pt) {
      if(!measurement.pixel_is_occluded(i,pt)) {
	if(!measurement.pixel_is_occluded(parent,pt)) {
	  pix        = measurement.pixel(i,pt);
	  cross_prod = F[i].e_ij().cross(pix);
	  scale      = 
	    cross_prod.dot(F[i]*measurement.pixel(parent,pt))/
	    cross_prod.dot(cross_prod);
	  measurement.pixel(i,pt) *= scale;
	} else {
	  // --- parent occluded so scale to 1.0
	  // --- for later separation
	  measurement.pixel(i,pt) /= measurement.pixel(i,pt)(2);
	}
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
/// Used by std::sort algorithm during occluded factorisation. Returns
/// true if column j1 of the measurement matrix is 'less than'
/// column j2 in a strict weak ordering.
///
/// The ordering is as follows:
///
/// 1) Order by positions of occluded pixels using Gray code ordering
///
/// Gray code ordering is used so that adjacent columns have the most
/// overlapping occlusions, and then the most overlapping unweighted pixels.
///
/// If (g_n...g_0) are the bits of a Gray code, then the binary
/// expansion of its position in order (b_n...b_0) can be calculated
/// using the recurrence relation:
///
/// b_{m-1} = b_m XOR g_m
///
/// with
///
/// b_n = g_n
///
///////////////////////////////////////////////////////////////////////////////
template<int M>
bool SturmTriggsSolver<M>::colCompare::operator ()(int j1, int j2) {
  int i;
  bool b1, b2;

  // --- row 0 is MSB
  b1 = false;
  b2 = false;
  for(i = 0; i<M; ++i) {
    // --- this is equivalent to b = b XOR g
    b1 = (b1 == measurement.pixel_is_occluded(i,j1));
    b2 = (b2 == measurement.pixel_is_occluded(i,j2));
    if(!b1 && b2) return(true);
    if(b1 && !b2) return(false);
  }


  // --- must be the equivalent
  return(false);
} 


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
template<int M>
void SturmTriggsSolver<M>::form_subspace(std::vector<int>::iterator cBegin,
					 MeasurementMatrix<M> &subspace) {
  int 				v,i,j,k;
  std::vector<int>::iterator 	cEnd = cBegin+4;
  std::vector<int>::iterator 	colIt;
  Eigen::Matrix<bool,M,4>	soloPixel; // single pixel on its own?
  Eigen::Matrix<bool,M,1>	occluded;
  Eigen::Matrix<int,M,1>	colPosition;	
  bool				updated;
  int				colCount;
  int				parent;

  colCount = 0;

  // --- First calculate occlusions
  // --- occluded     - true if row is occluded
  for(i = 0; i<M; ++i) {
    occluded(i) = false;
    for(colIt = cBegin; colIt != cEnd; ++colIt) {
      if(measurement.pixel_is_occluded(i,*colIt) && !occluded(i)) {
	occluded(i) = true;
	colCount += 3;
      }
    }    
  }

  // --- Now find positions of columns
  // --- with one pixel on its own.
  soloPixel.fill(false);
  for(v = 0; v<M; ++v) {
    i = tree_traversal(v);
    if(!occluded(i)) {
      j = 0;
      for(colIt = cBegin; colIt != cEnd; ++colIt) {
	if(spanning_tree(i) == -1 || occluded(spanning_tree(i))) {
	  soloPixel(i,j) = true;
	  colCount += 1;
	} else {
	  soloPixel(spanning_tree(i),j) = false;
	}
	++j;
      }
    }
  }

  // --- remove redundant solo-pixel expansions
  // --- (4 x solo pixels in one row)
  for(i = 0; i<M; ++i) {
    if(!occluded(i) &&
       soloPixel.row(i) == Eigen::Matrix<bool,1,4>(true,true,true,true)) {
      occluded(i) = true;
      colCount -= 1;
    }
  }
  
  // std::cout << "Occluded = " << occluded.transpose() << std::endl;
  // std::cout << "Solo pixels = " << std::endl << soloPixel << std::endl;
  // std::cout << "ColCount = " << colCount << std::endl;

  // --- now create subspace
  // --- column by column
  subspace.resize(3*M, colCount);
  subspace.setZero();
  j = 0;
  for(colIt = cBegin; colIt != cEnd; ++colIt) {
    for(v = 0; v<M; ++v) {
      i = tree_traversal(v);
      if(!occluded(i)) {
	parent = spanning_tree(i);
	if(parent == -1 || occluded(parent)) {
	  // --- start a new column
	  subspace.pixel(i,j) = measurement.pixel(i,*colIt);
	  colPosition(i) = j;
	  updated = true;
	  ++j;
	} else {
	  // --- add to an existing column
	  colPosition(i) = colPosition(parent);
	  subspace.pixel(i,colPosition(i)) = measurement.pixel(i,*colIt);
	  updated = true;
	}
      }
    }
  }

  //--- now add occlusion expansions
  for(i = 0; i<M; ++i) {
    if(occluded(i)) {
      subspace.template block<3,3>(3*i,j).setIdentity();
      j += 3;
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
/// Randomly permutes the elements of p
///////////////////////////////////////////////////////////////////////////////
template<int M>
void SturmTriggsSolver<M>::random_permute(std::vector<int> &p) {
  int i,j,d;
  

  for(i=0; i<p.size()-1; ++i) {
    j = i + 1 + (int)((p.size()-i-1)*(rand()*(1.0 - 1e-8)/RAND_MAX));
    d = p[j];
    p[j] = p[i];
    p[i] = d;
  }
}


///////////////////////////////////////////////////////////////////////////////
/// Takes a copy of the measurement matrix and builds a spanning tree of views
///////////////////////////////////////////////////////////////////////////////
template<int M>
void SturmTriggsSolver<M>::set_measurement_matrix(const MeasurementMatrix<M> &E) {
  AdjacencyMatrix<M> 		view_connections;   // graph of views
  Eigen::Matrix<bool,M,1>	visited;	    // has this view been visited
  int				v,i;
  int				treeRoot;
  bool				updated;

  measurement = E;

  create_connectivity_graph(view_connections);
  view_connections.shortest_spanning_tree(treeRoot, spanning_tree);
  std::cout << "Spanning tree = " << std::endl << spanning_tree << std::endl; 
  visited.fill(false);
  visited(treeRoot) = true;
  tree_traversal(0) = treeRoot;
  v = 1;
  do {
    updated = false;	
    for(i = 0; i<M; ++i) {
      if(!visited(i) && visited(spanning_tree(i))) {
	tree_traversal(v++) = i;
	visited(i) = true;
	updated = true;
      }
    }
  } while(updated);
}
