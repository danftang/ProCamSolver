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

  measurement = E;
  measurement.normalise(denormalisation);

  std::cout << "Normalised measurement norm =" << std::endl;
  std::cout << measurement.rowwise().norm() << std::endl;

  approx_scale();
  factorise_with_occlusions(motion, shape);

  std::cout << "First guess at motion matrix = " << std::endl;
  std::cout << motion << std::endl;
  std::cout << "shapenorm = " << std::endl;
  std::cout << shape.rowwise().norm() << std::endl;

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
  const double	singular_cutoff=0.001;	// cutoff value for rank calculation

  // --- calculate the best permutation of
  // --- columns of the measurement matrix
  // -------------------------------------
  for(i=0; i<permutation.size(); ++i) {
    permutation[i] = i;
  }
  random_permute(permutation);
  std::sort(permutation.begin(), permutation.end(), *this);

  std::cout << "Permuted measurement =" << std::endl;
  for(j = 0; j<permutation.size(); ++j) {
    for(i = 0; i<M; ++i) {
      std::cout << (measurement(3*i+2,permutation[j]) != 0) << " ";
    }
    std::cout << std::endl;
  }

  // --- calculate motion matrix
  // --- from subspace constraints
  // -----------------------------
  c = 0;	// current position in permutation
  nmRank = 0;
  while(c < permutation.size()-4) {

    std::cout << "Forming subspace from =" << std::endl;
    for(j = 0; j<4; ++j) {
      for(i = 0; i<M; ++i) {
	std::cout << (measurement(3*i+2,permutation[c+j]) != 0) << " ";
      }
      std::cout << std::endl;
    }    

    form_subspace(permutation.begin() + c, subspace);

    std::cout << "Subspace is =" << std::endl;
    for(j = 0; j<subspace.cols(); ++j) {
      for(i = 0; i<M; ++i) {
	std::cout << (subspace(3*i+2,j) != 0) << " ";
      }
      std::cout << std::endl;
    }    

    if(subspace.cols() < 3*M-1) {
      subspaceSvd.compute(subspace, Eigen::ComputeFullU);
      if(subspaceSvd.singularValues()(subspace.cols()-1) > singular_cutoff) {
	// -- found full rank subspace
	notMotion.new_columns(3*M-subspace.cols());
	notMotion.rightCols(3*M-subspace.cols()) = 
	  subspaceSvd.matrixU().rightCols(3*M-subspace.cols());
	c += 4;
	std::cout << "Added contraint to notmotion" << std::endl;
	std::cout << "Singular vals = " << notMotion.jacobiSvd().singularValues().transpose() << std::endl;
      } else {
	c += 2;
      }
    } else {
      c += 2;
    }
  }

  /******
  while(c < permutation.size()) {
    // --- create constraint
    for(j = 0; c<measurement.cols() && j<4; ++j) {
      subspace.delete_columns(0,subspace.cols());
      // --- add a column to subspace
      j0 = subspace.insert_column_and_expand(measurement.col(c));
      ++c;
      srank = subspace.rank(singular_cutoff);
      // --- subspace must be of full rank
      while(c<measurement.cols() && srank < subspace.cols() && srank < 3*M) {
	subspace.delete_columns(j0, subspace.cols()-j0);
	j0 = subspace.insert_column_and_expand(measurement.col(c));
	srank = subspace.rank(singular_cutoff);
	++c;
      }
    }

    // --- add nullspace of subspace to notMotion
    // --- if it is a meaningful constraint
    if(j == 4 && subspace.cols() < 3*M) {
      notMotion.new_columns(3*M-subspace.cols());
      notMotion.rightCols(3*M-subspace.cols()) = subspace.jacobiSvd(Eigen::ComputeFullU).matrixU().rightCols(3*M-subspace.cols());

      nmRank = notMotion.rank(singular_cutoff); // test
      //      std::cout << "nmRank = " << nmRank << std::endl;

    }
  }

  if(nmRank < 3*M-4) {
    std::cout << "Not motion singular values " << std::endl 
	      << notMotion.jacobiSvd().singularValues() << std::endl;
    throw("Too many occlusions to fill in the Measurement matrix.");
  }

  ****/  

  // --- motion is approximate nullspace of the complement
  motion = notMotion.jacobiSvd(Eigen::ComputeFullU).matrixU().rightCols(4);

  // --- remove 'unscalable' flags in measurement matrix
  for(i = 0; i<M; ++i) {
    for(j = 0; j<measurement.cols(); ++j) {
      if(measurement.pixel(i,j)(2) == 0.0 && std::signbit(measurement.pixel(i,j)(2)) != 0) {
	measurement.pixel(i,j)(2) = 1.0;
      }
    }
  }

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
  AdjacencyMatrix<M> 		view_connections;   // graph of views
  int				tree_root;          
  Eigen::Matrix<bool,M,1>	scaled;		    // is this view scaled
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
      F[v].eight_point_algorithm(measurement.view(v),
				 measurement.view(spanning_tree(v)));      
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
	for(pt = 0; pt<measurement.cols(); ++pt) {
	  if(!measurement.pixel_is_occluded(v,pt)) {
	    if(!measurement.pixel_is_occluded(parent,pt)) {
	      pix        = measurement.pixel(v,pt);
	      cross_prod = F[v].e_ij().cross(pix);
	      scale      = 
		cross_prod.dot(F[v]*measurement.pixel(parent,pt))/
		cross_prod.dot(cross_prod);
	      measurement.pixel(v,pt) *= scale;
	    } else {
	      // --- parent occluded
	      // --- so start scale at 1.0 again
	      // --- for later separation
	      measurement.pixel(v,pt) /= measurement.pixel(v,pt)(2);
	      // measurement.pixel(v,pt)(2) = -0.0;
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
bool SturmTriggsSolver<M>::operator ()(int j1, int j2) {
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
  int i,j,k;
  std::vector<int>::iterator 	cEnd = cBegin+4;
  std::vector<int>::iterator 	colIt;
  Eigen::Matrix<bool,M,4>	leafExpanded; // single pixel on its own?
  Eigen::Matrix<bool,M,1>	occluded;
  Eigen::Matrix<int,M,1>	colPosition;	
  bool				updated;
  int				colCount;

  // --- first calculate total number of columns we need
  // --- and map out expansions:
  // --- occluded     - true if row is occluded
  // --- leafExpanded - true if i'th row of jth col is expanded
  // ---                to one pixel on its own.
  
  colCount = 0;
  for(i = 0; i<M; ++i) {
    occluded(i) = false;
    j = 0;
    for(colIt = cBegin; colIt != cEnd; ++colIt) {
      leafExpanded(i,j) = false;
      if(measurement.pixel_is_occluded(i,*colIt)) {
	if(!occluded(i)) {
	  occluded(i) = true;
	  colCount += 3;
	}
      } else if(spanning_tree(i) == -1 || 
		measurement.pixel_is_occluded(spanning_tree(i),*colIt)) {
	leafExpanded(i,j) = true;
	colCount += 1;
      } else {
	leafExpanded(spanning_tree(i),j) = false;
      }
      ++j;
    }
  }
  // --- remove redundant leaf expansions
  // --- (4xleaves and leaves in occluded rows)
  for(i = 0; i<M; ++i) {
    if(occluded(i)) {
      for(j = 0; j<4; ++j) {
	if(leafExpanded(i,j)) {
	  colCount -= 1;
	}
      }
    } else {
      if(leafExpanded.row(i).sum() == 4) {
	occluded(i) = true;
	colCount -= 1;
      }
    }
  }
  
  std::cout << "Occluded = " << occluded.transpose() << std::endl;
  std::cout << "Leaf expanded = " << std::endl << leafExpanded << std::endl;
  std::cout << "ColCount = " << colCount << std::endl;

  // --- now create subspace
  // --- column by column
  subspace.resize(3*M, colCount);
  subspace.setZero();
  j = 0;
  k = 0;
  for(colIt = cBegin; colIt != cEnd; ++colIt) {
    colPosition.fill(-1);
    do {
      updated = false;
      for(i = 0; i<M; ++i) {
	std::cout << "Starting " << i << " " << *colIt << " " << j << std::endl;
	if(colPosition(i) == -1) {
	  if(!(occluded(i) && leafExpanded(i,k))) {
	    if(!measurement.pixel_is_occluded(i,*colIt) &&
	       (spanning_tree(i) == -1 ||
		measurement.pixel_is_occluded(spanning_tree(i),*colIt))) {
	      // --- start a new column
	      std::cout << "Starting new column" << std::endl;
	      subspace.pixel(i,j) = measurement.pixel(i,*colIt);
	      colPosition(i) = j;
	      updated = true;
	      ++j;
	    } else if(colPosition(spanning_tree(i)) != -1) {
	      // --- add to an existing column
	      std::cout << "Adding to column column" << std::endl;
	      colPosition(i) = colPosition(spanning_tree(i));
	      subspace.pixel(i,colPosition(i)) = measurement.pixel(i,*colIt);
	      updated = true;
	    }
	  }
	}
      }
    } while(updated);
    ++k;
  }

  std::cout << "Setting occlusions " << std::endl;
  //--- now add occlusion expansions
  for(i = 0; i<M; ++i) {
    if(occluded(i)) {
      std::cout << "Setting occlusion " << i << std::endl;
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
