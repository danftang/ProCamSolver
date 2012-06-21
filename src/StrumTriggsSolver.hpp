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
StrumTriggsSolver<M>::StrumTriggsSolver(MotionMatrix<M> &P,
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
void StrumTriggsSolver<M>::solve(const MeasurementMatrix<M> &E) {
  Eigen::JacobiSVD<typename MeasurementMatrix<M>::Base> 	svd;
  ImageTransform<M>				denormalisation;

  measurement = E;
  measurement.normalise(denormalisation);
  approx_scale();
  factorise_with_occlusions(motion, shape);
  measurement.scale_and_fill(motion, shape);
  svd.compute(measurement, Eigen::ComputeThinV | Eigen::ComputeThinU);
  motion = svd.matrixU().leftCols(4);
  motion *= svd.singularValues().topRows(4).asDiagonal();

  std::cout << "Strum-Triggs singular values are:" << std::endl;
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
void StrumTriggsSolver<M>::factorise_with_occlusions(MotionMatrix<M> &motion, 
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
  while(c<measurement.cols()) {
    // --- create constraint
    subspace.delete_columns(0,subspace.cols());
    for(j = 0; c<measurement.cols() && j<4; ++j) {
      // --- add a column to subspace
      j0 = subspace.insert_column_and_expand(measurement.col(c));
      ++c;
      srank = subspace.rank(singular_cutoff);
      while(c<measurement.cols() && srank < subspace.cols() && srank < 3*M) {
	subspace.delete_columns(j0, subspace.cols()-j0);
	j0 = subspace.insert_column_and_expand(measurement.col(c));
	++c;
      }
    }

    // --- add nullspace of constraint to notMotion
    if(j == 4 && subspace.cols() < 3*M) {
      notMotion.new_columns(3*M-subspace.cols());
      notMotion.rightCols(3*M-subspace.cols()) = subspace.jacobiSvd(Eigen::ComputeFullU).matrixU().rightCols(3*M-subspace.cols());
      nmRank = notMotion.rank(singular_cutoff);
    }
  }
  
  if(nmRank < 3*M-4) {
    throw("Too many occlusions to fill in the Measurement matrix.");
  }

  // --- motion is approximate nullspace of the complement
  motion = notMotion.jacobiSvd(Eigen::ComputeFullU).matrixU().rightCols(4);

  shape.solve(measurement, motion);
}


///////////////////////////////////////////////////////////////////////////////
/// Sets G to the graph whose vertices represent views and whose edge
/// weights are a measure of the expected error in the fundamental
/// matrix between views. G is in the form of an adjacency matrix. 
///////////////////////////////////////////////////////////////////////////////
template<int M>
void StrumTriggsSolver<M>::
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
void StrumTriggsSolver<M>::approx_scale() {
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
	  if(measurement.pixel(v,pt)(2) != 0.0) {
	    if(measurement.pixel(parent,pt)(2) != 0.0) {
	      pix = measurement.pixel(v,pt);
	      cross_prod = F[v].e_ij().cross(pix);
	      scale = cross_prod.dot(F[v]*measurement.pixel(parent,pt))/
		cross_prod.dot(cross_prod);
	      measurement.pixel(v,pt) *= scale;
	    } else {
	      // --- parent missing so all descendants are unscaled
	      measurement.pixel(v,pt)(2) = -0.0;
	    }
	  }
	}
	scaled(v) = true;
	updated = true;
      }
    }
  }
}


