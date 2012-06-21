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
#ifndef ADJACENCYMATRIX_H
#define ADJACENCYMATRIX_H

#include "stdincludes.h"

///////////////////////////////////////////////////////////////////////////////
/// Class to represent a directed graph as an adjacency matrix
///////////////////////////////////////////////////////////////////////////////
template<int M>
class AdjacencyMatrix : public Eigen::Matrix<double,M,M> {
public:

  void shortest_paths_tree(int, Eigen::Matrix<int,M,1> &,
			   Eigen::Matrix<double,M,1> &);
  void shortest_spanning_tree(int &, Eigen::Matrix<int,M,1> &);

};


//////////////////////////////////////////////////////////////////////////////
/// Calculates the shortest paths spanning tree whose root is 'root',
/// putting the resulting tree in 'predecessor' and distances in
/// 'cost'. Uses Dijkstra's algorithm.
///
/// \param root the node we want to be the root of the resulting tree
/// \param predecessor vector to hold the resulting tree. Each element
/// holds the index of its parent.
/// \param cost vector holding the cost to go from the root to this node.
//////////////////////////////////////////////////////////////////////////////
template<int M>
void AdjacencyMatrix<M>::shortest_paths_tree(int root, 
					  Eigen::Matrix<int,M,1> &predecessor,
					  Eigen::Matrix<double,M,1> &cost) {
  int i;
  int freeVertices;
  int closestNode;
  double dist;
  Eigen::Matrix<bool,M,1>     isFree;

  // --- initialise
  for(i=0; i<M; i++) {
    predecessor(i) = -1;
    cost(i) = std::numeric_limits<double>::infinity();
    isFree(i) = true;
  }
  cost(root)= 0;

  for(freeVertices = M; freeVertices > 0; --freeVertices) {
    // --- find closest node
    dist = std::numeric_limits<double>::infinity();
    for(i = 0; i<M; i++) {
      if(isFree(i) && cost(i) < dist) {
        dist = cost(i);
        closestNode = i;
      }
    }
    isFree[closestNode] = false;

    // --- relax neighbours of closestNode
    for(i=0; i<M; i++) {
      if(isFree[i] && cost[i] > dist+(*this)(closestNode,i)) {
        cost[i] = dist+(*this)(closestNode,i);
        predecessor[i] = closestNode;
      }
    }
  }
}


//////////////////////////////////////////////////////////////////////////////
/// Finds the tree that, among all spanning trees of *this, has the shortest
/// maximum depth, putting the resulting tree in 'predecessor' and the index of
/// the root node in 'root'.
///
/// \param root set to the root of the shortest tree \param
/// predecessor set to the resulting tree. Each element holds the
/// index of its parent.
//////////////////////////////////////////////////////////////////////////////
template<int M>
void AdjacencyMatrix<M>::shortest_spanning_tree(int &root, 
					Eigen::Matrix<int,M,1> &predecessor) {
  Eigen::Matrix<double,M,1> cost;
  Eigen::Matrix<int,M,1> rooted_tree;
  double min_cost;
  int r,i;

  min_cost = std::numeric_limits<double>::infinity();
  for(r = 0; r<M; ++r) {
    shortest_paths_tree(r, rooted_tree, cost);
    for(i=0; i<M; ++i) {
      if(cost(i) <= min_cost) {
	min_cost = cost(i);
	root = r;
	predecessor = rooted_tree;
      }
    }
  }
}

#endif
