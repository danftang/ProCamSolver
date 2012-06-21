#include <iostream>
#include <limits>
#include <Eigen/Core>
#include <vector>

using namespace Eigen;
//using namespace std;

void Dijkstra(const MatrixXd &adjMatrix, MatrixXd &result) {
  int i;
  int freeVertices;
  int closestNode;
  double dist;
  int numOfVertices = adjMatrix.cols();
  std::vector<int>	predecessor(numOfVertices);
  std::vector<double> 	cost(numOfVertices);
  std::vector<bool>	isFree(numOfVertices);

  // --- initialise
  for(i=0;i<numOfVertices;i++) {
    predecessor[i] = -1;
    cost[i] = std::numeric_limits<double>::infinity();
    isFree[i] = true;
  }
  cost[2]= 0;

  for(freeVertices = numOfVertices; freeVertices > 0; --freeVertices) {
    // --- find closest node
    dist = std::numeric_limits<double>::infinity();
    for(i = 0; i<numOfVertices; i++) {
      if(isFree[i] && cost[i] < dist) {
	dist = cost[i];
	closestNode = i;
      }
    }
    isFree[closestNode] = false;

    // --- relax neighbours of closestNode
    for(i=0;i<numOfVertices;i++) {
      if(isFree[i] && cost[i] > dist+adjMatrix(closestNode,i)) {
	cost[i] = dist+adjMatrix(closestNode,i);
	predecessor[i] = closestNode;
      }
    }
  }

  result.resize(numOfVertices, numOfVertices);
  result.fill(0.0);
  for(i = 0; i<numOfVertices; ++i) {
    std::cout << predecessor[i] << " -> " << i << std::endl;
    if(predecessor[i] != -1) {
      result(predecessor[i],i) = 1.0;
    }
  }
}



int main(){
  const int N = 10;
  MatrixXd myMatrix(N,N);
  MatrixXd mintree(N,N);
  int i,j;
  double cost;

  for(i = 0; i<N; ++i) {
    myMatrix(i,i) = std::numeric_limits<double>::infinity();
    for(j = 0; j<i; ++j) {
      cost = rand()*1.1/RAND_MAX + 0.5;
      if(cost > 1.0) cost = std::numeric_limits<double>::infinity();
      myMatrix(i,j) = myMatrix(j,i) = cost;
    }
  }

  Dijkstra(myMatrix, mintree);
  std::cout << myMatrix << std::endl;
  std::cout << "Min tree = " << std::endl;
  std::cout << mintree << std::endl;
}



