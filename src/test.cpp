#include <iostream>
#include <fstream>
#include "EquivalencePartition.h"

int main() {
  EquivalencePartition<int> myEq;

  myEq.equate(1,2);
  myEq.equate(1,3);
  myEq.equate(3,5);

  std::cout << myEq;

  return(0);
}
