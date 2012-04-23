#include <iostream>
#include "transformMatrix.h"
#include "transformMatrix.cpp"

int main() {
  coord 	t(3);
  coord		Q(3);
  coord		R(3);
  coord		pi(3);
  coord		QxtxQ(3);
  coord		txQ(3);
  double 	k;

  t[0] = 1.0;
  t[1] = 2.0;
  t[2] = 3.0;

  pi[0] = -0.5;
  pi[1] = -2.345;
  pi[2] = 1.2345;

  Q = pi;
  R = pi - t;

  transformMatrix tx(cross_prod, t[0], t[1], t[2]);
  transformMatrix Qx(cross_prod, Q[0], Q[1], Q[2]);


  //  txQ = tx*Q;
  QxtxQ = Qx*tx*Q;

  std::cout << "R0 = " << (R*Q).sum() << std::endl;
  std::cout << "R1 = " << (R*QxtxQ).sum() << std::endl;
  //std::cout << "R2 = " << (R*txQ).sum() << std::endl;

  k = ((t*Q).sum() - (t*QxtxQ).sum()*(R*Q).sum()/(R*QxtxQ).sum())/(Q*Q).sum();

  std::cout << "pi = " << k*Q << std::endl;

  return(0);
}
