#include <iostream>
#include <cstdlib>

#include "transformMatrix.h"
#include "transformMatrix.cpp"

int main() {
  const double PI = 3.1415926535897932;

  double		xOffset = 1.0;
  double		sceneOffset = 2.0;
  transformMatrix 	M(intrinsic, 1.0, 1.0, 0.0, 0.0);
  transformMatrix	M1(inverse_intrinsic, 1.0, 1.0, 0.0, 0.0);
  transformMatrix 	MP(intrinsic, 1.0, 1.2, 0.01, 0.02);
  transformMatrix	MP1(inverse_intrinsic, 1.0, 1.2, 0.01, 0.02);
  transformMatrix	MP1T(inverse_intrinsic, 1.0, 1.2, 0.01, 0.02);
  transformMatrix	R(yrotate,atan(xOffset/sceneOffset));
  transformMatrix	R1(yrotate,-atan(xOffset/sceneOffset));
  transformMatrix	S(cross_prod,xOffset,0.0,0.0);
  transformMatrix	F;
  coord		t(3);
  coord		Q(3);
  coord		q(3);
  coord		qp(3);
  coord		Q1(3);
  coord		Q2(3);
  int			i;

  std::cout << "# t=(" << xOffset << ",0,0)" << "M=(1,1,0,0) MP=(1,1.2,0.01,0.02) R=yrotate(atan(0.5))" << std::endl;

  std::cout.precision(15);
  MP1T.transpose();
  F = MP1T;
  F *= R;
  F *= S;
  F *= M1;

  t[0] = xOffset;
  t[1] = 0.0;
  t[2] = 0.0;

  for(i = 0; i<100; ++i) {
    Q[0] = rand()*1.0/RAND_MAX - 0.5;
    Q[1] = rand()*1.0/RAND_MAX - 0.5;
    Q[2] = rand()*1.0/RAND_MAX -0.5 + sceneOffset;
    
    q  = M*Q;
    qp = MP*(R*(Q - t));
    // Q1 = M1*q;
    // Q2 = t + R1*(MP1*qp);
    
    // std::cout << (qp*(MP1T*(R*(S*(M1*q))))).sum() << std::endl;
    //std::cout << (qp*(F*q)).sum() << std::endl;

    std::cout << std::scientific << "0 " << q[0]/q[2] << " " << q[1]/q[2]
	      << "\t0 " << qp[0]/qp[2] << " " << qp[1]/qp[2]
	      << "\t1.0" << std::endl;
  }

  // std::cout << F;
}
