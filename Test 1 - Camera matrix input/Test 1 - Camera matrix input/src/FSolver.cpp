#include <fstream>
#include <cmath>
#include "FSolver.h"
#include "numerical.h"


///////////////////////////////////////////////////////////////////////////////
//
// Calculate the fundamental matrix and set to F
//
///////////////////////////////////////////////////////////////////////////////
void FSolver::solve(transformMatrix &F) {
  Frprmn<FSolver>	mySolver(*this, 1e-6);	// solver
  VecDoub 		p(9);			// parameters for solver
  int			n;
  transformMatrix	M1(inverse_intrinsic,1.0,1.0,0.0,0.0);
  transformMatrix	MP1T(inverse_intrinsic,1.0,1.2,0.01,0.02);
  transformMatrix	S(cross_prod,1.0,0.0,0.0);
  transformMatrix	R(yrotate,atan(0.5));

  // --- initial guess
  // -----------------
  MP1T.transpose();
  //F = MP1T;
  //F *= R;
  //F *= S;
  //F *= M1;
  F = M1;

  for(n = 0; n<9; ++n) {
    p[n] = F[n/3][n%3];
  }
  p = mySolver.minimize(p);
  for(n = 0; n<9; ++n) {
    F[n/3][n%3] = p[n];
  }
}


///////////////////////////////////////////////////////////////////////////////
//
// Load set of correspondences from a file in format
// Index1 x y Index2 x' y' weight
//
///////////////////////////////////////////////////////////////////////////////
bool FSolver::load(const char *filename) {
  double 	weight;
  coord		cam_pixel(3);
  coord		pro_pixel(3);
  ifstream 	in(filename);

  if(!in) return(false);
  q.clear();
  qp.clear();
  while(in.peek() == '#') {
    in.ignore(1024,'\n');
  }
  while(in >> cam_pixel[2] >> cam_pixel[0] >> cam_pixel[1]
	   >> pro_pixel[2] >> pro_pixel[0] >> pro_pixel[1] >> weight) {
    cam_pixel[2] = 1.0;
    pro_pixel[2] = 1.0;
    qp.push_back(pro_pixel);
    q.push_back(cam_pixel);
    //    std::cout << q.back() << " " << qp.back() << std::endl;
  }
  e.resize(data_size());

  return(true);
}


///////////////////////////////////////////////////////////////////////////////
//
// returns the sum of squared errors
//
///////////////////////////////////////////////////////////////////////////////
Doub FSolver::operator()(VecDoub_I &p) {
  int 		n,i,j;
  double 	e2;

  e2 = 0.0;
  for(n=0; n<data_size(); ++n) {
    e[n] = 0.0;
    for(i=0; i<3; ++i) {
      for(j=0; j<3; ++j) {
	e[n] += qp[n][i] * p[3*i+j] * q[n][j];
      }
    }
    e2 += e[n]*e[n];
  }

  // --- add scale forcing
  // ---------------------
  e2 +=(p[2]*p[2] + p[5]*p[5] + p[8]*p[8] - 0.694722)*
    (p[2]*p[2] + p[5]*p[5] + p[8]*p[8] - 0.694722);

  return(e2);
}


///////////////////////////////////////////////////////////////////////////////
//
// calculates the rate of change of e^2 with each parameter in 'p' at 'p'
// returning the result in 'deriv'.
//
// assumes e is up to date with p
//
///////////////////////////////////////////////////////////////////////////////
void FSolver::df(VecDoub_I &p, VecDoub_O &deriv) {
  coord		etot(3);	// total error vector
  int		i,j,n;

  for(n=0; n<data_size(); ++n) {
    e[n] = 0.0;
    for(i=0; i<3; ++i) {
      for(j=0; j<3; ++j) {
	e[n] += qp[n][i] * p[3*i+j] * q[n][j];
      }
    }
  }


  for(i = 0; i<3; ++i) {
    for(j = 0; j<3; ++j) {
      deriv[3*i+j] = 0.0;
      for(n=0; n<data_size(); ++n) {
	deriv[3*i+j] += 2.0 * e[n] * qp[n][i] * q[n][j];
      }
    }
  }

  // --- add scale forcing
  // ---------------------
  double t;
  t = 4.0*(p[2]*p[2] + p[5]*p[5] + p[8]*p[8] - 0.694722);
  deriv[2] += t*p[2];
  deriv[5] += t*p[5];
  deriv[8] += t*p[8];
}


///////////////////////////////////////////////////////////////////////////////
//
// calculate the normalised vector pointing towards the epipole (intersection
// between the projection plane and the line connecting the projection centres
// of the two cameras) for a given fundamental matrix, F.
//
///////////////////////////////////////////////////////////////////////////////
coord FSolver::epipole(transformMatrix &F) {
  MatDoub 	FP;
  coord		e(3);
  int i,j;

  for(i=0; i<3; ++i) {
    for(j=0; j<3; ++j) {
      FP[i][j] = F[i][j];
    }
  }
  SVD	decomposition(FP);
  FP = decomposition.nullspace(0.0);
  e[0] = FP[0][0];
  e[1] = FP[1][0];
  e[2] = FP[2][0];
  return(e);
}
