#include <fstream>
#include <cmath>
#include "FSolver.h"


///////////////////////////////////////////////////////////////////////////////
//
// default place to put the residual after fitting.
//
///////////////////////////////////////////////////////////////////////////////
double FSolver::residual = 0.0;


///////////////////////////////////////////////////////////////////////////////
//
// Solve from correspondences in a file
//
///////////////////////////////////////////////////////////////////////////////
FSolver::FSolver(const char *filename, double &r) {
  load(filename);
  r = solve();
}


///////////////////////////////////////////////////////////////////////////////
//
// Solve from a vector of correspondences.
// if v1 == -1 then ignore camera IDs and solve for all correspondences
// if v1 != -1 then solve only for correspondences between v1 and v2
//
///////////////////////////////////////////////////////////////////////////////
FSolver::FSolver(std::vector<Correspondence> &C,int v1, int v2, double &r) {
  int 	n;
  coord	pixel(3);

  pixel[2] = 1.0;
  for(n = 0; n<C.size(); ++n) {
    if((v1 == -1) || 
       (C[n].i == v1 && C[n].j == v2) || 
       (C[n].i == v2 && C[n].j == v1)) {
      pixel[0] = C[n].xi;
      pixel[1] = C[n].yi;
      q.push_back(pixel);
      pixel[0] = C[n].xj;
      pixel[1] = C[n].yj;
      qp.push_back(pixel);
    }
  }
  e.resize(data_size());
  r = solve();
}

///////////////////////////////////////////////////////////////////////////////
//
// Calculate the fundamental matrix and set to F, returning the residual of
// the fit.
//
///////////////////////////////////////////////////////////////////////////////
double FSolver::solve() {
  Frprmn<FSolver>	mySolver(*this, 1e-7);	// solver
  VecDoub 		p(9);			// parameters for solver
  int			n;
  transformMatrix	M1(inverse_intrinsic,1.0,1.0,0.0,0.0);

  // --- initial guess
  // -----------------
  F() = M1;

  for(n = 0; n<9; ++n) {
    p[n] = F()[n/3][n%3];
  }
  p = mySolver.minimize(p);
  for(n = 0; n<9; ++n) {
    F()[n/3][n%3] = p[n];
  }
  return((*this)(p));
}


///////////////////////////////////////////////////////////////////////////////
//
// Load set of correspondences from a file
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

  //std::cout << "in operator(p)" << std::endl;

  //std::cout << "p = ";
  //for(n=0; n<p.size(); ++n) {
  //  std::cout << p[n] << " ";
  //}
  //std::cout << std::endl;

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
  //std::cout << "alpha = " << e2 << std::endl;

  // --- add scale forcing
  // ---------------------
  e2 +=(p[2]*p[2] + p[5]*p[5] + p[8]*p[8] - 0.694722)*
    (p[2]*p[2] + p[5]*p[5] + p[8]*p[8] - 0.694722);

  // --- add determinant forcing
  // ---------------------------
  /*****
  e2 += pow(p[0]*p[4]*p[8] + p[1]*p[5]*p[6] + p[2]*p[3]*p[7] -
    p[2]*p[4]*p[6] - p[1]*p[3]*p[8]  - p[0]*p[5]*p[7],2);
  *****/

  // std::cout << "done, e2 = " << e2 << std::endl;

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


  // std::cout << "in df" << std::endl;

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

  // --- add determinant forcing
  // ---------------------------
  /******
  t = 2.0*(p[0]*p[4]*p[8] + p[1]*p[5]*p[6] + p[2]*p[3]*p[7] -
    p[2]*p[4]*p[6] - p[1]*p[3]*p[8]  - p[0]*p[5]*p[7]);
  deriv[0] += t*(p[4]*p[8] - p[5]*p[7]);
  deriv[1] += t*(p[5]*p[6] - p[3]*p[8]);
  deriv[2] += t*(p[3]*p[7] - p[4]*p[6]);
  deriv[3] += t*(p[2]*p[7] - p[1]*p[8]);
  deriv[4] += t*(p[0]*p[8] - p[2]*p[6]);
  deriv[5] += t*(p[1]*p[6] - p[0]*p[7]);
  deriv[6] += t*(p[1]*p[5] - p[2]*p[4]);
  deriv[7] += t*(p[2]*p[3] - p[0]*p[5]);
  deriv[8] += t*(p[0]*p[4] - p[1]*p[3]);
  *****/

  /****
  for(n=0; n<deriv.size(); ++n) {
    std::cout << deriv[n] << " ";
  }
  std::cout << std::endl;
  ***/
}


///////////////////////////////////////////////////////////////////////////////
//
// calculate the normalised vector pointing towards the epipole (intersection
// between the projection plane and the line connecting the projection centres
// of the two cameras) for a given fundamental matrix, F.
//
///////////////////////////////////////////////////////////////////////////////
/********
coord FSolver::epipole(transformMatrix &F) {
  MatDoub 	FP(3,3);
  coord		e(3);
  int i,j;

  for(i=0; i<3; ++i) {
    for(j=0; j<3; ++j) {
      FP[i][j] = F[i][j];
    }
  }

  SVD	decomposition(FP);
  FP = decomposition.nullspace(1e-3);
  if(FP.ncols() > 0) {
    e[0] = FP[0][0];
    e[1] = FP[1][0];
    e[2] = FP[2][0];
  } else {
    throw("can't find epipole");
  }
  return(e);
}
********/
