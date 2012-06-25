#include <fstream>
#include "sceneSolver.h"



///////////////////////////////////////////////////////////////////////////////
//
// setup reasonable default conditions
//
///////////////////////////////////////////////////////////////////////////////
sceneParams::sceneParams() : t(3),
			     N(intrinsic, 1.0, 1.0), 
			     M(intrinsic, 1.0, 1.0) {
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
std::ostream &operator <<(std::ostream &out, sceneParams &p) {
  out << "t = (" 
      << p.t[0] << ", " 
      << p.t[1] << ", " 
      << p.t[2] << ")" << std::endl;
  out<<"   ("<< p.N[0][0] <<" "<< p.N[0][1] <<" "<< p.N[0][2] <<")"<<std::endl;
  out<<"N =("<< p.N[1][0] <<" "<< p.N[1][1] <<" "<< p.N[1][2] <<")"<<std::endl;
  out<<"   ("<< p.N[2][0] <<" "<< p.N[2][1] <<" "<< p.N[2][2] <<")"<<std::endl;

  out<<"   ("<< p.M[0][0] <<" "<< p.M[0][1] <<" "<< p.M[0][2] <<")"<<std::endl;
  out<<"M =("<< p.M[1][0] <<" "<< p.M[1][1] <<" "<< p.M[1][2] <<")"<<std::endl;
  out<<"   ("<< p.M[2][0] <<" "<< p.M[2][1] <<" "<< p.M[2][2] <<")"<<std::endl;

  return(out);
}

///////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////
void sceneSolver::solve() {
  Frprmn<sceneSolver>	mySolver(*this, 1e-15);	// solver
  VecDoub 		p(2*data_size()+16);	// parameters for solver
  int			i;

  // --- initial guess
  // -----------------
  params.t[0] = 1.0;
  params.t[1] = 0.0;
  params.t[2] = 0.0;
  for(i = 0; i<data_size(); ++i) {
    params.z[i]  = 1.0;
    params.zp[i] = 1.0;
  }
  
  packParams(p);
  p = mySolver.minimize(p);
  std::cout << "solved" << std::endl;
  unpackParams(p);
  std::cout << "unpacked params" << std::endl;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void sceneSolver::packParams(VecDoub &p) {
  int n;

  p[0] = params.t[0];
  p[1] = params.t[1];
  p[2] = params.t[2];
  for(n = 0; n<9; ++n) {
    p[n+3] = params.N[n/3][n%3];
  }
  p[12] = params.M[0][0];
  p[13] = params.M[1][1];
  p[14] = params.M[0][2];
  p[15] = params.M[1][2];
  for(n = 0; n<data_size(); ++n) {
    p[16+2*n] = params.zp[n];
    p[17+2*n] = params.z[n];
  }
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void sceneSolver::unpackParams(VecDoub &p) {
  int n;
  
  params.t[0] = p[0];
  params.t[1] = p[1];
  params.t[2] = p[2];
  for(n = 0; n<9; ++n) {
    params.N[n/3][n%3] = p[3+n];
  }
  params.M[0][0] = p[12];
  params.M[1][1] = p[13];
  params.M[0][2] = p[14];
  params.M[1][2] = p[15];
  for(n = 0; n<data_size(); ++n) {
    params.zp[n] = p[16+2*n]; 
    params.z[n] = p[17+2*n];
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// Load set of correspondences from a file
//
///////////////////////////////////////////////////////////////////////////////
bool sceneSolver::load(const char *filename) {
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
    q.push_back(pro_pixel);
    qp.push_back(cam_pixel);
  }
  params.z.resize(data_size());
  params.zp.resize(data_size());
  cam_pixel = 0.0;
  e.resize(data_size(),cam_pixel);


  return(true);
}


///////////////////////////////////////////////////////////////////////////////
//
// Calcuate e and etot from q, qp and params;
//
///////////////////////////////////////////////////////////////////////////////
void sceneSolver::calcErrors(VecDoub_I &p) {
  int n;
  transformMatrix	N;
  transformMatrix	M;
  coord			t(3);

  std::cout << "Calculating errors" << std::endl;

  t[0] = p[0];
  t[1] = p[1];
  t[2] = p[2];
  for(n = 0; n<9; ++n) {
    N[n/3][n%3] = p[3+n];
  }
  M[0][0] = p[12];
  M[1][1] = p[13];
  M[0][2] = p[14];
  M[1][2] = p[15];

  for(n = 0; n<data_size(); ++n) {
    e[n] = t + p[2*n+16]*(N*qp[n]) - p[2*n+17]*(M*q[n]);
  }
  std::cout << "Done" << std::endl;

}


///////////////////////////////////////////////////////////////////////////////
//
// returns the sum of squared errors
//
// e^2 = sum_n (t + z'_nNq'_n - z_nMq_n)^T . (t + z'_nNq'_n - z_nMq_n)
//	+ ((sum_n z'_n) - n)^2 + ((sum_n z_n) - n)^2
//
// the final two terms break the scale invariance
// and enforce unity average distance from camera to scene and average distance
// from projector to scene.
//
///////////////////////////////////////////////////////////////////////////////
Doub sceneSolver::operator()(VecDoub_I &p) {
  int 		n;
  double 	e2;

  std::cout << "in operator(p)" << std::endl;

  //for(n=0; n<p.size(); ++n) {
  //std::cout << p[n] << " ";
  //}
  //std::cout << std::endl;


  calcErrors(p);
  e2 = 0.0;
  z_bar = -data_size();
  zp_bar = -data_size();
  for(n=0; n<data_size(); ++n) {
    e2 += e[n][0]*e[n][0] + e[n][1]*e[n][1] + e[n][2]*e[n][2];
    z_bar += p[17 + 2*n];
    zp_bar += p[16 + 2*n];
  }

  std::cout << "z_bar = " << z_bar << std::endl;
  std::cout << "zp_bar = " << zp_bar << std::endl;

  e2 += z_bar*z_bar + zp_bar*zp_bar;

  std::cout << "done " << e2 << std::endl;

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
void sceneSolver::df(VecDoub_I &p, VecDoub_O &deriv) {
  coord		etot(3);	// total error vector
  int		i,j,n;


  std::cout << "in df" << std::endl;

  // --- de/dt
  // ---------
  etot = 0.0;
  for(n=0; n<data_size(); ++n) {
    etot += e[n];
  }
  deriv[0] = etot[0];
  deriv[1] = etot[1];
  deriv[2] = etot[2];

  // --- de/N_{ij}
  // -------------
  for(i = 0; i<3; ++i) {
    for(j = 0; j<3; ++j) {
      deriv[3+3*i+j] = 0.0; 
      for(n=0; n<data_size(); ++n) {
	deriv[3+3*i+j] += e[n][i] * qp[n][j] * p[16+2*n];
      }
    }
  }

  // --- de/M_{ij}
  // -------------
  deriv[12] = 0.0;
  deriv[13] = 0.0;
  deriv[14] = 0.0;
  deriv[15] = 0.0;
  for(n=0; n<data_size(); ++n) {
    deriv[12] -= e[n][0] * q[n][0] * p[17+2*n];
    deriv[13] -= e[n][1] * q[n][1] * p[17+2*n];
    deriv[14] -= e[n][0] * q[n][2] * p[17+2*n];
    deriv[15] -= e[n][1] * q[n][2] * p[17+2*n];
  }

  // --- de/z'_n
  // -----------
  for(n=0; n<data_size(); ++n) {
    deriv[16+2*n] = etot[0]*(p[3]*qp[n][0] + p[4]*qp[n][1] + p[5]*qp[n][2]) 
      +	etot[1]*(p[6]*qp[n][0] + p[7]*qp[n][1] + p[8]*qp[n][2]) 
      +	etot[2]*(p[9]*qp[n][0] + p[10]*qp[n][1] + p[11]*qp[n][2])
      + zp_bar;
  }

  // --- de/z_n
  // ----------
  for(n=0; n<data_size(); ++n) {
    deriv[17+2*n] = etot[0]*(p[12]*q[n][0] + p[14]*q[n][2])
      + etot[1]*(p[13]*q[n][1] + p[15]*q[n][2])
      + etot[2]*q[n][2]
      + z_bar;
  }

  for(n=0; n<deriv.size(); ++n) {
    std::cout << deriv[n] << " ";
  }
  std::cout << std::endl;
}
