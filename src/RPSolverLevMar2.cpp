#include "RPSolverLevMar2.h"
#include "nrlib/gaussj.h"
#include "nrlib/fitmrq.h"

void RPSolver::operator()(const Doub x, VecDoub_I &p, Doub &y, VecDoub_O &dy_dp) {
  int eqSet = (int)(x + 0.5);

  y = (*this)(eqSet, p);

  df(eqSet,p,dy_dp);
}


///////////////////////////////////////////////////////////////////////////////
//
// Default place to put the residual of the fit after solving
//
///////////////////////////////////////////////////////////////////////////////
double RPSolver::residual = 0.0;


///////////////////////////////////////////////////////////////////////////////
//
// set to false to fix radial distortion of all views to 0.0
//
///////////////////////////////////////////////////////////////////////////////
bool RPSolver::radial = false;


///////////////////////////////////////////////////////////////////////////////
//
// cor 			- pixel correlations
// r			- contains residual of fit on return
//
///////////////////////////////////////////////////////////////////////////////
RPSolver::RPSolver(std::vector<Correspondence> &cor, double &r) {
  init_views(cor);
  r = solve();
}


///////////////////////////////////////////////////////////////////////////////
//
// cor 			- pixel correlations
// fixedIntrinsics	- prior known intrinsic matrices for certain views
// r			- contains residual of fit on return
//
///////////////////////////////////////////////////////////////////////////////
RPSolver::RPSolver(std::vector<Correspondence> &cor,
		 std::map<int, transformMatrix> &fixedIntrinsics,
		 double &r) {
  init_intrinsics(fixedIntrinsics);
  init_views(cor);
  r = solve();
}


///////////////////////////////////////////////////////////////////////////////
//
// cor 			- pixel correlations
// initials 		- should contain initial guesses at translation and
//                        rotation of some (not necessarily all) views
// 			  and any prior known intrinsic matrices.
//			  If fx is non-zero in a view then the intrinsics
//			  for that view are held fixed.
// r			- contains residual of fit on return
//
///////////////////////////////////////////////////////////////////////////////
RPSolver::RPSolver(std::vector<Correspondence> &cor,
		 std::map<int, RadialViewProjectionMatrix> &initials,
		 double &r) {
  init_VPM(initials);
  init_views(cor);
  r = solve();
}



///////////////////////////////////////////////////////////////////////////////
//
// initialises the intrinsic matrices of *this and marks them as fixed
//
///////////////////////////////////////////////////////////////////////////////
void RPSolver::init_intrinsics(std::map<int, transformMatrix> &initials) {
  std::map<int, transformMatrix>::iterator it;
  coord	t_guess(3);
  t_guess[0] = 0.03; t_guess[1] = 0.02; t_guess[2] = 0.01; 
  
  for(it = initials.begin(); it != initials.end(); ++it) {
    (*this)[it->first].fx = it->second[0][0];
    (*this)[it->first].fy = it->second[1][1];
    (*this)[it->first].cx = it->second[0][2];
    (*this)[it->first].cy = it->second[1][2];
    (*this)[it->first].t  = t_guess * (double)it->first;
    (*this)[it->first].rot = t_guess*(double)it->first;
    paramMask.insert(it->first);
  }
}


///////////////////////////////////////////////////////////////////////////////
//
// initialises the view-projection matrices of this to allow initial guesses
// at the correct values. intrinsic matrices are taken to be fixed if fx
// is not equal to zero, otherwise they are variable.
//
///////////////////////////////////////////////////////////////////////////////
void RPSolver::init_VPM(std::map<int, RadialViewProjectionMatrix> &initials) {
  std::map<int, RadialViewProjectionMatrix>::iterator it;

  for(it = initials.begin(); it != initials.end(); ++it) {
    insert(*it);
    if(vpmat(it).fx != 0.0) paramMask.insert(it->first);
  }
}


///////////////////////////////////////////////////////////////////////////////
//
// Initialise the map from view IDs to positions in the parameter vector
//
///////////////////////////////////////////////////////////////////////////////
void RPSolver::init_views(std::vector<Correspondence> &cor) {
  iterator      view;
  int           i;
  coord         t_guess(3);
  GPixel        pixel1, pixel2;

  // --- create views for all correspondences
  t_guess[0] = 0.03; t_guess[1] = 0.02; t_guess[2] = 0.01; 
  for(i = 0; i< cor.size(); ++i) {
    // --- setup ViewProjectionMatrices
    if(find(cor[i].i) == end()) {
      // --- initial guess at translation
      (*this)[cor[i].i].t = t_guess*(double)cor[i].i;
      (*this)[cor[i].i].rot = t_guess*(double)cor[i].i;
    }
    if(find(cor[i].j) == end()) {
      (*this)[cor[i].j].t = t_guess*(double)cor[i].j;
      (*this)[cor[i].j].rot = t_guess*(double)cor[i].j;
    }

    // --- setup pixel equivalence partitioning
    pixel1.id = cor[i].i;
    pixel1.x = cor[i].xi;
    pixel1.y = cor[i].yi;
    pixel2.id = cor[i].j;
    pixel2.x = cor[i].xj;
    pixel2.y = cor[i].yj;
    C.equate(pixel1, pixel2);
    std::cout << cor[i].i << " " << cor[i].j << std::endl;
  }

  // view 0 has reference intrinsics
  paramMask.insert(0);
  (*this)[0].fx = 1.0;
  (*this)[0].fy = 1.0;
  (*this)[0].cx = 0.0;
  (*this)[0].cy = 0.0;
  
  i = 0;
  for(view = begin(); view != end(); ++view) {
    paramBase[id(view)] = i;
    if(id(view) > 0) {
      // view 0 has fixed r and t by definition
      i += 6;
    }
    /****
    if(id(view) == 1) {
      // view 1 has fixed scale by definition
      i += 5;
    }
    *****/
    if(has_variable_M(id(view))) {
      // intrinsic matrix is not fixed
      vpmat(view).fx = 1.0;
      i += 4;
    }
    if(radial) ++i;
  }
  pn = i;
}


///////////////////////////////////////////////////////////////////////////////
//
// Find the camera matrices which minimise the sum of squared errors
//
///////////////////////////////////////////////////////////////////////////////
double RPSolver::solve() {
  //  Frprmn<RPSolver>	mySolver(*this, 1e-6);
  VecDoub		p(pn);
  VecDoub		xx(C.size());
  VecDoub		yy(C.size());
  VecDoub		ss(C.size());
  int 			n,b;
  iterator		view;
  double		r;


  for(n=0; n<xx.size(); ++n) {
    xx[n] = n;
    yy[n] = 0.0;
    ss[n] = 1.0;
  }

  vpmatrices_to_params(p);

  Fitmrq<RPSolver> myFit(xx, yy, ss, p, *this, 1e-16);
  myFit.fit();

  params_to_vpmatrices(p);
  r = 0.0;

  // --- re-scale the translations
  // --- distance view1->view 2 = 1.0
  // --------------------------------
  double scale = 1.0/sqrt(((*this)[1].t * (*this)[1].t).sum());
  std::cout << "Scaling by " << scale << std::endl;
  for(view = begin(); view != end(); ++view) {
    vpmat(view).t = vpmat(view).t * scale;
  }

  return(r);
}


///////////////////////////////////////////////////////////////////////////////
//
// returns the error for a given set of parameters
//
///////////////////////////////////////////////////////////////////////////////
double RPSolver::operator()(int n, VecDoub_I &p) {
  double 	e2;
  coord		pi1(3);
  coord		pi2(3);
  coord 	pi_sum(3);
  coord 	pi2_sum(3);
  int		N;
  EquivalencePartition<GPixel>::member_iterator pixel1, pixel2;

  params_to_vpmatrices(p);
  pi_sum = 0.0;
  pi2_sum = 0.0;
  N = 0;
  for(pixel2= ++(C[n].begin()); pixel2 != C[n].end(); ++pixel2) {
    for(pixel1 = C[n].begin(); pixel1 != pixel2; ++pixel1) {
      closest_points(*pixel1, *pixel2, pi1, pi2);
      pi_sum += pi1;
      pi_sum += pi2;
      pi2_sum += pi1*pi1;
      pi2_sum += pi2*pi2;
      N += 2;
    }
  }
  e2 = (pi2_sum - pi_sum*pi_sum/N).sum()/((pi_sum*pi_sum).sum()/(N*N));

  /***
  double scale = 0.0;
  for(iterator view = begin(); view != end(); ++view) {
    scale += sqrt((vpmat(view).t * vpmat(view).t).sum());
  }
  scale -= size();
  e2 += scale*scale;
  ***/

  return(e2);
}


///////////////////////////////////////////////////////////////////////////////
//
// calculate derivative of error with each parameter using finite difference
// approximation
//
///////////////////////////////////////////////////////////////////////////////
void RPSolver::df(int eqSet, VecDoub_I &p, VecDoub_O &deriv) {
  double e_old;
  double e;
  double x_old;
  VecDoub ph = p;
  double dx;
  const double delta = 1e-8;

  // ---- derivative by finite difference
  // ------------------------------------
  e_old = (*this)(eqSet,p);
  // std::cout.precision(10);
  // std::cout << "original\n" << (*this)[2] << std::endl;
  for(int i=0; i<deriv.size(); ++i) {
    x_old = p[i];
    dx = delta*fabs(x_old);
    if(dx == 0.0) dx = delta;
    ph[i] = x_old + dx;
    dx = ph[i] - x_old;
    e = (*this)(eqSet,ph);
    ph[i] = x_old;
    deriv[i] = (e - e_old)/dx;
  }
  if(eqSet == 0)  {
    std::cout << e_old << std::endl;
    iterator view;
    for(view = begin(); view != end(); ++view) {
      std::cout << view->second << std::endl;
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
//
// copies the parameters in p to this map of view projection matrices
//
///////////////////////////////////////////////////////////////////////////////
void RPSolver::params_to_vpmatrices(VecDoub_I &p) {
  iterator	view;
  int 		i;

  i = 0;
  for(view = begin(); view != end(); ++view) {
    if(id(view) > 0) {
      // view 0 has fixed r and t by definition
      vpmat(view).rot[0] 	= p[i];
      vpmat(view).rot[1] 	= p[i+1];
      vpmat(view).rot[2] 	= p[i+2];
      vpmat(view).t[0] 		= p[i+3];
      vpmat(view).t[1] 		= p[i+4];
      vpmat(view).t[2] 		= p[i+5];
      i += 6;
    }
    /*****
    if(id(view) == 1) {
      // view 1 has fixed scale by definition
      vpmat(view).rot[0] 	= p[i];
      vpmat(view).rot[1] 	= p[i+1];
      vpmat(view).rot[2] 	= p[i+2];
      vpmat(view).t[0] 		= cos(p[i+4])*cos(p[i+3]);
      vpmat(view).t[1] 		= cos(p[i+4])*sin(p[i+3]);
      vpmat(view).t[2] 		= sin(p[1+4]);
      i += 5;
    }
    *****/
    if(has_variable_M(id(view))) {
      // intrinsic matrix is not fixed
      vpmat(view).fx 	= p[i];
      vpmat(view).fy 	= p[i+1];
      vpmat(view).cx 	= p[i+2];
      vpmat(view).cy 	= p[i+3];
      i += 4;
    }
    if(radial) {
      vpmat(view).d = p[i];
      ++i;
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
//
// copies the parameters in this map of view projection matrices to p
//
///////////////////////////////////////////////////////////////////////////////
void RPSolver::vpmatrices_to_params(VecDoub &p) {
  iterator	view;
  int 		i;

  i = 0;
  for(view = begin(); view != end(); ++view) {
    if(id(view) > 0) {
      // view 0 has fixed r and t by definition
      p[i]   = vpmat(view).rot[0];
      p[i+1] = vpmat(view).rot[1];
      p[i+2] = vpmat(view).rot[2];
      p[i+3] = vpmat(view).t[0];
      p[i+4] = vpmat(view).t[1];
      p[i+5] = vpmat(view).t[2];
      i += 6;
    }
    /******
    if(id(view) == 1) {
      // view 1 has fixed scale by definition
      p[i]   = vpmat(view).rot[0];
      p[i+1] = vpmat(view).rot[1];
      p[i+2] = vpmat(view).rot[2];
      p[i+3] = atan(vpmat(view).t[1]/vpmat(view).t[0]);
      p[i+4] = asin(vpmat(view).t[2]);
      i += 5;
    }
    *******/
    if(has_variable_M(id(view))) {
      // intrinsic matrix is not fixed
      p[i]   = vpmat(view).fx;
      p[i+1] = vpmat(view).fy;
      p[i+2] = vpmat(view).cx;
      p[i+3] = vpmat(view).cy;
      i += 4;
    }
    if(radial) {
      p[i] = vpmat(view).d;	// radial distortion
      ++i;
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
// 
// returns the Sampson reprojection error due to correspondence, c.
// R is the (pre-computed) radial fundamental matrix for the two
// views of the correspondence.
//
///////////////////////////////////////////////////////////////////////////////
void RPSolver::closest_points(const GPixel &i, const GPixel &j, coord &p, coord &q) {
  //std::cout << "solving " << i << " " << j << std::endl;
  
  RadialViewProjectionMatrix &	Pi((*this)[i.id]);
  RadialViewProjectionMatrix &	Pj((*this)[j.id]);	  
  coord pixi(3);// pixel i
  coord pixj(3);// pixel j
  coord I(3);	// direction of ray i
  coord J(3);	// direction of ray j
  coord t(3);	// tranlation between camera centres
  coord yhat(3);// basis vector in epipolar plane perp. to I
  double k;	// distance along I of closest point

  pixi[0] = i.x;
  pixi[1] = i.y;
  pixi[2] = 1.0;
  pixj[0] = j.x;
  pixj[1] = j.y;
  pixj[2] = 1.0;


  I = Pi.R1()*(Pi.M1()*pixi);
  J = Pj.R1()*(Pj.M1()*pixj);
  //I = P1[i.id]*pixi;
  //J = P1[j.id]*pixj;
  t = Pj.t - Pi.t;

  transformMatrix tx(cross_prod,t[0], t[1], t[2]);
  transformMatrix Ix(cross_prod,I[0], I[1], I[2]);
  transformMatrix Jx(cross_prod,J[0], J[1], J[2]);

  yhat = Ix*tx*I;
  k = ((t*I).sum() - (t*yhat).sum()*(J*I).sum()/(J*yhat).sum())/(I*I).sum();
  p = k*I + Pi.t;

  yhat = Jx*tx*J;
  k = (-(t*J).sum() + (t*yhat).sum()*(J*I).sum()/(I*yhat).sum())/(J*J).sum();
  q = k*J + Pj.t;

  //std::cout << p << " -> " << q << std::endl;
}


///////////////////////////////////////////////////////////////////////////////
//
// returns true iff the n'th view has a variable intrinsic matrix (i.e.
// the intrinsic matrix is not already known)
//
///////////////////////////////////////////////////////////////////////////////
bool RPSolver::has_variable_M(int n) {
  return(paramMask.find(n) == paramMask.end());
}
