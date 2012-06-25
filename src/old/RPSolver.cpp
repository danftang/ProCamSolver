#include "RPSolver.h"


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
  iterator 	view;
  int		i;
  coord		t_guess(3);
  GPixel	pixel1, pixel2;

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
    // std::cout << cor[i].i << " " << cor[i].j << std::endl;
  }

  // --- view 0 has reference intrinsics
  paramMask.insert(0);
  (*this)[0].fx = 1.0;
  (*this)[0].fy = 1.0;
  (*this)[0].cx = 0.0;
  (*this)[0].cy = 0.0;
  
  // --- calculate reference point and reference pixels
  EquivalencePartition<GPixel>::member_iterator pix;
  coord pix3(3);
  pix3[2] = 1.0;
  // find equivalence set that can be seen by all cameras
  for(i=0; i < C.size() && C[i].size() < size(); ++i) {}
  if(i == C.size())
    throw("Need a point that can be seen by all cameras to make an RPSolver");
  for(pix = C[i].begin(); pix != C[i].end(); ++pix) {
    pix3[0] = pix->x;
    pix3[1] = pix->y;
    qref[pix->id].resize(3);
    qref[pix->id] = pix3;
    // std::cout << "qref " << pix->id << " = " << qref[pix->id] << std::endl;
  }
  Qref.resize(3);
  Qref = qref[0];

  // --- calculate number of parameters
  i = 0;
  for(view = begin(); view != end(); ++view) {
    paramBase[id(view)] = i;
    if(id(view) > 0) {
      // view 0 has fixed r and t by definition
      i += 4;
    }
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
  Frprmn<RPSolver>	mySolver(*this, 1e-8);
  VecDoub		p(pn);
  int 			n,b;
  iterator		view;
  double		r;
  int			kicks;	 // number of times we've had to kick

  kicks = 0;
  do {
    vpmatrices_to_params(p);
    p = mySolver.minimize(p);
    r = (*this)(p);
    if(r > 1e-6 && ++kicks < 10) {
      for(view = ++begin(); view != end(); ++view) {
	std::cout << view->second << std::endl;
	view->second.t[0] += rand()*1e-6/RAND_MAX - 5e-7;
	view->second.t[1] += rand()*1e-6/RAND_MAX - 5e-7;
	view->second.t[2] += rand()*1e-6/RAND_MAX - 5e-7;
	view->second.rot[0] += rand()*1e-6/RAND_MAX - 5e-7;
	view->second.rot[1] += rand()*1e-6/RAND_MAX - 5e-7;
	view->second.rot[2] += rand()*1e-6/RAND_MAX - 5e-7;
	if(has_variable_M(view->first)) {
	  view->second.fx += rand()*1e-6/RAND_MAX - 5e-7;
	  view->second.fy += rand()*1e-6/RAND_MAX - 5e-7;
	  view->second.cx += rand()*1e-6/RAND_MAX - 5e-7;
	  view->second.cy += rand()*1e-6/RAND_MAX - 5e-7;
	}
      }
    }
    std::cout << "----" << std::endl;
  } while(r > 1e-6 & kicks < 10);
  params_to_vpmatrices(p);

  // --- re-scale the translations
  // --- distance view1->view 2 = 1.0
  // --------------------------------
  double scale = 1.0/sqrt(((*this)[1].t * (*this)[1].t).sum());
  std::cout << "scaling by " << scale << std::endl;
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
double RPSolver::operator()(VecDoub_I &p) {
  int		n,i;
  long double 	e2;
  long double   sigma;
  coord		pi1(3);
  coord		pi2(3);
  std::valarray<long double> 	pi_sum(3);
  std::valarray<long double> 	pi2_sum(3);
  int		N;
  EquivalencePartition<GPixel>::member_iterator pixel1, pixel2;

  params_to_vpmatrices(p);
  e2 = 0.0;
  for(n=0; n<C.size(); ++n) {
    // --- for each equivalence set of global pixels
    pi_sum = 0.0;
    pi2_sum = 0.0;
    N = 0;
    for(pixel2= ++(C[n].begin()); pixel2 != C[n].end(); ++pixel2) {
      for(pixel1 = C[n].begin(); pixel1 != pixel2; ++pixel1) {
	closest_points(*pixel1, *pixel2, pi1, pi2);
	for(i = 0; i<3; ++i) {
	  pi_sum[i] += pi1[i];
	  pi_sum[i] += pi2[i];
	  pi2_sum[i] += pi1[i]*pi1[i];
	  pi2_sum[i] += pi2[i]*pi2[i];
	}
	N += 2;
      }
    }
    pi_sum = pi_sum*pi_sum/N;
    sigma = (pi2_sum - pi_sum).sum()/(pi_sum.sum()/N);
    // std::cout << "Error of point " << n << " is " << sigma << std::endl;
    e2 += sigma;
  }

  return(e2);
}


///////////////////////////////////////////////////////////////////////////////
//
// calculate derivative of error with each parameter using finite difference
// approximation
//
///////////////////////////////////////////////////////////////////////////////
void RPSolver::df(VecDoub_I &p, VecDoub_O &deriv) {
  double e_old;
  double e;
  double x_old;
  VecDoub ph = p;
  double dx;
  const double delta = 1e-7;

  // ---- derivative by finite difference
  // ------------------------------------
  e_old = (*this)(p);
  std::cout.precision(10);
  // std::cout << "original\n" << (*this)[2] << std::endl;
  for(int i=0; i<deriv.size(); ++i) {
    x_old = p[i];
    dx = delta*fabs(x_old);
    if(dx == 0.0) dx = delta;
    ph[i] = x_old + dx;
    dx = ph[i] - x_old;
    e = (*this)(ph);
    ph[i] = x_old;
    deriv[i] = (e - e_old)/dx;
    //    std::cout << deriv[i] << " ";
  }
  //  std::cout << std::endl;
  std::cout << e_old << std::endl;
}


///////////////////////////////////////////////////////////////////////////////
//
// copies the parameters in p to this map of view projection matrices
//
///////////////////////////////////////////////////////////////////////////////
void RPSolver::params_to_vpmatrices(VecDoub_I &p) {
  iterator	view;
  int 		i;
  int		t_index;

  i = 0;
  for(view = begin(); view != end(); ++view) {
    t_index = -1;
    if(id(view) > 0) {
      // view 0 has fixed r and t by definition
      vpmat(view).rot[0] 	= p[i];//*100.0;
      vpmat(view).rot[1] 	= p[i+1];//*100.0;
      vpmat(view).rot[2] 	= p[i+2];//*100.0;
      t_index = i+3;
      i += 4;
    }
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
    P1[view->first] = view->second.R1()*view->second.M1();
    if(t_index > 0) {
      //      std::cout << "Calculating t for " << id(view) << std::endl;
      vpmat(view).t = Qref - (P1[id(view)]*qref[id(view)])*p[t_index];
      // std::cout << "t = " << vpmat(view).t << std::endl;
    }
  }

  /*******
  for(i=0; i<p.size(); ++i) {
    std::cout << p[i] << " ";
  }
  std::cout << std::endl;
  ******/
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
      p[i]   = vpmat(view).rot[0];///100.0;
      p[i+1] = vpmat(view).rot[1];///100.0;
      p[i+2] = vpmat(view).rot[2];///100.0;
      p[i+3] = (vpmat(view).M()*vpmat(view).R()*(Qref - vpmat(view).t))[2];
      i += 4;
    }
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
// returns the point on ray 'i' that is the closest to ray 'j'
//
///////////////////////////////////////////////////////////////////////////////
void RPSolver::closest_points(const GPixel &i, const GPixel &j, coord &p, coord &q) {
  //  std::cout << "solving " << i << " " << j << std::endl;
  
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


  // I = Pi.R1()*(Pi.M1()*pixi);
  // J = Pj.R1()*(Pj.M1()*pixj);
  I = P1[i.id]*pixi;
  J = P1[j.id]*pixj;
  t = Pj.t - Pi.t;

  transformMatrix tx(cross_prod,t[0], t[1], t[2]);
  transformMatrix Ix(cross_prod,I[0], I[1], I[2]);
  transformMatrix Jx(cross_prod,J[0], J[1], J[2]);

  //yhat = Ix*tx*I;
  //k = ((t*I).sum() - (t*yhat).sum()*(J*I).sum()/(J*yhat).sum())/(I*I).sum();
  yhat = Jx*Ix*J;
  k = (t*yhat).sum()/(I*yhat).sum();
  p = k*I + Pi.t;

  //yhat = Jx*tx*J;
  //k = (-(t*J).sum() + (t*yhat).sum()*(J*I).sum()/(I*yhat).sum())/(J*J).sum();
  yhat = Ix*Jx*I;
  k = -(t*yhat).sum()/(J*yhat).sum();
  q = k*J + Pj.t;

  // std::cout << p << " -> " << q << std::endl;
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
