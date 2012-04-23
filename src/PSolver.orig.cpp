#include "PSolver.h"
#include "FSolver.h"


///////////////////////////////////////////////////////////////////////////////
//
// Measure of error introduced by the weighting of the ViewRelations 
// (i.e. the Fundamental matrices) when calculating the error function
//
///////////////////////////////////////////////////////////////////////////////
const double PSolver::k_scale = 16.0;


///////////////////////////////////////////////////////////////////////////////
//
// Default place to put the residual of the fit after solving
//
///////////////////////////////////////////////////////////////////////////////
double PSolver::residual = 0.0;


///////////////////////////////////////////////////////////////////////////////
//
// cor 			- pixel correlations
//
///////////////////////////////////////////////////////////////////////////////
PSolver::PSolver(std::vector<Correspondence> &cor, double &r) {
  init_V(cor);
  r = solve();
}


///////////////////////////////////////////////////////////////////////////////
//
// cor 			- pixel correlations
// fixedIntrinsics 	- prior known intrinsic matrices for certain views
//
///////////////////////////////////////////////////////////////////////////////
PSolver::PSolver(std::vector<Correspondence> &cor,
		 std::map<ViewID, transformMatrix> &fixedIntrinsics,
		 double &r) {

  init_intrinsics(fixedIntrinsics);
  init_V(cor);
  r = solve();
}


///////////////////////////////////////////////////////////////////////////////
//
// VR 			- Fundamental matrices from pairs of views
// fixedIntrinsics	- Information we already have about the
//			  camera/projector intrinsics
//
///////////////////////////////////////////////////////////////////////////////
PSolver::PSolver(std::vector<ViewRelation> &VR,
		 std::map<int, transformMatrix> &fixedIntrinsics,
		 double &r) {
  init_intrinsics(fixedIntrinsics);
  init_V(VR);
  r = solve();
}


///////////////////////////////////////////////////////////////////////////////
//
// Copies data to this for the fixed intrinsics
//
///////////////////////////////////////////////////////////////////////////////
void PSolver::init_intrinsics(std::map<ViewID, transformMatrix> &fixedM) {
  std::map<int, transformMatrix>::iterator it;

  for(it = fixedM.begin(); it != fixedM.end(); ++it) {
    (*this)[it->first].fx = it->second[0][0];
    (*this)[it->first].fy = it->second[1][1];
    (*this)[it->first].cx = it->second[0][2];
    (*this)[it->first].cy = it->second[1][2];
  }
}


///////////////////////////////////////////////////////////////////////////////
//
// initialise the vector of ViewRelations and set up paramBase and
// paramMask
//
///////////////////////////////////////////////////////////////////////////////
void PSolver::init_V(std::vector<ViewRelation> &VR) {

  V.resize(VR.size());
  for(int n=0; n<VR.size(); ++n) {
    V[n] = VR[n];
    V[n].w = 1.0/k_scale;
  }

  // --- create entries in this for all views
  // ----------------------------------------
  for(int n = 0; n<V.size(); ++n) {
    (*this)[V[n].i];
    (*this)[V[n].j];
  }

  // --- set up indeces into the parameter vector
  // --- for each view
  // --------------------------------------------
  iterator	view;
  pn = 0;
  for(view = begin(); view != end(); ++view) {
    paramBase[id(view)] = pn;
    if(id(view) == 0) {
      // reference view has fixed rotation/translation
      if(vpmat(view).fx == 0.0) {
	vpmat(view).fx = 1.0;
	paramMask[id(view)] = true;
	pn += 4;
      } else {
	paramMask[id(view)] = false;
      }
    } else if(vpmat(view).fx == 0.0) {
      vpmat(view).fx = 1.0;
      paramMask[id(view)] = true;
      pn += 10;
    } else {
      paramMask[id(view)] = false;
      pn += 6;
    }
  }

  // --- add parameters for scale of each fundamental matrix
  // --- apart from the first.
  // -------------------------------------------------------
  pn += V.size() - 1;
}


///////////////////////////////////////////////////////////////////////////////
//
// create the fundamental matrices for the set of correspondences
//
///////////////////////////////////////////////////////////////////////////////
void PSolver::init_V(std::vector<Correspondence> &cor) {
  ViewRelation 					view;
  std::map<ViewID, ViewProjectionMatrix>::iterator vi, vj;	
  std::vector<ViewRelation>			VR;

  // --- create entries in this for all views
  // ----------------------------------------
  for(int n = 0; n<cor.size(); ++n) {
    (*this)[cor[n].i];
    (*this)[cor[n].j];
  }

  for(vj = begin(); vj != end(); ++vj) {
    for(vi = begin(); vi != vj; ++vi) {
      view.i = vi->first;
      view.j = vj->first;
      view.F = FSolver(cor, view.i, view.j);
      VR.push_back(view);
    }
  }

  init_V(VR);
}


///////////////////////////////////////////////////////////////////////////////
//
// Find the camera matrices which minimise the sum of squared errors
//
///////////////////////////////////////////////////////////////////////////////
double PSolver::solve() {
  Frprmn<PSolver>	mySolver(*this, 1e-7);
  VecDoub		p(pn);
  int 			n,b;
  iterator		view;

  // --- setup initial conditions
  // ----------------------------
  for(view = begin(); view != end(); ++view) {
    vpmat(view).t[0] = 0.011*id(view);
    vpmat(view).t[1] = 0.012*id(view);
    vpmat(view).t[2] = 0.013*id(view);
  }

  vpmatrices_to_params(p);
  p = mySolver.minimize(p);
  params_to_vpmatrices(p);
  return((*this)(p));
}


///////////////////////////////////////////////////////////////////////////////
//
// returns the error for a given set of parameters
//
///////////////////////////////////////////////////////////////////////////////
double PSolver::operator()(VecDoub_I &p) {
  int	n;
  double e2;

  params_to_vpmatrices(p);
  e2 = 0.0;
  for(n=0; n<V.size(); ++n) {
    e2 += view_error(V[n]);
  }
  return(e2);
}


///////////////////////////////////////////////////////////////////////////////
//
// calculate derivative of error with each parameter
//
///////////////////////////////////////////////////////////////////////////////
void PSolver::df(VecDoub_I &p, VecDoub_O &deriv) {
  int i;
  
  for(i = 0; i<deriv.size(); ++i) {
    deriv[i] = 0.0;
  }
  
  params_to_vpmatrices(p);
  for(i=0; i<V.size(); ++i) {
    view_deriv(V[i], deriv, deriv.size()-i);
  }
}


///////////////////////////////////////////////////////////////////////////////
//
// copies the parameters in p to this map of view projection matrices
//
///////////////////////////////////////////////////////////////////////////////
void PSolver::params_to_vpmatrices(VecDoub_I &p) {
  iterator	view;
  int 		i;

  for(view = begin(); view != end(); ++view) {
    i = paramBase[id(view)];
    if(id(view) != 0) {
      // view 0 has fixed r and t by definition
      vpmat(view).rot[0] 	= p[i];
      vpmat(view).rot[1] 	= p[i+1];
      vpmat(view).rot[2] 	= p[i+2];
      vpmat(view).t[0] 		= p[i+3];
      vpmat(view).t[1] 		= p[i+4];
      vpmat(view).t[2] 		= p[i+5];
      i += 6;
    }
    if(paramMask[id(view)]) {
      // intrinsic matrix is not fixed
      vpmat(view).fx 	= p[i];
      vpmat(view).fy 	= p[i+1];
      vpmat(view).cx 	= p[i+2];
      vpmat(view).cy 	= p[i+3];
    }
  }
  // --- copy weights of fundamental matrices
  for(i = 1; i<V.size(); ++i) {
    V[i].w = p[p.size()-i];
  }
}


///////////////////////////////////////////////////////////////////////////////
//
// copies the parameters in this map of view projection matrices to p
//
///////////////////////////////////////////////////////////////////////////////
void PSolver::vpmatrices_to_params(VecDoub &p) {
  iterator	view;
  int 		i;

  for(view = begin(); view != end(); ++view) {
    i = paramBase[id(view)];
    if(id(view) != 0) {
      // view 0 has fixed r and t by definition
      p[i]   = vpmat(view).rot[0];
      p[i+1] = vpmat(view).rot[1];
      p[i+2] = vpmat(view).rot[2];
      p[i+3] = vpmat(view).t[0];
      p[i+4] = vpmat(view).t[1];
      p[i+5] = vpmat(view).t[2];
      i += 6;
    }
    if(paramMask[id(view)]) {
      // intrinsic matrix is not fixed
      p[i]   = vpmat(view).fx;
      p[i+1] = vpmat(view).fy;
      p[i+2] = vpmat(view).cx;
      p[i+3] = vpmat(view).cy;
    }
  }  
  // --- copy weights of fundamental matrices
  for(i = 1; i<V.size(); ++i) {
    p[p.size()-i] = V[i].w;
  }
}


///////////////////////////////////////////////////////////////////////////////
// 
// returns the sum of squared errors due to view relation, vr.
//
///////////////////////////////////////////////////////////////////////////////
double PSolver::view_error(WeightedViewRelation &vr) {
  ViewProjectionMatrix & pi((*this)[vr.i]);
  ViewProjectionMatrix & pj((*this)[vr.j]);
  Matrix<double>	E;
  double 		e2;
  
  E = 
    pj.M1T() * pj.R() * (pj.T() - pi.T()) * pi.R1() * pi.M1() - 
    (vr.F * (k_scale * vr.w));
  e2 = elemental_dot_prod(E,E);
  return(e2);
}


///////////////////////////////////////////////////////////////////////////////
//
// Add contribution of view relation vr to the derivative of the error, deriv.
// n identifies the element of deriv that corresponds to the derivative wrt
// the weight of the fundamental matrix of this view.
//
///////////////////////////////////////////////////////////////////////////////
void PSolver::view_deriv(WeightedViewRelation &vr, VecDoub_O &deriv, int n) {
  ViewProjectionMatrix & pi((*this)[vr.i]);
  ViewProjectionMatrix & pj((*this)[vr.j]);
  transformMatrix 	E;
  transformMatrix	a;
  transformMatrix	b;
  int			ibase = paramBase[vr.i];
  int			jbase = paramBase[vr.j];

  E = 
    pj.M1T() * pj.R() * (pj.T() - pi.T()) * pi.R1() * pi.M1() - 
    (vr.F * (k_scale * vr.w));

  if(vr.i != 0) {
    // --- i rotation
    // --------------
    a = pj.M1T() * pj.R() * (pj.T() - pi.T());
    b = pi.M1();
    deriv[ibase]     += 2.0*elemental_dot_prod(E, a * pi.dR1_drx() * b);
    deriv[ibase + 1] += 2.0*elemental_dot_prod(E, a * pi.dR1_dry() * b);
    deriv[ibase + 2] += 2.0*elemental_dot_prod(E, a * pi.dR1_drz() * b);

    // --- i translation
    // -----------------
    a = pj.M1T() * pj.R();
    b = pi.R1() * pi.M1();
    deriv[ibase + 3] += -2.0*elemental_dot_prod(E, a * pi.dT_dtx() * b);
    deriv[ibase + 4] += -2.0*elemental_dot_prod(E, a * pi.dT_dty() * b);
    deriv[ibase + 5] += -2.0*elemental_dot_prod(E, a * pi.dT_dtz() * b);

    ibase += 6;
  } else {
    a = pj.M1T() * pj.R(); // setup for j translation
    b = pi.R1() * pi.M1();
  }
  if(vr.j != 0) {
    // --- j translation
    // -----------------
    // a & b already set
    deriv[jbase + 3] += 2.0*elemental_dot_prod(E, a * pj.dT_dtx() * b);
    deriv[jbase + 4] += 2.0*elemental_dot_prod(E, a * pj.dT_dty() * b);
    deriv[jbase + 5] += 2.0*elemental_dot_prod(E, a * pj.dT_dtz() * b);

    // --- j rotation
    // --------------
    a = pj.M1T();
    b = (pj.T() - pi.T()) * pi.R1() * pi.M1();
    deriv[jbase]     += 2.0*elemental_dot_prod(E, a * pj.dR_drx() * b);
    deriv[jbase + 1] += 2.0*elemental_dot_prod(E, a * pj.dR_dry() * b);
    deriv[jbase + 2] += 2.0*elemental_dot_prod(E, a * pj.dR_drz() * b);
    
    jbase += 6;
  }

  // --- i intrinsics
  // ----------------
  if(paramMask[vr.i]) {
    a = pj.M1T() * pj.R() * (pj.T() - pi.T()) * pi.R1();
    deriv[ibase]     += 2.0*elemental_dot_prod(E, a * pi.dM1_dfx());
    deriv[ibase + 1] += 2.0*elemental_dot_prod(E, a * pi.dM1_dfy());
    deriv[ibase + 2] += 2.0*elemental_dot_prod(E, a * pi.dM1_dcx());
    deriv[ibase + 3] += 2.0*elemental_dot_prod(E, a * pi.dM1_dcy());
  }

  // --- j intrinsics
  // ----------------
  if(paramMask[vr.j]) {
    b = pj.R() * (pj.T() - pi.T()) * pi.R1() * pi.M1();
    deriv[jbase]     += 2.0*elemental_dot_prod(E, pj.dM1_dfx().transpose() *b);
    deriv[jbase + 1] += 2.0*elemental_dot_prod(E, pj.dM1_dfy().transpose() *b);
    deriv[jbase + 2] += 2.0*elemental_dot_prod(E, pj.dM1_dcx().transpose() *b);
    deriv[jbase + 3] += 2.0*elemental_dot_prod(E, pj.dM1_dcy().transpose() *b);
  }

  // --- weight of fundamental matrix
  if(n < deriv.size()) deriv[n] = -2.0*k_scale*elemental_dot_prod(E,vr.F);

  /****
  std::cout << "Derivative after F" << vr.i << vr.j << std::endl;
  for(int i=0; i<deriv.size(); ++i) {
    std::cout << deriv[i] << " ";
  }
  std::cout << std::endl;
  ****/
}


///////////////////////////////////////////////////////////////////////////////
//
// returns sum_i sum_j M_{ij}N_{ij}
//
///////////////////////////////////////////////////////////////////////////////
double PSolver::elemental_dot_prod(const Matrix<double> &M, 
				   const Matrix<double> &N) {
  int i,j;
  double result;

  result = 0.0;
   for(i = 0; i<M.isize(); ++i) {
     for(j = 0; j<M.jsize(); ++j) {
       result += M[i][j]*N[i][j];
     }
   }
   return(result);
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
PSolver::WeightedViewRelation &PSolver::WeightedViewRelation::
operator =(const ViewRelation &vr) {
  ViewRelation::operator =(vr);
  w = 1.0;
  return(*this);
}
