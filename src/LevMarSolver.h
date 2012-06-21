///////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2012 Daniel Tang.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing,
//   software distributed under the License is distributed on an "AS
//   IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
//   express or implied.  See the License for the specific language
//   governing permissions and limitations under the License.
//
///////////////////////////////////////////////////////////////////////////////
#ifndef LEVMARSOLVER_H
#define LEVMARSOLVER_H

#include "stdincludes.h"

namespace DifferenceScheme {
  enum {
    Forward,
    Central
  };
};
///////////////////////////////////////////////////////////////////////////////
/// This is a wrapper for Eigen's Levenberg-Marquardt solver. Inherit
/// from this class to add the levmar_solve() method to your
/// class. The parent class should provide an
///
/// int operator()(const VectorXd &params, VectorXd &err)
///
/// method which calculates the error vector for a given set of parameters.
/// The return value should be zero.
///
/// You may also optionally provide a
///
/// int df(const VectorXd &params, MatrixXd &jac)
///
/// method that calculates the Jacobian for a given set of
/// parameters. If this is not present, then the Jacobian will be
/// calculated by numerical differentiation, but this will result in more
/// calls to operator().
///
/// \params DERIVED the type of the parent class
/// \params DIFFSCHEME the finite difference scheme to use when df() is not
/// present in the inheriting. Can be 'DifferenceScheme::Forward' (default) or
/// 'DifferenceScheme::Central'
///
///////////////////////////////////////////////////////////////////////////////
template<class DERIVED, int DIFFSCHEME = DifferenceScheme::Forward>
class LevMarSolver {
public:

  DERIVED &derived()	 		{return(*static_cast<DERIVED *>(this));}
  const DERIVED &derived() const	{return(*static_cast<const DERIVED *>(this));}

  template<class S> int levmar_solve(Eigen::Matrix<S,Eigen::Dynamic,1> &, int);
  template<class S> int operator ()(const Eigen::Matrix<S,Eigen::Dynamic,1> &,
				          Eigen::Matrix<S,Eigen::Dynamic,1> &);
  template<class S> int df(const Eigen::Matrix<S,Eigen::Dynamic,1> &,
			         Eigen::Matrix<S,Eigen::Dynamic,Eigen::Dynamic> &);
  int 			values() const {return(valsize);}

protected:
  int valsize;
};


template<class DERIVED, int DIFFSCHEME>
template<class S>
int LevMarSolver<DERIVED,DIFFSCHEME>::
operator ()(const Eigen::Matrix<S,Eigen::Dynamic,1> &params,
	    Eigen::Matrix<S,Eigen::Dynamic,1> &vals) {
    return(derived()(params, vals));
}


///////////////////////////////////////////////////////////////////////////////
/// Implementation of a numerical differentiation calculation of the
/// Jacobian. Called by the solver. Overload this in order to provide
/// your own Jacobian calculation.
///
/// \param p the parameter values.
/// \param jac matrix to store the Jacobian in.
///
///////////////////////////////////////////////////////////////////////////////
template<class DERIVED, int DIFFSCHEME>
template<class S>
int LevMarSolver<DERIVED, DIFFSCHEME>::
df(const Eigen::Matrix<S,Eigen::Dynamic,1> &p,
         Eigen::Matrix<S,Eigen::Dynamic,Eigen::Dynamic> &jac) {
  Eigen::Matrix<S,Eigen::Dynamic,1> e(jac.rows());
  Eigen::Matrix<S,Eigen::Dynamic,1> e2(jac.rows());
  Eigen::Matrix<S,Eigen::Dynamic,1> ph(p);
  double x_old;
  double dx;
  const double delta  = std::sqrt(std::numeric_limits<S>::epsilon())*4.0;
  const double dx_min = delta;
  int	i;

  switch(DIFFSCHEME) {
  case DifferenceScheme::Forward:
    (*this)(ph, e2);
    for(i=0; i<p.size(); ++i) {
      x_old = p(i);
      dx = delta*fabs(x_old);
      if(dx < dx_min) dx = dx_min;
      ph(i) = x_old + dx;
      dx = ph(i) - x_old;
      (*this)(ph, e);
      ph(i) = x_old;
      jac.col(i) = (e - e2)/dx;
    }
    return(0);

  case DifferenceScheme::Central:
    for(i=0; i<p.size(); ++i) {
      x_old = p(i);
      dx = delta*fabs(x_old);
      if(dx < dx_min) dx = dx_min;
      ph(i) = x_old - dx;
      (*this)(ph, e2);
      ph(i) = x_old + dx;
      (*this)(ph, e);
      dx = ph(i) - (x_old - dx); // rounding error trick?
      ph(i) = x_old;
      jac.col(i) = (e - e2)/dx;
    }
    return(0);
  otherwise:
    ;
  }
  return(-1);
}


///////////////////////////////////////////////////////////////////////////////
/// Call this from your inheriting class to solve with a given start state
///
/// \param params the state to start the solver in.
/// \param n the number of values in the error vector.
///////////////////////////////////////////////////////////////////////////////
template<class DERIVED, int DIFFSCHEME>
template<class S>
int LevMarSolver<DERIVED,DIFFSCHEME>::
levmar_solve(Eigen::Matrix<S,Eigen::Dynamic,1> &params, int n) {
  int info;

  valsize = n;
//  Eigen::LevenbergMarquardt<LevMarSolver, S> levMar(*this);
  Eigen::LevenbergMarquardt<DERIVED, S> levMar(derived());
  info = levMar.minimize(params);
  return(info);
}


#endif
