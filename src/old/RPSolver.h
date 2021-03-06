///////////////////////////////////////////////////////////////////////////////
//
// This class represents a set of cameras/projectors in a single scene. We
// suppose that we have a set of pixel correspondences between some pairs
// of cameras and possibly know some of the intrinsic matrices of the 
// cameras, and wish to know the view-projection matrices of each camera
// and the relative rotation and translation of the cameras with respect
// to some reference.
//
// This class can be constructed with a vector of pixel correspondences.
// Within the correspondences there should be one 3D point that is visible
// by all cameras. 
// If some or all of the cameras have known intrinsic matrices, these can be
// specified in an optional map from viewIDs to intrinsic matrices.
// These will then be held fixed during the fit.
//
// An initial guess at the translation/rotation of the cameras can be supplied
// by constructing with a map from viewIDs to RadialViewProjectionMatrices
// in this case, set the x focal length to zero if you DON'T want the
// camera intrinsics to be fixed.
//
// The camera with viewID of 0 is taken to be the reference view whose
// rotation and translation is held fixed during the fit.
//
// On construction, a fit is performed to find the view-projection matrices
// of all cameras from the set of pixel correspondences.
// 
//
//
// The solution is scaled so that the distance from view 0 to view 1
// is equal to 1.0.
//
///////////////////////////////////////////////////////////////////////////////
#pragma once

#include <valarray>
#include <vector>
#include <map>
#include <set>
#include "numerical.h"
#include "RadialViewProjectionMatrix.h"
#include "Correspondence.h"
#include "EquivalencePartition.h"
#include "GPixel.h"

class RPSolver : public std::map<int,RadialViewProjectionMatrix> {
public:
  RPSolver(std::vector<Correspondence> &,
	   double & =residual);

  RPSolver(std::vector<Correspondence> &, 
	  std::map<int, RadialViewProjectionMatrix> &,
	  double & =residual);

  RPSolver(std::vector<Correspondence> &, 
	  std::map<int, transformMatrix> &,
	  double & =residual);

  double	operator()(VecDoub_I &);
  void		df(VecDoub_I &, VecDoub_O &);

protected:
  void		init_intrinsics(std::map<int, transformMatrix> &);
  void		init_VPM(std::map<int, RadialViewProjectionMatrix> &);
  void		init_views(std::vector<Correspondence> &);
  double	solve();
  void		closest_points(const GPixel &, const GPixel &, coord &, coord &);
  void		params_to_vpmatrices(VecDoub_I &);
  void		vpmatrices_to_params(VecDoub &);
  bool		has_variable_M(int);

  static int	id(iterator it) 			{return(it->first);}
  static RadialViewProjectionMatrix &
                vpmat(iterator it)			{return(it->second);}

protected:
  EquivalencePartition<GPixel> 	C;   	   // the pixel correspondences to fit

  std::map<int,Matrix<double> >	P1;	   // inverse projection matrices
  std::map<int,coord>		qref;	   // pixels of reference point
  coord				Qref;	   // 3D reference point

  std::map<int,int>		paramBase; // indeces into the parameter vector
  std::set<int>			paramMask; // which views have fixed intrinsics
  int				pn;	   // parameter count
  static double			residual;  // residual of fit
  static  bool			radial;	   // include radial distortion
};


