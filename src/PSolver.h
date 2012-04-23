///////////////////////////////////////////////////////////////////////////////
//
// This class represents a set of cameras/projectors in a single scene. We
// suppose that we know the fundamental matrices between some pairs of cameras
// and possibly know some of the intrinsic matrices of the cameras, and wish
// to know the view-projection matrices of each camera and the relative
// rotation and translation of the cameras with respect to some reference.
//
// This can be constructed with a vector of pixel correspondences or a vector
// of ViewRelations that specify the fundamental matrices between cameras.
// If some or all of the cameras have known intrinsic matrices, these can be
// specified in an optional map from viewIDs to intrinsic matrices.
// These will then be held fixed during the fit.
//
// The camera with viewID of 0 is taken to be the reference view whose
// rotation and translation is held fixed during the fit.
//
// On construction, a fit is performed to find the view-projection matrices
// of all cameras from the set of fundamental matrices.
// We assume that each view-projection matrix is of the form P=MR[I|-t]
// where
//     (fx 0  cx )
// M = (0  fy cy )
//     (0  0  1  )
// R = rotation matrix
// I = identity
// t = translation of centre of projection
//
// The fit finds the minimum of the sum of the squares of the errors in the
// fundamental matrices, compared to the supplied fundamental matrices (if
// constructed with a vector of ViewRelations) or compared to the fundamental
// matrices created by using FSolver to solve for the pixel correspondences
// (if constructed with a vector of Correspondences).
//
///////////////////////////////////////////////////////////////////////////////
#ifndef PSOLVER_H
#define PSOLVER_H

#include <valarray>
#include <vector>
#include <map>
#include "numerical.h"
#include "ViewRelation.h"
#include "ViewProjectionMatrix.h"
#include "Correspondence.h"

class PSolver : public std::map<int,ViewProjectionMatrix> {
public:
  typedef int ViewID;

  class WeightedViewRelation : public ViewRelation {
  public:
    double w;
    WeightedViewRelation &operator =(const ViewRelation &);
  };

public:
  PSolver(std::vector<Correspondence> &, double & =residual);
  PSolver(std::vector<Correspondence> &, 
	  std::map<ViewID, transformMatrix> &,
	  double & =residual);
  PSolver(std::vector<ViewRelation> &,
	  std::map<ViewID, transformMatrix> &,
	  double & =residual);

  double	operator()(VecDoub_I &);
  void		df(VecDoub_I &, VecDoub_O &);

protected:
  void		init_V(std::vector<ViewRelation> &);
  void		init_V(std::vector<Correspondence> &);
  void		init_intrinsics(std::map<ViewID, transformMatrix> &);
  double	solve();
  double	view_error(WeightedViewRelation &);
  void		view_deriv(WeightedViewRelation &, VecDoub_O &, int);
  void		params_to_vpmatrices(VecDoub_I &);
  void		vpmatrices_to_params(VecDoub &);
  double	elemental_dot_prod(const Matrix<double> &, 
				   const Matrix<double> &);
  static int	id(iterator it) 			{return(it->first);}
  static ViewProjectionMatrix &	vpmat(iterator it)	{return(it->second);}

protected:
  std::vector<WeightedViewRelation>	V;
  std::map<int,int>		paramBase; // indeces into the parameter vector
  std::map<int,bool>		paramMask; // which views have fixed intrinsics
  int				pn;	   // parameter count
  static double			residual;  // residual of fit
  static const double		k_scale;   // Sensitivity of error to scaling
};


#endif
