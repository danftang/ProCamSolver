///////////////////////////////////////////////////////////////////////////////
//
// Solver
//
///////////////////////////////////////////////////////////////////////////////
#ifndef SCENESOLVER_H
#define SCENESOLVER_H

#include <vector>
#include <ostream>
#include "transformMatrix.h"
#include "nrlib/conjgrad.h"

class sceneParams {
public:
  sceneParams();

  transformMatrix	N;	// rotation and intrinsics of camera
  transformMatrix	M;	// intrinsics of projector
  coord			t;	// translation of camera from projector
  std::valarray<double>	z;	// projector-z of correspondence points
  std::valarray<double>	zp;	// camera-z of correspondence points
  
};
std::ostream &operator <<(std::ostream &, sceneParams &);


class sceneSolver {
public:

  bool		load(const char *);		// load correspondences
  int		data_size();
  Doub 		operator ()(VecDoub_I &);   	// calculate squared error
  void		df(VecDoub_I &, VecDoub_O &); 	// calculate rates of change
  void		solve();
  
  void		packParams(VecDoub &);
  void		unpackParams(VecDoub &);
  sceneParams &	Params() 			{return(params);}

protected:
  void		calcErrors(VecDoub_I &);	// update 'e' and 'etot'

protected:
  std::vector<coord>	q;	// homogeneous projector pixel correspondences
  std::vector<coord>	qp;	// homogeneous camera pixel correspondences
  std::valarray<coord>	e;	// error vector for each correspondence
  double		z_bar;  // (sum_n z_n) - n
  double		zp_bar;	// (sum_n z'_n) - n

  sceneParams		params;	// parameters of the scene

};


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline int sceneSolver::data_size() {
  return(q.size());
}

#endif
