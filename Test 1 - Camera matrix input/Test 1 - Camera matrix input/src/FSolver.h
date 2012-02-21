///////////////////////////////////////////////////////////////////////////////
//
// Solver for calculateing the fundamental matrix from a set of correspondences
//
///////////////////////////////////////////////////////////////////////////////
#ifndef FSOLVER_H
#define FSOLVER_H

#include <vector>
#include <ostream>
#include "transformMatrix.h"
#include "nrlib/nr3.h"

class FSolver {
public:

  bool		load(const char *);		// load correspondences
  int		data_size();
  Doub 		operator ()(VecDoub_I &);   	// calculate squared error
  void		df(VecDoub_I &, VecDoub_O &); 	// calculate rates of change
  void		solve(transformMatrix &);
  static coord	epipole(transformMatrix &);

protected:
  void		calcErrors(VecDoub_I &);	// update 'e' and 'etot'

protected:
  std::vector<coord>	q;	// projector pixel of the correspondences
  std::vector<coord>	qp;	// camera pixel of the correspondences
  std::valarray<double>	e;	// error in n'th correspondence
};


///////////////////////////////////////////////////////////////////////////////
//
// Returns the number of correspondences 
//
///////////////////////////////////////////////////////////////////////////////
inline int FSolver::data_size() {
  return(q.size());
}

#endif
