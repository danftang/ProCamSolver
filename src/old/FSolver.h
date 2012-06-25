///////////////////////////////////////////////////////////////////////////////
//
// This class fits a fundamental matrix to a set of pixel correspondences
//
///////////////////////////////////////////////////////////////////////////////
#ifndef FSOLVER_H
#define FSOLVER_H

#include <vector>
#include <ostream>
#include "numerical.h"
#include "transformMatrix.h"
#include "Correspondence.h"

class FSolver : public transformMatrix {
public:
  FSolver(const char *, double & =residual);	// solve from file
  FSolver(std::vector<Correspondence> &,
	  int =-1, int =0, double & =residual); // solve correspondences

  bool		load(const char *);		// load correspondences
  int		data_size();
  transformMatrix &F() 				{return(*this);}

  Doub 		operator ()(VecDoub_I &);   	// calculate squared error
  void		df(VecDoub_I &, VecDoub_O &); 	// calculate rates of change

protected:
  double	solve();
  void		calcErrors(VecDoub_I &);	// update 'e' and 'etot'

protected:
  std::vector<coord>	q;	// projector pixel of the correspondences
  std::vector<coord>	qp;	// camera pixel of the correspondences
  std::valarray<double>	e;	// error in n'th correspondence
  static double		residual;// residual of the fit
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
