#include "LevMarSolver.h"

using namespace Eigen;

class MySolver : public LevMarSolver<MySolver> {
public:

  int 		operator()(const VectorXd &, VectorXd &);
  //int		df(const VectorXd &, MatrixXd &);
};


int MySolver::operator ()(const VectorXd &x, VectorXd &y) {
  int i;
  for(i = 0; i<x.size(); ++i) {
    y(i) = x(i)*x(i) - i + x(0);
  }
  std::cout << "Err = " << y.norm() << std::endl;
  return(0);
}

/****
int MySolver::df(const VectorXd &x, MatrixXd &jac) {
  int i,j;

  LevMarSolver<MySolver>::df(x,jac);

  std::cout << "jac finite diff = " << std::endl;
  std::cout << jac << std::endl;

  for(i = 0; i<jac.rows(); ++i) {
    for(j = 0; j<jac.cols(); ++j) {
      if(i == j) {
	jac(i,j) = 2*x(i);
      } else {
	jac(i,j) = 0.0;
      }
      if(j == 0) {
	jac(i,j) += 1.0;
      }
    }
  }

  std::cout << "jac analytic = " << std::endl;
  std::cout << jac << std::endl;

  return(0);
}
****/

int main() {
  MySolver 		mysolve;
  VectorXd		params;
  int			info;

  params.resize(10);
  params.fill(1.0);
  info = mysolve.levmar_solve(params, 10);

  std::cout << "info = " << info << std::endl;
  std::cout << "params = " << std::endl
  	    << params << std::endl;

  return(0);
}
