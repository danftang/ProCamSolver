#include <Eigen/Core>
#include <unsupported/Eigen/NonLinearOptimization>

using namespace Eigen;

template<class DERIVED>
class LevMarFunctor {
public:
  typedef double		 		Scalar;
  typedef Matrix<Scalar,Eigen::Dynamic,1> 	InputType;
  typedef Matrix<Scalar,Eigen::Dynamic,1>	ValueType;
  typedef MatrixXd 				JacobianType;

  enum {
    InputsAtCompileTime = 6,//DERIVED::InputsAtCompileTime,
    ValuesAtCompileTime = 6 //DERIVED::ValuesAtCompileTime
  };

  DERIVED &derived()    	{return(*static_cast<DERIVED *>(this));}
  const DERIVED &derived() const{return(*static_cast<const DERIVED *>(this));}
  int inputs() const 		{return(InputsAtCompileTime);}
  int values() const 		{return(ValuesAtCompileTime);}

  int operator ()(const InputType &params, ValueType &vals) const {
    return(derived().operator()(params,vals));
  }

};


template<class DERIVED>
class LevMarSolver : public NumericalDiff<LevMarFunctor<DERIVED> > {
public:
  typedef typename LevMarFunctor<DERIVED>::Scalar	Scalar;
  
  Scalar solve(Matrix<Scalar,Eigen::Dynamic,1> &);

};


template<class DERIVED>
typename LevMarSolver<DERIVED>::Scalar LevMarSolver<DERIVED>::solve(Matrix<Scalar,Eigen::Dynamic,1> &params) {
    LevenbergMarquardt<LevMarSolver<DERIVED> > levMar(*this);
    std::cout << "Info = " << levMar.minimize(params) << std::endl;
    return(0.0);
}


class MySolver : public LevMarSolver<MySolver> {
public:
  typedef double	Scalar;

  int 		operator()(const InputType &, ValueType &) const;
  Scalar	solver();

  InputType 	params;
};

double MySolver::solver() {
  return(solve(params));
}

int MySolver::operator ()(const InputType &x, ValueType &y) const {
  int i;
  for(i = 0; i<values(); ++i) {
    y(i) = x(i)*x(i) - i;
  }
  return(0);
}


int main() {
  MySolver 		mysolve;

  mysolve.params.resize(6);
  mysolve.params.setOnes();
  mysolve.solver();

  std::cout << "params = " << std::endl
  	    << mysolve.params << std::endl;

  return(0);
}
