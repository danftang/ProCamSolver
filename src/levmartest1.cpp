#include <Eigen/Core>
#include <unsupported/Eigen/NonLinearOptimization>

using namespace Eigen;

template<int M>
class LMParams : public VectorXd {
public:
  typedef double 	Scalar;
  typedef VectorXd 	InputType;
  typedef VectorXd 	ValueType;
  typedef MatrixXd 	JacobianType;

  enum {
    InputsAtCompileTime = M,
    ValuesAtCompileTime = M
  };

  LMParams() : VectorXd(InputsAtCompileTime), diff(*this) {}

  int inputs() const { return(InputsAtCompileTime); }
  int values() const { return(ValuesAtCompileTime); }

  int operator ()(const InputType &, ValueType &) const;

};


template<int M>
int LMParams<M>::operator ()(const InputType &x, ValueType &y) const {
  int i;
  for(i = 0; i<values(); ++i) {
    y(i) = x(i)*x(i) - i;
  }
  return(0);
}


int main() {
  LMParams<6> 		myParams;
  int 			info;

  
  NumericalDiff<LMParams<6> > numDiff(myParams);
  LevenbergMarquardt<NumericalDiff<LMParams<6> > > levMar(numDiff);
  //LevenbergMarquardt<LMParams<6> > levMar(myParams);

  myParams.fill(1.0);
  info = levMar.minimize(myParams);

  std::cout << "Info = " << info << std::endl;
  std::cout << "params = " << std::endl
  	    << myParams << std::endl << std::endl
	    << myParams.asDiagonal()*myParams << std::endl;	
  return(0);
}
