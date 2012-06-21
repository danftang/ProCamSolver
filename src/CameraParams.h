///////////////////////////////////////////////////////////////////////////////
//
// This class represents a set of cameras/projectors in a single scene. We
// suppose that we have a set of pixel correspondences between some pairs
// of cameras and possibly know some of the intrinsic matrices of the 
// cameras, and wish to know the view-projection matrices of each camera
// and the relative rotation and translation of the cameras with respect
// to some reference.
//
///////////////////////////////////////////////////////////////////////////////
#ifndef CAMERAPARAMS_H
#define CAMERAPARAMS_H

#include "MotionMatrix.h"
#include "MeasurementMatrix.h"
#include "ImageTransform.h"

template<int M, int F>
class CameraParams {
public:
  typedef double 	       		Scalar;
  typedef Eigen::VectorXd      		InputType;
  typedef Eigen::VectorXd      		ValueType;
  typedef Eigen::MatrixXd      		JacobianType;
  typedef Eigen::Matrix3d      		Matrix3d;
  typedef Eigen::Block<VectorXd,3,1>	Translation;
  typedef Eigen::Block<VectorXd,2,1>	PrincipalPoint;

  enum {
    InputsAtCompileTime = M,
    ValuesAtCompileTime = M
  };


  CameraParams(const MotionMatrix &);

  double			solve_levmar(const MeasurementMatrix &);
  RowsBlock			translation(int);
  PrincipalPoint		principal_pt(int);
  double &			focal_len(int);
  const Eigen::Transform2d	intrinsics(int);
  Eigen::AngleAxisd		rotation(int);

  int 		inputs() const { return(InputsAtCompileTime); }
  int 		values() const { return(ValuesAtCompileTime); }
  int 		operator ()(const InputType &, ValueType &) const;

protected:
  void		params_to_matrices();

protected:
  InputType		variableParams;
  InputType		fixedParams;
  MotionMatrix<M>	P;
  ImageTransform<M>  	Pinv;
};


#endif
