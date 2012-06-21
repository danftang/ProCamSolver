///////////////////////////////////////////////////////////////////////////////
#ifndef VIEWPROJECTIONMATRIX_H
#define VIEWPROJECTIONMATRIX_H

#include "stdincludes.h"
#include "TransformMatrix.h"

class ViewProjectionMatrix {
public:
  ViewProjectionMatrix();

  
  Eigen::Matrix3d	R() const;	// rotation matrix
  Eigen::Matrix3d	R1() const;	// inverse rotation matrix
  Eigen::Matrix3d	M() const;	// projection matrix
  Eigen::Matrix3d	M1() const;	// inverse projection matrix
  Eigen::Matrix3d	M1T() const;	// transpose inverse projection matrix
  Eigen::Matrix3d	T() const;	// translation cross-product
  Eigen::Matrix<double,3,4> P() const; // Camera matrix

  /*****
  Eigen::Matrix3d	dR_drx(); // derivatives
  Eigen::Matrix3d	dR_dry();
  Eigen::Matrix3d	dR_drz();
  Eigen::Matrix3d	dR1_drx();
  Eigen::Matrix3d	dR1_dry();
  Eigen::Matrix3d	dR1_drz();
  Eigen::Matrix3d	dM1_dfx();
  Eigen::Matrix3d	dM1_dfy();
  Eigen::Matrix3d	dM1_dcx();
  Eigen::Matrix3d	dM1_dcy();
  Eigen::Matrix3d	dT_dtx();
  Eigen::Matrix3d	dT_dty();
  Eigen::Matrix3d	dT_dtz();
  ***/

  coord		rot;	// rotation vector
  coord		t;	// translation vector
  double	cx;	// origin of sensor
  double	cy;
  double	fx;	// x-focal length
  double	fy;	// y-focal length

};




///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline ViewProjectionMatrix::ViewProjectionMatrix() : 
  rot(3), 
  t(3), 
  fx(0.0),
  fy(1.0),
  cx(0.0),
  cy(0.0)
{
  rot(0) = 0.0;
  rot(1) = 0.0;
  rot(2) = 0.0;
  t(0) = 0.0;
  t(1) = 0.0;
  t(2) = 0.0;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//inline ViewProjectionMatrix &
//ViewProjectionMatrix::operator =(const Eigen::Matrix3d &P) {
//  Eigen::Matrix3d::operator=(P);
//  return(*this);
//}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline Eigen::Matrix3d ViewProjectionMatrix::R() const {
  return(TransformMatrix(xrotate,rot(0)) *
	 TransformMatrix(yrotate,rot(1)) *
	 TransformMatrix(zrotate,rot(2)));

}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline Eigen::Matrix3d ViewProjectionMatrix::R1() const {
  return(TransformMatrix(zrotate,-rot(2)) *
	 TransformMatrix(yrotate,-rot(1)) *
	 TransformMatrix(xrotate,-rot(0)));
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline Eigen::Matrix3d ViewProjectionMatrix::M() const {
  return(TransformMatrix(intrinsic,fx,fy,cx,cy));
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline Eigen::Matrix3d ViewProjectionMatrix::M1() const {
  return(TransformMatrix(inverse_intrinsic,fx,fy,cx,cy));
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline Eigen::Matrix3d ViewProjectionMatrix::M1T() const {
  return(TransformMatrix(transpose_inverse_intrinsic,fx,fy,cx,cy));
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline Eigen::Matrix3d ViewProjectionMatrix::T() const {
  return(TransformMatrix(cross_prod,t(0),t(1),t(2)));
}


/*****************
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline Eigen::Matrix3d ViewProjectionMatrix::dR_drx() {
  const double PI = 3.1415926535897932;
  return(TransformMatrix(xrotate,rot(0) + PI/2.0) *
	 TransformMatrix(yrotate,rot(1)) *
	 TransformMatrix(zrotate,rot(2)));
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline Eigen::Matrix3d ViewProjectionMatrix::dR_dry() {
  const double PI = 3.1415926535897932;
  return(TransformMatrix(xrotate,rot(0)) *
	 TransformMatrix(yrotate,rot(1) + PI/2.0) *
	 TransformMatrix(zrotate,rot(2)));
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline Eigen::Matrix3d ViewProjectionMatrix::dR_drz() {
  const double PI = 3.1415926535897932;
  return(TransformMatrix(xrotate,rot(0)) *
	 TransformMatrix(yrotate,rot(1)) *
	 TransformMatrix(zrotate,rot(2) + PI/2.0));
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline Eigen::Matrix3d ViewProjectionMatrix::dR1_drx() {
  const double PI = 3.1415926535897932;
  return(TransformMatrix(zrotate,-rot(2)) *
	 TransformMatrix(yrotate,-rot(1)) *
	 TransformMatrix(xrotate,-rot(0) - PI/2.0));
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline Eigen::Matrix3d ViewProjectionMatrix::dR1_dry() {
  const double PI = 3.1415926535897932;
  return(TransformMatrix(zrotate,-rot(2)) *
	 TransformMatrix(yrotate,-rot(1) - PI/2.0) *
	 TransformMatrix(xrotate,-rot(0)));
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline Eigen::Matrix3d ViewProjectionMatrix::dR1_drz() {
  const double PI = 3.1415926535897932;
  return(TransformMatrix(zrotate,-rot(2) - PI/2.0) *
	 TransformMatrix(yrotate,-rot(1)) *
	 TransformMatrix(xrotate,-rot(0)));
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline Eigen::Matrix3d ViewProjectionMatrix::dM1_dfx() {
  Eigen::Matrix3d result;
  result.fill(0.0);
  result(0,0) = -1.0/(fx*fx);
  return(result);
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline Eigen::Matrix3d ViewProjectionMatrix::dM1_dfy() {
  Eigen::Matrix3d result;
  result.fill(0.0);
  result(1,1) = -1.0/(fy*fy);
  return(result);
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline Eigen::Matrix3d ViewProjectionMatrix::dM1_dcx() {
  Eigen::Matrix3d result;
  result.fill(0.0);
  result(0,2) = -1.0/fx;
  return(result);
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline Eigen::Matrix3d ViewProjectionMatrix::dM1_dcy() {
  Eigen::Matrix3d result;
  result.fill(0.0);
  result(1,2) = -1.0/fy;
  return(result);
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline Eigen::Matrix3d ViewProjectionMatrix::dT_dtx() {
  Eigen::Matrix3d result;
  result.fill(0.0);
  result(1,2) = -1.0;
  result(2,1) = 1.0;
  return(result);
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline Eigen::Matrix3d ViewProjectionMatrix::dT_dty() {
  Eigen::Matrix3d result;
  result.fill(0.0);
  result(0,2) = 1.0;
  result(2,0) = -1.0;
  return(result);
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline Eigen::Matrix3d ViewProjectionMatrix::dT_dtz() {
  Eigen::Matrix3d result;
  result.fill(0.0);
  result(0,1) = -1.0;
  result(1,0) = 1.0;
  return(result);
}
***************/

///////////////////////////////////////////////////////////////////////////////
/// Returns the camera matrix of this
///////////////////////////////////////////////////////////////////////////////
inline Eigen::Matrix<double,3,4> ViewProjectionMatrix::P() const {
  Eigen::Matrix<double,3,4>	Projection;
  Projection.topLeftCorner(3,3) = Eigen::Matrix3d::Identity();
  Projection.col(3) = -t;
  Projection = M() * R() * Projection;
  return(Projection);
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline std::ostream &operator <<(std::ostream &out, ViewProjectionMatrix &v) {
  out << "rotation = " 
      << v.rot(0) << ", " << v.rot(1) << ", " << v.rot(2) << std::endl;
  out << "translation = "
      << v.t(0) << ", " << v.t(1) << ", " << v.t(2) << std::endl;
  out << "intrinsics = " 
      << v.fx << ", " << v.fy << " | " << v.cx << ", " << v.cy << std::endl;
  return(out);
}


#endif
