///////////////////////////////////////////////////////////////////////////////
//
// Class to represent a view projection matrix, with added radial distortion
// coefficient, d.
//
// For a description of the distortion matrix see Yamazaki et. al (2011)
// "Simultaneous self-calibration of a projector and a camera using structured
// light", Computer vision and pattern recognition workshops, 2011, 60-67 
// DOI:10.1109/CVPRW.2011.5981781
//
///////////////////////////////////////////////////////////////////////////////

#ifndef RADIALVIEWPROJECTIONMATRIX_H
#define RADIALVIEWPROJECTIONMATRIX_H

#include "ViewProjectionMatrix.h"

class RadialViewProjectionMatrix : public ViewProjectionMatrix {
public:
  RadialViewProjectionMatrix();

  Matrix<double>	D();		// radial distortion matrix
  Matrix<double>	DT();		// transpose of D()

  Matrix<double>	dD_dd();	// change of D() with distortion coeff
  Matrix<double>	dD_dcx();	// change of D() with x-centre
  Matrix<double>	dD_dcy();	// change of D() with y-centre

  double	d;	// radial distortion coefficient
};


//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
inline RadialViewProjectionMatrix::RadialViewProjectionMatrix() :
  ViewProjectionMatrix(),
  d(0.0)
{  
}


//////////////////////////////////////////////////////////////////////////////
//
// The radial distortion operator
// 
//////////////////////////////////////////////////////////////////////////////
inline Matrix<double> RadialViewProjectionMatrix::D() {
  return(transformMatrix(radial_distortion,cx,cy,d));
}

//////////////////////////////////////////////////////////////////////////////
inline Matrix<double> RadialViewProjectionMatrix::DT() {
  return(transformMatrix(transpose_radial_distortion,cx,cy,d));
}


//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
inline Matrix<double> RadialViewProjectionMatrix::dD_dd() {
  Matrix<double>	r(3,4);

  r[0][0] = cx; 	r[0][1] = -2.0*cx*cx;
  r[0][2] = -2.0*cx*cy; r[0][3] = cx*(cx*cx + cy*cy);
  r[1][0] = cy; 	r[1][1] = -2.0*cx;
  r[1][2] = -2.0*cy*cy; r[1][3] = cy*(cx*cx + cy*cy);
  r[2][0] = 1.0; 	r[2][1] = -2.0*cx;
  r[2][2] = -2.0*cy; 	r[2][3] = cx*cx + cy*cy;
  return(r);
}

//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
inline Matrix<double> RadialViewProjectionMatrix::dD_dcx() {
  Matrix<double>	r(3,4);

  r[0][0] = d; 		r[0][1] = -4.0*d*cx;
  r[0][2] = -2.0*d*cy; 	r[0][3] = d*(3.0*cx*cx + cy*cy);
  r[1][0] = 0.0; 	r[1][1] = -2.0*d;
  r[1][2] = 0.0; 	r[1][3] = 2.0*d*cy*cx;;
  r[2][0] = 0.0; 	r[2][1] = -2.0*d;
  r[2][2] = 0.0; 	r[2][3] = 2.0*cx*d;
  return(r);
}

//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
inline Matrix<double> RadialViewProjectionMatrix::dD_dcy() {
  Matrix<double>	r(3,4);

  r[0][0] = 0.0;	r[0][1] = 0.0;
  r[0][2] = -2.0*d*cx; 	r[0][3] = 2.0*d*cx*cy;
  r[1][0] = d;	 	r[1][1] = 0.0;
  r[1][2] = 04.0*d*cy; 	r[1][3] = d*(cx*cx + 3.0*cy*cy);
  r[2][0] = 0.0; 	r[2][1] = 0.0;
  r[2][2] = -2.0*d; 	r[2][3] = 2.0*cy*d;
  return(r);
}


//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
inline std::ostream &operator <<(std::ostream &out, RadialViewProjectionMatrix &vp) {
  out << (ViewProjectionMatrix &)vp;
  out << "radial distortion = " << vp.d << std::endl;
  return(out);
}


#endif
