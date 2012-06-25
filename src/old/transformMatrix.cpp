#include "transformMatrix.h"

///////////////////////////////////////////////////////////////////////////////
//
// Initialise matrix to various useful matrices:
//
// Intrinsic matrix 	- transformMatrix(intrinsic,fx,fy,cx,cy)
// Inverse intrinsic	- transformMatrix(inverse_intrinsic,fx,fy,cx,cy)
// Rotate about x-axis	- transformMatrix(xrotate,angle)
// Rotate about y-axis	- transformMatrix(yrotate,angle)
// Rotate about z-axis	- transformMatrix(zrotate,angle)
// Cross product 	- transformMatrix(cross_prod,x,y,z)
// Radial distortion    - transformMatrix(radial_distortion, cx, cy, d)
//
///////////////////////////////////////////////////////////////////////////////
transformMatrix::transformMatrix(matrixType t, 
				double a, double b, 
				double c, double d) : Matrix<double>(3,3) {
  int i,j;

  (*this)[0][1] = 0.0;
  (*this)[0][2] = 0.0;
  (*this)[1][0] = 0.0;
  (*this)[1][2] = 0.0;
  (*this)[2][0] = 0.0;
  (*this)[2][1] = 0.0;
  
  switch(t) {
  case intrinsic:
    (*this)[0][0] = a;
    (*this)[1][1] = b;
    (*this)[0][2] = c;
    (*this)[1][2] = d;
    (*this)[2][2] = 1.0;
    break;
  case inverse_intrinsic:
    (*this)[0][0] = 1.0/a;
    (*this)[1][1] = 1.0/b;
    (*this)[0][2] = -c/a;
    (*this)[1][2] = -d/b;
    (*this)[2][2] = 1.0;
    break;
  case transpose_inverse_intrinsic:
    (*this)[0][0] = 1.0/a;
    (*this)[1][1] = 1.0/b;
    (*this)[2][0] = -c/a;
    (*this)[2][1] = -d/b;
    (*this)[2][2] = 1.0;
    break;
  case xrotate:
    (*this)[2][2] = (*this)[1][1] = cos(a);
    (*this)[2][1] = sin(a);
    (*this)[1][2] = -(*this)[2][1];
    (*this)[0][0] = 1.0;
    break;
  case yrotate:
    (*this)[2][2] = (*this)[0][0] = cos(a);
    (*this)[0][2] = sin(a);
    (*this)[2][0] = -(*this)[0][2];
    (*this)[1][1] = 1.0;
  break;
  case zrotate:
    (*this)[1][1] = (*this)[0][0] = cos(a);
    (*this)[1][0] = sin(a);
    (*this)[0][1] = -(*this)[1][0];
    (*this)[2][2] = 1.0;
    break;
  case cross_prod:
    (*this)[0][0] = 0.0; (*this)[0][1] = -c;  (*this)[0][2] = b;
    (*this)[1][0] = c;   (*this)[1][1] = 0.0; (*this)[1][2] = -a;
    (*this)[2][0] = -b;  (*this)[2][1] = a;   (*this)[2][2] = 0.0;
    break;
  case zero:
    (*this)[0][0] = 0.0;
    (*this)[1][1] = 0.0;
    (*this)[2][2] = 0.0;
    break;
  case camera:
    resize(3,4);
    break;
  case radial_distortion:
    resize(3,4);
    (*this)[0][0] = c*a; 	   (*this)[0][1] = 1.0 - 2.0*c*a*a; 
    (*this)[0][2] = -2.0*c*a*b;    (*this)[0][3] = c*a*(a*a+b*b);
    (*this)[1][0] = c*b; 	   (*this)[1][1] = -2.0*c*a; 
    (*this)[1][2] = 1.0-2.0*c*b*b; (*this)[1][3] = c*b*(a*a+b*b);
    (*this)[2][0] = c;	 	   (*this)[2][1] = -2.0*c*a; 
    (*this)[2][2] = -2.0*c*b;      (*this)[2][3] = 1.0 + c*(a*a+b*b);
    break;
  case transpose_radial_distortion:
    resize(4,3);
    (*this)[0][0] = c*a; 	   (*this)[1][0] = 1.0 - 2.0*c*a*a; 
    (*this)[2][0] = -2.0*c*a*b;    (*this)[3][0] = c*a*(a*a+b*b);
    (*this)[0][1] = c*b; 	   (*this)[1][1] = -2.0*c*a; 
    (*this)[2][1] = 1.0-2.0*c*b*b; (*this)[3][1] = c*b*(a*a+b*b);
    (*this)[0][2] = c;	 	   (*this)[1][2] = -2.0*c*a; 
    (*this)[2][2] = -2.0*c*b;      (*this)[3][2] = 1.0 + c*(a*a+b*b);
    break;
  case identity:
  default:
    (*this)[0][0] = 1.0;
    (*this)[1][1] = 1.0;
    (*this)[2][2] = 1.0;
  }
}


///////////////////////////////////////////////////////////////////////////////
//
// Debugging: prints out std::valarray<double>s
//
///////////////////////////////////////////////////////////////////////////////
std::ostream & operator <<(std::ostream &out, const std::valarray<double> &v) {
  out << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")";
  return(out);
}


