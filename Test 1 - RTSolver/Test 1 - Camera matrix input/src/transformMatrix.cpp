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
//
///////////////////////////////////////////////////////////////////////////////
transformMatrix::transformMatrix(matrixType t, 
				double a, double b, 
				double c, double d) : Matrix<double>(3,3) {
  int i,j;

  (*this)[0][0] = 1.0;
  (*this)[0][1] = 0.0;
  (*this)[0][2] = 0.0;
  (*this)[1][0] = 0.0;
  (*this)[1][1] = 1.0;
  (*this)[1][2] = 0.0;
  (*this)[2][0] = 0.0;
  (*this)[2][1] = 0.0;
  (*this)[2][2] = 1.0;

  switch(t) {
  case intrinsic:
    (*this)[0][0] = a;
    (*this)[1][1] = b;
    (*this)[0][2] = c;
    (*this)[1][2] = d;
    break;
  case inverse_intrinsic:
    (*this)[0][0] = 1.0/a;
    (*this)[1][1] = 1.0/b;
    (*this)[0][2] = -c/a;
    (*this)[1][2] = -d/b;
    break;
  case xrotate:
    (*this)[2][2] = (*this)[1][1] = cos(a);
    (*this)[2][1] = sin(a);
    (*this)[1][2] = -(*this)[2][1];
    break;
  case yrotate:
    (*this)[2][2] = (*this)[0][0] = cos(a);
    (*this)[0][2] = sin(a);
    (*this)[2][0] = -(*this)[0][2];
    break;
  case zrotate:
    (*this)[1][1] = (*this)[0][0] = cos(a);
    (*this)[1][0] = sin(a);
    (*this)[0][1] = -(*this)[1][0];
    break;
  case cross_prod:
    (*this)[0][0] = 0.0; (*this)[0][1] = -c;  (*this)[0][2] = b;
    (*this)[1][0] = c;   (*this)[1][1] = 0.0; (*this)[1][2] = -a;
    (*this)[2][0] = -b;  (*this)[2][1] = a;   (*this)[2][2] = 0.0;
    break;
  case camera:
    resize(3,4);
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


