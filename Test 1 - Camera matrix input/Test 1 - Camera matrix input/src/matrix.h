///////////////////////////////////////////////////////////////////////////////
//
// Simple class to store a rectangular matrix (x_00...x_nm) as a 
// 1D array of doubles in the form
// (x_00, x_10...x_n0, x_01, x_11...x_nn)
//
// (c) Daniel Tang 2009
//
///////////////////////////////////////////////////////////////////////////////

#ifndef POLYFIT_MATRIX_H
#define POLYFIT_MATRIX_H

#include <cmath>
#include <new>
#include <iostream>

template<class DATATYPE>
class Matrix {
public:
  Matrix(int isize, int jsize);
  Matrix();
  ~Matrix();

  int 		jsize() const		       	{return(cols);}
  int		isize() const 			{return(SIZE/cols);}
  int 		size() const		       	{return(SIZE);}
  void		resize(int, int);
  DATATYPE *	array() 			{return(x);}
  const DATATYPE *array() const			{return(x);}
  void		clear();
  void		transpose();

  DATATYPE *		operator[](int i) 		{return(x+i*cols);}
  const DATATYPE *	operator[](int i) const 	{return(x+i*cols);}
  const DATATYPE &	operator()(int i, int j) const	{return(x[i*cols+j]);}
  std::valarray<DATATYPE> operator *(const std::valarray<DATATYPE> &);
  Matrix  		operator *(const Matrix &);
  Matrix  		operator -(const Matrix &);
  Matrix  		operator +(const Matrix &);
  Matrix  		operator |(const std::valarray<DATATYPE> &);
  Matrix &     		operator *=(Matrix &);
  Matrix &     		operator =(const Matrix &);
  Matrix &     		operator =(double);

protected:
  DATATYPE * 	x;
  int 		cols;
  int		SIZE;
};


///////////////////////////////////////////////////////////////////////////////
//
// Constructors/Deconstructors
//
///////////////////////////////////////////////////////////////////////////////

template<class DATATYPE>
Matrix<DATATYPE>::Matrix() : cols(0), SIZE(0) {
}


template<class DATATYPE>
Matrix<DATATYPE>::Matrix(int isize, int jsize) : cols(jsize) {
  SIZE = isize*jsize;
  x = new double [SIZE];
}


template<class DATATYPE>
Matrix<DATATYPE>::~Matrix() {
  delete [] x;
}


///////////////////////////////////////////////////////////////////////////////
//
// resizes this matrix, invalidating its contents.
//
///////////////////////////////////////////////////////////////////////////////
template<class DATATYPE>
void Matrix<DATATYPE>::resize(int isize, int jsize) {
  delete [] x;
  cols = jsize;
  SIZE = isize*jsize;
  x = new double [SIZE];
}


///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
template<class DATATYPE>
void Matrix<DATATYPE>::clear() {
  (*this) = 0;
}


///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
template<class DATATYPE>
std::valarray<DATATYPE> Matrix<DATATYPE>::operator *(const std::valarray<DATATYPE> &v) {
  std::valarray<DATATYPE>  result(isize());
  int           	   i,j;

  if(jsize() != v.size()) throw("Multiplying matrix with incompatible vector");

  for(i = 0; i<isize(); ++i) {
    result[i] = 0.0;
    for(j = 0; j<jsize(); ++j) {
      result[i] += (*this)[i][j]*v[j];
    }
  }
  return(result);
}


///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
template<class DATATYPE>
Matrix<DATATYPE> Matrix<DATATYPE>::operator *(const Matrix<DATATYPE> &M) {
  Matrix<DATATYPE>	result(isize(),M.jsize());
  int           	i,j,k;
  
  if(jsize() != M.isize()) throw("Multiplying incompatible matrices");

  for(i = 0; i<isize(); ++i) {
    for(j = 0; j<M.jsize(); ++j) {
      result[i][j] = 0.0;
      for(k = 0; k<jsize(); ++k) {
	result[i][j] += (*this)[i][k]*M[k][j];
      }
    }
  }
  return(result);
}


///////////////////////////////////////////////////////////////////////////////
//
// N = [M|v]
//
// add a column, v, to M and return the result
//
///////////////////////////////////////////////////////////////////////////////
template<class DATATYPE>
Matrix<DATATYPE> Matrix<DATATYPE>::operator |(const std::valarray<DATATYPE> &v) {
  Matrix<DATATYPE>	result(isize(),jsize()+1);
  int           	i,j;
  
  for(i = 0; i<isize(); ++i) {
    for(j = 0; j<jsize(); ++j) {
      result[i][j] = (*this)[i][j];
    }
    result[i][j] = v[i];
  }
  return(result);
}


///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
template<class DATATYPE>
Matrix<DATATYPE> Matrix<DATATYPE>::operator +(const Matrix<DATATYPE> &M) {
  Matrix<DATATYPE>	result(isize(),jsize());
  int           	i;
  
  for(i = 0; i<size(); ++i) {
      result.x[i] = x[i] + M.x[i];
  }
  return(result);
}


///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
template<class DATATYPE>
Matrix<DATATYPE> Matrix<DATATYPE>::operator -(const Matrix<DATATYPE> &M) {
  Matrix<DATATYPE>	result(isize(),jsize());
  int           	i;
  
  for(i = 0; i<size(); ++i) {
      result.x[i] = x[i] - M.x[i];
  }
  return(result);
}


///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
template<class DATATYPE>
Matrix<DATATYPE> &Matrix<DATATYPE>::operator *=(Matrix<DATATYPE> &N) {
  int i,j,k;
  std::valarray<DATATYPE> t(jsize());

  if(jsize() != N.isize()) throw("Multiplying incompatible matrices");

  for(i = 0; i<isize(); ++i) {
    for(j=0; j<jsize(); ++j) {
      t[j] = 0.0;
      for(k=0; k<jsize(); ++k) {
	t[j] += (*this)[i][k]*N[k][j];
      }
    }
    for(j=0; j<3; ++j) {
      (*this)[i][j]  = t[j];
    }
  }
  return(*this);
}


///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
template<class DATATYPE>
Matrix<DATATYPE> &Matrix<DATATYPE>::operator =(const Matrix<DATATYPE> &N) {
  int i;

  if(SIZE != N.SIZE) {
    resize(N.isize(), N.jsize());
  } else {
    cols = N.cols;
  }

  for(i = 0; i<SIZE; ++i) {
    x[i] = N.x[i];
  }
}


///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
template<class DATATYPE>
Matrix<DATATYPE> &Matrix<DATATYPE>::operator =(double c) {
  int i;

  for(i = 0; i<SIZE; ++i) {
    x[i] = c;
  }
}

///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
template<class DATATYPE>
void Matrix<DATATYPE>::transpose() {
  DATATYPE tmp;
  int i,j;

  for(i = 0; i<isize(); ++i) {
    for(j = i+1; j<jsize(); ++j) {
      tmp = (*this)[i][j];
      (*this)[i][j] = (*this)[j][i];
      (*this)[j][i] = tmp;
    }
  }
}



///////////////////////////////////////////////////////////////////////////////
//
// Prints the martix to 'out'
//
///////////////////////////////////////////////////////////////////////////////
template<class DATATYPE>
std::ostream &operator <<(std::ostream &out, const Matrix<DATATYPE> &M) {
  int i,j;

  out.precision(6);

  for(i=0; i<M.isize(); ++i) {
    out << "( ";
    for(j=0; j<M.jsize(); ++j) {
      out << std::scientific << M[i][j] << "\t";
    }
    out << ")" << std::endl;
  }
  return(out);
}


#endif
