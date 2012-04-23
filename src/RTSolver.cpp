#include "RTSolver.h"
#include "numerical.h"

const double RTSolver::PI = 3.1415926535897932;

///////////////////////////////////////////////////////////////////////////////
//
// Constructor/solver
// f  = the fundamental matrix between two views
// m  = the intrinsic matrix of the reference camera
// mp = the intrinsic matrix of the second camera (M')
// 
// on construction, 'this' is set to be equal to the camera matrix of
// the second view.
//
// If F is exact, we have that 
// F = M'^{-1}RTM^{-1}
// so, if we let r be the residual matrix:
// r = M'^{-1}RTM^{-1} - F
// where
// M'^{-1T} 	= transpose of the inverse intrinsic matrix of the 2nd camera
// R		= rotation matrix (from reference to second view)
// T		= translation cross-product matrix (TQ = t x Q)
// M^{-1}	= inverse intrinsic matrix of the reference camera
// F		= the fundamental matrix
//
// then R and T are the unknowns with 6 degrees of freedom.
// These are calculated by minimising the sum of the squres of the elements
// of r. The camera matrix pair is then given by
//
// P = M[I|0] and P'=M'R[I|-t]
//
// So P' is the camera matrix we require.
//
///////////////////////////////////////////////////////////////////////////////
RTSolver::RTSolver(Matrix<double> &f,
		   Matrix<double> &m,
		   Matrix<double> &mp) : 
  transformMatrix(camera),
  N(inverse_intrinsic, m[0][0], m[1][1], m[0][2], m[1][2]),
  NPT(inverse_intrinsic, mp[0][0], mp[1][1], mp[0][2], mp[1][2])
 {
   transformMatrix	I(identity);
   Frprmn<RTSolver>	mySolver(*this, 1e-8);	// solver
   VecDoub 		p(6); 			// parameters for solver
   coord		t(3);			// translation vector

   F = f;
   NPT.transpose();

   // --- initial guess
   // -----------------
   p[0] = 1.2;
   p[1] = 0.0;
   p[2] = 0.0;
   p[3] = 0.0;
   p[4] = 0.0;
   p[5] = 0.0;

   // --- solve
   // ---------
   p = mySolver.minimize(p);

   // --- extract solution
   // --------------------
   t[0] = -p[0];
   t[1] = -p[1];
   t[2] = -p[2];
   R = transformMatrix(xrotate,p[3]) *
     transformMatrix(yrotate,p[4]) *
     transformMatrix(zrotate,p[5]);

   P() = mp*R*(I|t);
}


///////////////////////////////////////////////////////////////////////////////
//
// return sum of squares of errors, where 
//
// e^2 = (M'^{-1T}RTM^{-1} - F).(M'^{-1T}RTM^{-1} - F)
//
// where the dot product means element-wise multiplication and summing.
//
///////////////////////////////////////////////////////////////////////////////
double RTSolver::operator()(VecDoub_I &p) {
   transformMatrix T(cross_prod,p[0],p[1],p[2]);	// translation matrix
   double e2;
   int i,j;

   R = transformMatrix(xrotate,p[3]) *
     transformMatrix(yrotate,p[4]) *
     transformMatrix(zrotate,p[5]);
   T = NPT*R*T*N - F;
   e2 = elemental_dot_prod(T,T);
   return(e2);
}


///////////////////////////////////////////////////////////////////////////////
//
// calculate rate of change of error with parameters, de2/dp, at p and
// put result in 'deriv'.
//
// Uses the identity d(sum_n(y_n^2))/dx = sum_n(2y_n * dy_n/dx)
//
// and the identities d(sin(t))/dt = sin(t + PI/2) and 
// d(cos(t))/dt = cos(t + PI/2)
//
///////////////////////////////////////////////////////////////////////////////
void RTSolver::df(VecDoub_I &p, VecDoub_O &deriv) {
   transformMatrix T(cross_prod,p[0],p[1],p[2]);	// translation matrix
   transformMatrix e;
   transformMatrix NPT_R;
   transformMatrix T_N;
   transformMatrix de_dp;

   R = 	
     transformMatrix(xrotate,p[3]) *
     transformMatrix(yrotate,p[4]) *
     transformMatrix(zrotate,p[5]);
   NPT_R = NPT*R;
   T_N   = T*N;
   e = NPT_R*T_N - F;

   de_dp = 0.0; de_dp[1][2] = -1.0; de_dp[2][1] = 1.0;
   deriv[0] = 2.0*elemental_dot_prod(e, NPT_R*de_dp*N);

   de_dp = 0.0; de_dp[0][2] = 1.0; de_dp[2][0] = -1.0;
   deriv[1] = 2.0*elemental_dot_prod(e, NPT_R*de_dp*N);

   de_dp = 0.0; de_dp[0][1] = -1.0; de_dp[1][0] = 1.0;
   deriv[2] = 2.0*elemental_dot_prod(e, NPT_R*de_dp*N);

   de_dp = 
     transformMatrix(xrotate,p[3] + PI/2) * 
     transformMatrix(yrotate,p[4]) *
     transformMatrix(zrotate,p[5]);
   deriv[3] = 2.0*elemental_dot_prod(e, NPT*de_dp*T_N);

   de_dp = 
     transformMatrix(xrotate,p[3]) * 
     transformMatrix(yrotate,p[4] + PI/2) *
     transformMatrix(zrotate,p[5]);
   deriv[4] = 2.0*elemental_dot_prod(e, NPT*de_dp*T_N);

   de_dp = 
     transformMatrix(xrotate,p[3]) * 
     transformMatrix(yrotate,p[4]) *
     transformMatrix(zrotate,p[5] + PI/2);
   deriv[5] = 2.0*elemental_dot_prod(e, NPT*de_dp*T_N);
}


///////////////////////////////////////////////////////////////////////////////
//
// Returns sum_i sum_j M_{ij}N_{ij}
//
///////////////////////////////////////////////////////////////////////////////
double RTSolver::elemental_dot_prod(const Matrix<double> &M, 
				    const Matrix<double> &N) {
  int i,j;
  double result;

  result = 0.0;
   for(i = 0; i<M.isize(); ++i) {
     for(j = 0; j<M.jsize(); ++j) {
       result += M[i][j]*N[i][j];
     }
   }
   return(result);
}

