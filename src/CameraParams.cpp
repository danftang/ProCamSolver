#include "CameraParams.h"

///////////////////////////////////////////////////////////////////////////////
//
// returns the error for a given set of parameters
//
///////////////////////////////////////////////////////////////////////////////
double RPSolver::operator()(const InputType &params, ValueType &errs) {
  double 	e2;
  coord		pi1(3);
  coord		pi2(3);
  coord		point3D(4);
  int		N;
  int		i;
  int		n; // 3D point id
  int		d; // dimension
  int		cam; // camera id
  EquivalencePartition<GPixel>::member_iterator pixel1, pixel2;

  // std::cout << "Called f for m = " << m << std::endl;

  params_to_vpmatrices(p);
  e2  = 0.0;
  n   = m%C.size();
  cam = (m/C.size())%size();
  d   = m/(C.size()*size());
  N = 0;
  point3D = 0.0;
  for(pixel2= ++(C[n].begin()); pixel2 != C[n].end(); ++pixel2) {
    for(pixel1 = C[n].begin(); pixel1 != pixel2; ++pixel1) {
      closest_points(*pixel1, *pixel2, pi1, pi2);
      for(i = 0; i<3; ++i) {	
	point3D[i] += pi1[i] + pi2[i];
      }
      N += 2;
    }
  }
  point3D /= N;
  point3D[3] = 1.0;

  // --- calculate reprojection error for camera m
  for(pixel1= C[n].begin(); pixel1 != C[n].end(); ++pixel1) {
    // calc reprojected pixel
    pi1 = P[pixel1->id]*point3D;
    if(pixel1->id == cam) {
      if(d == 0) {
	e2 = pixel1->x - pi1[0]/pi1[2];
      } else {
	e2 = pixel1->y - pi1[1]/pi1[2];
      }
    }
  }

  return(e2);
}


///////////////////////////////////////////////////////////////////////////////
/// Returns the rotation of the v'th view as an AngleAxis
///////////////////////////////////////////////////////////////////////////////
template<int M, int F>
Eigen::AngleAxisd CameraParams<M,F>::rotation(int v) {
  if(v == 0) return(Eigen::AngleAxisd(0,Vector3d::UnitX()));
  v -= 1;
  return(Eigen::AngleAxisd(variableParams(3*v),Vector3d::UnitX())*
	 Eigen::AngleAxisd(variableParams(3*v+1),Vector3d::UnitY())*
	 Eigen::AngleAxisd(VariableParams(3*v+2),Vector3d::UnitZ()));
}


///////////////////////////////////////////////////////////////////////////////
/// Returns the camera intrinsics of the v'th view
///////////////////////////////////////////////////////////////////////////////
template<int M, int F>
inline const Eigen::Transform2d CameraParams<M,F>::intrinsics(int v) {
  return(Translation2d(principal_pt(v)) * Scaling(focal_len(v)));
}


///////////////////////////////////////////////////////////////////////////////
/// returns the principal point of the v'th view.
///////////////////////////////////////////////////////////////////////////////
template<int M, int F>
inline CameraParams::PrincipalPoint CameraParams<M,F>::principal_pt(int v) {
  if(v < F) return(fixedParams.middleRows<2>(3*v));
  return(variableParams.middleRows<2>(6*(M-1) + 3*(v-F)));
}

///////////////////////////////////////////////////////////////////////////////
/// returns the focal length of the v'th view.
///////////////////////////////////////////////////////////////////////////////
template<int M, int F>
inline double & CameraParams<M,F>::focal_len(int v) {
  if(v < F) return(fixedParams(3*v+2));
  return(variableParams(6*(M-1) + 3*(v-F) + 2));
}



///////////////////////////////////////////////////////////////////////////////
//
// copies the parameters in p to this map of view projection matrices
//
///////////////////////////////////////////////////////////////////////////////
template<int M, int F>
void CameraParams<M,F>::params_to_matrices() {
  int			v;
  ImageTransform<M>  	R;
  ImageTransfrom<M>  	R1;
  

  for(v = 0; v < M; ++v) {
    Pinv.view(v) = rotation(v).inverse() * intrinsics(v).inverse();
    P.M(v) = intrinsics(v) * rotation(v);
    P.T(v) = P.M(v)*translation(v);
  }
}


///////////////////////////////////////////////////////////////////////////////
// 
// returns the Sampson reprojection error due to correspondence, c.
// R is the (pre-computed) radial fundamental matrix for the two
// views of the correspondence.
//
///////////////////////////////////////////////////////////////////////////////
void RPSolver::closest_points(const GPixel &i, const GPixel &j, coord &p, coord &q) {
  //std::cout << "solving " << i << " " << j << std::endl;
  
  RadialViewProjectionMatrix &	Pi((*this)[i.id]);
  RadialViewProjectionMatrix &	Pj((*this)[j.id]);	  
  coord pixi(3);// pixel i
  coord pixj(3);// pixel j
  coord I(3);	// direction of ray i
  coord J(3);	// direction of ray j
  coord t(3);	// tranlation between camera centres
  coord yhat(3);// basis vector in epipolar plane perp. to I
  double k;	// distance along I of closest point

  pixi[0] = i.x;
  pixi[1] = i.y;
  pixi[2] = 1.0;
  pixj[0] = j.x;
  pixj[1] = j.y;
  pixj[2] = 1.0;


  //I = Pi.R1()*(Pi.M1()*pixi);
  //J = Pj.R1()*(Pj.M1()*pixj);
  I = P1[i.id]*pixi;
  J = P1[j.id]*pixj;
  t = Pj.t - Pi.t;

  //transformMatrix tx(cross_prod,t[0], t[1], t[2]);
  transformMatrix Ix(cross_prod,I[0], I[1], I[2]);
  transformMatrix Jx(cross_prod,J[0], J[1], J[2]);

  //yhat = Ix*tx*I;
  //k = ((t*I).sum() - (t*yhat).sum()*(J*I).sum()/(J*yhat).sum())/(I*I).sum();
  yhat = Jx*Ix*J;
  k = (t*yhat).sum()/(I*yhat).sum();

  p = k*I + Pi.t;

  //yhat = Jx*tx*J;
  //k = (-(t*J).sum() + (t*yhat).sum()*(J*I).sum()/(I*yhat).sum())/(J*J).sum();
  yhat = Ix*Jx*I;
  k = -(t*yhat).sum()/(J*yhat).sum();
  q = k*J + Pj.t;

  //std::cout << p << " -> " << q << std::endl;
}


