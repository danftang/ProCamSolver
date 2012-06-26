///////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2012 Daniel Tang.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing,
//   software distributed under the License is distributed on an "AS
//   IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
//   express or implied.  See the License for the specific language
//   governing permissions and limitations under the License.
//
///////////////////////////////////////////////////////////////////////////////
#include "MeasurementMatrix.h"
#include "BundleAdjuster.h"
#include "SturmTriggsSolver.h"

//////////////////////////////////////////////////////////////////////////////
/// returns a refrence to the v'th camera matrix
//////////////////////////////////////////////////////////////////////////////
template<int V>
typename MotionMatrix<V>::Block34 MotionMatrix<V>::view(int v) {
  return(Base::template middleRows<3>(v*3));
}

template<int V>
typename MotionMatrix<V>::constBlock34 MotionMatrix<V>::view(int v) const {
  return(Base::template middleRows<3>(v*3));
}


//////////////////////////////////////////////////////////////////////////////
/// returns a refrence to the first 3 columns of the v'th camera matrix
//////////////////////////////////////////////////////////////////////////////
template<int V>
typename MotionMatrix<V>::Block33 MotionMatrix<V>::M(int v) {
  return(Base::template block<3,3>(v*3,0));
}

template<int V>
typename MotionMatrix<V>::constBlock33 MotionMatrix<V>::M(int v) const {
  return(Base::template block<3,3>(v*3,0));
}


//////////////////////////////////////////////////////////////////////////////
/// returns a refrence to the last column of the v'th camera matrix
//////////////////////////////////////////////////////////////////////////////
template<int V>
typename MotionMatrix<V>::Block31 MotionMatrix<V>::T(int v) {
  return(Base::template block<3,1>(v*3,3));
}

template<int V>
typename MotionMatrix<V>::constBlock31 MotionMatrix<V>::T(int v) const {
  return(Base::template block<3,1>(v*3,3));
}


//////////////////////////////////////////////////////////////////////////////
/// Calculate the motion matrix from the given measurement matrix, holding
/// some or none of the camera intrinsics fixed. No initial guess of this
/// motion matrix is necessary.
///
/// \param E The measurement matrix containing the pixel
/// correspondences to be fitted. The measurement matrix may contain
/// occlusions and need not be scaled.
///
/// \param fixedIntrinsics A transformation matrix containing the
/// intrinsics of the cameras that have already been calibrated (if
/// any). If this parameter is missing, all intrinsics are variable.
///
/// The following method is used:
///
/// 1) Calculate projective depths from the Fundamental matrices
/// (Sturm and Triggs 1996)
///
/// 2) Approximate the shape/motion matrices using method of Martinec
/// and Pajdla (2002)
///
/// 3) Use the approximate motion matrix to fill in occlusions and
/// unknown projective depths
///
/// 4) Take Singular Value Decomposition to factorise into
/// shape/motion matrices
///
//////////////////////////////////////////////////////////////////////////////
template<int V>
void MotionMatrix<V>::svd_solve(const MeasurementMatrix<V> &measurement) {

  ShapeMatrix		shape;
  SturmTriggsSolver<V>	stSolver(*this, shape);

  stSolver.solve(measurement);

  /******
  MeasurementMatrix<V> 	scaledMeasurement(measurement);
  ImageTransform<V>	denormalisation;
  ShapeMatrix		shape;

  scaledMeasurement.normalise(denormalisation);
  scaledMeasurement.approx_scale();
  scaledMeasurement.factorise_with_occlusions(*this, shape);
  scaledMeasurement.scale_and_fill(*this, shape);

  std::cout << "Filled measurement = " << std::endl;
  std::cout << scaledMeasurement << std::endl;

  scaledMeasurement.svd_factorise(*this, shape);
  denormalisation.apply_to(*this);
  *****/
}


//////////////////////////////////////////////////////////////////////////////
/// Minimises reprojection error using bundle adjustment on Euclidean
/// parameters (i.e. camera intrinsics/rotation/translation).
//////////////////////////////////////////////////////////////////////////////
template<int V>
double MotionMatrix<V>::bundle_adjust(int f, const MeasurementMatrix<V> &measurement) {
  int			info;
  Eigen::VectorXd	err;
  MeasurementMatrix<V> 	scaledMeasurement(measurement);
  ImageTransform<V>	denormalisation;
  BundleAdjuster<V>	ba_solver(scaledMeasurement, *this);

  // scaledMeasurement.normalise(denormalisation, true);
  // denormalisation.apply_inverse_to(*this);
  ba_solver.levmar_solve(f);
  // denormalisation.apply_to(*this);

  reprojection_err(measurement, err);
  return(err.norm());
}


//////////////////////////////////////////////////////////////////////////////
///
/// Calculates the projective transform that lifts this to a Euclidean
/// form, i.e. so that each camera matrix can be expressed in the form
/// K[R|t] where K is a camera intrinsic matrix, R is a rotation
/// matrix and t is a translation. If the intrinsic matrices of some
/// of the views are known, these can be supplied in
/// fixedIntrinsics. The first 'F' views are then held fixed at the
/// supplied values.
///
/// The algorithm solves for the absolute dual quadric, on the
/// assumption of zero skew, aspect ratio of 1 (i.e. pixels are
/// square) and a principal point on the y axis. If aspect ratio is
/// known and not 1, the camera matrices in *this can be left
/// multiplied by [1 0 0] [0 a 0] [0 0 1] before euclidean lifting.
/// and is described in Han and Kanade (2000).
///
/// \param nCameras The first 'nCameras' views are assumed to be cameras
///        and so have principal points close to the origin
///
/// \param fixedIntrinsics The intrinsic matrices of pre-calibrated views.
///
//////////////////////////////////////////////////////////////////////////////
template<int V>
template<int F>
Eigen::Matrix4d MotionMatrix<V>::diac_euclidean_lift(int nCameras, const ImageTransform<F> &fixedIntrinsics) {
  int			i,j,k;
  Eigen::Matrix4d	H;
  Eigen::Matrix3d	diac;

  // H = [A|B]

  // --- first compute B, assuming origin
  // --- is at the first view so if P_0 = [R|r]
  // --- and B = [b|1] then Rb = -r
  // ------------------------------
  Eigen::Matrix3d 	R;
  Eigen::Vector3d	r;
 
  R = view(0).leftCols(3);
  r = view(0).col(3);
  H.col(3).topRows(3) = -R.inverse()*r;
  H(3,3) = 1.0;

  // --- now do linear least squares on constraints:
  // ---  	[ a_i  0   0   ]
  // --- w_i = 	[ 0   a_i  0   ] = PQP^T
  // ---	[ 0    0   s_i ]
  // --- and s_1 = 1.0;
  Eigen::Matrix<double, Eigen::Dynamic, 10> coeffs; // coefficients of Q
  Eigen::Matrix<double, Eigen::Dynamic, 1>  y;	// vals of fit coeffs*x=y
  Eigen::Matrix<double, 10, 1>		x;	// result of least square fit
  Eigen::Matrix4d			Q;	// result repacked as matrix

  // --- construct matrix of coefficients of elements of Q
  if(nCameras < F) nCameras = F;
  coeffs.resize(4*F + 3*(nCameras-F) + V + 1,10);
  y.resize(4*F + 3*(nCameras-F) + V + 1);
  y.fill(0.0);
  i = 0;
  for(j = 0; j < V; ++j) {
    if(j<F) {
      // --- constrains for known intrinsics:
      // --- i   -> w_00 - (fx^2 + px^2)w_22 = 0
      // --- i+1 -> w_01 - (pxpy)w_22 = 0
      // --- i+2 -> w_02 - (px)w_22 = 0
      // --- i+3 -> w_11 - (fy^2 + py^2)w_22 = 0
      // --- i+4 -> w_12 - (py)w_22 = 0
      quadratic_prod(view(j).row(2), view(j).row(2), coeffs.row(i+5));
      quadratic_prod(view(j).row(0), view(j).row(0), coeffs.row(i));
      quadratic_prod(view(j).row(0), view(j).row(1), coeffs.row(i+1));
      quadratic_prod(view(j).row(0), view(j).row(2), coeffs.row(i+2));
      quadratic_prod(view(j).row(1), view(j).row(1), coeffs.row(i+3));
      quadratic_prod(view(j).row(1), view(j).row(2), coeffs.row(i+4));
      diac = fixedIntrinsics.view(j) * fixedIntrinsics.view(j).transpose();
      diac /= diac(2,2);
      coeffs.row(i)   -= diac(0,0) * coeffs.row(i+5);
      coeffs.row(i+1) -= diac(0,1) * coeffs.row(i+5);
      coeffs.row(i+2) -= diac(0,2) * coeffs.row(i+5);
      coeffs.row(i+3) -= diac(1,1) * coeffs.row(i+5);
      coeffs.row(i+4) -= diac(1,2) * coeffs.row(i+5);
      i += 5;
    } else {
      // --- constrain for aspect skew + pxpy = 0 
      quadratic_prod(view(j).row(0), view(j).row(1), coeffs.row(i));
      ++i;
      if(j<nCameras) {
	// --- constraints for cameras:
	// --- i    -> x-principal point = 0
	// --- i+1  -> y-principal point = 0
	// --- i+2  -> aspect ratio = 1
	quadratic_prod(view(j).row(0), view(j).row(2), coeffs.row(i));
	quadratic_prod(view(j).row(1), view(j).row(2), coeffs.row(i+1));
	quadratic_prod(view(j).row(0), view(j).row(0), coeffs.row(i+2));
	quadratic_prod(view(j).row(1), view(j).row(1), coeffs.row(i+3)); // tmp
	coeffs.row(i+2) -= coeffs.row(i+3);
	i += 3;
      }
    }
  }
  // --- constraint on scale
  quadratic_prod(view(0).row(2), view(0).row(2), coeffs.row(i));
  y(i) = 1.0;

  // --- solve least squares and reform Q

  x = coeffs.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(y);

  std::cout << "Euclidean lift rms error on 1st pass = " << std::endl;
  std::cout << (coeffs*x - y).norm() << std::endl;

  k = 0;
  for(i = 0; i<4; ++i) {
    for(j = i; j<4; ++j){ 
      Q(i,j) = Q(j,i) = x(k++);
    }
    Q(i,i) *= 2.0;
  }

  // --- re-weight the coefficients and solve again
  double weight;
  i = 0;
  for(j=0; j<V; ++j) {
    weight = 1.0/(view(j) * Q * view(j).transpose())(2,2);
    if(j<F) {
      coeffs.middleRows<5>(i) *= weight;
      i += 5;
    }
    else {
      coeffs.row(i) *= weight;
      ++i;
      if(j<nCameras) {
	coeffs.middleRows<3>(i) *= weight;
	i += 3;
      }
    }
  }

  x = coeffs.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(y);

  std::cout << "Euclidean lift rms error on 2nd pass = " << std::endl;
  std::cout << (coeffs*x - y).norm() << std::endl;

  k = 0;
  for(i = 0; i<4; ++i) {
    for(j = i; j<4; ++j){ 
      Q(i,j) = Q(j,i) = x(k++);
    }
    Q(i,i) *= 2.0;
  }

  // --- extract A from Q
  // --------------------
  Eigen::JacobiSVD<Eigen::Matrix4d> svdQ(Q, Eigen::ComputeFullU);
  H.leftCols(3) = svdQ.matrixU().leftCols(3) * 
    svdQ.singularValues().topRows(3).cwiseSqrt().asDiagonal();

  // --- set reference rotation
  // --- so that P_0 has zero
  // --- rotation: MM^T = [M|t]Q[M|t]^T = KK^T
  // --- PH = K[R|0] so K^-1PH = [R|0]
  // --- so PH[R^-1|0] = K[I|0]
  // ---      [ 0  |1]
  // -----------------------------------------
  Eigen::Matrix3d	K1;   // inverse intrinsics of view 0
  Eigen::Matrix4d	Rot; // rotation
  Eigen::Matrix3d	MMT; // MM^T where P=[M|t]

  MMT = view(0)*Q*view(0).transpose();
  K1 = MMT.inverse().llt().matrixU();

  Rot.fill(0.0);
  Rot.topLeftCorner(3,3) = (K1*(view(0)*H).leftCols(3)).inverse();
  Rot(3,3) = 1.0;

  H *= Rot;
  
  return(H);
}

template<int V>
Eigen::Matrix4d MotionMatrix<V>::diac_euclidean_lift(int nCameras) {
  ImageTransform<0> dummy;
  return(diac_euclidean_lift(nCameras, dummy));
}


//////////////////////////////////////////////////////////////////////////////
/// Calculates the plane at infinity
//////////////////////////////////////////////////////////////////////////////
template<int V>
Eigen::Vector3d MotionMatrix<V>::plane_at_infinity(int nCameras) {
  int				i,j,k;
  Eigen::Matrix<double,4,3>	A;
  Eigen::Vector3d		p; // Plane at infinity

  // --- now do linear least squares on constraints:
  // ---  	[ a_i  0   0   ]
  // --- w_i = 	[ 0   a_i  0   ] = PQP^T
  // ---	[ 0    0   s_i ]
  // --- and s_1 = 1.0;
  Eigen::Matrix<double, Eigen::Dynamic, 10> coeffs; // coefficients of Q
  Eigen::Matrix<double, Eigen::Dynamic, 1>  y;	// vals of fit coeffs*x=y
  Eigen::Matrix<double, 10, 1>		x;	// result of least square fit
  Eigen::Matrix4d			Q;	// result repacked as matrix

  // --- construct matrix of coefficients of elements of Q
  coeffs.resize(3*nCameras + V + 1,10);
  y.resize(3*nCameras + V + 1);
  y.fill(0.0);
  i = 0;
  j = 0;
  for(j = 0; j< nCameras; ++j) {
    // --- equations enforcing zero principal point
    // --- i+1  -> x-principal point = 0
    // --- i+2  -> y-principal point = 0
    quadratic_prod(view(j).row(0), view(j).row(2), coeffs.row(i));
    quadratic_prod(view(j).row(1), view(j).row(2), coeffs.row(i+1));
    quadratic_prod(view(j).row(0), view(j).row(0), coeffs.row(i+2));
    quadratic_prod(view(j).row(1), view(j).row(1), coeffs.row(i+3)); // tmp
    coeffs.row(i+2) -= coeffs.row(i+3);
    i += 3;

  }
  for(j = 0; j<V; ++j) {
    // --- equation enforcing aspect skew + pxpy = 0 
    quadratic_prod(view(j).row(0), view(j).row(1), coeffs.row(i));
    ++i;
  }
  // --- equation to set scale
  quadratic_prod(view(0).row(2), view(0).row(2), coeffs.row(i));
  y(i) = 1.0;

  // --- solve least squares and reform Q

  x = coeffs.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(y);

  std::cout << "coeffs*x = " << std::endl;
  std::cout << coeffs*x << std::endl;

  k = 0;
  for(i = 0; i<4; ++i) {
    for(j = i; j<4; ++j){ 
      Q(i,j) = Q(j,i) = x(k++);
    }
    Q(i,i) *= 2.0;
  }

  // --- extract A from Q
  // --------------------
  Eigen::JacobiSVD<Eigen::Matrix4d> svdQ(Q, Eigen::ComputeFullU);
  A = svdQ.matrixU().leftCols(3) * 
    svdQ.singularValues().topRows(3).cwiseSqrt().asDiagonal();
  p.transpose() = -A.row(3) * A.topRows(3).inverse();
  
  return(p);
}


//////////////////////////////////////////////////////////////////////////////
/// Lift given plane at infinity and assuming zero skew and unit aspect
/// ratio. Solves for the image of the absolute conic.
//////////////////////////////////////////////////////////////////////////////
template<int V>
Eigen::Matrix4d MotionMatrix<V>::euclidean_lift(const Eigen::Vector3d &p) {
  int 					i,j,k;
  Eigen::Matrix<double, 2*V + 1, 6> 	coeffs;
  Eigen::Matrix<double, 2*V + 1, 1> 	y;
  Eigen::Matrix<double, 6, 1> 		x;
  Eigen::Matrix3d			HT;  // inverse transpose of the 
                                             // infinite homography
  Eigen::Matrix3d			iac; // image of the absolute conic
  Eigen::Matrix3d			K;   // camera intrinsics
  Eigen::Matrix4d			H;   // transform 

  y.fill(0.0);
  i = 0;
  for(j = 0; j<V; ++j) {
    // --- If P = [A|a], H^-T = (A - ap^T)^-T
    // --------------------------------------
    HT.transpose() = (view(j).leftCols(3) - 
		      view(j).col(3)*p.transpose()).inverse();
    if(j == 0) {
      // --- equation to set scale
      quadratic_prod(HT.row(0), HT.row(2), coeffs.row(i)); // zero skew
      quadratic_prod(HT.row(1), HT.row(2), coeffs.row(i+1)); // aspect ratio=1
      quadratic_prod(HT.row(2), HT.row(2), coeffs.row(i+2)); // aspect ratio=1
      coeffs.row(i) += coeffs.row(i+1);
      coeffs.row(i) += coeffs.row(i+2);
      y(i) = 1.0;
      ++i;
    }
    quadratic_prod(HT.row(0), HT.row(0), coeffs.row(i)); // aspect ratio=1
    quadratic_prod(HT.row(1), HT.row(1), coeffs.row(i+1));
    coeffs.row(i) -= coeffs.row(i+1);
    quadratic_prod(HT.row(0), HT.row(2), coeffs.row(i+1)); // zero skew
    i += 2;
  }

  x = coeffs.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(y);

  std::cout << "coeffs*x = " << std::endl;
  std::cout << coeffs*x << std::endl;

  k = 0;
  for(i = 0; i<3; ++i) {
    for(j = i; j<3; ++j){ 
      iac(i,j) = iac(j,i) = x(k++);
    }
    iac(i,i) *= 2.0;
  }

  K = iac.llt().matrixU();
  H.topLeftCorner(3,3) = K.inverse();
  H.leftCols(3).row(3) = -p.transpose()*H.topLeftCorner(3,3);

  // --- set origin at view 0
  H.col(3).topRows(3) = -view(0).leftCols(3).inverse()*view(0).col(3);
  H(3,3) = 1.0;

  return(H);
}

//////////////////////////////////////////////////////////////////////////////
/// Extracts the intrinsic / rotation matrices from view 'v of
/// *this. Works by calculating $MM^T$ where the camera matrix P =
/// [M|T] and M = sKR (see Han and Kanade, 2000) where s is a scaling
/// factor, K is an upper triangular (i.e. intrinsics) matrix with 1.0
/// in the lower-right element and R is a rotation matrix
///
/// \param K affine Transform to be filled with the camera intrinsics
/// \param R AngleAxis to be filled with the camera rotation
///
/// \return the scaling factor, s.
//////////////////////////////////////////////////////////////////////////////
template<int VIEWS>
double MotionMatrix<VIEWS>::KR_decompose(int i, Eigen::Affine2d &K,
					        Eigen::AngleAxisd &R) const {
  Eigen::Matrix3d 	MMT;		// M*M^T
  Eigen::Matrix3d 	RMat;		// Rotation in matrix form
  double		scale;		// scaling factor
  Eigen::Affine2d 	K1;		// inverse intrinsics

  // --- Do Cholesky decomposition of
  // --- (MM^T)^-1 = K^-TR^-TR^-1K^-1
  // ---           = K^-TK^-1 = LL^T
  // --- then    R = K^-1M
  // --------------------------------
  MMT   = M(i) * M(i).transpose();
  K1.matrix()  = MMT.inverse().llt().matrixU();
  RMat  = K1.matrix() * M(i);
  scale = 1.0/K1(2,2);
  K1.matrix() *= scale;
  K     = K1.inverse();

  /************
  // --- Do QR decomposition where
  // --- if X.fliprows^T = QR then
  // --- X = R^T.fliprows.flipcols Q^T.fliprows
  // ------------------------------------------
  K.matrix().triangularView<Eigen::Upper>() = 
    M(i).colwise().reverse().transpose().
    householderQr().matrixQR().
    transpose().colwise().reverse().rowwise().reverse();
  K(1,0) = 0.0;
  // --- ensure +ve focal lengths
  int j;
  for(j = 0; j<3; ++j) {
    if(K(j,j) < 0.0) {
      K.matrix().col(j) *= -1.0;
    }
  }
  scale = K(2,2);
  K.matrix() /= scale;
  RMat = (K.inverse() * M(i))/scale;
  **********/

  // --- remove reflections and
  // --- extract rotation angles
  if(RMat.determinant() < 0.0) {
    RMat *= -1.0;
    scale *= -1.0;
  }
  R = RMat;

  return(scale);
}


//////////////////////////////////////////////////////////////////////////////
/// Splits *this into intrinsic / rotation matrices. *this should previously
/// have been transformed to Euclidean. Works by calculating $MM^T$ where
/// the camera matrix P = [M|T] and M = KR (see Han and Kanade, 2000)
///
/// \param K matrix to be filled with the camera intrinsics
/// \param R matrix to be filled with the camera rotation
///
/// \return A measure of the extent that *this deviates from the
/// Euclidean. This is calculated as the sum over all views of the
/// absolute values of the skews, which should be zero for perfectly
/// Euclidean cameras.
//////////////////////////////////////////////////////////////////////////////
template<int VIEWS>
double MotionMatrix<VIEWS>::KR_decompose(ImageTransform<VIEWS> &K,
					 ImageTransform<VIEWS> &R) const {
  int 			i;
  Eigen::Affine2d	Kn;
  Eigen::AngleAxisd	Rn;

  for(i = 0; i<VIEWS; ++i) {
    KR_decompose(i,Kn,Rn);
    K.view(i) = Kn.matrix();
    R.view(i) = Rn.toRotationMatrix();
  }
  return(0.0);
}


//////////////////////////////////////////////////////////////////////////////
/// Calculates the reprojection error of *this compared to E
//////////////////////////////////////////////////////////////////////////////
template<int V>
void MotionMatrix<V>::reprojection_err(
				       const MeasurementMatrix<V> &E,
				       Eigen::VectorXd &err) {
  int 				i,j,n;
  Eigen::Matrix<double, 3*V, 1>	reproj;
  ShapeMatrix 			shape;
  bool				resizeErr;

  if(err.size() == 0) {
    err.resize(2*V*E.cols());
    resizeErr = true;
  } else {
    resizeErr = false;
  }

  shape.solve(E, *this);
  n = 0;
  for(j=0; j<E.cols(); ++j) {
    reproj = *this * shape.col(j);
    for(i = 0; i<V; ++i) {
      if(!E.pixel_is_occluded(i,j)) {
	err(n++) = reproj(3*i)  /reproj(3*i+2)-E.pixel(i,j)(0)/E.pixel(i,j)(2);
	err(n++) = reproj(3*i+1)/reproj(3*i+2)-E.pixel(i,j)(1)/E.pixel(i,j)(2);
      }
    }
  }

  if(resizeErr) {
    err.conservativeResize(n);
  } else {
    while(n < err.size()) {
      err(n++) = 0.0;
    }
  }
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
template<int V>
template<class D>
MotionMatrix<V> &MotionMatrix<V>::operator=(const Eigen::MatrixBase<D> &o) {
  Base::operator=(o);
  return(*this);
}


//////////////////////////////////////////////////////////////////////////////
/// Sets result to the elements of v*w^T
//////////////////////////////////////////////////////////////////////////////
template<int V>
void MotionMatrix<V>::quadratic_prod(const Eigen::Vector4d &v,
				     const Eigen::Vector4d &w,
				     typename Eigen::Matrix<double,Eigen::Dynamic,10>::RowXpr result) {
  int i,j,k;

  k = 0;
  for(i = 0; i<4; ++i) {
    for(j = i; j<4; ++j) {
      result(k++) = v(i)*w(j) + v(j)*w(i);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
template<int VIEWS>
void MotionMatrix<VIEWS>::quadratic_prod(const Eigen::Vector3d &v,
					 const Eigen::Vector3d &w,
		 typename Eigen::Matrix<double,2*VIEWS+1,6>::RowXpr result) {
  int i,j,k;

  k = 0;
  for(i = 0; i<3; ++i) {
    for(j = i; j<3; ++j) {
      result(k++) = v(i)*w(j) + v(j)*w(i);
    }
  }
}
