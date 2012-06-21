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
#include "ImageTransform.h"
#include "MotionMatrix.h"
#include "MeasurementMatrix.h"


///////////////////////////////////////////////////////////////////////////////
/// Construct with an unweighted measurement matrix and a motion
/// matrix that is to hold the solution.
///////////////////////////////////////////////////////////////////////////////
template<int M>
BundleAdjuster<M>::BundleAdjuster(const MeasurementMatrix<M> &E,
				  MotionMatrix<M> &P) :
  measurement(E),
  motion(P)
{
}

///////////////////////////////////////////////////////////////////////////////
/// Calculates the reprojection error for a given set of
/// parameters. Called by the Lev-Mar solving algorithm.
///
/// \param params The parameters holding the camera intrinsics/extrinsics
/// \param errs The repreojection errors
/// \returns Zero, not clear in the Eigen documentation what this
/// should be, but looking at the code, this should always be zero.
///////////////////////////////////////////////////////////////////////////////
template<int M>
int BundleAdjuster<M>::operator()(const Eigen::VectorXd &params, 
			 Eigen::VectorXd &errs) {

  params_to_matrices(params);
  motion.reprojection_err(measurement, errs);
  // --- add error for scale
  errs(errs.size()-1) = motion.col(3).norm() - M;

  //std::cout << "params are = " << std::endl;
  //std::cout << params << std::endl;
  //std::cout << "Size error = " << errs(errs.size()-1) << std::endl;
  //std::cout << "Reprojection error = " << std::endl;
  //std::cout << errs << std::endl;
  std::cout << "Err = " << errs.norm()/(2.0*M*measurement.cols()*0.92) << std::endl;
  return(0);
}


///////////////////////////////////////////////////////////////////////////////
/// Converts a rotation stored as 3 angles into an AngleAxis.
///
/// If the rotation angles are (theta, psi, phi) and the AngleAxis is
/// (a,v0,v1,v2) then
///
/// a = theta
/// v0 = cos(psi)cos(phi)
/// v1 = sin(psi)cos(phi)
/// v2 = sin(phi)
///
///////////////////////////////////////////////////////////////////////////////
template<int M>
Eigen::AngleAxisd BundleAdjuster<M>::rotation_to_angleaxis(constBlock3 rot) {
  double cosphi = cos(rot(2));
  return(Eigen::AngleAxisd(rot(0),
	  Eigen::Vector3d(cos(rot(1))*cosphi,
			  sin(rot(1))*cosphi,
			  sin(rot(2)))));
}


///////////////////////////////////////////////////////////////////////////////
/// Converts an AngleAxis into three angles
///
/// If the rotation angles are (theta, psi, phi) and the AngleAxis is
/// (a,v0,v1,v2) then
///
/// a = theta
/// v0 = cos(psi)cos(phi)
/// v1 = sin(psi)cos(phi)
/// v2 = sin(phi)
///
///////////////////////////////////////////////////////////////////////////////
template<int M>
Eigen::Vector3d BundleAdjuster<M>::angleaxis_to_rotation(const Eigen::AngleAxisd &A){
  return(Eigen::Vector3d(A.angle(), 
			 std::atan(A.axis()(1)/A.axis()(0)),
			 std::atan(A.axis()(2)/sqrt(A.axis()(0)*A.axis()(0) + 
						    A.axis()(1)*A.axis()(1)))
			 )
	 );
}


///////////////////////////////////////////////////////////////////////////////
/// returns the principal point of the v'th view in 'params'.
///////////////////////////////////////////////////////////////////////////////
template<int M>
inline typename BundleAdjuster<M>::Block2
BundleAdjuster<M>::principal_pt(int v, Eigen::VectorXd &params) {
  if(v < fixedIntrinsics.size()) 
    throw("Tried to vary principal point of a fixed intrinsic in BundleAdjuster");
  return(params.middleRows<2>(6*(M-1) + 3*(v-fixedIntrinsics.size())));
}

template<int M>
inline typename BundleAdjuster<M>::constBlock2
BundleAdjuster<M>::principal_pt(int v, const Eigen::VectorXd &params) {
  if(v < fixedIntrinsics.size()) 
    throw("Tried to read principal point of a fixed intrinsic in BundleAdjuster");
  return(params.middleRows<2>(6*(M-1) + 3*(v-fixedIntrinsics.size())));
}


///////////////////////////////////////////////////////////////////////////////
/// returns the focal length of the v'th view in 'params'.
///////////////////////////////////////////////////////////////////////////////
template<int M>
inline double &
BundleAdjuster<M>::focal_len(int v, Eigen::VectorXd &params) {
  if(v < fixedIntrinsics.size()) 
    throw("Tried to vary focal length of a fixed intrinsic in BundleAdjuster");
  return(params(6*(M-1) + 3*(v-fixedIntrinsics.size()) + 2));
}
template<int M>
inline const double &
BundleAdjuster<M>::focal_len(int v, const Eigen::VectorXd &params) {
  if(v < fixedIntrinsics.size()) 
    throw("Tried to read focal length of a fixed intrinsic in BundleAdjuster");
  return(params(6*(M-1) + 3*(v-fixedIntrinsics.size()) + 2));
}


///////////////////////////////////////////////////////////////////////////////
/// returns the translation of the v'th view in 'params'.
///////////////////////////////////////////////////////////////////////////////
template<int M>
inline typename BundleAdjuster<M>::Block3
BundleAdjuster<M>::translation(int v, Eigen::VectorXd &params) {
  return(params.middleRows<3>(6*(v-1)));
}

template<int M>
inline typename BundleAdjuster<M>::constBlock3 
BundleAdjuster<M>::translation(int v, const Eigen::VectorXd &params) {
  return(params.middleRows<3>(6*(v-1)));
}


///////////////////////////////////////////////////////////////////////////////
/// returns the rotation of the v'th view in 'params'.
///////////////////////////////////////////////////////////////////////////////
template<int M>
inline typename BundleAdjuster<M>::Block3
BundleAdjuster<M>::rotation(int v, Eigen::VectorXd &params) {
  return(params.middleRows<3>(6*(v-1)+3));
}

template<int M>
inline typename BundleAdjuster<M>::constBlock3
BundleAdjuster<M>::rotation(int v, const Eigen::VectorXd &params) {
  return(params.middleRows<3>(6*(v-1)+3));
}


///////////////////////////////////////////////////////////////////////////////
// copies the parameters in 'params' to the motion matrix
///////////////////////////////////////////////////////////////////////////////
template<int M>
void BundleAdjuster<M>::params_to_matrices(const Eigen::VectorXd &params) {
  int v;

  for(v = 0; v < M; ++v) {
    if(v < fixedIntrinsics.size()) {
      motion.M(v) = fixedIntrinsics[v].matrix();
    } else {
      motion.M(v) = (Eigen::Translation2d(principal_pt(v,params)) *
	            Eigen::Scaling(focal_len(v,params))).matrix();
    }
    if(v > 0) {
      motion.M(v) *= rotation_to_angleaxis(rotation(v,params)).toRotationMatrix();
      motion.T(v) = translation(v,params);
    }
  }
  motion.T(0).setZero();
}


///////////////////////////////////////////////////////////////////////////////
/// Copies the motion matrix to the internal parameters. Assumes T(0)
/// is at the origin.
///////////////////////////////////////////////////////////////////////////////
template<int M>
void BundleAdjuster<M>::matrices_to_params() {
  int			v;
  Eigen::Affine2d 	K; // intrinsics
  Eigen::AngleAxisd	R; // camera rotation
  double		translation_norm;
  double		Mscale;

  translation_norm = 0.0;
  for(v = 0; v<M; ++v) {
    Mscale = motion.KR_decompose(v,K,R);
    if(v < fixedIntrinsics.size()) {
      fixedIntrinsics[v] = K;
    } else {
      principal_pt(v, variableParams) = K.translation();
      focal_len(v, variableParams)    = (K(0,0) + K(1,1))*0.5;
    }
    if(v > 0) {
      rotation(v, variableParams)    = angleaxis_to_rotation(R);
      translation(v, variableParams) = motion.T(v) / Mscale;
      translation_norm += translation(v, variableParams).squaredNorm();
    }
  }

  // --- normalise scale of translations
  translation_norm = sqrt(translation_norm);
  for(v = 1; v<M; ++v) {
    // translation(v, variableParams) -= motion.T(0);
    translation(v, variableParams) *= M/translation_norm;
  }
}


///////////////////////////////////////////////////////////////////////////////
/// Uses the Levenberg-Marquardt algorithm to solve for the
/// measurement matrix supplied at construction, holding the
/// intrinsics of the first 'nFixedIntrinsics' views fixed. The result
/// is placed into the motion matrix supplied at construction.
///
/// \returns The exit status of the LevMar solver: 
///
/// ImproperInputParameters = 0,
/// RelativeReductionTooSmall = 1,
/// RelativeErrorTooSmall = 2,
/// RelativeErrorAndReductionTooSmall = 3,
/// CosinusTooSmall = 4,
/// TooManyFunctionEvaluation = 5,
/// FtolTooSmall = 6,
/// XtolTooSmall = 7,
/// GtolTooSmall = 8
///
///////////////////////////////////////////////////////////////////////////////
template<int M>
int BundleAdjuster<M>::levmar_solve(int nFixedIntrinsics) {
  int info;
  int nPixels;
  int i,j;

  // --- count number of un-occluded pixels
  nPixels = 0;
  for(i = 0; i<M; ++i) {
    for(j = 0; j<measurement.cols(); ++j) {
      if(!measurement.pixel_is_occluded(i,j)) {
	++nPixels;
      }
    }
  }
  
  // --- setup variables for solve
  fixedIntrinsics.resize(nFixedIntrinsics);
  variableParams.resize(6*(M-1) + 3*(M-nFixedIntrinsics));
  matrices_to_params();

  Eigen::VectorXd tmperr;
  params_to_matrices(variableParams);
  motion.reprojection_err(measurement, tmperr);
  std::cout << "Parameterised reprojection err = " << std::endl;
  std::cout << tmperr.norm()/(2.0*nPixels) << std::endl;

  // --- do solve
  std::cout << "Starting solve with params = " << std::endl;
  std::cout << variableParams << std::endl;
  info = levmar_solve(variableParams, 2*nPixels + 1);
  params_to_matrices(variableParams);
  std::cout << "Finished levmar with status = " << info << std::endl;

  motion.reprojection_err(measurement, tmperr);
  std::cout << "Final reprojection err = " << std::endl;
  std::cout << tmperr.norm()/(2.0*nPixels) << std::endl;

  return(info);
}

