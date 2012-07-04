#ProCamSolver 0.1

This is a C++ library for reconstructing a 3D scene from a number of 2D images of the scene. Using the library, the 3D coordinates of points in the scene can be calculated as well as the position/orientation/intrinsics of each camera. It is assumed that a set of pixel correspondences has already been calculated (that is, we have pairs of pixels - one from each of two different images - that correspond to the same 3D point).

If there are enough 2D images of the same 3D scene then a set of pixel correspondences define a 3D reconstruction of the scene (the 'shape') and the position, orientation and intrinsics of each camera (the 'motion') up to a projective transofrmation. This library provides various methods of finding these and a convenient platform for developing new methods. It also provides methods of finding the projective transofrmation that leads to a Euclidean (zero skew) set of cameras.

**API documentation can be found [here](http://danftang.github.com/ProCamSolver/).**

_The code is pretty much working at the moment but **isn't finished**, I'll be tidying it up, tesing it, adding more documentation and more functionality over the next few weeks_

# Introduction

## Quick mathematical background

For the mathematically minded, the problem can be expressed most succinctly in terms of three matrices: The 'shape matrix', the 'motion matrix' and the 'measurement matrix'. The shape matrix represents the 3D coordinates of the points in the scene in a 4xN matrix whose columns are of the form (x_n,y_n,z_n,1)^T, where N is the number of 3D points. The motion matrix consists of the camera projection matrices of each view placed one on top of the other in a 3Mx4 matrix, where M is the number of views. The measurement matrix is a representation of the pixel correspondences in a 3MxN matrix where each column represents to the correspondences of one 3D point, with the 2D pixel coordinates of each view placed one on top of the other. The pixel coordinates are of the form s(x,y,1)^T where s is a scalar called the 'projective depth'.

In this form, these matrices have the simple relationship:

E = P * S

where

E = measurement matrix

P = motion matrix

S = shape matrix

The difficulty is that we don't know P, S or the projective depths of E. However, because P and S have only 4 columns/rows respectively, we know that E must be a rank 4 matrix, and this is often enough to give a unique solution.up to a projective transformation.

## Overview of methodology

1. Use the __Measurement Matrix__ to find minimum __Spanning Tree__ of views and create __Fundamental Matrices__ for branches of the tree. Use Martinec Pajda's SVD method to approximate the __Motion Matrix__.
2. Use the __Motion Matrix__ from step 1 to fill in occlusions in the __Measurement Matrix__. Perform standard SVD on the __Measurement Matrix__ to get the second level approximation of the __Motion Matrix__.
2b. _Euclidean Lift_ the __Motion Matrix__.
3. Use _Levenberg-Marquardt_ to minimise the __Reprojection Error__ which gives the final __Motion Matrix__.
4. Use _KR Decomposition_ to recover __Intrinsics__ and __Extrinsics__.

## Usage

1. Fill the _Measurement Matrix_ (e.g. with the `CorrespenceSet` class, the first set of indices are for cameras, whilst remaining indices are projectors).
2. Use `MotionMatrix::svd_solve(const MeasurementMatrix &)` to gain the first level approximation by minimising the _Frobenius Norm_. This step performs 4 things (as described [here](http://danftang.github.com/ProCamSolver/classMotionMatrix.html#5a07d5b1459cb42e191a3892c9c122ea)). This stage gives you a solution with an arbitrary projective transformation.
3. Use `MotionMatrix::diac_euclidean_lift(int N)`, where N is the number of views which are cameras (i.e. not projectors), this transforms our matrices in the _Motion Matrix_ into Euclidean transforms. Optionally this can also be performed with `MotionMatrix::diac_euclidean_lift(int N, ImageTransform<F> &)` where we input F known intrinsics (i.e. for F cameras).
4. Perform Bundle Adjustment using `MotionMatrix::bundle_adjust(int f, const MeasurementMatrix &)`, this is a levmar fit to minimise _Reprojection Error_ where f is the number of cameras. __N.B. code also available to perform this stage with known intrinsics__.
5. (Optional) To create a 3D scan, create a `ShapeMatrix` and use `ShapeMatrix::solve(MeasurementMatrix &, MotionMatrix &)`.
6. Use `MotionMatrix::KR_decompose(ImageTransform<V>& K, ImageTransform<V>& R)` to get __Intrinsics__ (`K`) and __Extrinsics__ (`R`) where `V` is view count, and use `MotionMatrix::T(i)` to get the translation in projection space, use `K.view(i).inverse() * R.view(i).inverse() * MotionMatrix::T(i)` to recover the world translation of the view.

# Further work

## RANSAC implementation

Implement RANSAC into FundamentalMatrix class which will solve against subsets of the MeasurementMatrix, and then mask the MeasurementMatrix (by deletion of datapoints) for further calculation downstream.

RANSAC is best performed on step 1 (Fundamental Matrix) of [here](http://danftang.github.com/ProCamSolver/classMotionMatrix.html#5a07d5b1459cb42e191a3892c9c122ea). Svoboda et al (A Convenient Multi-Camera Self-Calibration for Virtual Environment) also perform this on the other steps

## Better sparse matrix 

With large datasets, we should be selecting a good subset to work with which is an involved process. A strategy would need to be made to find the best way of summarising the data without losing important features.

## Standalone app

Since we template against number of views, we could create version of the classes for 3-10 views for example (3 views is the minimum requirement).