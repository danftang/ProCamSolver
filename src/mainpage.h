/// \mainpage ProCamSolver
/// \author Dan Tang
///
/// This is a C++ library for reconstructing a 3D scene from a number of 2D images of the scene. Using the library, the 3D coordinates of points in the scene can be calculated as well as the position/orientation/intrinsics of each camera. It is assumed that a set of pixel correspondences has already been calculated (that is, we have pairs of pixels - one from each of two different images - that correspond to the same 3D point). This is useful for many applications, for example structure from motion, augmented reality and projection mapping.
///
///The library consists of three main classes: ShapeMatrix, MotionMatrix and MeasurementMatrix.
///
///A MeasurementMatrix can be thought of a container of pixel correspondences (i.e. pixels in different images that correspond to the same 3D point). Each column of this matrix corresponds to a single 3D point for which we have corresponding image points, and each block of 3 rows corresponds to a single viewpoint. So, suppose we have information about N 3D points seen from M viewpoints. If the \f$j^{th}\f$ 3D point is seen from the \f$i^{th}\f$ image as pixel \f$(x_{ij},y{ij},1)^T\f$, then the MeasurementMatrix would look like:
///\f[
///E = \left(
///\begin{array}{cccc}
///x_{00} 	& x_{01} & \cdots & x_{0N}\\
///y_{00} 	& y_{01} & \cdots & y_{0N}\\
///1	& 1      & \cdots & 1     \\
///\vdots  & \ddots &        & \vdots \\
///\vdots  &        & \ddots & \vdots \\
///x_{M0} 	& x_{M1} & \cdots & x_{MN}\\
///y_{M0} 	& y_{M1} & \cdots & y_{MN}\\
///1	& 1      & \cdots & 1     \\
///\end{array}
///\right)
///\f]
///
/// In general, a 3D point will only be visible from some of the viewpoints, other viewpoints being occluded by objects placed between the 3D point and the camera. When this is the case, occluded pixels are given the value \f$(0,0,0)^T\f$.
///
/// A MotionMatrix is just the camera matrices from each view placed one on top of the other:
///\f[
///P = \left(
///\begin{array}{c}
///P_0 \\
///P_1 \\
///\vdots \\
///P_M
///\end{array}
///\right)
///\f]
///
/// A ShapeMatrix contains the reconstructed 3D coordinates of the scene, one coordinate per column, so if the \f$n^{th}\f$ 3D point is \f$(x_n, y_n, z_n, 1)^T\f$ then the corresponding shape matrix is:
///\f[
///S = \left(
///\begin{array}{ccc}
///x_0 & \cdots & x_N\\
///y_0 & \cdots & y_N\\
///z_0 & \cdots & z_N\\
///1   & \cdots & 1  \\
///\end{array}
///\right)
///\f]
///
///The best way to get correspondences into a MeasurementMatrix is to construct it with a CorrespondenceSet. Alternatively use the load, add_correspondence, or pixel methods. Once constructed, the motion matrix can be calculated by constructing an empty MotionMatrix and using the svd_solve method. If extra accuracy is required, follow this with a call to the bundle_adjust method to find the least squares minimum reprojection error.
///
///Given a measurement matrix and a motion matrix, the shape matrix can easily be calculated by constructing an empty ShapeMatrix and calling the solve method.
///
