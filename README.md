##ProCamSolver 0.1

This is a C++ library for reconstructing a 3D scene from a number of 2D images of the scene and calculating the position/orientation/intrinsics of each image. It is assumed that a set of pixel correspondences has already been calculated (that is, we have pairs of pixels, one from each of two different images, that point to the same 3D point - although the coordinates of the 3D point and the position/orientation/intrinsics of the images need not be known).

If there are enough 2D images of the same 3D scene then a set of pixel correspondences can be split into a 3D reconstruction of the scene (the 'shape') and the position, orientation and intrinsics of each image (the 'motion'). This library provides various methods of doing this and a convenient platform for developing new methods.

_The code is pretty much working at the moment but *isn't finished*, I'll be tidying it up, tesing it, adding more documentation and more functionality over the next few weeks_

# Mathematical background

For the mathematically minded, the problem can be expressed most succinctly in terms of three matrices: The 'shape matrix', the 'motion matrix' and the 'measurement matrix'. The shape matrix represents the 3D coordinates of the points in the scene in a 4xN matrix whose columns are of the form (x_n,y_n,z_n,1)^T, where N is the number of 3D points. The motion matrix consists of the camera projection matrices of each view placed one on top of the other in a 3Mx4 matrix, where M is the number of views. The measurement matrix is a representation of the pixel correspondences in a 3MxN matrix where each column corresponds to the correspondences of one 3D point, with the 2D pixel coordinates of each view placed one on top of the other. The pixel coordinates are of the form s(x,y,1)^T where s is a scalar called the 'projective depth'.

In this form, these matrices have the simple relationship:

E = PX

where

E = measurement matrix
P = motion matrix
X = shape matrix
