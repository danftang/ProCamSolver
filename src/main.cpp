#include <iostream>
#include <cstdlib>

#include "RadialViewProjectionMatrix.h"
#include "SceneSynth.h"
#include "RPSolver.h"
//#include "RPSolverLevMar2.h"
#include "PSolver.h"
#include "CorrespondenceSet.h"

int main() {
  std::map<int, ViewProjectionMatrix>	cameras;
  std::map<int, RadialViewProjectionMatrix>	rcameras;
  std::map<int, RadialViewProjectionMatrix>	myFit;
  //std::map<int, ViewProjectionMatrix>	myFit;
  std::map<int, transformMatrix>    	intrinsics;
  std::vector<Correspondence>		correspondences;
  RadialViewProjectionMatrix		view;
  double				residual;
  int 					i;

  // --- Set up a scene with 3 cameras
  // ---------------------------------
  view.fx = 1.0;
  view.fy = 1.0;
  view.cx = 0.0;
  view.cy = 0.0;
  view.rot[0] = 0.0;
  view.rot[1] = 0.0;
  view.rot[2] = 0.0;
  view.t[0] = 0.0;
  view.t[1] = 0.0;
  view.t[2] = 0.0;
  cameras[0] = view;
  rcameras[0] = view;
  intrinsics[0] = transformMatrix(intrinsic, view.fx,view.fy,view.cx,view.cy);

  view.fx = 1.0;
  view.fy = 1.2;
  view.cx = 0.01;
  view.cy = 0.02;
  view.rot[0] = 0.0;
  view.rot[1] = 0.1;//46364;
  view.rot[2] = 0.0;
  view.t[0] = 1.0;
  view.t[1] = 0.0;
  view.t[2] = 0.0;
  cameras[1] = view;
  rcameras[1] = view;
  intrinsics[1] = transformMatrix(intrinsic, view.fx,view.fy,view.cx,view.cy);

  view.fx = 1.0123;
  view.fy = 1.0345;
  view.cx = 0.0678;
  view.cy = 0.0789;
  view.rot[0] = 0.0;
  view.rot[1] = 0.2;
  view.rot[2] = 0.0;
  view.t[0] = 2.0;
  view.t[1] = 1.0;
  view.t[2] = 0.0;
  cameras[2] = view;
  rcameras[2] = view;
  intrinsics[2] = transformMatrix(intrinsic, view.fx,view.fy,view.cx,view.cy);

  view.fx = 1.0123;
  view.fy = 1.0345;
  view.cx = 0.0678;
  view.cy = 0.0789;
  view.rot[0] = 0.1;
  view.rot[1] = 0.1;
  view.rot[2] = 0.1;
  view.t[0] = 3.0;
  view.t[1] = -1.0;
  view.t[2] = 0.1;
  cameras[3] = view;
  rcameras[3] = view;
  intrinsics[3] = transformMatrix(intrinsic, view.fx,view.fy,view.cx,view.cy);


  view.fx = 1.0123;
  view.fy = 1.0345;
  view.cx = 0.0678;
  view.cy = 0.0789;
  view.rot[0] = 0.0;
  view.rot[1] = 0.0;
  view.rot[2] = 0.0;
  view.t[0] = 4.0;
  view.t[1] = 0.0;
  view.t[2] = 0.0;
  cameras[4] = view;
  rcameras[4] = view;
  intrinsics[4] = transformMatrix(intrinsic, view.fx,view.fy,view.cx,view.cy);

  /****

  view.fx = 1.0123;
  view.fy = 1.0345;
  view.cx = 0.0678;
  view.cy = 0.0789;
  view.rot[0] = 0.0;
  view.rot[1] = 0.0;
  view.rot[2] = 0.0;
  view.t[0] = 5.0;
  view.t[1] = 0.0;
  view.t[2] = 0.0;
  cameras[5] = view;
  ***/

  // --- synthesise a vector of pixel correspondences
  // ------------------------------------------------
  SceneSynth 	myScene(cameras);
  myScene.generate_correspondences(correspondences, cameras, 20);

  // --- now solve for the correspondences,
  // --- given the intrinsics from 2 of the
  // --- cameras, returning the residual
  // --------------------------------------
  try{ 
    std::cout << "Fitting..." << std::endl;

    //CorrespondenceSet myCors;
    //myCors.load("correspondences.dat");
    //myFit = RPSolver(myCors, residual);

    //rcameras[2].t[0] += 0.1;
    //myFit = RPSolver(correspondences, rcameras, residual);
    myFit = RPSolver(correspondences, intrinsics, residual);

  } catch(const char *err) {
    std::cerr << "Caught error " << err << std::endl;
  }
  // --- Print out the results of the fit
  // ------------------------------------
  std::cout << "residual = " << residual << std::endl;
  for(i=0; i<myFit.size(); ++i) {
    std::cout << "Camera " << i << std::endl;
    std::cout << "--------" << std::endl;
    std::cout << myFit[i];
    std::cout << std::endl;
  }

  return(0);
}
