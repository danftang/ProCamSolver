#include <cstdlib>
#include "SceneSynth.h"


///////////////////////////////////////////////////////////////////////////////
///
/// Constructs a scene synthesiser containing a set of cameras with the
/// given CameraMatrices
///
/// \param cameras set of cameras to install in this synthetic scene.
///
///////////////////////////////////////////////////////////////////////////////
SceneSynth::SceneSynth(std::map<int, ViewProjectionMatrix> &cameras, double n) {
  std::map<int, ViewProjectionMatrix>::iterator vi, vj;
  ViewRelation	view;

  noise = n;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void SceneSynth::generate_correspondences(CorrespondenceSet &cor,
				    std::map<int,ViewProjectionMatrix> &cams,
				    int N) {
  const double  sceneOffset = 2.0;
  int 		i;
  coord		Q(3);		// 3D coord
  coord		qi(3);		// pixel in camera i
  coord		qj(3);		// pixel in camera j
  Correspondence C;
  std::map<int, ViewProjectionMatrix>::iterator vi, vj;
  
  for(i = 0; i<N; ++i) {
    Q[0] = rand()*1.0/RAND_MAX - 0.5;
    Q[1] = rand()*1.0/RAND_MAX - 0.5;
    Q[2] = rand()*1.0/RAND_MAX - 0.5 + sceneOffset;
    for(vj = cams.begin(); vj != cams.end(); ++vj) {
      for(vi = cams.begin(); vi != vj; ++vi) {
	qi = vi->second.M()*(vi->second.R()*(Q - vi->second.t));
	qj = vj->second.M()*(vj->second.R()*(Q - vj->second.t));
	C.i = vi->first;
	C.j = vj->first;
	C.xi = (qi[0]/qi[2]) + (rand()*(2.0*noise)/RAND_MAX) - noise;
	C.yi = (qi[1]/qi[2]) + (rand()*(2.0*noise)/RAND_MAX) - noise;
	C.xj = (qj[0]/qj[2]) + (rand()*(2.0*noise)/RAND_MAX) - noise;
	C.yj = (qj[1]/qj[2]) + (rand()*(2.0*noise)/RAND_MAX) - noise;
	cor.push_back(C);
      }
    }
  }
}
