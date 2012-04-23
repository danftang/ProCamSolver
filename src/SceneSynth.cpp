#include <cstdlib>
#include "SceneSynth.h"


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SceneSynth::SceneSynth(std::map<int, ViewProjectionMatrix> &cameras) {
  std::map<int, ViewProjectionMatrix>::iterator vi, vj;
  ViewRelation	view;

  for(vj = cameras.begin(); vj != cameras.end(); ++vj) {
    for(vi = cameras.begin(); vi != vj; ++vi) {
      view.i = vi->first;
      view.j = vj->first;
      view.F = F(vi->second, vj->second);
      push_back(view);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
Matrix<double> SceneSynth::F(ViewProjectionMatrix &pi, 
			     ViewProjectionMatrix &pj) {
  return(pj.M1T() * pj.R() * (pj.T() - pi.T()) * pi.R1() * pi.M1());
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void SceneSynth::generate_correspondences(std::vector<Correspondence> &cor,
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
	C.xi = qi[0]/qi[2];
	C.yi = qi[1]/qi[2];
	C.xj = qj[0]/qj[2];
	C.yj = qj[1]/qj[2];
	//std::cout << "Creating correspondence " << C << std::endl;
	//std::cout << "qj F qi = " << (qj*(F(vi->second, vj->second)*qi)).sum() << std::endl;
	cor.push_back(C);
      }
    }
  }
}
