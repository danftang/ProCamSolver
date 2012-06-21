///////////////////////////////////////////////////////////////////////////////
//
// Class for synthesising fundamental matrices and pixel correlations
// from scenes with multiple cameras
//
///////////////////////////////////////////////////////////////////////////////
#ifndef SCENESYNTH_H
#define SCENESYNTH_H

#include <vector>
#include <map>
#include "ViewRelation.h"
#include "ViewProjectionMatrix.h"
#include "Correspondence.h"
#include "CorrespondenceSet.h"

class SceneSynth {
public:
  SceneSynth(std::map<int, ViewProjectionMatrix> &);

  void			generate_correspondences(CorrespondenceSet &,
					std::map<int, ViewProjectionMatrix> &,
					int);
  double	noise;
};


#endif
