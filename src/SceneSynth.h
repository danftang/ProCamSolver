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
#include "Matrix.h"
#include "Correspondence.h"

class SceneSynth : public std::vector<ViewRelation> {
public:
  SceneSynth(std::map<int, ViewProjectionMatrix> &);

  Matrix<double>	F(ViewProjectionMatrix &, ViewProjectionMatrix &);
  void			generate_correspondences(std::vector<Correspondence> &,
					std::map<int, ViewProjectionMatrix> &,
					int);
};


#endif
