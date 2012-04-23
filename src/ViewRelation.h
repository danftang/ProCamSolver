///////////////////////////////////////////////////////////////////////////////
//
// Class to represent the fundamental matrix between two views so that
// is p_i is a pixel from view-i and p_j is a pixel from view-j then
//
// p_j^T F p_i = 0
//
///////////////////////////////////////////////////////////////////////////////
#ifndef VIEWRELATION_H
#define VIEWRELATION_H

#include "transformMatrix.h"

class ViewRelation {
public:
  transformMatrix	F;	// fundamental matrix between views
  int			i;	// ID of the first view
  int			j;	// ID of the second view
};

#endif
