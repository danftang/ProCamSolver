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
#ifndef CORRESPONDENCESET_H
#define CORRESPONDENCESET_H

#include "EquivalencePartition.h"
#include "Correspondence.h"

///////////////////////////////////////////////////////////////////////////////
/// Class to represent a set of correspondences. Create an instance of
/// this class and use add_correspondece to fill it. The instance can
/// then be used to construct a MeasurementMatrix. If you have a large
/// number (>1000) of correspondences in no particular order this is a
/// more efficient way of creating a MeasurementMatrix than using
/// MeasurementMatrix::add_correspondence
///////////////////////////////////////////////////////////////////////////////
class CorrespondenceSet : public EquivalencePartition<GPixel> {
public:
  int	add_correspondence(const Correspondence &);
};

///////////////////////////////////////////////////////////////////////////////
/// Adds a pair of corresponding pixels to this set
///////////////////////////////////////////////////////////////////////////////
inline int CorrespondenceSet::add_correspondence(const Correspondence &corr) {
  return(equate(corr.i, corr.j));
}

#endif
