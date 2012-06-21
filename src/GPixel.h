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
#ifndef GPIXEL_H
#define GPIXEL_H

///////////////////////////////////////////////////////////////////////////////
/// Class for global identification of pixels accross a number of views.
/// A GPixel identifies the view and the pixel within that view.
///////////////////////////////////////////////////////////////////////////////
class GPixel {
public:
  bool operator <(const GPixel &) const;
  void rescale(const double &, const double &);

  int	 id;	// first camera's ID
  double x;	// (x,y) of pixel in first view
  double y;
};


///////////////////////////////////////////////////////////////////////////////
/// rescales this pixel by 'sx' and 'sy'
///////////////////////////////////////////////////////////////////////////////
inline void GPixel::rescale(const double &sx, const double &sy) {
  x *= sx;
  y *= sy;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline bool GPixel::operator <(const GPixel &o) const {
  if(id < o.id) return(true);
  if(id > o.id) return(false);
  if(x < o.x) return(true);
  if(x > o.x) return(false);
  if(y < o.y) return(true);
  return(false);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline std::ostream &operator <<(std::ostream &out, const GPixel &c) {
  out << "(" << c.id << " " 
      << c.x << " "
      << c.y << ")";
  return(out);
}

#endif
