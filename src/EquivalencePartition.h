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
#ifndef EQUIVALENCEPARTITION_H
#define EQUIVALENCEPARTITION_H

///////////////////////////////////////////////////////////////////////////////
///
/// A class for representing a set of objects that are partitioned into
/// equivalence classes, such that all members of a partition are equivalent
///
/// The class enforces transitivity, i.e. A=B & B=C implies A=C
///
///////////////////////////////////////////////////////////////////////////////
template<class T>
class EquivalencePartition : public std::vector<std::set<T> > {
public:
  typedef typename std::set<T>		 	EquivalenceClass;
  using   std::vector<std::set<T> >::size;
  using   std::vector<std::set<T> >::begin;
  
  int	equate(const T &, const T &);	// equate two objects
  int	equivalenceClass(const T &) const; // find the equiv. class of an obj
  int	newClass();			// create a new equivalence class
  bool	insert(int, const T &);		// insert obj into given class
  int	merge(int, int);		// merge two equivalence classes
};

///////////////////////////////////////////////////////////////////////////////
///
/// returns the ID of the equivalence class that x belongs to,
/// or -1 if x is not in this partitioning.
///
///////////////////////////////////////////////////////////////////////////////
template<class T>
int EquivalencePartition<T>::equivalenceClass(const T &x) const {
  int i;

  for(i = 0; i<size(); ++i) {
    if((*this)[i].find(x) != (*this)[i].end()) {
      return(i);
    }
  }
  return(-1);
}


///////////////////////////////////////////////////////////////////////////////
///
/// create a new equaivalence class and return its id.
///
///////////////////////////////////////////////////////////////////////////////
template<class T>
int EquivalencePartition<T>::newClass() {
  std::set<T> emptySet;
  push_back(emptySet);
  return(size()-1);
}


///////////////////////////////////////////////////////////////////////////////
///
/// insert x into equivalence class i. Returns true if x is new.
///
///////////////////////////////////////////////////////////////////////////////
template<class T>
bool EquivalencePartition<T>::insert(int i, const T &x) {
  return(((*this)[i].insert(x)).second);
}


///////////////////////////////////////////////////////////////////////////////
///
/// Equates x and y. Both may or may not already exist in the partitioning.
/// returns the id of the class that they are in.
///
///////////////////////////////////////////////////////////////////////////////
template<class T>
int EquivalencePartition<T>::equate(const T &x, const T &y) {
  int xclass = equivalenceClass(x);
  int yclass = equivalenceClass(y);

  if(xclass < 0) {
    if(yclass < 0) {			// neither exist yet
      xclass = newClass();
      insert(xclass, x);
      insert(xclass, y);      
    } else {				// y exists, x doesn't
      insert(yclass, x);
      return(yclass);
    }
  } else {
    if(yclass < 0) {			// x exists, y doesn't
      insert(xclass, y);
    } else {				// both exist
      xclass = merge(xclass, yclass);
    }
  }
  return(xclass);
}


///////////////////////////////////////////////////////////////////////////////
///
/// merges two equivqalence classes, returning the id of the merged class
///
///////////////////////////////////////////////////////////////////////////////
template<class T>
int EquivalencePartition<T>::merge(int i, int j) {
  if(i > j) {
    int t;
    t = i; i = j; j = t;
  }
  if(i == j) return(i);

  (*this)[i].insert((*this)[j].begin(), (*this)[j].end());
  erase(begin() + j);
  return(i);
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
template<class T>
std::ostream &operator <<(std::ostream &out, const EquivalencePartition<T> &E){
  typename EquivalencePartition<T>::const_iterator it;
  typename EquivalencePartition<T>::EquivalenceClass::const_iterator cit;

  for(it = E.begin(); it != E.end(); ++it) {
    out << '{';
    for(cit = it->begin(); cit != it->end(); ++cit) {
      out << " " << *cit; 
    }
    out << " }" << std::endl;
  }
  return(out);
}

#endif
