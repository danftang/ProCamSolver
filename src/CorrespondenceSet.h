#include "stdincludes.h"
#include "Correspondence.h"

///////////////////////////////////////////////////////////////////////////////
/// Class to represent a set of 
//////////////////////////////////////////////////////////////////////////////
class CorrespondenceSet : public std::vector<Correspondence> {
public:
  class RawCorrespondence {
  public:
    uint32_t i;
    char d;
    double xi;
    double yi;
    uint32_t j;
    char k;
    double xj;
    double yj;
    double w;
  };

  void 	load(const char *);
  void	synthesize_correspondences(std::map<int,ViewProjectionMatrix> &, int, double);

};

inline void CorrespondenceSet::load(const char *filename) {
  std::ifstream myFile(filename, std::ios::binary);
  int i;
  Correspondence c;
  RawCorrespondence rc;

  myFile.read((char *)&i,sizeof(uint32_t));
  i = 10000;
  while(i >0 && myFile) {    
    myFile.read((char *)&rc, sizeof(RawCorrespondence));
    c.i = rc.i;
    c.xi = rc.xi;
    c.yi = rc.yi;
    c.j = rc.j;
    c.xj = rc.xj;
    c.yj = rc.yj;
    c.w  = rc.w;
    push_back(c);
    //    std::cout << c << std::endl;
    --i;
  }
  myFile.close();
}
