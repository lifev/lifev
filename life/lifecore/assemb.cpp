#include "assemb.hpp"

//
void compute_vec(Real constant,ElemVec& elvec,const CurrentFE& fe,int iblock){
  int i,ig;
  Tab1dView vec = elvec.block(iblock);
  Real s;
  for(i=0;i<fe.nbNode;i++){
    s = 0;
    for(ig=0;ig<fe.nbQuadPt;ig++) s += fe.phi(i,ig)*fe.weightDet(ig);
    vec(i) += constant*s;    
  }
}
