#include "mesh_util_base.hpp"
GetCoordComponent::GetCoordComponent() :comp(-1) {}

GetCoordComponent::GetCoordComponent(int i):comp(i) {}

void GetCoordComponent::operator()(Real const x, Real const y, Real const z, Real ret[3])const
{
  switch(comp){
  case(0):
    ret[0]=x;
    ret[1]=0.0;
    ret[2]=0.0;
    break;
  case(1):
    ret[0]=0.0;
    ret[1]=y;
    ret[2]=0.0;
    break;
  case(2):
    ret[0]=0.0;
    ret[1]=0.0;
    ret[2]=z;
    break;
  default:
    ret[0]=x;
    ret[1]=y;
    ret[2]=z;
  }
}

void GetOnes::operator()(Real const x, Real const y, Real const z, Real ret[3])
  const{
  ret[0]=1.0;
  ret[1]=1.0;
  ret[2]=1.0;
}
