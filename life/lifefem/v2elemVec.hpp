/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#ifndef V2ELEMVEC_H
#define V2ELEMVEC_H
#include "lifeV.hpp"
#include "elemVec.hpp"
#include "currentFE.hpp"

namespace LifeV
{
//======================================================================
template<typename DOF, typename Vector, typename ElemVec>
void
vec2elemVec(Vector& V,ElemVec& elvec,const CurrentFE& fe, const DOF& dof,
	    int iblock=0,int nb=1)
{
  UInt totdof = dof.numTotalDof();
  int i;
  UInt ig;
  UInt eleId = fe.currentId();
  for(int ib = iblock ; ib <iblock+nb ; ib++){
    Tab1dView vec=elvec.block(ib);
    for(i=0 ; i<fe.nbNode ; i++){
      ig = dof.localToGlobal(eleId,i+1) - 1+ib*totdof;
      vec(i) = V[ig];
    }
  }
};

//======================================================================
/*
  obsolete ?
  template<typename Vector,typename FE1,typename GeoMap,typename DOF>
  void assemble(Vector& V,ElemVec& elvec,const FE1& fe1,const GeoMap& geo,
  const DOF& dof)
  {
  if(elvec.nBlockRow()!=1){
  cout << "assemble for vector elem vec not yet implemented\n";
  exit(1);
  }
  Tab1dView vec=elvec.block(0);
  UInt i;
  UInt ig,eleID=geo.currentID();
  for(i=0 ; i<FE1::nbNode ; i++){
  ig = dof.localToGlobal(eleID,i+1) - 1;
  V(ig) += vec(i);
  }
  };
*/
}
#endif

