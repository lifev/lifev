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
#ifndef ELEMVEC_H
#define ELEMVEC_H

#include <vector>

#include "lifeV.hpp"
#include <iomanip> 
#include "tab.hpp"
#include "currentFE.hpp"

class ElemVec
    :
    public Tab1d
{
    int _nBlockRow; // number of block rows
    std::vector<int> _nRow; // _nRow[i]=nb of rows in the i-th block row 
    std::vector<int> _firstRow;//_firstRow[i]=index of first row of i-th block row
public:
    typedef Tab1d super;

    ElemVec(int nNode1,int nbr1);
    ElemVec(int nNode1,int nbr1,
	    int nNode2,int nbr2);
    ElemVec(int nNode1,int nbr1,
	    int nNode2,int nbr2,
	    int nNode3,int nbr3);

    ElemVec& operator=( super const& __v )
	{
	    if ( this == &__v )
		return *this;
	    super& __super = (super&)*this;
	    __super = super(*this);
	    return *this;
	}

    Tab1d& vec(){return *this;};
    const Tab1d& vec() const {return *this;};
    int nBlockRow()const{return _nBlockRow;}
    Tab1dView block(int i)
	{
	    return (*this)(SubArray(_nRow[i],_firstRow[i]));
	}
    //  inline void zero(){_vec=Tab1d(_nBlockRow*_vec.N(),0.0);};
    inline void zero() 
	{ 
	    super& __super = (super&)*this;
	    __super = super(this->N(),0.0);
	}
    void showMe(ostream& c=cout);
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

#endif

