#ifndef ELEMVEC_H
#define ELEMVEC_H

#include <vector>

#include "lifeV.hpp"
#include <iomanip> 
#include "tab.hpp"
#include "currentFE.hpp"

class ElemVec
{
  Tab1d _vec; // the array
  int _nBlockRow; // number of block rows
  vector<int> _nRow; // _nRow[i]=nb of rows in the i-th block row 
  vector<int> _firstRow;//_firstRow[i]=index of first row of i-th block row
public:
  ElemVec(int nNode1,int nbr1);
  ElemVec(int nNode1,int nbr1,
	  int nNode2,int nbr2);
  ElemVec(int nNode1,int nbr1,
	  int nNode2,int nbr2,
	  int nNode3,int nbr3);
  inline Tab1d& vec(){return _vec;};
  inline const Tab1d& vec() const {return _vec;};
  inline int nBlockRow()const{return _nBlockRow;}
  inline Tab1dView block(int i){
    return _vec(SubArray(_nRow[i],_firstRow[i]));
  }
  //  inline void zero(){_vec=Tab1d(_nBlockRow*_vec.N(),0.0);};
  inline void zero(){_vec=Tab1d(_vec.N(),0.0);};//Alex 11/01/02: _vec.N() accounts for the number of block, doesn't it ?
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

// $Id: elemVec.hpp,v 1.1 2004-02-08 09:09:24 prudhomm Exp $
