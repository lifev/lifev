#include "refFE.hpp"
# include <set> 


RefFE::RefFE(string _name,int _type,ReferenceShapes _shape,int _nbDofPerVertex,
	     int _nbDofPerEdge,int _nbDofPerFace,int _nbDofPerVolume,
	     int _nbDof,int _nbCoor,const Fct* phi,const Fct* dPhi,
	     const Fct* d2Phi,const Real* _refCoor,
	     const SetOfQuadRule& sqr,PatternType _patternType,const RefFE* bdRefFE):
  RefEle(_name,_shape,_nbDof,_nbCoor,phi,dPhi,d2Phi,_refCoor,sqr),
  LocalDofPattern(_nbDof,_nbDofPerVertex,_nbDofPerEdge,_nbDofPerFace,_nbDofPerVolume,_patternType),
 _boundaryFE(bdRefFE), type(_type)
{
  CONSTRUCTOR("RefFE");
}

RefFE::~RefFE()
{
  DESTRUCTOR("RefFE");
}

ostream& operator << (ostream& f,const RefFE& fe)
{
  f << "-------------------------\n";
  f << "Reference Finite Element: " << fe.name << endl;
  f << "-------------------------\n";
  f << "*** Shape : " << fe.shape << endl;
  f << "*** Local coordinate of the nodes :\n";
  for(int i=0;i<fe.nbDof;i++){
    f << fe.xi(i) << " " << fe.eta(i) << " " << fe.zeta(i) << endl;
  }
  f << "*** Pattern :\n";
  for(int i=0;i<fe.nbPattern();i++) f << "(" << fe.patternFirst(i) << "," << fe.patternSecond(i) << ") \n";
  for(int k=0;k<fe._sqr->nbQuadRule;k++){
    const QuadRule& qr = fe._sqr->quadRule(k);
    f << "\n*** Quadrature rule : " << qr.name << endl;
    for(int ig=0;ig<qr.nbQuadPt;ig++){
      f << "    - Quadrature point : " << ig << endl;
      //      f << "     number and values of basis functions = " << fe.phiQuadPt(ig,qr) << endl;
      for(int i=0;i<fe.nbDof;i++){
	f << "      Basif fct " << i << endl;
	f << "         Value = " << fe.phi(i,ig,qr) << endl;
	f << "         Derivatives = " ;
	for(int icoor=0;icoor<fe.nbCoor;icoor++) f << " " << fe.dPhi(i,icoor,ig,qr);
	f << endl;
	f << "         Second derivatives = " ;
	for(int icoor=0;icoor<fe.nbCoor;icoor++)
	  for(int jcoor=0;jcoor<fe.nbCoor;jcoor++)
	  f << " " << fe.d2Phi(i,icoor,jcoor,ig,qr);
	f << endl;
      }
    }
  }
  return f;
}


