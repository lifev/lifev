#include "refHybridFE.hpp"
#include <set> 


RefHybridFE::RefHybridFE(const UInt& nbdfe, const StaticBdFE* bdfelist,
			 string _name, int _type, ReferenceShapes _shape,
			 int _nbDofPerVertex,int _nbDofPerEdge,
			 int _nbDofPerFace,int _nbDofPerVolume,
			 int _nbDof,int _nbCoor,const Real* refCoor,
			 PatternType _patternType):
  LocalDofPattern(_nbDof,_nbDofPerVertex,_nbDofPerEdge,_nbDofPerFace,_nbDofPerVolume,_patternType),
  _nBdFE(nbdfe), _bdfeList(bdfelist), _refCoor(refCoor), 
  name(_name), type(_type), shape(_shape),
  nbDof(_nbDof),nbCoor(_nbCoor)
{
  CONSTRUCTOR("RefHybridFE");
  
  //! simple consistency test: (to be removed some day)
  Int nbdofbdfe =0;
  for(UInt nf = 0 ; nf < _nBdFE ;  nf ++) nbdofbdfe += _bdfeList[ nf ].nbNode;
  ASSERT_PRE( nbdofbdfe == nbDof , \
	      "There should be as many dof in the refHybrid element as the sum of all dof in the boundary elements.");
}

RefHybridFE::~RefHybridFE()
{
  DESTRUCTOR("RefHybridFE");
}

//! extracting a BdFE from the boundary elements list. //to be checked
const StaticBdFE& RefHybridFE::operator[](const Index_t& i) const 
{
  ASSERT_BD( i < (Index_t) _nBdFE );
  return _bdfeList[i];
}

void RefHybridFE::check() const 
{
  Real sumphi,sumdphi;
  int nbdofbdfe=0;
  cout << "*** Check " << name << endl;
  for(UInt nf = 0 ; nf < _nBdFE ;  nf ++){
    const StaticBdFE& bdfe = _bdfeList[ nf ];
    const QuadRule& qr = bdfe.qr;
    cout << endl << "    " << bdfe.refFE.name << ", whose identity is " << bdfe.currentId() << endl;
    cout << "    " << qr.name << endl;
    nbdofbdfe += bdfe.nbNode;

    for(int i = 0 ; i <  bdfe.nbNode ; i ++){
      sumphi=sumdphi=0.;      
      for(int ig = 0 ; ig < qr.nbQuadPt ; ig ++){
	sumphi +=  bdfe.phi(i,ig) *  bdfe.weightMeas(ig);
 	for(int icoor = 0 ; icoor < bdfe.nbCoor ; icoor ++) 
	  sumdphi += bdfe.dPhiRef(i,icoor,ig) * bdfe.weightMeas(ig);
      }
    cout << "      integral_Face phi_i                                   = " << sumphi << endl;
    cout << "      sum_{icoor} integral dphi_i/dx_icoor                  = " << sumdphi << endl;
    }
  }
  cout << "     number of dof (should be : "<< nbDof << " )             = " << nbdofbdfe << endl;
  //! simple consistency test:
  ASSERT_PRE( nbdofbdfe == nbDof , \
    "There should be as many dof in the refHybrid element as the sum of all dof in the boundary elements.");
  cout << endl;
}

ostream& operator << (ostream& f,const RefHybridFE& fe)
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
  return f;
}
