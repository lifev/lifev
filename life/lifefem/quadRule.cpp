#include "quadRule.hpp"

QuadRule::QuadRule(const QuadPoint* pt,int _id,string _name,
		   ReferenceShapes _shape,int _nbQuadPt,int _degOfExact):
  _pt(pt),shape(_shape),id(_id),name(_name),
  nbQuadPt(_nbQuadPt),degOfExact(_degOfExact)
{
  CONSTRUCTOR("QuadRule");
}

QuadRule::~QuadRule(){
  DESTRUCTOR("QuadRule");
}

ostream& operator << (ostream& c,const QuadRule& qr)
{
  c << " name: " << qr.name << endl;
  c << " shape:" << (int)qr.shape << endl;
  c << " id: " << qr.id << endl;
  c << " nbQuadPt: " << qr.nbQuadPt << endl;
  c << " Points: \n";
  for(int i=0;i<qr.nbQuadPt;i++) c << qr._pt[i] << endl;
  return c;
}

SetOfQuadRule::SetOfQuadRule(const QuadRule* qr,int _nb)
  :_qr(qr),nbQuadRule(_nb)
{
  CONSTRUCTOR("SetOfQuadRule");
  _totalNbQuadPoint=0;
  _maxIdQuadRule=0;
  for(int i=0;i<nbQuadRule;i++){
    if(qr[i].shape != qr[0].shape ){
      cout << "Error in File : " << __FILE__ << " Line : " << __LINE__ << endl;
      ERROR_MSG("All quadrature rules of a set of quadrature rules \n should have the same geometric reference shape");
    }
    _totalNbQuadPoint += qr[i].nbQuadPt;
    if (qr[i].id > _maxIdQuadRule) _maxIdQuadRule = qr[i].id;
  }
  _posQuadRule = new int[ _maxIdQuadRule+1 ]; // don't forget the +1 !
  for(int i=0;i<nbQuadRule;i++){
    _posQuadRule[ qr[i].id ] = i;
  }
};

SetOfQuadRule::~SetOfQuadRule()
{
  DESTRUCTOR("SetOfQuadRule");
}

