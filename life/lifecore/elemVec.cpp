#include "elemVec.hpp"


ElemVec::ElemVec(int nNode1,int nbr1):
  _vec(nNode1*nbr1)
{
  _nBlockRow = nbr1;
  _nRow.resize(_nBlockRow);
  _firstRow.resize(_nBlockRow);
  int first=0,n;
  for(n=0;n<nbr1;n++){
    _nRow[n] = nNode1;
    _firstRow[n] = first;
    first += nNode1;
  }
}

ElemVec::ElemVec(int nNode1,int nbr1,
		 int nNode2,int nbr2):
  _vec(nNode1*nbr1+nNode2*nbr2)
{
  _nBlockRow = nbr1+nbr2;
  _nRow.resize(_nBlockRow);
  _firstRow.resize(_nBlockRow);
  int first=0,n;
  for(n=0;n<nbr1;n++){
    _nRow[n] = nNode1;
    _firstRow[n] = first;
    first += nNode1;
  }
  for(n=nbr1;n<nbr1+nbr2;n++){
    _nRow[n] = nNode2;
    _firstRow[n] = first;
    first += nNode2;
  }
}

ElemVec::ElemVec(int nNode1,int nbr1,
		 int nNode2,int nbr2,
		 int nNode3,int nbr3):
  _vec(nNode1*nbr1+nNode2*nbr2+nNode3*nbr3)
{
  _nBlockRow = nbr1+nbr2+nbr3;
  _nRow.resize(_nBlockRow);
  _firstRow.resize(_nBlockRow);
  int first=0,n;
  for(n=0;n<nbr1;n++){
    _nRow[n] = nNode1;
    _firstRow[n] = first;
    first += nNode1;
  }
  for(n=nbr1;n<nbr1+nbr2;n++){
    _nRow[n] = nNode2;
    _firstRow[n] = first;
    first += nNode2;
  }
  for(n=nbr1+nbr2;n<nbr1+nbr2+nbr3;n++){
    _nRow[n] = nNode3;
    _firstRow[n] = first;
    first += nNode3;
  }
}


void ElemVec::showMe(ostream& c)
{
  for(int i=0;i<_nBlockRow;i++)
    c << "Block (" << i << "), " << block(i) << endl;
}
