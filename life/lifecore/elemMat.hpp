#ifndef _ELEMMAT_H_INCLUDED
#define _ELEMMAT_H_INCLUDED
#include <vector>
/*
#include "comprow_double.hpp"   // SparseLib++ matrix for "assemble_SparseLib"
*/
#include "tab.hpp"

class ElemMat
{
  Tab2d _mat; // the array
  UInt _nBlockRow; // number of block rows
  UInt _nBlockCol; // number of block columns
  vector<UInt> _nRow; // _nRow[i]=nb of rows in the i-th block row 
  vector<UInt> _firstRow;//_firstRow[i]=index of first row of i-th block row
  vector<UInt> _nCol; // _nCol[i]=nb of col in the i-th block col 
  vector<UInt> _firstCol;//_firstCol[i]=index of first col of i-th block col
public:
  ~ElemMat();
  ElemMat(UInt nNode1,UInt nbr1,UInt nbc1);// constructor for 1 finite element
  ElemMat(UInt nNode1,UInt nbr1,UInt nbc1, 
	  UInt nNode2,UInt nbr2,UInt nbc2);// constructor for 2 finite elements
  ElemMat(UInt nNode1,UInt nbr1,UInt nbc1,
	  UInt nNode2,UInt nbr2,UInt nbc2,
	  UInt nNode3,UInt nbr3,UInt nbc3);// constructor for 3 finite elements
  inline Tab2d& mat(){return _mat;}
  inline UInt nBlockRow()const{return _nBlockRow;}
  inline UInt nBlockCol()const{return _nBlockCol;}
  
  inline Tab2dView block(UInt i,UInt j){
    return _mat(SubArray(_nRow[i],_firstRow[i]),
		SubArray(_nCol[j],_firstCol[j]));}
  inline void zero(){_mat = 0;};
  void showMe(ostream& c=cout);
};

#endif


