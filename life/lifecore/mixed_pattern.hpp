/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003 LifeV Team
  
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
/*----------------------------------------------------------------------*
|          
| $Header: /cvsroot/lifev/lifev/life/lifecore/Attic/mixed_pattern.hpp,v 1.2 2004-03-30 12:14:06 prudhomm Exp $ 
|
|
| #Version  0.0 Experimental   19/6/00. Luca Formaggia                  |
|
| #Purposes Defines Patterns for block  matrices (Stokes Problem)
|
|  This release require a pattern class with the following methods:
|
|  UInt nRows(), UInt nCols(); UInt nNz(); bool empty(); Uint nbNeighbours(UInt ), UInt neighbour(UInt , UInt )
|  void Showme(ostream &), void spy(), buildPattern<DOF>(DOF &)
*----------------------------------------------------------------------*/
#ifndef _MIXED_PATTERNS_HH
#define _MIXED_PATTERNS_HH
#ifndef _LIFEV_HH_ 
//more correct version
typedef size_t UInt;
//original version
typedef unsigned int UInt;
#define INLINE inline
#endif
#include<set>
#include<algorithm>
#include "bareItems.hpp"
  

/* The mixed pattern class is a very general class able to held
  multiblock matrix patterns.  The block numbering starts from
  (0,0). Each block contains a pointer to a Pattern class and a couple
  of offsets, which indicate how the LOCAL block row/cols numbering
  has to be increased to get the GLOBAL numbering (i.e. the numbering
  associated to to the global matrix).  */

template <UInt BROWS, UInt BCOLS, typename PATTERN=CSRPatt>
class
MixedPattern : public PatternDefs
{
public:
  
  MixedPattern();
  ~MixedPattern();
  
  pair<Diff_t,Diff_t> nBlocks() const; // Number of blocks (rows and columns)
  
  UInt nRows(Diff_t const m; Diff_t const n) const; // Number of rows in block (m,n)
  UInt nCols(Diff_t const m; Diff_t const n) const; 
  UInt nNz(Diff_t m; Diff_t n) const; // Non zeros on Block (m,n)
  UInt nRows() const; // Global number of rows
  UInt nCols() const; 
  UInt nNz() const; // Non zeros in global matrix

  // Neigbours at block level. (local ID numbering)
  UInt nbNeighbours(Diff_t const m, Diff_t ID n,ID const d) const ;
  ID neighbour(Diff_t const m, Diff_t const n,ID const i, ID const d) const;
  void neighbours(Diff_t const m, Diff_t const n, ID const d,Container & neighs) const;

  // Neigbours at global level (global ID numbering)
  UInt nbNeighbours(ID const d_g) const ;
  ID neighbour(ID const i_g, ID const d_g) const;
  void neighbours(ID const d_g, Container & neighs) const;

  PATTERN * block_ptr(Diff_t const m, Diff_t const n); // Pointer to a block
  
  pair<Diff_t,Diff_t> blockOffset(Diff_t const m, Diff_t const n) const;// The row/col offsets of the block
  pair<Diff_t,Diff_t> locateElBlock(Index_t const i_g, Index_t const j_g) const; // 
  //  Give the block corresponding to the  GLOBAL matrix index (i_g,j_g) Returns
  // (BROWS,BCOLS) if element not found
  pair<Diff_t,Diff_t> locateDofBlock(ID const di_g, ID const dj_g) const; 
  // Give the block correponding to a THE GLOBAL DOF (di_g,dj_g) 
  // Returns (BROWS,BCOLS) if dofs not found
  pair<ID,ID> localNumber(Diff_t const m, Diff_t const n, ID const i_g, ID const j_g) const; 
  pair<ID,ID> globalNumber(Diff_t const m, Diff_t const n, ID const i, ID const j) const; 
  // local/global numbering in the block, given global numbering. It can be applied to indices and IDs. 
  // I rely on implicit conversion ID->Index_t if Index_t != ID 
  // If that does not wark do a template function pair<A,A>localNumber<T>(Diff_t,Diff_t,T,T) and the
  // necessary specialisations.

  pair<Diff_t,bool> locateDof(Diff_t const m, Diff_t const n, ID const di, ID const dj) // Locate DOF in block (local numbering) 
    pair<DIff_t,bool> locateIndex(Diff_t const m, Diff_t const n, Index_t const i, Index_t const j) // Locate element (i,j) in block (local numbering)
  
  // I can set the size of a block without linking the block to a pattern (I may have a matrix of ZEROS!)
  // To that purpose I use the following
  void buildZeroPattern(Diff_t const m, Diff_t const n,UInt const nrows, UInt const ncols);

  // This just links the block to an existing pattern
  void linkBlockToPattern(Diff_t const m, Diff_t const n,PATTERN & pattern);
   
  void linkBlocks(Diff_t const m_from, Diff_t const n_from, Diff_t const m_to, Diff_t const n_to); // if two blocks have teh same pattern
  
  void deleteBlock(Diff_t const m, Diff_t const n); // Use with extreme care
  
  
  // Checking and probing..
  INLINE bool isZero(Diff_t const m, Diff_t const n) const; // This is true if a block is set but not linked to a pattern (zero matrix)
  INLINE bool isSet(Diff_t const m, Diff_t const n) const;
  void ShowMe(ostream & c=cout, bool verbose=false) const;
  //bool check(bool verbose=false, ostream & c=cout) const;
  //void spy() const;
  bool consistentOffset(bool verbose=false, ostream & c=cout) const; // Tests if offset is consistent with that of a global matrix
  // If vercose is set it also wites on ostream

protected:
  
private:
  void resetOffset(Diff_t const m=0, Diff_t const n=0);
  

  PATTERN * _blocks[BROWS][BCOLS];
  bool _linked[BROWS][BCOLS]; // indicates whether a block is actually a link to another block
  UInt _rowoff[BROWS][BCOLS]; // rows offset
  UInt _coloff[BROWS][BCOLS];// cols offset
  UInt _nrows[BROWS][BCOLS]; // rows in block
  UInt _ncols[BROWS][BCOLS]; // cols in block
  
};



/*****************************************************************************************/
//                               IMPLEMENTATIONS
/*****************************************************************************************/

/*****************************************************************************************/
/*                              Mixed Pattern                */
/*****************************************************************************************/

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
MixedPattern<BROWS,BCOLS,PATTERN>::
MixedPattern()
{
  for (UInt i=0; i<BROWS; i++){
    for (UInt j=0; j<BCOLS; j++){
      _rowoff[i][j]=0;
      _coloff[i][j]=0;
      _linked[i][j]=false;
      _blocks[i][j]=0;
      _nrows[i][j]=0;
      _ncols[i][j]=0;
    }
  }
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
MixedPattern<BROWS,BCOLS,PATTERN>::
~MixedPattern()
{
  for (UInt i=0; i<BROWS; i++){
    for (UInt j=0; j<BCOLS; j++){
      if (_blocks[i][j] != 0 && !_linked[i][j]){
	delete _blocks[i][j]; // blocks are actually deleted only if not linked
	_blocks[i][j]=0;
      }
    }
  }
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
INLINE
pair<UInt,UInt>
MixedPattern<BROWS,BCOLS,PATTERN>::nBlocks() const{
  return make_pair(BROWS,BCOLS);
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
INLINE
UInt
MixedPattern<BROWS,BCOLS,PATTERN>::nRows(Diff_t m; Diff_t n) const{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BROWS && n>=0, " Invalid block column address");
  return _nrows[m][n];
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
INLINE
UInt
MixedPattern<BROWS,BCOLS,PATTERN>::nCols(Diff_t m; Diff_t n) const{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BROWS && n>=0, " Invalid block column address");
  return _ncols[m][n];
}


template
<UInt BROWS, UInt BCOLS, typename PATTERN>
INLINE
UInt
MixedPattern<BROWS,BCOLS,PATTERN>::nNz(Diff_t m; Diff_t n) const{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BROWS && n>=0, " Invalid block column address");

  return _blocks[m][n] ? _blocks[i][j]->nNz():0;
}


template
<UInt BROWS, UInt BCOLS, typename PATTERN>
INLINE
UInt
MixedPattern<BROWS,BCOLS,PATTERN>::
nRows() const{
  // I assume that the blocks are consistent
  return _rowoff[BROWS-1][BCOLS-1] + _nrows[BROWS-1][BCOLS-1];
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
INLINE
UInt
MixedPattern<BROWS,BCOLS,PATTERN>::
nCols() const{
  // I assume that the blocks are consistent
  return _rowoff[BROWS-1][BCOLS-1] + _ncols[BROWS-1][BCOLS-1];
}


template
<UInt BROWS, UInt BCOLS, typename PATTERN>
UInt
MixedPattern<BROWS,BCOLS,PATTERN>::nNz() const{
  UInt _n=0;
  for (UInt m=0; m<BROWS; m++){
    for (UInt n=0; n<BCOLS; n++){
      _n+=nNz(m,n);
    }
  }
  return _n;
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
UInt
MixedPattern<BROWS,BCOLS,PATTERN>::nbNeighbours(Diff_t const m,  Diff_t const n,ID const d) const{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BROWS && n>=0, " Invalid block column address");
  return (_blocks[m][n] ? _blocks[m][n]->nbNeighbours(d):0);
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
ID
MixedPattern<BROWS,BCOLS,PATTERN>::neighbour(Diff_t const m, Diff_t const n,ID const i, ID const d) const
{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BROWS && n>=0, " Invalid block column address");
  ASSERT_PRE(_blocks[m][n] !=0 , "Cannot access an empty block") ;
  return _blocks[m][n]->neighbour(i,d);
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
INLINE
void
MixedPattern<BROWS,BCOLS,PATTERN>::neighbours(Diff_t const m, Diff_t const n,ID const d, Container & neighs) const
{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BROWS && n>=0, " Invalid block column address");
  ASSERT_PRE(_blocks[m][n] !=0 , "Cannot access an empty block") ;
  return _blocks[m][n]->neighbours(d, neighs);
}


template
<UInt BROWS, UInt BCOLS, typename PATTERN>
UInt
MixedPattern<BROWS,BCOLS,PATTERN>::nbNeighbours(ID const d_g) const{
  // This is MORE complicated
  //Locate line corresponding to dof d_g in block
  pair<UInt,UInt> _b=locateDofBlock(d_g,1); // Remember that dofs are ALWAYS numbered from 1
  // line local nomber
  UInt m=b.first;
  UInt _d=d_g-_rowoff[m][0];
  UInt _n=0;
  // Sweep columns
  for (UInt n=0; n<BCOLS; n++){
    _n+=(_blocks[m][n] ? _blocks[m][n]->nbNeighbours(_d):0);
  }
  return _n;
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
ID
MixedPattern<BROWS,BCOLS,PATTERN>::neighbour(ID const i_g, ID const d_g) const //ith neighbour of d (Global!)
{
  // This is MUCH MORE complicated
  //Locate line in blocks
  pair<UInt,UInt> _b=locateBlock(d_g,0);
  ASSERT_PRE(_b<BROWS, "Invalid dof entry d_g");
  // Get local numbering of dof
  UInt m=b.first;
  UInt _d=d_g-_rowoff[m][0];
  
  // In which block does i_g fall?
  // I remember that all high level interfaces use the convention of having numbering from 1
  // So i_g and d_g have always numbering from 1

  UInt _first=0;
  UInt _last=0;
  UInt _i=-1;
  UInt n=0;
  // Sweep columns
  for (UInt n=0; n<BCOLS; n++){
    _last=_first+nbNeighbours(m,n,_d);
    if(i_g >_first && i_g <= _last){
      _i=i_g-_first;// local numbering 
      break;
    }
    _first=_last;
  }
  ASSERT_PRE(_i !=1, " Invalid i_g") ;

  return _block[m][n]->neighbour(_i,_d) + _coloff[m][n];
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
void
MixedPattern<BROWS,BCOLS,PATTERN>::neighbours(ID const d_g,Container & neighs) const //ith neighbour of d (Global!)
{
  // This is MUCH MORE complicated
  //Locate line in blocks
  pair<UInt,UInt> _b=locateBlock(d_g,0);
  ASSERT_PRE(_b<BROWS, "Invalid dof entry d_g");
  // Get local numbering of dof
  UInt m=b.first;
  UInt _d=d_g-_rowoff[m][0];
  
  // In which block does i_g fall?
  // I remember that all high level interfaces use the convention of having numbering from 1
  // So i_g and d_g have always numbering from 1
  neighs.clear();
  neighs.resize(nbNeighbours(d_g));
  if neighs.empty() return;
  Container piece;
  // Sweep columns
  for (UInt n=0; n<BCOLS; n++){
    if (_blocks[m][n] != 0 ){
      _blocks[m][n]->neighbours(_d,piece);
      for(Container::iterator ip=piece.begin();ip=!piece.end();++ip) *ip+=_coloff[m][n];
      neighs.insert(neighs.end(),piece.begin(),piece.end());
    }
  }
}


template
<UInt BROWS, UInt BCOLS, typename PATTERN>
INLINE
PATTERN *
MixedPattern<BROWS,BCOLS,PATTERN>::block_ptr(Diff_t const m, Diff_t const n){
  // VERY DANGEROUS
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BROWS && n>=0, " Invalid block column address");
  return _blocks[m][n]; // Beware of null pointers!!
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
INLINE
pair<Diff_t,Diff_t>
MixedPattern<BROWS,BCOLS,PATTERN>::blockOffset(UInt const m, UInt const n) const{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BROWS && n>=0, " Invalid block column address");
  return make_pair(_rowoff[m][n],_coloff[m][n]);
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
pair<Diff_t,Diff_t>
MixedPattern<BROWS,BCOLS,PATTERN>::
locateElBlock(Index_t const i_g, Index_t const j_g) const
{
  UInt m=0;
  UInt n=0;
  UInt _itest=_i2o(i_g);
  UInt _jtest=_i2o(j_g);
  
  while(m<BROWS &&(_itest <_rowoff[m][0] || _itest >=_rowoff[m][0]+ _nrows[m][n])) ++m;
  while(n<BROWS &&(_jtest <_coloff[0][n] || _jtest >=_coloff[0][n]+ _ncols[m][n])) ++n;
  return make_pair(m,n);
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
pair<Diff_t,Diff_t>
MixedPattern<BROWS,BCOLS,PATTERN>::
locateDofBlock(ID const i_g, ID const j_g) const
{
  UInt m=0;
  UInt n=0;
  UInt _itest=_d2o(i_g);
  UInt _jtest=_d2o(j_g);
  
  while(m<BROWS &&(_itest <_rowoff[m][0] || _itest >=_rowoff[m][0]+ _nrows[m][n])) ++m;
  while(n<BROWS &&(_jtest <_coloff[0][n] || _jtest >=_coloff[0][n]+ _ncols[m][n])) ++n;
  return make_pair(m,n);
}


template
<UInt BROWS, UInt BCOLS, typename PATTERN>
INLINE
pair<ID,ID>
MixedPattern<BROWS,BCOLS,PATTERN>::localNumber(Diff_t const m, Diff_t const n,ID const i_g, ID const j_g) const{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BROWS && n>=0, " Invalid block column address");

  return make_pair(i_g-_rowoff[m][n],j_g-_coloff[m][n]);
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
INLINE
pair<ID,ID>
MixedPattern<BROWS,BCOLS,PATTERN>::globalNumber(Diff_t const m, Diff_t const n,ID const i, ID const j) const{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BROWS && n>=0, " Invalid block column address");
  return make_pair(i+_rowoff[m][n],j+_coloff[m][n]);
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
pair<Diff_t,bool>
MixedPattern<BROWS,BCOLS,PATTERN>::locateIndex(Diff_t const m, Diff_t const n,Index_t const i, Index_t const j) const{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BROWS && n>=0, " Invalid block column address");
  if (_blocks[m][n]==0)return make_pair(0,false);
  return _blocks[m][n]->locateIndex(i,j);
}


<UInt BROWS, UInt BCOLS, typename PATTERN>
pair<Diff_t,bool>
MixedPattern<BROWS,BCOLS,PATTERN>::locateDof(Diff_t const m, Diff_t const n,ID const i, ID const j) const{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BROWS && n>=0, " Invalid block column address");
  if (_blocks[m][n]==0)return make_pair(0,false);
  return _blocks[m][n]->locateDof(i,j);
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
void
MixedPattern<BROWS,BCOLS,PATTERN>::
buildZeroPattern(Diff_t const m, Diff_t const n,UInt const nrows, UInt const ncols){
  _nrows[m][n]=nrows;
  _ncols[m][n]=ncols;
  resetOffset(m,n);  
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
void
MixedPattern<BROWS,BCOLS,PATTERN>::
linkBlockToPattern(UInt const m, UInt const n,PATTERN & pattern){
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BROWS && n>=0, " Invalid block column address");
  if(_blocks[m][n]!=0 && !_linked[m][n]) delete _blocks[m][n];
  
  _blocks[m][n]=&pattern;
  _nrows[m][n]=pattern.nRows();
  _ncols[m][n]=pattern.nCols();
  _linked[m][n]=true;
  resetOffset(m,n);  
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
void
MixedPattern<BROWS,BCOLS,PATTERN>::
linkBlocks(Diff_t const m_from, Diff_t const n_from, Diff_t const m_to, Diff_t const n_to);{
  // ASSERT_PRE if it is ok!
  ASSERT_PRE( m_from<BROWS && m_from>=0, " Invalid block row address");
  ASSERT_PRE( n_from<BROWS && n_from>=0, " Invalid block column address");
  ASSERT_PRE( m_to<BROWS && m_to>=0, " Invalid block row address");
  ASSERT_PRE( n_to<BROWS && n_to>=0, " Invalid block column address");

  _blocks[m_to][n_to]=_blocks[m_from][n_from];
  _nrows[m_to][n_to]=_nrows[m_from][n_from];
  _ncols[m_to][n_to]=_ncols[m_from][n_from];
  _linked[m_to][n_to]=true;
  resetOffset(m,n);  
}


template
<UInt BROWS, UInt BCOLS, typename PATTERN>
void
MixedPattern<BROWS,BCOLS,PATTERN>::
deleteBlock(Diff_t const m, Diff_t const n);{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BROWS && n>=0, " Invalid block column address");
  if(_blocks[m][n]!=0&& ! linked[m][n]) delete _blocks[m][n];
  _blocks[m][n]=0;
  _nrows[m][n]=0;
  _ncols[m][n]=0;
  _linked[m][n]=false;
  resetOffset(m,n);
}



template
<UInt BROWS, UInt BCOLS, typename PATTERN>
INLINE
bool
MixedPattern<BROWS,BCOLS,PATTERN>::isZero(Diff_t const m, DIff_t const n) const{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BROWS && n>=0, " Invalid block column address");
  return _blocks[m][n]=0 && _nrows[m][n]+_ncols[m][n] != 0;
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
INLINE
bool
MixedPattern<BROWS,BCOLS,PATTERN>::isSet(Diff_t const m, Diff_t const n) const{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BROWS && n>=0, " Invalid block column address");
  return _blocks[m][n]!=0;
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
void
MixedPattern<BROWS,BCOLS,PATTERN>::resetOffset(Diff_t const m, Diff_t const n) // Puts offset right!
{
  UInt _m=m;
  UInt _n=n;
  
  while( _m<BROWS && _n <BCOLS){
    
    for (UInt i=_m+1; i<BROWS; i++){
      _rowoff[i][_n]=_rowoff[i-1][_n]+_nrows[i-1][_n];
    }
    
    for (UInt i=_n+1; i<BCOLS; i++){
      _coloff[_m][i]=_coloff[_m][i]+_ncols[_m][i-1];
    }
    _m++;
    _n++;
  }
}


template
<UInt BROWS, UInt BCOLS, typename PATTERN>
bool
MixedPattern<BROWS,BCOLS,PATTERN>::consistentOffset(bool verbose, ostream & c) const{
  bool ok=true;
  
  for (Diff_t m=1; m<BROWS; m++){
    for (Diff_t n=1; n<BCOLS; n++){
      ok=ok && _rowoff[m][n]=_rowoff[m][n-1] &&_coloff[m][n]=_coloff[m-1][n];
    }
  }
  for (Diff_t m=1; m<BROWS; m++){
    for (Diff_t n=1; n<BCOLS; n++){
      ok=ok && _rowoff[m][n]-_rowoff[m-1][n] +_nrows[m-1][n]==0 &&
	_coloff[m][n]-_coloff[m][n-1]-_ncols[m][n-1]==0;
    }
  }

  if (verbose && ! ok){
    c<< " Error checking mixed pattern"<<endl;
    for (Diff_t m=1; m<BROWS; m++){
      for (Diff_t n=1; n<BCOLS; n++){
	if(_rowoff[m][n]!=_rowoff[m][n-1] || _coloff[m][n]!=_coloff[m-1][n]){
	  c<<"Inconsistent Offsets of Blocks:"<<endl;
	  c<<"("<<m<<", "<<n<<") Row offset= "<<_rowoff[m][n]<<
	   <<"("<<m<<", "<<n-1<<") Row offset= "<<_rowoff[m][n-1]<<endl;
	  c<<"("<<m<<", "<<n<<") Col offset= "<<_coloff[m][n]<<
	   <<"("<<m-1<<", "<<n<<") Col offset= "<<_coloff[m-1][n]<<endl;
	}
      }
    }
    for (Diff_t m=1; m<BROWS; m++){
      for (Diff_t n=1; n<BCOLS; n++){
	if(_rowoff[m][n]-_rowoff[m-1][n] +_nrows[m-1][n]!=0 ||
	   _coloff[m][n]-_coloff[m][n-1]-_ncols[m][n-1]!=0){
	  c<<"Inconsistent Number of cols/rows:"<<endl;
	  c<<"("<<m<<", "<<n<<") Row offset= "<<_rowoff[m][n];
	  c<<"("<<m-1<<", "<<n<<") Row offset= "<<_rowoff[m-1][n]<<endl;
	  c<<"("<<m-1<<", "<<n<<") N Rows= "<<_nrows[m-1][n]<<endl;
	  c<<"("<<m<<", "<<n<<") Col offset= "<<_coloff[m][n];
	  c<<"("<<m<<", "<<n-1<<") Col offset= "<<_coloff[m][n-1]<<endl;
	  c<<"("<<m<<", "<<n-1<<") N Cols= "<<_ncols[m][n-1]<<endl;
	}
      }
    }
  }
  return ok;
}


template
<UInt BROWS, UInt BCOLS, typename PATTERN>
bool
MixedPattern<BROWS,BCOLS,PATTERN>::showMe(bool verbose, ostream & c) const{
  c<<endl;
  c<<" ******************* MIXED PATTERN *******************"<<endl;
  c<<BROWS<<"X"<<BCOLS<<" Blocks"<<endl;
  c<<"Block    Row Off   Cols OFF    Nrows     Ncols     Set      Linked  "<<endl;
  for (Diff_t m=1; m<BROWS; m++){
    for (Diff_t n=1; n<BCOLS; n++){
      c<<"("<<m<<", "<<n<<") "<<_rowoff[m][n]<<" "<<_coloff[m][n]<<" "<<_nrows[m][n]<<" "<<_ncols[m][n]<<" ";
      c<<isSet(m,n)<<" "<<_linked[m][n]<<endl;
    }
  }
  if(verbose){
    for (Diff_t m=1; m<BROWS; m++){
      for (Diff_t n=1; n<BCOLS; n++){
	if (_blocks[m][n] != 0 )_blocks[m][n]->showMe(false,c); 
      }
    }
  }
}
#endif
