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
#include "pattern.hpp"
#include <fstream>

namespace LifeV
{
using namespace std;

//////////////////////////////////////////////////////////////
// BasePattern
//////////////////////////////////////////

BasePattern::BasePattern():
  _nnz(_defaultsize), _nrows(_defaultsize), _ncols(_defaultsize), _filled(false), _diagfirst(false)
{}

BasePattern::BasePattern(UInt ex_nnz, UInt ex_nrow, UInt ex_ncol):
  _nnz(ex_nnz), _nrows(ex_nrow), _ncols(ex_ncol), _filled(false),_diagfirst(false)
{}

// ALAIN : this constructor should exist as well for other matrix classes
BasePattern::BasePattern(CSRPatt const &RightHandCSRP):
  _nnz(RightHandCSRP.nNz()), _nrows(RightHandCSRP.nRows()),
  _ncols(RightHandCSRP.nCols()), _filled(!RightHandCSRP.isEmpty()),
  _diagfirst(RightHandCSRP.diagFirst())
{}

void
BasePattern::showMe(bool const verbose,ostream & out) const{
  out <<"==== BASE Pattern Class (BasePattern)===="<<endl;
  out << "Number of rows :" << nRows() <<endl;
  out << "Number of cols :" << nCols() <<endl;
  out << "Non zero terms :" << nNz() <<endl;
  out << "Empty          :" << isEmpty() <<endl;
}

////////////////////////////////////////////////////////////////////////
//
// C S R Pattern
//
////////////////////////////////////////////////////////////////////////
CSRPatt::CSRPatt()
{
    // nothing to do here
}

CSRPatt::CSRPatt(UInt ex_nnz, UInt ex_nrow, UInt ex_ncol):
  BasePattern(ex_nnz,ex_nrow,ex_ncol)
{
  _ia.reserve(ex_nrow+1);
  _ja.reserve(ex_nnz);
}

CSRPatt::CSRPatt(UInt ex_nnz, UInt ex_nrow, UInt ex_ncol, const vector<Index_t> &ex_ia, const vector<Index_t> &ex_ja):
  BasePattern(ex_nnz,ex_nrow,ex_ncol),_ia(ex_ia),_ja(ex_ja)
{
  ASSERT_PRE(ex_ia.size()==ex_nrow+1 && ex_ja.size()==ex_nnz,"Error in CSR Pattern Life V"); // check on the compatibility of the external data
  _filled = !ex_ia.empty(); // Test if the given containers are empty!
  if(_filled)_diagfirst=_i2o(_ja[_i2o(_ia[_nrows-1])])==_nrows-1; //stupid test
}

CSRPatt::CSRPatt(const CSRPatt &RightHandCSRP):BasePattern(RightHandCSRP),
					       _ia(RightHandCSRP.ia()),_ja(RightHandCSRP.ja()){}


CSRPatt::CSRPatt(MSRPatt const& msrPatt)
    : BasePattern(msrPatt.nNz(), msrPatt.nRows(), msrPatt.nCols()) {
    Container const& bindx = msrPatt.give_bindx();
    _ja.resize(_nnz);
    _ia.resize(_nrows+1);
    _jaT.resize(_nnz);
    UInt nOffset = 0;

    _ia[0] = 0;
    for(UInt iRow=0; iRow<_nrows; ++iRow) {
        UInt rowLength = bindx[iRow+1] - bindx[iRow];
        //std::cout << nOffset << "/" << _ja.size() << std::endl;
        //std::cout << iRow << "," << bindx[bindx[iRow]] << std::endl;
        _ja[nOffset++] = iRow;
        for(UInt ii=0; ii<rowLength; ++ii) {
            //std::cout << nOffset << "/" << _ja.size() << std::endl;
            _ja[nOffset++] = bindx[ii+bindx[iRow]-PatternOffset]-PatternOffset;
            //std::cout << iRow << "," << bindx[ii+bindx[iRow]] << std::endl;
        }
        _ia[iRow+1] = rowLength+1+_ia[iRow];
    }

    // fill _jaT

    Container iaT(_ncols+1, 0);
    // in a first step, fill iaT[k+1] with the number of entries in column k
    for(UInt i=0; i<_nnz; ++i) {
        iaT[_ja[i]+1]++;
    }
    // then update iaT to store absolute positions of in _jaT
    for(UInt iCol=2; iCol<_ncols+1; ++iCol) {
        iaT[iCol] = iaT[iCol]+iaT[iCol-1];
    }

    // fill diagonal entries
    // update iaT[iCol] to hold next position in _jaT for column iCol
    for(UInt iCol=0; iCol<_ncols; ++iCol) {
        if (iaT[iCol] >= _nnz) {
            std::cout << iaT[iCol] << ">=" << _nnz << std::endl;
        }
        _jaT[iaT[iCol]++] = _ia[iCol];
    }

    for(UInt i=0; i<_nnz; ++i) {
        if ( i != _ia[_ja[i]] ) {
            _jaT[iaT[_ja[i]]++] = i;
        }
    }

    // take PatternOffset into account
    if (PatternOffset != 0) {
        for (Container::iterator ip=_ia.begin(); ip != _ia.end(); ++ip) {
            *ip += PatternOffset;
        }
        for (Container::iterator ip=_ja.begin(); ip != _ja.end(); ++ip) {
            *ip += PatternOffset;
        }
        for (Container::iterator ip=_jaT.begin(); ip != _jaT.end(); ++ip) {
            *ip += PatternOffset;
        }
   }
    _filled = true;
    _diagfirst = true;
}

CSRPatt& CSRPatt::operator= (const CSRPatt& RhCsr)
{
  if (&RhCsr != this)
    {
      _nnz = RhCsr.nNz();
      _nrows = RhCsr.nRows();
      _ncols = RhCsr.nCols();
      _ia = RhCsr.ia();
      _ja = RhCsr.ja();
      _filled= !  RhCsr.isEmpty();
      _diagfirst=RhCsr.diagFirst();
    }
  return *this;
}

pair<PatternDefs::Diff_t,bool>
CSRPatt::locate_pattern(Index_t const i,Index_t const j) const {
  ASSERT_BD(i>=PatternOffset && i<static_cast<Index_t>(_nrows)+PatternOffset);
  ASSERT_BD(j>=PatternOffset && j<static_cast<Index_t>(_ncols)+PatternOffset);
  if (! _filled ) {
    return make_pair(0,false);
  }
  Diff_t _off=_row_off(i); // offset of point associated to index i
  Container::const_iterator start=_ja.begin()+_off; // the start of row i

  // This part is here to account for the possibility of having the diagonal at the first place (without ASSUMING it)
  if (*start == j){
    return  make_pair(_off,true);
  }
  Container::const_iterator finish=_ja.begin()+_row_off(i+1); //the real end (remember STL convention for ranges!)
  Container::const_iterator current=start+1;
  current = search_binary(current, finish, j); // search in the rest (the first position has been already checked)
  // difference of pointers should  return distance, which is of integral type, like Diff_t
  return make_pair(current-_ja.begin(),current !=finish);
}

void  CSRPatt::showMe(bool verbose, ostream& c) const
{
  BasePattern::showMe(verbose,c);

  typedef vector<Index_t>::iterator found;
  int i_first;
  string pare="[";
  c << "**************************" << endl;
  c << "     CSR Matrix Pattern   " << endl;
  c << endl;
  if (verbose){

    c << pare;
    for (Diff_t i_index=0; i_index < static_cast<Diff_t>(_nrows); i_index++)
      {
	c << pare;
	pare = " [";
	i_first=_i2o(_ia[i_index])+1; // In _ia[i_index] there is the diagonal entry
	UInt jj=0;
	for(Diff_t j=0;j<static_cast<Diff_t>(_ncols);j++)
	  {
	    if (j==i_index) c << " * ";
	    else {
	      if (j==_i2o(_ja[i_first+jj])){
		c << " * "; jj++;}
	      else c << " 0 ";
	    }
	  }
	if (i_index==static_cast<Diff_t>(_nrows-1))
	  c << " ]] " << endl;
	else
	  c << " ]  " << endl;
      }
  }
  c << "**************************" << endl;
  return;
}

void  CSRPatt::spy(string  const &filename) const
{
  // Purpose: Matlab dumping and spy
  string nome=filename, uti=" , ";
  //
  // check on the file name
  //
  unsigned int i=filename.find(".");

  if (i<=0) nome=filename+".m";
  else {
    if (i!=filename.size()-2  || filename[i+1]!='m')
      {cerr << "Wrong file name ";
      nome=filename+".m";}
  }

  ofstream file_out(nome.c_str());


  file_out << "S = [ ";
  for (UInt i=0;i<_nrows;++i){
    for (Index_t ii=_ia[i]-PatternOffset;ii<_ia[i+1]-PatternOffset;++ii)
      file_out << i+1 << uti << _ja[ii]+1-PatternOffset << uti << "1.0" << endl; /* */
  }
  file_out << "];" << endl;

  file_out << "I=S(:,1); J=S(:,2); S=S(:,3); A=sparse(I,J,S); spy(A);"<<endl;
}

// column-concatenation of two CSR block patterns

CSRPatt colUnify(CSRPatt const &patt1, CSRPatt const &patt2)
{

    cout << "colUnify: OBSOLETE function ==> use CSRPatt constructor or MixedPattern class" << endl;

    typedef PatternDefs::Container::const_iterator ConstIter;
    typedef PatternDefs::Container::iterator Iter;

    UInt nnz1=patt1._nnz, nnz2=patt2._nnz;
    UInt nrows1=patt1._nrows, nrows2=patt2._nrows;
    UInt ncols1=patt1._ncols, ncols2=patt2._ncols;

    ASSERT(nrows1==nrows2, "pattern1 and pattern2 must have same nrows");


    UInt nnz= nnz1 + nnz2;
    UInt nrows= nrows1;
    UInt ncols= ncols1 + ncols2;

    UInt end = 0;

    bool diag = patt1._diagfirst;

    if (diag)
        if (nrows > ncols1)
        {
            pair<UInt,bool> where;
            if (ncols2 <= (nrows - ncols1))
                end = ncols2;
            else
                end = nrows - ncols1;
            for (UInt i = 0; i < end; ++i)
            {
                where = patt2.locate_index(i+ncols1,i);
                if (!where.second)
                {
                    diag = false;
                    break;
                }
            }
        }

    CSRPatt ans(nnz, nrows, ncols);

    ans._ja.resize(nnz,PatternOffset);

    Iter ja_start = ans._ja.begin();

    if (!diag && patt1._diagfirst)
    {
        for (UInt i=0; i < nrows; i++)
        {
            // construction of ia
            ans._ia.push_back(patt1._ia[i]+patt2._ia[i]-PatternOffset);

            // construction of ja
            ConstIter ja1_start = patt1._ja.begin() + patt1._i2o(patt1._ia[i]);
            ConstIter ja2_start = patt2._ja.begin() + patt2._i2o(patt2._ia[i]);
            ConstIter ja1_end = patt1._ja.begin() + patt1._i2o(patt1._ia[i+1]);
            ConstIter ja2_end = patt2._ja.begin() + patt2._i2o(patt2._ia[i+1]);

            Iter ja_start1 = ja_start + ans._i2o(ans._ia[i]);
            Iter ja_end1 = ja_start1 + (patt1._ia[i+1]-patt1._ia[i]);

            // copy of the first block
            partial_sort_copy(ja1_start, ja1_end, ja_start1, ja_end1);
            //the starting place of the second block
            Iter ja_start2 = ja_end1;
            // the end of the second block
            Iter ja_end2 = ja_start2 + (patt2._ia[i+1]-patt2._ia[i]);
            // copy of the second block
            partial_sort_copy(ja2_start, ja2_end, ja_start2, ja_end2);

            // adds ncols1 to the second block
            for (Iter ip=ja_start2; ip != ja_end2; ++ip)
                *ip+=ncols1;
        }
    }
    else
    {
        for (UInt i=0; i < nrows; i++)
        {
            // construction of ia
            ans._ia.push_back(patt1._ia[i]+patt2._ia[i]-PatternOffset);

            // construction of ja
            ConstIter ja1_start = patt1._ja.begin() + patt1._i2o(patt1._ia[i]);
            ConstIter ja2_start = patt2._ja.begin() + patt2._i2o(patt2._ia[i]);
            ConstIter ja1_end = patt1._ja.begin() + patt1._i2o(patt1._ia[i+1]);
            ConstIter ja2_end = patt2._ja.begin() + patt2._i2o(patt2._ia[i+1]);

            Iter ja_start1 = ja_start + ans._i2o(ans._ia[i]);

            // copy of the first block
            copy(ja1_start, ja1_end, ja_start1);
            //the starting place of the second block
            Iter ja_start2 = ja_start1 + (patt1._ia[i+1]-patt1._ia[i]);
            // the end of the second block
            Iter ja_end = ja_start2 + (patt2._ia[i+1]-patt2._ia[i]);
            // copy of the second block
            partial_sort_copy(ja2_start, ja2_end, ja_start2, ja_end);

            // adds ncols1 to the second block
            for (Iter ip=ja_start2; ip != ja_end; ++ip)
                *ip+=ncols1;
        }
    }

    //ultimate _ia
    ans._ia.push_back(patt1._ia[nrows]+patt2._ia[nrows]-PatternOffset);
    ans._filled = true;
    ans._diagfirst = diag;

    if (!diag)
        return ans;

    for (UInt i = 0; i < end; ++i)
    {
        UInt row = i + ncols1;

        Iter ja_start1 = ja_start + ans._i2o(ans._ia[row]);
        Iter ja_end1 = ja_start + ans._i2o(ans._ia[row+1]);

        // find diagonal element
        Iter ja_diag = find(ja_start1, ja_end1, row);

        // reorder
        rotate(ja_start1, ja_diag, ja_diag + 1);
    }

    return ans;
}


// column-concatenation of one CSR block patterns and ncolZero null columns
// zero on the right
CSRPatt colUnify(CSRPatt const &patt1, UInt const ncolZero)
{
  cout << "colUnify: OBSOLETE function ==> use CSRPatt constructor or MixedPattern class" << endl;
  UInt nnz=patt1._nnz;
  UInt nrows=patt1._nrows;
  UInt ncols= patt1._ncols + ncolZero;

  // nothing changes but the number of columns of the matrix
  CSRPatt ans(nnz, nrows, ncols, patt1._ia, patt1._ja);

  ans._filled = true;
  ans._diagfirst = patt1._diagfirst;
  return ans;
}

//zero on the left
CSRPatt colUnify(UInt const ncolZero, CSRPatt const &patt1)
{
  cout << "colunify: OBSOLETE function ==> use CSRPatt constructor or MixedPattern class" << endl;
  typedef PatternDefs::Container::iterator Iter;
  UInt nnz=patt1._nnz;
  UInt nrows=patt1._nrows;
  UInt ncols1=patt1._ncols;
  UInt ncols= ncols1 + ncolZero;
  CSRPatt ans(nnz, nrows, ncols, patt1._ia, patt1._ja);

  // adds ncolZero to the block
  for (Iter ip=ans._ja.begin(); ip != ans._ja.end(); ++ip)
    *ip+=ncolZero;

  ans._filled=true;
  ans._diagfirst = false;
  if (!patt1._diagfirst)
    return ans;

  // sort the col indexes
  for (UInt i=0; i < nrows; i++)
    {
      Iter ja_start = ans._ja.begin() + ans._i2o(patt1._ia[i]);
      Iter ja_end   = ans._ja.begin() + ans._i2o(patt1._ia[i+1]);
      sort(ja_start, ja_end);
    }
  return ans;
}

// rows-concatenation of two CSR block patterns
CSRPatt rowUnify(CSRPatt const &patt1, CSRPatt const &patt2)
{
    cout << "rowUnify: OBSOLETE function ==> use CSRPatt constructor or MixedPattern class" << endl;
    typedef PatternDefs::Container::const_iterator ConstIter;
    typedef PatternDefs::Container::iterator Iter;

    UInt nnz1=patt1._nnz, nnz2=patt2._nnz;
    UInt nrows1=patt1._nrows, nrows2=patt2._nrows;
    UInt ncols1=patt1._ncols, ncols2=patt2._ncols;

    ASSERT(ncols1==ncols2, "pattern1 and pattern2 must have same ncols");

    UInt nnz= nnz1 + nnz2;
    UInt nrows= nrows1 + nrows2;
    UInt ncols= ncols1;

    CSRPatt ans(nnz, nrows, ncols);
    ans._ia.resize(nrows+1,0);
    ans._ja.resize(nnz,0);

    UInt end = 0;

    bool diag = patt1._diagfirst;

    if (diag)
        if (ncols > nrows1)
        {
            pair<UInt,bool> where;
            if (nrows2 <= (ncols - nrows1))
                end = nrows2;
            else
                end = ncols - nrows1;
            for (UInt i = 0; i < end; ++i)
            {
                where = patt2.locate_index(i,i+nrows1);
                if (!where.second)
                {
                    diag = false;
                    break;
                }
            }
        }

    //construction of ia
    ConstIter ia1_start = patt1._ia.begin();
    ConstIter ia2_start = patt2._ia.begin();
    ConstIter ia1_end = patt1._ia.end()-1;
    ConstIter ia2_end = patt2._ia.end();

    // copy of the first block
    Iter ia_start1 = ans._ia.begin();
    copy(ia1_start, ia1_end, ia_start1);
    // copy of the second block
    Iter ia_start2 = ans._ia.begin() + nrows1;
    copy(ia2_start, ia2_end, ia_start2);
    // adds nnz1 to the second block
    Iter ia_end = ans._ia.end();
    for (Iter ip=ia_start2; ip != ia_end; ++ip)
        *ip+=nnz1;

    //construction of ja
    ConstIter ja1_start = patt1._ja.begin();
    ConstIter ja1_end = patt1._ja.end();

    //copy of the first block
    Iter ja_start = ans._ja.begin();

    if (!diag && patt1._diagfirst)
    {
        for (UInt i = 0; i < nrows1; ++i)
        {
            //construction of ja
            ConstIter ja1_start1 = ja1_start + patt1._i2o(patt1._ia[i]);
            ConstIter ja1_end1 = ja1_start  + patt1._i2o(patt1._ia[i+1]);

            Iter ja_start1 = ja_start + ans._i2o(ans._ia[i]);
            Iter ja_end1 = ja_start + ans._i2o(ans._ia[i+1]);

            //copy of the first block
            partial_sort_copy(ja1_start1, ja1_end1, ja_start1, ja_end1);
        }
    }
    else
        //copy of the first block
        copy(ja1_start, ja1_end, ja_start);

    //construction of ja
    ConstIter ja2_start = patt2._ja.begin();
    ConstIter ja2_end = patt2._ja.end();

    //copy of the second block
    if (patt2._diagfirst)
    {
        for (UInt i=0; i < nrows2; i++)
        {
            UInt row= nrows1+i;
            ConstIter ja2_start2 = ja2_start + patt2._i2o(patt2._ia[i]);
            ConstIter ja2_end2 = ja2_start + patt2._i2o(patt2._ia[i+1]);

            //the starting place of the second block
            Iter ja_start2 =  ja_start + ans._i2o(ans._ia[row]);
            // the end of the second block
            Iter ja_end2 = ja_start + ans._i2o(ans._ia[row+1]);

            // copy of the second block
            partial_sort_copy(ja2_start2, ja2_end2, ja_start2, ja_end2);
        }
    }
    else
        copy(ja2_start, ja2_end, ja_start + ans._ia[nrows1]);


    ans._filled=true;
    ans._diagfirst = diag;

    if (!diag)
        return ans;

    for (UInt i = 0; i < end; ++i)
    {
        UInt row = i + nrows1;

        Iter ja_start1 = ans._ja.begin() + ans._i2o(ans._ia[row]);
        Iter ja_end1 = ans._ja.begin() + ans._i2o(ans._ia[row+1]);

        // find diagonal element
        Iter ja_diag = find(ja_start1, ja_end1, row);

        // reorder
        rotate(ja_start1, ja_diag, ja_diag + 1);
    }

    return ans;
}

// row-concatenation of one CSR block patterns and nrowZero null rows
// zero on the below
CSRPatt rowUnify(CSRPatt const &patt1, UInt const nrowZero)
{
  cout << "rowUnify: OBSOLETE function ==> use CSRPatt constructor or MixedPattern class" << endl;
  UInt nnz=patt1._nnz;
  UInt nrows=patt1._nrows + nrowZero;
  UInt ncols=patt1._ncols;

  // nothing changes but the number of rows of the matrix
  CSRPatt ans(nnz, nrows, ncols);

  ans._ia.resize(nrows+1,nnz+PatternOffset); // points at the end !

  //copy of patt1._ia
  for (UInt i=0; i < patt1._ia.size(); i++)
    ans._ia[i] = patt1._ia[i];
  //copy of patt1._ja
  for (UInt i=0; i < patt1._ja.size(); i++)
    ans._ja.push_back(patt1._ja[i]);

  ans._filled = true;
  ans._diagfirst = patt1._diagfirst;
  return ans;
}

// zero on the top
CSRPatt rowUnify(UInt const nrowZero, CSRPatt const &patt1)
{
  cout << "rowUnify: OBSOLETE function ==> use CSRPatt constructor or MixedPattern class" << endl;
  typedef PatternDefs::Container::iterator Iter;
  UInt nnz=patt1._nnz;
  UInt nrows=patt1._nrows + nrowZero;
  UInt ncols=patt1._ncols;

  // nothing changes but the number of rows of the matrix
  CSRPatt ans(nnz, nrows, ncols);
  ans._ia.resize(nrows+1,PatternOffset); // points at the begining !

  //copy of patt1._ia
  for (UInt i=0; i < patt1._ia.size(); i++)
    ans._ia[i+patt1._nrows] = patt1._ia[i];
  //copy of patt1._ja
  for (UInt i=0; i < patt1._ja.size(); i++)
    ans._ja.push_back(patt1._ja[i]);

  ans._filled=true;
  ans._diagfirst=false;

  if (!patt1._diagfirst)
    return ans;

  // sort of col indexes
  for (UInt i=nrowZero; i < nrows; i++)
    {
      Iter ja_start = ans._ja.begin() + ans._i2o(patt1._ia[i]);
      Iter ja_end   = ans._ja.begin() + ans._i2o(patt1._ia[i+1]);
      sort(ja_start, ja_end);
    }

  return ans;
}

// construction of a block diagonal matrix
CSRPatt diagblockMatrix(CSRPatt const &patt, UInt const nblock)
{
    typedef PatternDefs::Container::const_iterator ConstIter;
    typedef PatternDefs::Container::iterator Iter;

    ASSERT(patt._nrows == patt._ncols, "Matrix must be square!");

    UInt nrowsblock = patt._nrows; // Numero di righe di un blocco
    UInt nnz = nblock*patt._nnz;
    UInt nrows = nblock*nrowsblock;
    UInt ncols = nblock*nrowsblock;

    CSRPatt ans(nnz, nrows, ncols);
    ans._ia.resize(nrows+1,0);
    ans._ja.resize(nnz,0);

    ConstIter ja1_start = patt._ja.begin();

    Iter ja_start = ans._ja.begin();

    UInt i, nrow = 0, brow = 0;

    for (UInt block = 1; block <= nblock; ++block)
    {
        for (i = 0; i < nrowsblock; ++i)
        {
            // current row number of global matrix
            nrow = i + nrowsblock*(block-1);

            ans._ia[nrow] = patt._ia[i] + brow - PatternOffset;

            Iter ja_start1 = ja_start + ans._i2o(ans._ia[nrow]);

            // construction of ja
            ConstIter ja1_start1 = ja1_start + patt._i2o(patt._ia[i]);
            ConstIter ja1_end1 = ja1_start + patt._i2o(patt._ia[i+1]);

            // copy of ja
            copy(ja1_start1, ja1_end1, ja_start1);

            Iter ja_end1 = ja_start1 + patt._ia[i+1] - patt._ia[i];

            for (Iter ip = ja_start1; ip != ja_end1; ++ip)
                *ip += nrowsblock*(block-1);
        }
        // row number of global matrix where block number block+1 begins
        brow = ans._ia[nrow] + patt._ia[i] - patt._ia[i-1] - PatternOffset;
    }

    ans._ia[nrow+1] = nnz;
    ans._filled=true;
    ans._diagfirst = patt._diagfirst;
    return ans;
}

////////////////////////////////////////////////////////////////////////
//
// V B R Pattern
//
////////////////////////////////////////////////////////////////////////

VBRPatt::VBRPatt(UInt ex_nnz, UInt ex_nrow, UInt ex_ncol):
  CSRPatt(ex_nnz, ex_nrow, ex_ncol)
{
  _indx.reserve(ex_nnz+1);
  _rpntr.reserve(ex_nrow+1);
  _cpntr.reserve(ex_ncol+1);
}

VBRPatt::VBRPatt(UInt ex_nnz, UInt ex_nrow, UInt ex_ncol, const
		 vector<Index_t> &ex_ia, const vector<Index_t> &ex_ja,
		 const vector<Index_t> &ex_indx, const vector<Index_t>
		 &ex_rpntr, const vector<Index_t> &ex_cpntr):
  CSRPatt(ex_nnz, ex_nrow, ex_ncol, ex_ia,ex_ja), _indx(ex_indx),
  _rpntr(ex_rpntr),_cpntr(ex_cpntr)
{
  // check on the compatibility of the external data
  ASSERT_PRE((ex_indx.size()==ex_nnz+1 && ex_rpntr.size()==ex_nrow+1
	      && ex_cpntr.size()==ex_ncol+1), "Error in VBR Pattern Life V");
}
//Alain (23/10/02)
// Copy constructor of CSR might be not necessary to call ?
// I will test it soon.
VBRPatt::VBRPatt(const VBRPatt &RightHandVBRP)://CSRPatt(RightHandVBRP),
  _indx(RightHandVBRP.indx()), _rpntr(RightHandVBRP.rpntr()),
  _cpntr(RightHandVBRP.cpntr())
{
  //CSRPatt::CSRPatt(RightHandVBRP);
}

bool VBRPatt::buildPattern(UInt const blockSize)
{
  bool test=(blockSize>0);
  if (!test) return test;

  _indx.resize(_nnz+1);
  _rpntr.resize(_nrows+1);
  _cpntr.resize(_ncols+1);

  // implementation for a fixed size block.
  // Remark: _rpntr=_cpntr if _nrows=_ncols
  for (UInt i=0; i<_nnz+1; ++i) _indx[i]= i*blockSize*blockSize;
  for (UInt i=0; i<_nrows+1; ++i) _rpntr[i]=i*blockSize;

  for (UInt i=0; i<_ncols+1; ++i) _cpntr[i]=i*blockSize;

  // taking into account of PatternOffset
  if ( PatternOffset != 0){
    for (Container::iterator ip=_indx.begin(); ip != _indx.end(); ++ip)
      *ip+=PatternOffset;
    for (Container::iterator ip=_rpntr.begin(); ip != _rpntr.end(); ++ip)
      *ip+=PatternOffset;
    for (Container::iterator ip=_cpntr.begin(); ip != _cpntr.end(); ++ip)
      *ip+=PatternOffset;
  }
  return true;
}
VBRPatt &VBRPatt::operator= (const VBRPatt& RhVbr){
  if (&RhVbr != this)
    {
      this->CSRPatt::operator=(RhVbr);
      _indx = RhVbr.indx();
      _rpntr = RhVbr.rpntr();
      _cpntr = RhVbr.cpntr();
    }
  return *this;
}


void  VBRPatt::showMe(bool verbose, ostream& c) const
{
  BasePattern::showMe(verbose,c);

  UInt blsize=_rpntr[1]-_rpntr[0]; // block size
  int i_first;
  string pare="[";
  c << "**************************" << endl;
  c << "     VBR Matrix Pattern   " << endl;
  c << endl;
  if (verbose){

    c << pare;
    for (Diff_t i_index=0; i_index < static_cast<Diff_t>(_nrows); i_index++)
      for (UInt jb=0;jb<blsize;jb++){
	c << pare;
	pare = " [";
	i_first=_i2o(_ia[i_index])+1; // In _ia[i_index] there is the
	                              // diagonal entry
	UInt jj=0;
	for(Diff_t j=0;j<static_cast<Diff_t>(_ncols);j++)
	  {
	    if (j==i_index) for (UInt ib=0;ib<blsize;ib++) c << " * ";
	    else {
	      if (j==_i2o(_ja[i_first+jj])){
		for (UInt ib=0;ib<blsize;ib++) c << " * ";
		jj++;
	      }
	      else for (UInt ib=0;ib<blsize;ib++) c << " 0 ";
	    }
	  }
	if (i_index==static_cast<Diff_t>(_nrows-1))
	  c << " ]] " << endl;
	else
	  c << " ]  " << endl;
	c<< endl;
      }
  }
  c << "**************************" << endl;
  return;
}

void VBRPatt::spy(string  const &filename) const
{
    // Purpose: Matlab dumping and spy
    string nome=filename, uti=" , ";
    UInt nblocrow=_nrows, blocsize=_rpntr[1]-_rpntr[0];
    //
    // check on the file name
    //
    int i=filename.find(".");

    if (i<=0) nome=filename+".m";
    else {
        if (i!=(int)filename.size()-2  || filename[i+1]!='m')
        {cerr << "Wrong file name ";
        nome=filename+".m";}
    }

    ofstream file_out(nome.c_str());
    ASSERT(file_out,"Error: Output Matrix (Values) file cannot be open");


    file_out << "S = [ ";
    for (UInt irb=0;irb<nblocrow;++irb){
        for (UInt ic=_ia[irb];ic<_ia[irb+1];++ic)
            for (UInt i=0;i<blocsize;++i)
                for (UInt j=0;j<blocsize;++j)
                    file_out << irb*blocsize+i-PatternOffset+1 << uti <<
                        _ja[ic]*blocsize+j-PatternOffset+1  << uti << "1.0" << endl;
    }
    file_out << "];" << endl;
    file_out << "I=S(:,1); J=S(:,2); S=S(:,3); A=sparse(I,J,S); spy(A);"<<endl;
}

////////////////////////////////////////////////////////////////////////
//
// C S R Symmetric Pattern
//
////////////////////////////////////////////////////////////////////////
CSRPattSymm::CSRPattSymm()
{
    // nothing to do here
}

CSRPattSymm::CSRPattSymm(UInt ex_nnz, UInt ex_nrow, UInt ex_ncol):
  BasePattern(ex_nnz,ex_nrow,ex_ncol)
{
  _ia.reserve(ex_nrow+1);
  _ja.reserve((ex_nnz+ex_nrow)/2); // Notice that in the symmetric pattern....
}

CSRPattSymm::CSRPattSymm(UInt ex_nnz, UInt ex_nrow, UInt ex_ncol,  const vector<Index_t> &ex_ia, const vector<Index_t> &ex_ja):
  BasePattern(ex_nnz,ex_nrow,ex_ncol),_ia(ex_ia),_ja(ex_ja)
{
  ASSERT_PRE(ex_ia.size()==ex_nrow+1,"Error in CSR Pattern Life V"); // check on the compatibility of the external data
  _filled=ex_ia.size()>0;
  _diagfirst=true;
}

CSRPattSymm::CSRPattSymm(const CSRPattSymm &RightHandCSRP):BasePattern(RightHandCSRP),
							   _ia(RightHandCSRP.ia()),_ja(RightHandCSRP.ja())
{
    // nothing to do here
}




bool CSRPattSymm::isThere(Index_t i,Index_t j) const
{
  if (i<=j)
    {
      ASSERT_BD(i>=PatternOffset && i<static_cast<Index_t>(_nrows)+PatternOffset);
      ASSERT_BD(j>=PatternOffset && j<static_cast<Index_t>(_ncols)+PatternOffset);

      Container::const_iterator start=_ja.begin()+_row_off(i);
      Container::const_iterator finish=_ja.begin()+_row_off(i+1);
      return binary_search(start, finish, j);
    }
  else return isThere(j,i);
}

pair<PatternDefs::Diff_t,bool>
CSRPattSymm::locate_pattern(Index_t const i,Index_t const j) const {
  ASSERT_BD(i>=PatternOffset && i<static_cast<Index_t>(_nrows)+PatternOffset);
  ASSERT_BD(j>=PatternOffset && j<static_cast<Index_t>(_ncols)+PatternOffset);

  if (! _filled ) {
    return make_pair(0,false);
  }
  if (i<=j){
    Diff_t _off=_row_off(i);
    Container::const_iterator finish=_ja.begin()+ _row_off(i+1);
    Container::const_iterator current = search_binary(_ja.begin()+ _off,finish, j); // search with a binary search
    return make_pair(current-_ja.begin(),current !=finish);
  } else {
    return locate_pattern(j,i);
  }
}


CSRPattSymm& CSRPattSymm::operator= (const CSRPattSymm& RhCsr){
  if (&RhCsr != this)
    {
      _nnz = RhCsr.nNz();
      _nrows = RhCsr.nRows();
      _ncols = RhCsr.nCols();
      _ia = RhCsr.ia();
      _ja = RhCsr.ja();
      _filled= !RhCsr.isEmpty();
      _diagfirst=RhCsr.diagFirst();
    }
  return *this;
}

UInt
CSRPattSymm::nbNeighbours(ID const d) const
{
  ASSERT_BD(d>0 && d<=_nrows);
  ASSERT_PRE(_filled, "Cannot access an empty pattern");

  UInt counter=_ia[d]-_ia[d-1]; // upper diagonal
  Index_t _ind=_d2i(d);
  for(UInt i=0; i<d-1;++i){
    if(binary_search(_ja.begin()+_row_off(i) +1,_ja.begin()+_row_off(i+1),_ind)) ++counter;
  }
  return counter;
}

ID
CSRPattSymm::neighbour(ID const n, ID const d) const{
  ASSERT_BD(d>0 && d<=_nrows);
  ASSERT_BD(n>0);
  ASSERT_PRE(_filled, "Cannot access an empty pattern");

  if (n==1) return d; // Diagonal is the first

  Index_t _ind=_d2i(d);
  Container::const_iterator start;
  Container::const_iterator finish;
  UInt counter=1;

  for(UInt i=0; i<_d2o(d);++i){
    start=_ja.begin()+(_i2o(_ia[i]) +1); // no need to search diag
    finish=_ja.begin()+_i2o(_ia[i+1]);
    if(binary_search(start, finish, _ind)) ++counter;
    if (counter==n) return i+1;// Number from 1
  }
  counter=n-counter;

  if (static_cast<Index_t>(counter)<_ia[d]-_ia[d-1]){
    return _i2d(_ja[_i2o(_ia[_d2o(d)]) +counter]);
  }else{
    return 0; // ERROR
  }
}

void
CSRPattSymm::neighbours(ID const d, Container & neighs) const
{
  ASSERT_BD(d>0 && d<=_nrows);
  ASSERT_PRE(_filled, "Cannot access an empty pattern");

  Container::const_iterator finish;
  Container::const_iterator start;
  //neighs.clear();
  neighs.reserve(nbNeighbours(d));
  neighs.push_back(d);// Diagonal first
  Index_t _row=_d2i(d);

  for(UInt i=0; i<d-1;++i){
    start=_ja.begin()+_i2o(_ia[i]) +1; // no need to search diag
    finish=_ja.begin()+_i2o(_ia[i+1]);
    if(binary_search(start, finish, _row))neighs.push_back(i+1);
  }
  for(start=_ja.begin()+_row_off(_row)+1;start!=_ja.begin()+_row_off(_row+1);++start)neighs.push_back(_i2d(*start));
}

void  CSRPattSymm::showMe(bool verbose, ostream& c) const{
  BasePattern::showMe(verbose,c);

  typedef vector<Index_t>::iterator found;
  int i_first;
  string pare="[";
  c << "********************************" << endl;
  c << "  CSR Matrix Symmetric Pattern   " << endl;
  c << endl;
  c << pare;
  if(verbose){
    for (UInt i_index=0; i_index < _nrows; i_index++){
      c << pare;
      pare = " [";
      i_first=_ia[i_index]+1 - PatternOffset; // In _ia[i_index] there is the diagonal entry
      for(unsigned int j=0;j<_ncols;j++){
	if (j == i_index) c << " * ";
	else {
	  if (isThere(i_index,j)==true){
	    c << " * ";}
	  else
	    c << " 0 ";
	}
      }
      if (i_index==_nrows-1)
	c << " ]] " << endl;
      else
	c << " ]  " << endl;
    }
  }
  c << "********************************" << endl;
  return;
}

void  CSRPattSymm::spy(string  const &filename) const{
  // Purpose: Matlab dumping and spy
  string nome=filename, uti=" , ";
  //
  // check on the file name
  //
  unsigned int i=filename.find(".");


  if (i<=0) nome=filename+".m";
  else {
    if (i!=filename.size()-2  || filename[i+1]!='m'){
      cerr << "Wrong file name " << i << endl;
      nome=filename+".m";}
  }

  ofstream file_out(nome.c_str());


  file_out << "S = [ ";
  for (Diff_t i=0;i<_nrows;++i){
    for (Index_t ii=_ia[i]-PatternOffset;ii<_ia[i+1]-PatternOffset;++ii){
      file_out << i+1 << uti << _ja[ii]+1-PatternOffset << uti << "1.0" << endl;
      if (i!=_i2o(_ja[ii])) file_out <<  _i2o(_ja[ii]+1)<< uti << i+1 << uti << "1.0" << endl; }
  }
  file_out << "];" << endl;

  file_out << "I=S(:,1); J=S(:,2); S=S(:,3); A=sparse(I,J,S); spy(A);"<<endl;
}
/////////////////////////////
////////////////////////////////////////////////////////////////////////
//
// M S R Pattern
//
////////////////////////////////////////////////////////////////////////
MSRPatt::MSRPatt()
{
    // nothing to do here
}

MSRPatt::MSRPatt(UInt ex_nnz, UInt ex_nrow, UInt ex_ncol):
  BasePattern(ex_nnz,ex_nrow,ex_ncol)
{
  _bindx.reserve(ex_nnz+1);
  _diagfirst=true;// default for MSR
  _ybind.reserve(ex_nnz-ex_nrow);
}
MSRPatt::MSRPatt(UInt ex_nnz, UInt ex_nrow, UInt ex_ncol, const vector<Index_t> &ex_bindx, const vector<Index_t> &ex_ybind):
  BasePattern(ex_nnz,ex_nrow,ex_ncol),_bindx(ex_bindx),_ybind(ex_ybind)
{
    ASSERT_PRE(ex_bindx.size()==ex_nnz+1,"Compatibility error in MSR Pattern "); // check on the compatibility of the external data
    _filled=!ex_bindx.empty();
    _diagfirst=true; // default for MSR
}

MSRPatt::MSRPatt(const MSRPatt &RightHandMSRP):BasePattern(RightHandMSRP),
					       _bindx(RightHandMSRP.bindx()),_ybind(RightHandMSRP.ybind())
{
    // nothing to do here
}

// TODO: It would be useful to write some convertors, not only as constructors
//       but also as functions

// version assuming that the CSR pattern has all its diagonal terms non null
MSRPatt::MSRPatt(const CSRPatt &RightHandCSRP):BasePattern(RightHandCSRP),
					       _bindx(RightHandCSRP.nNz()+1)
{
    PatternDefs::Container::const_iterator ia=RightHandCSRP.give_ia().begin();
    PatternDefs::Container::const_iterator ja=RightHandCSRP.give_ja().begin();

    Diff_t jj=0;
    _bindx[0] = _nrows+1;
    for (UInt i=1; i<=_nrows; ++i)
    {
        Index_t ifirst = *(ia+ i-1);
        Index_t ilast = *(ia + i);

        _bindx[i] = _bindx[i-1] + ilast - ifirst -1; // Pay attention to the circumstance that only off-diagonal elements are counted

        for (Index_t j=ifirst; j < ilast;++j)
        {
            Index_t jjj = *(ja + _i2o(j));
            if ((i-1)!=_i2o(jjj))
                _bindx[_nrows+1+jj++] = jjj;
        }
    }

    //construction of ybind
    _ybind.resize(RightHandCSRP.nNz()-_nrows,0);
    //temporary bindx
    PatternDefs::Container bindAux(_bindx);

    jj=0;
    for (UInt i=0; i<_nrows; ++i)
    {
        Index_t ifirst= _bindx[i];
        Index_t ilast= _bindx[i+1];
        for (Index_t j=ifirst; j < ilast;++j)
            _ybind[jj++]= bindAux[_bindx[j]]++;
    }

    _diagfirst=true;// default for MSR
    return;
}


MSRPatt & MSRPatt::operator= (const MSRPatt& RhMsr){
  if (&RhMsr != this)
    {
      _nnz = RhMsr.nNz();
      _nrows= RhMsr.nRows();
      _ncols= RhMsr.nCols();
      _bindx = RhMsr.bindx();
      _ybind = RhMsr.ybind();
      _filled= !  RhMsr.isEmpty();
      _diagfirst=true;// default for MSR

    }
  return *this;
}

void MSRPatt::neighbours(ID const d, Container & neigh) const
{
  ASSERT_BD(d>0 && d<=_nrows);
  ASSERT_PRE(_filled, "Cannot access an empty pattern");
  //neigh.clear();
  neigh.reserve(nbNeighbours(d));
  neigh.push_back(d); // diagonal, which is NOT explicitely stored.
  for(Container::const_iterator start1=_bindx.begin()+_i2o(_bindx[d-1]);
      start1!=_bindx.begin()+_i2o(_bindx[d]);++start1)neigh.push_back(_i2d(*start1));
}
// locate function for CSR Pattern. It returns a pair. First member
// is the position in the array correponding to (i,j), the second
// is a bool telling if that position exists (i.e. if i,j is in the
//pattern). If the boolean value is false the first member is meaningless!
pair<PatternDefs::Diff_t,bool>
MSRPatt::locate_pattern(Index_t const i,Index_t const j)const
{
  if (i==j)
    { return make_pair(_i2o(i),true);} // In MSR Format the diagonal entries are alway part of the pattern
  else
    {
      Container::const_iterator start=_bindx.begin()+_row_off(i); // the real start
      Container::const_iterator finish=_bindx.begin()+_row_off(i+1); //the real end (remember STL convention for ranges!)
      Container::const_iterator current= search_binary(start, finish, j); // the off-diagonal terms have been ordered
      // difference of pointers should  return distance, which is of integral type
      return make_pair(current-_bindx.begin(),current !=finish);
    }
}

void  MSRPatt::showMe(bool verbose, ostream& c) const{
  unsigned int i_first;
  string pare="[";
  BasePattern::showMe(verbose, c);
  cout << "**************************" << endl;
  cout << "     MSR Matrix Pattern   " << endl;
  cout << endl;
  if (verbose){
    cout << pare;
    for (UInt i_index=0; i_index < _nrows; i_index++){
      cout << pare;
      pare = " [";
      i_first=_bindx[i_index];
      UInt jj=0;
      for(unsigned int j=0;j<_ncols;j++){
	if (j==i_index)
	  cout << " * ";
	else{
	  if (j==_i2o(_bindx[i_first+jj])){
	    c << " * "; jj++;}
	  else
	    c << " 0 ";
	}
      }
      if (i_index==_nrows-1)
	cout << " ]]; " << endl;
      else
	cout << " ]  " << endl;
    }
  }
  return;
}

void  MSRPatt::spy(string  const &filename)const
{
  // Purpose: Matlab dumping and spy
  string nome=filename, uti=" , ";
  //
  // check on the file name
  //
  unsigned int i=filename.find(".");

  if (i<=0) nome=filename+".m";
  else {
    if (i!=filename.size()-2  || filename[i+1]!='m'){
      cerr << "Wrong file name ";
      nome=filename+".m";}
  }

  ofstream file_out(nome.c_str());


  file_out << "S = [ ";
  for (UInt i=0;i<_nrows;++i){
    if (i < _ncols)
      file_out << i+1 << uti << i+1 << uti << "1.0" << endl;
    for (Index_t ii=_bindx[i];ii<_bindx[i+1];++ii)
      file_out << i+1 << uti << _bindx[ii]+1-PatternOffset << uti << "1.0" << endl;
  }
  file_out << "];" << endl;

  file_out << "I=S(:,1); J=S(:,2); S=S(:,3); A=sparse(I,J,S); spy(A);"<<endl;
}

// Construction of diagonal block matrix. Done by A. Gilardi.
// Alain (nov. 2002), update of ybind as well !
void diagblockMatrix(MSRPatt &ans, MSRPatt const &patt, UInt const nblock)
{
  ASSERT(patt._nrows == patt._ncols, "Matrix must be square!");

  UInt nrowsblock = patt._nrows; // Number of rows for one block
  UInt nnz = nblock*patt._nnz;
  UInt nrows = nblock*nrowsblock;
  UInt ncols = nblock*nrowsblock;

  ans._nnz = nnz;
  ans._nrows = nrows;
  ans._ncols = ncols;
  ans._bindx.resize(nnz+1,0);
  ans._ybind.resize(nnz-nrows);

  UInt pos = nrows + 1;
  UInt row,coloffset, valoffset;
  UInt block, i, ii;

  for (block = 0; block < nblock; ++block)
    {
      for (i = 0; i < nrowsblock; ++i)
	{
	  row = i + block*nrowsblock;
	  ans._bindx[row] = pos;

	  coloffset = block*patt._ncols;
	  // offset for ybind
	  valoffset = nrows-nrowsblock+block*(patt._nnz-nrowsblock);

	  for (ii = patt._bindx[i]; ii < patt._bindx[i+1]; ++ii)
	    {
	      ans._bindx[pos] = patt._bindx[ii] + coloffset;
	      // construction of ybind
	      ans._ybind[pos-nrows-1]= patt._ybind[ii-nrowsblock-1]
		+ valoffset;
	      ++pos;
	    }
	}
    }

  ans._bindx[nrows] = nnz+1;

  ans._filled=true;
  ans._diagfirst = true; //default for MSR
}
}
