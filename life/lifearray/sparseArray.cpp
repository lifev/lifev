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
#include "sparseArray.hpp"

namespace LifeV
{
//Realizza l'inversa di una matrice diagonale e la moltiplica per un'altra matrice
void MultInvDiag(const vector<Real> &Diag,
		 const CSRMatr<CSRPatt,Real> &Mat, CSRMatr<CSRPatt,Real> &ans)
{
  ASSERT(find(Diag.begin(),Diag.end(),0) == Diag.end(), "La matrice diagonale deve essere invertibile");

  ASSERT(Diag.size() == Mat._Patt->nRows(), "Le matrice non sono compatibili per poter calcolare il loro prodotto");

  // Esecuzione del prodotto
  for (UInt row = 0; row < Mat._Patt->nRows(); ++row)
    for (UInt pos = Mat._Patt->ia()[row]; pos < Mat._Patt->ia()[row+1]; ++pos)
      ans._value[pos] = (Mat._value[pos] / Diag[row]);

  return;
}

//version for Datatype=Tab2d
CSRMatr<CSRPatt,Tab2d>::
CSRMatr(const CSRPatt &ex_pattern, UInt const nr, UInt const nc)
{
  _Patt = &ex_pattern;
  Tab2d mzero(nr,nc);
  mzero= 0.;
  _value.resize(ex_pattern.nNz(), mzero);
}

// version for block matrices
VectorBlock
CSRMatr<CSRPatt,Tab2d>::
trans_mult(const VectorBlock &v)
{
  UInt nblockr=_Patt->nRows(); // for square matrices...
  int blsize=_value[0].N(); // for square block

  ASSERT(nblockr==v.size(),"Error in Matrix Vector product");
  VectorBlock ans(nblockr,blsize);
  ans=0.0;

  for (UInt ir=0+OFFSET;ir<nblockr+OFFSET;++ir)
    {
      for (UInt ii=_Patt->ia()[ir]-OFFSET;ii<_Patt->ia()[ir+1]-OFFSET;++ii)
	{
	  Tab2dView tr_block(_value[ii].t()); // transpose-block
	  ans.numBlock(_Patt->ja()[ii]-OFFSET)+=(tr_block*v.numBlock(ir));
	};
    };
  return ans;
}

// version for block matrices
VectorBlock CSRMatr<CSRPatt,Tab2d>::
operator*(const VectorBlock &v) const
{
  UInt nblockr=_Patt->nRows();
  UInt nblockc=_Patt->nCols();
  int blsize=_value[0].N(); // for square block

  ASSERT(nblockc==v.size(),"Error in Matrix Vector product");
  VectorBlock ans(nblockr,blsize);
  ans=0.0;

  for (UInt ir=0+OFFSET;ir<nblockr+OFFSET;++ir)
    {
      for (UInt ii=_Patt->give_ia()[ir]-OFFSET;ii<_Patt->give_ia()[ir+1]-OFFSET;++ii)
	ans.numBlock(ir)+=(_value[ii]*v.numBlock(_Patt->give_ja()[ii]-OFFSET));
    };
  return ans;
}

// the case of block matrices with Tab2d block type.
void CSRMatr<CSRPatt,Tab2d>::
spy(string  const &filename)
{
  // Purpose: Matlab dumping and spy
  string nome=filename, uti=" , ";
  UInt nrows=_Patt->nRows();
  Container ia=_Patt->ia(), ja=_Patt->ja();
  int blsize=_value[0].N(); // for square block
  //
  // check on the file name
  //
  UInt i=filename.find(".");

  if (i<=0) nome=filename+".m";
  else {
    if (i!=filename.size()-2  || filename[i+1]!='m')
      {cerr << "Wrong file name ";
      nome=filename+".m";}
  }

  ofstream file_out(nome.c_str());
  ASSERT(file_out,"Error: Output Matrix (Values) file cannot be open");


  file_out << "S = [ ";
  for (UInt rblock=0; rblock<nrows; ++rblock){
    for (UInt ii=ia[rblock]-OFFSET; ii<ia[rblock+1]-OFFSET; ++ii){
      UInt cblock=ja[ii]-OFFSET;
      for (int rloc=0; rloc<blsize; ++rloc)
	for (int cloc=0; cloc<blsize; ++cloc)
	  file_out << rblock*blsize + rloc +1 << uti << cblock*blsize
	    + cloc +1 << uti << _value[ii](rloc,cloc) << endl; /* */
    }
  }
  file_out << "];" << endl;

  file_out << "I=S(:,1); J=S(:,2); S=S(:,3); A=sparse(I,J,S); spy(A);"<<endl;
};
//version without using static (I think it is better)
// Modified by A. Gilardi. 03/02.
void colUnify(CSRMatr<CSRPatt,double> &ans, const CSRMatr<CSRPatt,double> &Mat1,
	      const CSRMatr<CSRPatt,double> &Mat2)
{
  typedef vector<double>::const_iterator ConstIter;
  typedef vector<double>::iterator Iter;
  typedef PatternDefs::Container::const_iterator PaConstIter;
  typedef PatternDefs::Container::iterator PaIter;

  // the pattern is assumed having been constructed first
  ASSERT(ans._Patt->nNz()==Mat1._Patt->nNz()+Mat2._Patt->nNz(),
	 "pattern concatenation must be done first");

  // concatenation of the matrix values
  vector<double> val(ans._Patt->nNz(),0.0);

  UInt end = 0;

  if (Mat1._Patt->nRows() > Mat1._Patt->nCols())
    {
      if (Mat2._Patt->nCols() <= (Mat1._Patt->nRows() - Mat1._Patt->nCols()))
	end = Mat2._Patt->nCols() + Mat1._Patt->nCols();
      else
	end = Mat1._Patt->nRows();
    }

  if (!Mat2._Patt->diagFirst())
    {
      for (UInt i = 0; i < ans._Patt->nRows(); ++i)
	{
	  ConstIter value1_start = Mat1._value.begin() + Mat1._Patt->give_ia()[i]-OFFSET;
	  ConstIter value2_start = Mat2._value.begin() + Mat2._Patt->give_ia()[i]-OFFSET;
	  ConstIter value1_end = Mat1._value.begin() + Mat1._Patt->give_ia()[i+1]-OFFSET;
	  ConstIter value2_end = Mat2._value.begin() + Mat2._Patt->give_ia()[i+1]-OFFSET;

	  Iter value_start1 = val.begin() + ans._Patt->give_ia()[i];

	  // copy of the first block
	  copy(value1_start, value1_end, value_start1);

	  //the starting place of the second block
	  Iter value_start2 = value_start1 + (Mat1._Patt->give_ia()[i+1]-Mat1._Patt->give_ia()[i]);

	  // copy of the second block
	  copy(value2_start, value2_end, value_start2);
	};

    }
  else if (Mat1._Patt->diagFirst() == ans._Patt->diagFirst())
    {
      for (UInt i = 0; i < ans._Patt->nRows(); ++i)
	{
	  ConstIter value1_start = Mat1._value.begin() + Mat1._Patt->give_ia()[i]-OFFSET;
	  ConstIter value2_start = Mat2._value.begin() + Mat2._Patt->give_ia()[i]-OFFSET;
	  ConstIter value1_end = Mat1._value.begin() + Mat1._Patt->give_ia()[i+1]-OFFSET;
	  ConstIter value2_end = Mat2._value.begin() + Mat2._Patt->give_ia()[i+1]-OFFSET;

	  PatternDefs::Container ja = Mat2._Patt->ja();
	  PaConstIter ja2_start = ja.begin() + Mat2._Patt->_i2o(Mat2._Patt->give_ia()[i]);
	  PaConstIter ja2_end = ja.begin() + Mat2._Patt->_i2o(Mat2._Patt->give_ia()[i+1]);

	  //the starting place of the first block
	  Iter value_start1 = val.begin() + ans._Patt->_i2o(ans._Patt->give_ia()[i]);

	  // copy of the first block
	  copy(value1_start, value1_end, value_start1);

	  //the starting place of the second block
	  Iter value_start2_part1 = value_start1 + Mat1._Patt->give_ia()[i+1] - Mat1._Patt->give_ia()[i];

	  PaConstIter ja_pos_diag = find_if(ja2_start, ja2_end, bind2nd(greater<UInt>(),*ja2_start));

	  ConstIter value2_pos_diag = value2_start + (ja_pos_diag - ja2_start);

	  // copy of the second block
	  Iter value_start2_part2 = rotate_copy(value2_start, value2_start + 1, value2_pos_diag, value_start2_part1);
	  copy(value2_pos_diag, value2_end, value_start2_part2);
	};
    }
  else
    {
      for (UInt i = 0; i < ans._Patt->nRows(); i++)
	{
	  ConstIter value1_start = Mat1._value.begin() + Mat1._Patt->give_ia()[i]-OFFSET;
	  ConstIter value2_start = Mat2._value.begin() + Mat2._Patt->give_ia()[i]-OFFSET;
	  ConstIter value1_end = Mat1._value.begin() + Mat1._Patt->give_ia()[i+1]-OFFSET;
	  ConstIter value2_end = Mat2._value.begin() + Mat2._Patt->give_ia()[i+1]-OFFSET;

	  PatternDefs::Container ja = Mat1._Patt->ja();
	  PaConstIter ja1_start = ja.begin() + Mat1._Patt->_i2o(Mat1._Patt->give_ia()[i]);
	  PaConstIter ja1_end = ja.begin() + Mat1._Patt->_i2o(Mat1._Patt->give_ia()[i+1]);

	  //the starting place of the first block
	  Iter value_start1_part1 = val.begin() + ans._Patt->_i2o(ans._Patt->give_ia()[i]);

	  PaConstIter ja_pos_diag = find_if(ja1_start, ja1_end, bind2nd(greater<UInt>(),*ja1_start));

	  ConstIter value_pos_diag = value1_start + (ja_pos_diag - ja1_start);

	  // copy of the first block
	  Iter value_start1_part2 = rotate_copy(value1_start, value1_start + 1, value_pos_diag, value_start1_part1);
	  copy(value_pos_diag, value1_end, value_start1_part2);

	  ja = Mat2._Patt->ja();
	  PaConstIter ja2_start = ja.begin() + Mat1._Patt->_i2o(Mat2._Patt->give_ia()[i]);
	  PaConstIter ja2_end = ja.begin() + Mat1._Patt->_i2o(Mat2._Patt->give_ia()[i+1]);

	  //the starting place of the second block
	  Iter value_start2_part1 = value_start1_part1 + Mat1._Patt->give_ia()[i+1] - Mat1._Patt->give_ia()[i];

	  ja_pos_diag = find_if(ja2_start, ja2_end, bind2nd(greater<UInt>(),*ja2_start));

	  value_pos_diag = value2_start + (ja_pos_diag - ja2_start);

	  // copy of the second block
	  Iter value_start2_part2 = rotate_copy(value2_start, value2_start + 1, value_pos_diag, value_start2_part1);
	  copy(value_pos_diag, value2_end, value_start2_part2);
	};
    }

  if (!ans._Patt->diagFirst() || !end)
    {
      ans._value.resize(val.size(),0.0);
      ans._value = val;
      return;
    }

  for (UInt i = Mat1._Patt->nCols(); i < end; ++i)
    {
      Iter value_start = val.begin() + ans._Patt->give_ia()[i] - OFFSET;

      PatternDefs::Container ja = ans._Patt->ja();
      PaIter ja_start = ja.begin() + ans._Patt->_i2o(ans._Patt->give_ia()[i]);
      PaIter ja_end = ja.begin() + ans._Patt->_i2o(ans._Patt->give_ia()[i+1]);

      PaIter ja_pos_diag = find_if(ja_start, ja_end, bind2nd(greater<UInt>(),*ja_start));

      Iter value_pos_diag = value_start + (ja_pos_diag - ja_start) - 1;

      // ordinamento
      rotate(value_start, value_pos_diag, value_pos_diag + 1);

    };

  ans._value.resize(val.size(),0.0);
  ans._value = val;
}
// row-concatenation of two blocks of CSRMatr
// Modified by A. Gilardi. 03/02.
void rowUnify(CSRMatr<CSRPatt,double> &ans, const CSRMatr<CSRPatt,double> &Mat1,
	      const CSRMatr<CSRPatt,double> &Mat2)
{
  typedef vector<double>::const_iterator ConstIter;
  typedef vector<double>::iterator Iter;
  typedef PatternDefs::Container::const_iterator PaConstIter;
  typedef PatternDefs::Container::iterator PaIter;

  // the pattern is assumed having been constructed first
  ASSERT(ans._Patt->nNz()==Mat1._Patt->nNz()+Mat2._Patt->nNz(),"pattern concatenation must be done first");

  // concatenation of the matrix values
  vector<double> val(ans._Patt->nNz(),0.0);

  bool diag = ans._Patt->diagFirst();

  UInt end = 0;

  if (diag)
    if (Mat1._Patt->nRows() < Mat1._Patt->nCols())
      {
	if (Mat2._Patt->nRows() <= (Mat1._Patt->nCols() - Mat1._Patt->nRows()))
	  end = Mat2._Patt->nRows() + Mat1._Patt->nRows();
	else
	  end = Mat1._Patt->nCols();
      }

  if (!diag && Mat1._Patt->diagFirst())
    {
      for (UInt i = 0; i < Mat1._Patt->nRows(); i++)
	{
	  ConstIter value1_start = Mat1._value.begin() + Mat1._Patt->give_ia()[i] - OFFSET;
	  ConstIter value1_end = Mat1._value.begin() + Mat1._Patt->give_ia()[i+1] - OFFSET;

	  PatternDefs::Container ja = Mat1._Patt->ja();
	  PaConstIter ja1_start = ja.begin() + Mat1._Patt->_i2o(Mat1._Patt->give_ia()[i]);
	  PaConstIter ja1_end = ja.begin() + Mat1._Patt->_i2o(Mat1._Patt->give_ia()[i+1]);

	  //the starting place of the first block
	  Iter value_start1_part1 = val.begin() + ans._Patt->_i2o(ans._Patt->give_ia()[i]);

	  PaConstIter ja_pos_diag = find_if(ja1_start, ja1_end, bind2nd(greater<UInt>(),*ja1_start));

	  ConstIter value_pos_diag = value1_start + (ja_pos_diag - ja1_start);

	  // copy of the first block
	  Iter value_start1_part2 = rotate_copy(value1_start, value1_start + 1, value_pos_diag, value_start1_part1);
	  copy(value_pos_diag, value1_end, value_start1_part2);
	};
    }
  else
    {
      ConstIter value1_start = Mat1._value.begin();
      ConstIter value1_end = Mat1._value.end();

      //the starting place of the first block
      Iter value_start1 = val.begin();

      // copy of the first block
      copy(value1_start, value1_end, value_start1);
    }

  //the starting place of the first block
  Iter value_start2 = val.begin() + Mat1._Patt->nNz();

  if (Mat2._Patt->diagFirst())
    {
      for (UInt i = 0; i < Mat2._Patt->nRows(); i++)
	{
	  ConstIter value2_start = Mat2._value.begin() + Mat2._Patt->give_ia()[i] - OFFSET;
	  ConstIter value2_end = Mat2._value.begin() + Mat2._Patt->give_ia()[i+1] - OFFSET;

	  PatternDefs::Container ja = Mat2._Patt->ja();
	  PaConstIter ja2_start = ja.begin() + Mat2._Patt->_i2o(Mat2._Patt->give_ia()[i]);
	  PaConstIter ja2_end = ja.begin() + Mat2._Patt->_i2o(Mat2._Patt->give_ia()[i+1]);

	  Iter value_start2_part1 = value_start2 + Mat2._Patt->give_ia()[i] - OFFSET;

	  PaConstIter ja_pos_diag = find_if(ja2_start, ja2_end, bind2nd(greater<UInt>(),*ja2_start));

	  ConstIter value_pos_diag = value2_start + (ja_pos_diag - ja2_start);

	  // copy of the second block
	  Iter value_start2_part2 = rotate_copy(value2_start, value2_start + 1, value_pos_diag, value_start2_part1);
	  copy(value_pos_diag, value2_end, value_start2_part2);
	};
    }
  else
    {
      ConstIter value2_start = Mat2._value.begin();
      ConstIter value2_end = Mat2._value.end();

      // copy of the second block
      copy(value2_start, value2_end, value_start2);
    }

  if (!diag)
    {
      ans._value.resize(val.size(),0.0);
      ans._value = val;
      return;
    }

  for (UInt i = Mat1._Patt->nRows(); i < end; ++i)
    {
      Iter value_start = val.begin() + ans._Patt->give_ia()[i] - OFFSET;

      PatternDefs::Container ja = ans._Patt->ja();
      PaIter ja_start = ja.begin() + ans._Patt->_i2o(ans._Patt->give_ia()[i]);
      PaIter ja_end = ja.begin() + ans._Patt->_i2o(ans._Patt->give_ia()[i+1]);

      PaIter ja_pos_diag = find_if(ja_start, ja_end, bind2nd(greater<UInt>(),*ja_start));

      Iter value_pos_diag = value_start + (ja_pos_diag - ja_start) - 1;

      // ordinamento
      rotate(value_start, value_pos_diag, value_pos_diag + 1);
    };

  ans._value.resize(val.size(),0.0);
  ans._value = val;
}

//Inversion of a diagonal matrix and multiplication with a sparse matrix
//Miguel: 4/2003, the last version was too expensive.
void MultInvDiag(const vector<Real> &Diag, const MSRMatr<Real> &Mat, MSRMatr<Real> &ans)
{
  ASSERT(find(Diag.begin(),Diag.end(),0) == Diag.end(), "La matrice diagonale deve essere invertibile");

  ASSERT(Diag.size() == Mat._Patt->nRows(), "Le matrice non sono compatibili per poter calcolare il loro prodotto");
  UInt begin, end;
  Real diag;
  UInt nRows= Mat._Patt->nRows();
  vector<Index_t>::const_iterator start = Mat._Patt->give_bindx().begin();

  for (UInt row = 0; row < nRows; ++row)
    {
      begin = *(start+row);
      end   = *(start+row+1);
      diag  = 1.0/Diag[row];
      ans._value[row] = Mat._value[row]*diag;
      for (UInt pos = begin; pos < end; ++pos)
	ans._value[pos] = Mat._value[pos]*diag;
    }
}

double nihil(double val){
  return 0.;
}

//-----------------------------------------------------------------------
// DiagPreconditioner
//-----------------------------------------------------------------------

//for CSR or MSR normal pattern
DiagPreconditioner<Vector>::
DiagPreconditioner(const CSRMatr<CSRPatt,double> &M)
{
  double loc_val=0.0;
  UInt M_size=M.Patt()->nRows();
  Vector v_zero(M_size);
  v_zero=0.;

  _diag=v_zero;
  for (UInt i=0; i < M_size; i++)
    {
      loc_val = M.get_value(i,i);
      ASSERT(loc_val != 0. , "zero detected in diagonal preconditioner");
      _diag(i) = 1./loc_val ;
    };
}
DiagPreconditioner<Vector>::
DiagPreconditioner(const MSRMatr<double> &M)
{
  double loc_val=0.0;
  UInt M_size=M.Patt()->nRows();
  Vector v_zero(M_size);
  v_zero=0.;

  _diag=v_zero;
  for (UInt i=0; i < M_size; i++)
    {
      loc_val=M.get_value(i,i);
      ASSERT(loc_val != 0. , "zero detected in diagonal preconditioner");
      _diag(i) = 1./loc_val ;
    };
}
//for VBR pattern
DiagPreconditioner<Vector>::
DiagPreconditioner(const VBRMatr<double> &M)
{
  double loc_val=0.0;
  UInt Nblocks=M.Patt()->nRows();
  UInt blockSize = M.Patt()->rpntr()[1]-M.Patt()->rpntr()[0];
  UInt M_size=Nblocks*blockSize;
  Vector v_zero(M_size);
  v_zero=0.;

  _diag=v_zero;
  for (UInt i=0; i < M_size; i++)
    {
      loc_val=M.get_value(i,i);
      ASSERT(loc_val != 0. , "zero detected in diagonal preconditioner");
      _diag(i) = 1./loc_val ;
    };
}

//for CSR block pattern
DiagPreconditioner<VectorBlock>::
DiagPreconditioner(const CSRMatr<CSRPatt,Tab2d> &M)
{
  UInt Nblocks=M.Patt()->nRows();
  int blockSize =M.value()[0].N();
  Tab2d loc_val(blockSize,blockSize);
  loc_val=0.;
  VectorBlock v_zero(Nblocks,blockSize);
  v_zero=0.;

  _diag=v_zero;
  for (UInt i=0; i < Nblocks; i++)
    {
      loc_val=M.get_value(i,i);
      for (int j=0; j < blockSize; j++)
	{
	  ASSERT(loc_val(j,j) != 0. ,
		 "zero detected in diagonal preconditioner");
	  _diag.numBlock(i)(j)=1./loc_val(j,j);
	};
    };
}
//solve the diagonal system
Vector
DiagPreconditioner<Vector>::
solve(const Vector &x) const
{
  Vector y(x.size());

  for (UInt i=0; i<x.size(); i++)
    y(i) = x(i) * diag(i);

  return y;
}
VectorBlock
DiagPreconditioner<VectorBlock>::
solve(const VectorBlock &x) const
{
  VectorBlock y(x.size(),x.numBlock(0).N());

  for (UInt i=0; i<x.size(); i++)
    for (int j=0; j<x.numBlock(0).N(); j++)
      y.numBlock(i)(j) = x.numBlock(i)(j) * diagBlock(i)(j);

  return y;
}
//-----------------------------------------------------------------------
// DiagPreconditioner
//-----------------------------------------------------------------------

//for CSR or MSR normal pattern
IDPreconditioner<Vector>::
IDPreconditioner(const CSRMatr<CSRPatt,double> &M)
{
  UInt M_size=M.Patt()->nRows();
  Vector v_id(M_size);
  v_id=1.;
  _diag=v_id;
}

IDPreconditioner<Vector>::
IDPreconditioner(const MSRMatr<double> &M)
{
  UInt M_size=M.Patt()->nRows();
  Vector v_id(M_size);
  v_id=1.;
  _diag=v_id;
}

//for VBR pattern
IDPreconditioner<Vector>::
IDPreconditioner(const VBRMatr<double> &M)
{
  UInt Nblocks=M.Patt()->nRows();
  UInt blockSize = M.Patt()->rpntr()[1]-M.Patt()->rpntr()[0];
  UInt M_size=Nblocks*blockSize;
  Vector v_id(M_size);
  v_id=1.;
  _diag=v_id;
}
//for CSR block pattern
IDPreconditioner<VectorBlock>::
IDPreconditioner(const CSRMatr<CSRPatt,Tab2d> &M)
{
  int Nblocks=M.Patt()->nRows();
  int blockSize =M.value()[0].N();
  VectorBlock v_id(Nblocks,blockSize);
  v_id=1.;
  _diag=v_id;
}

//solve the diagonal system
Vector
IDPreconditioner<Vector>::
solve(const Vector &x) const
{return x;}

VectorBlock
IDPreconditioner<VectorBlock>::
solve(const VectorBlock &x) const
{return x;}

}
