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
/* --------------------------------------------------------------------------*
/                                                                            /
/      ...                                                                   /
/                                                                            /
/                                                                            /
/ MATRIX VALUES                                                              /
/                                                                            /
/ 21/06/2000 Alessandro Veneziani                                            /
/                                                                            /
/ #Purpose: Provides the basic classes for sparse matrix handling            /
/ including common formats (CSR,MSR),                                        /
/ assembling in the LifeV environment and utilities                          /
/                                                                            /
/ #Note: It will be suitably modified and linked with                        /
/ yAMLib (by Bertolazzi, Formaggia, Manzini, Veneziani)                      /
/ and its sparse environment                                                 /
/                                                                            /
/ #attempt to define the VBR matrix class in order to handle vectorial       /
/  problems. 26/10/01, Alain Gauthier.                                       /
/                                                                            /
/ #modification of CSR matrix: _value can be a vector of matrices, each      /
/  matrix being a RNM class. Attempt of 15/11/01.                            /
/                                                                            /
/ #class used in IML++: construction of a DiagPreconditioner class.          /
/  21/11/01.                                                                 /
/                                                                            /
/ #definition of MixedMatr class: holds the values of block matrix           /
/  stored under a MixedPattern format. 30/04/02.                             /
/                                                                            /
/---------------------------------------------------------------------------*/
#ifndef _SPARSE_ARRAY_HH
#define _SPARSE_ARRAY_HH
#ifndef OFFSET
#define OFFSET 0 // for the Fortran vs C numbering
#endif
#include "lifeV.hpp"
#ifndef _LIFEV_HH_
//more correct version
typedef size_t  UInt;
typedef vector<UInt>::iterator UIIter;
#endif
#include <fstream>
#include<set>
#include<algorithm>
#include<string>
#include<utility>
#include "pattern.hpp"
#ifndef _VEC_UNKNOWN_HH
#include "vecUnknown.hpp"
#endif

namespace LifeV
{
//more correct version
typedef vector<INDEX_T> Container;
typedef Container::iterator ContIter;
//original version
//typedef vector<UInt> Container;
//typedef vector<UInt>::iterator ContIter;
// Alain : I am a bit lost in all typedef...

using namespace std;

////////////////////////////////////////////////////////////////
//
// CSR Format
//
///////////////////////////////////////////////////////////////
//
// The matrix is: a reference to a Pattern and a values vector
// Each of them is a template: the Pattern could be symmetric or not
// The values: obviously !
//
// So I can handle non symmetric matrices with a symmetric pattern
//

template<typename PatternType, typename DataType>
class CSRMatr
{
public:
  CSRMatr(); // default constructor : NULL pattern
  //
  // Note that the constructors MUST be based on an existing pattern
  //
  CSRMatr(const PatternType &ex_pattern);
  CSRMatr(const PatternType &ex_pattern, const vector<DataType> &ex_value);
  CSRMatr(const CSRMatr<PatternType,DataType> &RightHandCSR);
  const PatternType * Patt() const {return _Patt;};
  vector<DataType> & value()  {return _value;};
  DataType * giveRawCSR_value() {return &(_values.front());}

  CSRMatr& operator= (const CSRMatr<PatternType,DataType> &RhCsr  );// Warning: the two matrices will point to the same pattern
  void set_mat(UInt where, DataType loc_val);
  void set_mat(UInt row, UInt col, DataType loc_val);
  void set_mat_inc(UInt row, UInt col, DataType loc_val);
  DataType get_value(UInt i, UInt j)
  { return _value[_Patt->locate_index(i,j).first];};

  void ShowMe();
  void spy(string  const &filename);

private:
  vector<DataType> _value;
  const PatternType *_Patt; // I want to link the values to a pattern, NOT to change the pattern itself  (which is const)
  //     static const DataType _DefaultValue = 0;
};

////////////////////////////////////////////////////////////////
//
// CSR Format SPECIALIZATION FOR Usual CSR
//
///////////////////////////////////////////////////////////////
//

template<typename DataType>
class CSRMatr<CSRPatt,DataType>
{
public:
  CSRMatr(); //!< default constructor : NULL pattern
  //
  // Note that the constructors MUST be based on an existing pattern
  //
  CSRMatr(const CSRPatt &ex_pattern);
  //! Version for DataType=Tab2d
  CSRMatr(const CSRPatt &ex_pattern, UInt const nr, UInt const nc);
  CSRMatr(const CSRPatt &ex_pattern, const vector<DataType> &ex_value);
  CSRMatr(const CSRMatr<CSRPatt,DataType> &RightHandCSR);
  const CSRPatt * Patt() const {return _Patt;};
  const vector<DataType> & value() const {return _value;};
  vector<DataType> & value()  {return _value;};
  DataType * giveRawCSR_value() {return &(_value.front());}

  CSRMatr& operator= (const CSRMatr<CSRPatt,DataType> &RhCsr  );
  // Warning: the two matrices will point to the same pattern

  //! Matrix-vector product
  Vector operator*(const Vector &v) const;
  //! Version for block matrices
  VectorBlock operator*(const VectorBlock &v) const;
  //! Version for C pointer vector. BEWARE: no check on bounds done !
  template <typename DataT>
  friend void operMatVec(DataT * const mv,
			 const CSRMatr<CSRPatt,DataT> &Mat,
			 const DataT *v);

  //! necessary for IML++ library: transpose-matrix by vector product
  Vector trans_mult(const Vector &v) const;
  //! version for block matrices
  VectorBlock trans_mult(const VectorBlock &v);

  //! assign a row to zero. Remark, zero might be defined for any DataType
  void zero_row(UInt const row);

  //! diagonalize a row r and rhs for Dirichlet BC treatment
  void diagonalize(UInt const r, DataType const coeff, Vector &b,
		   DataType datum);
  //! diagonalize the row r of the matrix only (cf diagonalize for taking into
  //! account the rhs of the linear system as well)
  void diagonalize_row(UInt const r, DataType const coeff);

  //! determine the lumped matrix for P1
  vector<DataType> MassDiagP1() const;
  //! build the inverse of a diagonal matrix and multiply it by another matrix
  friend void MultInvDiag(const vector<Real> &Diag,
			  const CSRMatr<CSRPatt,Real> &Mat,
			  CSRMatr<CSRPatt,Real> &ans);

  void set_mat(UInt row, UInt col, DataType loc_val);
  void set_mat_inc(UInt row, UInt col, DataType loc_val);
  DataType get_value(UInt i, UInt j)
  { return _value[_Patt->locate_index(i,j).first];};
  const DataType get_value(UInt i, UInt j) const
  { return _value[_Patt->locate_index(i,j).first];};

  void ShowMe();
  void spy(string  const &filename);

  //!column-concatenation of two blocks of CSRMatr
  /*  friend CSRMatr<CSRPatt,double>
      colUnify(const CSRMatr<CSRPatt,double> &Mat1,
      const CSRMatr<CSRPatt,double> &Mat2);*/
  //version without using static : the one to keep (13/12/01)
  friend void
  colUnify(CSRMatr<CSRPatt,double> &ans, const CSRMatr<CSRPatt,double> &Mat1,
	   const CSRMatr<CSRPatt,double> &Mat2);

  //!row-concatenation of two blocks of CSRMatr (without using static)
  friend void
    rowUnify(CSRMatr<CSRPatt,double> &ans, const CSRMatr<CSRPatt,double> &Mat1,
	     const CSRMatr<CSRPatt,double> &Mat2);

  //!set to zero the row trD and the corresponding column=row for D
  template<typename VectorType>
    friend void
    zero_row_col(UInt const row, CSRMatr<CSRPatt,double> &trD,
		 CSRMatr<CSRPatt,double> &D, VectorType &bp,
		 DataType const datum);

  //! set to zero the matrix;
  void zeros();

 private:
  vector<DataType> _value;
  const CSRPatt *_Patt; // I want to link the values to a pattern, NOT to change the pattern itself  (which is const)
  //  static const DataType _DefaultValue = 0;
};

////////////////////////////////////////////////////////////////
//
// VBR Format based on CSR pattern for SQUARE blocks !!!
//
///////////////////////////////////////////////////////////////

template<typename DataType>
class VBRMatr
{
public:
  VBRMatr():_Patt(0){}; // default constructor : NULL pattern
  //
  // Note that the constructors MUST be based on an existing pattern
  //
  VBRMatr(const VBRPatt &ex_pattern);
  VBRMatr(const VBRPatt& ex_pattern, const vector<DataType> &ex_value);
  VBRMatr(const VBRMatr<DataType> &RightHandVBR):
    _Patt(RightHandVBR.Patt()), _value(RightHandVBR.value()){}

  const VBRPatt *Patt() const {return _Patt;}
  vector<DataType> value() const {return _value;}
  DataType* giveRaw_value() {return &(_value.front());} // give the
  // value vector in Raw format (suitable for C)

  VBRMatr& operator=(const VBRMatr<DataType> &RhVbr);// Warning: the
  // two matrices will point to the same pattern

  vector<DataType>  operator*(const vector<DataType> &v) const; //Matrix-vector product
  //Matrix vector product for IML++ (uses the Vector class)
  Vector operator*(const Vector &v) const;
  //necessary for IML++ library: transpose-matrix by vector product
  Vector trans_mult(const Vector &v) const;

  VBRMatr operator*(const DataType num);
  VBRMatr & operator*=(const DataType num); //Real*Matrix product
  VBRMatr & flop(const DataType num, VBRMatr<DataType> & M, VBRMatr<DataType> & A); //Matrix=num*M+A
  VBRMatr & flop(const DataType num, VBRMatr<DataType> & M); //Matrix=num*M+Matrix
  void set_mat(UInt where, DataType loc_val);
  void set_mat(UInt row, UInt col, DataType loc_val);
  void set_mat_inc(UInt row, UInt col, DataType loc_val);
  // version for block operation
  void set_mat(UInt row, UInt col, vector<DataType> &loc_block);
  void set_mat_inc(UInt row, UInt col, vector<DataType> &loc_block);

  DataType& get_value(UInt i, UInt j)
  { return _value[_Patt->locate_index(i,j).first];};
  const DataType& get_value(UInt i, UInt j) const
  { return _value[_Patt->locate_index(i,j).first];};
  //  void ShowMe() const;
  void spy(string  const &filename);
  //  void diagonalize_row ( UInt const r, DataType const coeff);
  //  void diagonalize ( UInt const r, DataType const coeff, vector<DataType> &b, DataType datum);

private:
  vector<DataType> _value;
  const VBRPatt *_Patt; // I want to link the values to a pattern, NOT
  // changing the pattern itself  (which is const)
  //  static const DataType _DefaultValue = 0;
};

////////////////////////////////////////////////////////////////
//
// MSR Format
//
///////////////////////////////////////////////////////////////

template<typename DataType>
class MSRMatr
{
public:
  MSRMatr(); // default constructor : NULL pattern
  //
  // Note that the constructors MUST be based on an existing pattern
  //
  MSRMatr(const MSRPatt &ex_pattern);
  // Alain : constructor from CSR pattern
  MSRMatr(const CSRPatt &ex_pattern);
  MSRMatr(const MSRPatt* ex_pattern, const vector<DataType> &ex_value);
  MSRMatr(const MSRMatr<DataType> &RightHandMSR);
  MSRMatr(const MSRPatt &ex_pattern,const CSRMatr<CSRPatt,DataType> &RightHandCSR);
  // CSR_values -> MSR_values
  /* to convert from CSR to MSR:
     MSRPatt nameMSRP(CSRPatt nameCSRP);
     MSRMatr nameMSR_matrix(nameMSRP,nameCSR_matrix);*/

  // Alain : for conversion, the pattern must have been already converted
  // REMARK: I haven't managed to do it with the template DataType
  friend void CSRmat2MSRmat(MSRMatr<double> &MSRmat,
			    CSRMatr<CSRPatt,double> const &CSRmat);

  const MSRPatt * Patt() const {return _Patt;};

  const vector<DataType> & value() const {return _value;};
  vector<DataType> & value()  {return _value;};

  DataType * giveRaw_value() {return &(_value.front());} // give the value vector in Raw format (suitable for C)
  DataType const * giveRaw_value() const {return &(_value.front());} // give the value vector in Raw format (suitable for C)

  //! determine the lumped matrix for P1
  vector<DataType> MassDiagP1() const;
  //! build the inverse of a diagonal matrix and multiply it by another matrix
  friend void MultInvDiag(const vector<Real> &Diag,
			  const MSRMatr<Real> &Mat, MSRMatr<Real> &ans);

  //! give the diagonal of an MSR matrix
  vector<DataType> giveDiag() const;

  MSRMatr & operator= (const MSRMatr<DataType> &RhMsr);// Warning: the two matrices will point to the same pattern
  void set_mat(UInt where, DataType loc_val);
  void set_mat(UInt row, UInt col, DataType loc_val);
  void set_mat_inc(UInt row, UInt col, DataType loc_val);
  DataType get_value(UInt i, UInt j)
  { return _value[_Patt->locate_index(i,j).first];};
  const DataType get_value(UInt i, UInt j) const
  { return _value[_Patt->locate_index(i,j).first];};
  void ShowMe();
  void spy(string  const &filename);
  void diagonalize_row ( UInt const r, DataType const coeff);
  void diagonalize ( UInt const r, DataType const coeff, vector<DataType> &b, DataType datum);
  //Version for the type Vector
  void diagonalize( UInt const r, DataType const coeff, Vector &b, DataType datum);
  vector<DataType>  operator* (const vector<DataType> &v) const; //Matrix-vector product
  //Matrix-vector product for the class Vector (useful for IML++)
  Vector operator*(const Vector &v) const;
  //! Version for C pointer vector. BEWARE: no check on bounds done !
  template <typename DataT>
  friend void operMatVec(DataT * const mv,
			 const MSRMatr<DataT> &Mat,
			 const DataT *v);
  //necessary for IML++ library: transpose-matrix by vector product
  Vector trans_mult(const Vector &v) const;

  MSRMatr operator*(const DataType num);
  MSRMatr & operator*=(const DataType num); //Real*Matrix product
  MSRMatr & flop(const DataType num, MSRMatr<DataType> & M, MSRMatr<DataType> & A); //Matrix=num*M+A
  MSRMatr & flop(const DataType num, MSRMatr<DataType> & M); //Matrix=num*M+Matrix
  //    MSRMatr & operator* (const DataType num,const MSRMatr<DataType> &M); //Real*Matrix product

  //! set to zero the matrix;
  void zeros();

private:
  vector<DataType> _value;
  const MSRPatt *_Patt; // I want to link the values to a pattern, NOT to change the pattern itself  (which is const)
  //     static const DataType _DefaultValue = 0;
};

////////////////////////////////////////////////////////////////
//
// MixedMatr Format
//
///////////////////////////////////////////////////////////////
// Alain Gauthier. 05/02.
/*! \class MixedMatr
    Contains the values for MixedPattern block patterns. It is useful
    for the management of any block matrices, for instance a matrix of
    a vectorial operator.
    AZTEC library can be used with this type of matrix format through
    a user defined matrix-vector product.
 */
template<UInt BRows, UInt BCols, typename PatternType, typename DataType>
class MixedMatr
{
public:
  //! Default constructor.
  MixedMatr();
  //! Constructor from an existing MixedPattern pattern.
  //! Version for MSR format.
  MixedMatr(const MixedPattern<BRows, BCols, MSRPatt> &ex_pattern);
  //! Version for CSR format.
  MixedMatr(const MixedPattern<BRows, BCols, CSRPatt> &ex_pattern);
  //! gives the pattern pointer.
  const MixedPattern<BRows, BCols, PatternType> * Patt() const {return _Patt;}
  //! points to the values vector array.
  vector<DataType> * bValues() const {return _bValues;}
  //! The container of block (i,j).
  //! BEWARE: no check is done on the array size.
  const vector<DataType> & bValues(UInt i,UInt j) const {return _bValues[i][j];}
  vector<DataType> & bValues(UInt i,UInt j) {return _bValues[i][j];}
  //! give the block (i,j) value vector in Raw format (suitable for C)
  DataType * giveRaw_value(UInt i, UInt j)
  {return &(_bValues[i][j].front());}
  //! give the block (i,j) value vector in Raw format (suitable for C)
  DataType const * giveRaw_value(UInt i, UInt j) const
  {return &(_bValues[i][j].front());}

  //! Determines the lumped diagonal of P1 mass matrix.
  vector<DataType> MassDiagP1() const;
  //! Inverts the diagonal matrix Diag and mulitply it by the matrix Mat.
  template<UInt BR, UInt BC, typename PattType>
  friend void MultInvDiag(const vector<Real> &Diag,
			  const MixedMatr<BR,BC,PattType,Real> &Mat,
			  MixedMatr<BR,BC,PattType,Real> &ans);
  //! Gives the diagonal of a block matrix.
  vector<DataType> giveDiag() const;

  //! Warning: the two matrices will point to the same pattern.
  MixedMatr & operator= (const MixedMatr<BRows,BCols,PatternType,DataType>
			 &RhMtr);
  //! Assigns the matrix to loc_val at the place of index where in the
  //! block (ib,jb).
  void set_mat(UInt ib, UInt jb, UInt where, DataType loc_val);
  //! Assigns the matrix element (row,col) to loc_val.
  void set_mat(UInt row, UInt col, DataType loc_val);
  //! Assigns the matrix element (row,col) of block (ib,jb) to loc_val.
  void set_mat(UInt ib, UInt jb, UInt row, UInt col, DataType loc_val);
  //! Adds loc_val to the matrix element (row,col).
  void set_mat_inc(UInt row, UInt col, DataType loc_val);
  //! Adds loc_val to the matrix element (row,col) of block (ib,jb).
  void set_mat_inc(UInt ib, UInt jb, UInt row, UInt col, DataType loc_val);
  //! Returns the matrix element (i,j) value.
  DataType get_value(UInt i, UInt j);
  const DataType get_value(UInt i, UInt j) const;
  //! Returns the matrix element (i,j) value of block (ib,jb).
  DataType get_value(UInt ib, UInt jb, UInt i, UInt j);
  const DataType get_value(UInt ib, UInt jb, UInt i, UInt j) const;
  //! Shows the matrix (only the pattern here).
  void ShowMe()
  {return _Patt->showMe(true);}
  //! Matrix visualization a la matlab.
  void spy(string  const &filename);
  //! Assigns matrix diagonal element (r,r) to coeff, other elts
  //! of row r to zero.
  void diagonalize_row( UInt const r, DataType const coeff);
  //! assign a row to zero. Remark, zero might be defined for any DataType
  void zero_row(UInt const row);
  //! Assigns matrix diagonal element (r,r) to coeff, other elts
  //! of row r to zero, and vector b element b(r) to coeff*datum.
  template< typename VectorType >
    void diagonalize( UInt const row, DataType const coeff, VectorType &b,
		      DataType datum);

  //! Matrix-vector product.
  vector<DataType>  operator*(const vector<DataType> &v) const;
  //! Version for type Vector.
  Vector operator*(const Vector &v) const;
  //! Version for C pointer vector. BEWARE: no check on bounds done !
  template<typename DataT, UInt BR, UInt BC, typename PatternT>
  friend void operMatVec(DataT * const mv,
			 const MixedMatr<BR, BC, PatternT,DataT> &Mat,
			 const DataT *v);
  //! necessary for IML++ library: transpose-matrix by vector product.
  Vector trans_mult(const Vector &v) const;

  //! set to zero the matrix;
  void zeros();

private:
  //! A pointer to the pattern.
  const MixedPattern<BRows, BCols, PatternType> * _Patt;
  //! An array of values vector: there is one vector of values for
  //! each block.
  vector<DataType> _bValues[BRows][BCols];

};

////////////////////////////////////////////////////////////////
//
// MixedMatr Format SPECIALIZATION FOR MSR
//
///////////////////////////////////////////////////////////////
// Alain Gauthier. 11/02.
/*! \class MixedMatr
    Contains the values for MixedPattern block patterns. It is useful
    for the management of any block matrices, for instance a matrix of
    a vectorial operator.
    AZTEC library can be used with this type of matrix format through
    a user defined matrix-vector product.

    This a more efficient implementation in the case of the use of MSR pattern

 */
template<UInt BRows, UInt BCols>
class MixedMatr<BRows, BCols, MSRPatt, double>
{
 public:
  //! Default constructor.
  MixedMatr();
  //! Constructor from an existing MixedPattern pattern.
  //! Version for MSR format.
  MixedMatr(const MixedPattern<BRows, BCols, MSRPatt> &ex_pattern);
  //! Version for CSR format.
  MixedMatr(const MixedPattern<BRows, BCols, CSRPatt> &ex_pattern);
  //! gives the pattern pointer.
  const MixedPattern<BRows, BCols, MSRPatt> * Patt() const {return _Patt;}
  //! points to the values vector array.
  vector<double> * bValues() const {return _bValues;}
  //! The container of block (i,j).
  //! BEWARE: no check is done on the array size.
  const vector<double> & bValues(UInt i,UInt j) const {return _bValues[i][j];}
  vector<double> & bValues(UInt i,UInt j) {return _bValues[i][j];}
  //! give the block (i,j) value vector in Raw format (suitable for C)
  double * giveRaw_value(UInt i, UInt j)
  {return &(_bValues[i][j].front());}
  //! give the block (i,j) value vector in Raw format (suitable for C)
  double const * giveRaw_value(UInt i, UInt j) const
  {return &(_bValues[i][j].front());}

  //! Determines the lumped diagonal of P1 mass matrix.
  vector<double> MassDiagP1() const;
  //! Inverts the diagonal matrix Diag and mulitply it by the matrix Mat.
  template<UInt BR, UInt BC>
  friend void MultInvDiag(const vector<Real> &Diag,
			  const MixedMatr<BR,BC,MSRPatt,Real> &Mat,
			  MixedMatr<BR,BC,MSRPatt,Real> &ans);
  //! Gives the diagonal of a block matrix.
  vector<double> giveDiag() const;

  //! Warning: the two matrices will point to the same pattern.
  MixedMatr & operator= (const MixedMatr<BRows,BCols,MSRPatt,double>
			 &RhMtr);
  //! Assigns the matrix to loc_val at the place of index where in the
  //! block (ib,jb).
  void set_mat(UInt ib, UInt jb, UInt where, double loc_val);
  //! Assigns the matrix element (row,col) to loc_val.
  void set_mat(UInt row, UInt col, double loc_val);
  //! Assigns the matrix element (row,col) of block (ib,jb) to loc_val.
  void set_mat(UInt ib, UInt jb, UInt row, UInt col, double loc_val);
  //! Adds loc_val to the matrix element (row,col).
  void set_mat_inc(UInt row, UInt col, double loc_val);
  //! Adds loc_val to the matrix element (row,col) of block (ib,jb).
  void set_mat_inc(UInt ib, UInt jb, UInt row, UInt col, double loc_val);
  //! Returns the matrix element (i,j) value.
  double get_value(UInt i, UInt j);
  const double get_value(UInt i, UInt j) const;
  //! Returns the matrix element (i,j) value of block (ib,jb).
  double get_value(UInt ib, UInt jb, UInt i, UInt j);
  const double get_value(UInt ib, UInt jb, UInt i, UInt j) const;
  //! Shows the matrix (only the pattern here).
  void ShowMe()
  {return _Patt->showMe(true);}
  //! Matrix visualization a la matlab.
  void spy(string  const &filename);
  //! Assigns matrix diagonal element (r,r) to coeff, other elts
  //! of row r to zero.
  void diagonalize_row( UInt const r, double const coeff);
  //! assign a row to zero. Remark, zero might be defined for any double
  void zero_row(UInt const row);
  //! Assigns matrix diagonal element (r,r) to coeff, other elts
  //! of row r to zero, and vector b element b(r) to coeff*datum.
  template< typename VectorType >
    void diagonalize( UInt const row, double const coeff, VectorType &b,
		      double datum);

  //! Matrix-vector product.
  vector<double>  operator*(const vector<double> &v) const;
  //! Version for type Vector.
  Vector operator*(const Vector &v) const;
  //! Version for C pointer vector. BEWARE: no check on bounds done !
  template<UInt BR, UInt BC>
  friend void operMatVec(double * const mv,
			 const MixedMatr<BR, BC, MSRPatt, double> &Mat,
			 const double *v);
  //! necessary for IML++ library: transpose-matrix by vector product.
  Vector trans_mult(const Vector &v) const;
  //! set to zero the matrix;
  void zeros();

 private:
  //! A pointer to the pattern.
  const MixedPattern<BRows, BCols, MSRPatt> * _Patt;
  //! An array of values vector: there is one vector of values for
  //! each block.
  vector<double> _bValues[BRows][BCols];

};

////////////////////////////////////////////////////////////////
//
// MixedMatr Format SPECIALIZATION FOR CSR
//
///////////////////////////////////////////////////////////////
// Alain Gauthier. 05/02.
/*! \class MixedMatr
    Contains the values for MixedPattern block patterns. It is useful
    for the management of any block matrices, for instance a matrix of
    a vectorial operator.
    AZTEC library can be used with this type of matrix format through
    a user defined matrix-vector product.

    This a more efficient implementation in the case of the use of CSR pattern

 */
template<UInt BRows, UInt BCols>
class MixedMatr<BRows, BCols, CSRPatt, double>
{
 public:
  //! Default constructor.
  MixedMatr();
  //! Constructor from an existing MixedPattern pattern.
  //! Version for CSR format.
  MixedMatr(const MixedPattern<BRows, BCols, CSRPatt> &ex_pattern);
  //! gives the pattern pointer.
  const MixedPattern<BRows, BCols, CSRPatt> * Patt() const {return _Patt;}
  //! points to the values vector array.
  vector<double> * bValues() const {return _bValues;}
  //! The container of block (i,j).
  //! BEWARE: no check is done on the array size.
  const vector<double> & bValues(UInt i,UInt j) const {return _bValues[i][j];}
  vector<double> & bValues(UInt i,UInt j) {return _bValues[i][j];}
  //! give the block (i,j) value vector in Raw format (suitable for C)
  double * giveRaw_value(UInt i, UInt j)
  {return &(_bValues[i][j].front());}
  //! give the block (i,j) value vector in Raw format (suitable for C)
  double const * giveRaw_value(UInt i, UInt j) const
  {return &(_bValues[i][j].front());}

  //! Determines the lumped diagonal of P1 mass matrix.
  vector<double> MassDiagP1() const;
  //! Inverts the diagonal matrix Diag and mulitply it by the matrix Mat.
  template<UInt BR, UInt BC>
  friend void MultInvDiag(const vector<Real> &Diag,
			  const MixedMatr<BR,BC,CSRPatt,Real> &Mat,
			  MixedMatr<BR,BC,CSRPatt,Real> &ans);

  //! Gives the diagonal of a block matrix.
  vector<double> giveDiag() const;

  //! Warning: the two matrices will point to the same pattern.
  MixedMatr & operator= (const MixedMatr<BRows,BCols,CSRPatt,double>
			 &RhMtr);
  //! Assigns the matrix to loc_val at the place of index where in the
  //! block (ib,jb).
  void set_mat(UInt ib, UInt jb, UInt where, double loc_val);
  //! Assigns the matrix element (row,col) to loc_val.
  void set_mat(UInt row, UInt col, double loc_val);
  //! Assigns the matrix element (row,col) of block (ib,jb) to loc_val.
  void set_mat(UInt ib, UInt jb, UInt row, UInt col, double loc_val);
  //! Adds loc_val to the matrix element (row,col).
  void set_mat_inc(UInt row, UInt col, double loc_val);
  //! Adds loc_val to the matrix element (row,col) of block (ib,jb).
  void set_mat_inc(UInt ib, UInt jb, UInt row, UInt col, double loc_val);
  //! Returns the matrix element (i,j) value.
  double get_value(UInt i, UInt j);
  const double get_value(UInt i, UInt j) const;
  //! Returns the matrix element (i,j) value of block (ib,jb).
  double get_value(UInt ib, UInt jb, UInt i, UInt j);
  const double get_value(UInt ib, UInt jb, UInt i, UInt j) const;
  //! Shows the matrix (only the pattern here).
  void ShowMe()
  {return _Patt->showMe(true);}
  //! Matrix visualization a la matlab.
  void spy(string  const &filename);
  //! Assigns matrix diagonal element (r,r) to coeff, other elts
  //! of row r to zero.
  void diagonalize_row( UInt const r, double const coeff);
  //! assign a row to zero. Remark, zero might be defined for any double
  void zero_row(UInt const row);
  //! Assigns matrix diagonal element (r,r) to coeff, other elts
  //! of row r to zero, and vector b element b(r) to coeff*datum.
  template< typename VectorType >
    void diagonalize( UInt const row, double const coeff, VectorType &b,
		      double datum);

  //!set to zero the row trD and the corresponding column=row for D
  template<UInt BR, UInt BC, typename VectorType, typename DataType>
    friend void
    zero_row_col(UInt const row, MixedMatr<BR,BC,CSRPatt,Real> &trD,
		 MixedMatr<BC,BR,CSRPatt,Real> &D, VectorType &bp,
		 DataType const datum);

  //! Matrix-vector product.
  vector<double>  operator*(const vector<double> &v) const;
  //! Version for type Vector.
  Vector operator*(const Vector &v) const;
  //! Version for C pointer vector. BEWARE: no check on bounds done !
  template<UInt BR, UInt BC>
  friend void operMatVec(double * const mv,
			 const MixedMatr<BR, BC, CSRPatt, double> &Mat,
			 const double *v);
  //! necessary for IML++ library: transpose-matrix by vector product.
  Vector trans_mult(const Vector &v) const;
  //! set to zero the matrix;
  void zeros();

private:
  //! A pointer to the pattern.
  const MixedPattern<BRows, BCols, CSRPatt> * _Patt;
  //! An array of values vector: there is one vector of values for
  //! each block.
  vector<double> _bValues[BRows][BCols];

};

////////////////////////////////////////////////////////////////
//
// DiagPreconditioner class for IML++ dense preconditioner matrix
//
///////////////////////////////////////////////////////////////
/* !\class DiagPreconditioner
   class useful for IML++ dense preconditioner matrix
*/
template<typename VectorType>
class DiagPreconditioner
{
  VectorType _diag;
public:
  DiagPreconditioner(){};
  //!for CSR or MSR normal pattern
  DiagPreconditioner(const CSRMatr<CSRPatt,double> &M);
  DiagPreconditioner(const MSRMatr<double> &M);
  //!for VBR pattern
  DiagPreconditioner(const VBRMatr<double> &M);
  //!for CSR block pattern
  DiagPreconditioner(const CSRMatr<CSRPatt,Tab2d> &M);

  Vector solve(const Vector &x) const;
  VectorBlock solve(const VectorBlock &x) const;
  VectorType trans_solve(const VectorType &x) const
  {return solve(x);}

  const double & diag(UInt i) const {return _diag(i);}
  double & diag(UInt i) {return _diag(i);}
  const Tab1d & diagBlock(UInt i) const {return _diag.numBlock(i);}
  Tab1d & diagBlock(UInt i) {return _diag.numBlock(i);}
};

////////////////////////////////////////////////////////////////
//
// IDPreconditioner class for IML++ dense preconditioner matrix
// Here the preconditioner is simply the identity matrix
//
///////////////////////////////////////////////////////////////
/*! \class IDPreconditioner
  class useful for IML++ dense preconditioner matrix.
  Here the preconditioner is simply the identity matrix (no preconditioning).
 */
template<typename VectorType>
class IDPreconditioner
{
  VectorType _diag;
public:
  IDPreconditioner(){};
  //! for CSR or MSR normal pattern
  IDPreconditioner(const CSRMatr<CSRPatt,double> &M);
  IDPreconditioner(const MSRMatr<double> &M);
  //! for VBR pattern
  IDPreconditioner(const VBRMatr<double> &M);
  //! for CSR block pattern
  IDPreconditioner(const CSRMatr<CSRPatt,Tab2d> &M);

  Vector solve(const Vector &x) const;
  VectorBlock solve(const VectorBlock &x) const;
  VectorType trans_solve(const VectorType &x) const
  {return solve(x);}

  const double & diag(UInt i) const {return _diag(i);}
  double & diag(UInt i) {return _diag(i);}
  const Tab1d & diagBlock(UInt i) const {return _diag.numBlock(i);}
  Tab1d & diagBlock(UInt i) {return _diag.numBlock(i);}
};


//-------------------------------------------------------------------------------------------------------
// CSR - VALUES
//------------------------------------------------------------------------------------------------------
template<typename PatternType, typename DataType>
CSRMatr<PatternType,DataType>::CSRMatr():_Patt(0){};

template<typename PatternType,typename DataType>
CSRMatr<PatternType,DataType>::
CSRMatr(const PatternType &ex_pattern)
{
  _Patt = &ex_pattern;
  _value.reserve(ex_pattern.nNz()); // no initialization of values !
}


template<typename PatternType, typename DataType>
CSRMatr<PatternType,DataType>::
CSRMatr(const PatternType &ex_pattern, const vector<DataType> &ex_value)
{
  _value.reserve(ex_pattern.nNz()); // in case of block matrix, there is
  // no default constructor for class KNM
  _Patt = &ex_pattern;
  ASSERT( _Patt->nNz() == ex_value.size(),
	  "Error in CSR Matrix Values LifeV");
  // Warning: if PatternType = CSRPattSymm => _ja.size() != _nnz.===> remember: _ja.size() = (_nnz+_nrows)/2
  _value = ex_value;
}

template<typename PatternType, typename DataType>
CSRMatr<PatternType,DataType>::
CSRMatr(const CSRMatr<PatternType,DataType> &RightHandCSR):
  _value(RightHandCSR.value()),_Patt(RightHandCSR.Patt()) {};

template<typename PatternType, typename DataType>
CSRMatr<PatternType,DataType>&
CSRMatr<PatternType,DataType>::operator= (const CSRMatr<PatternType,DataType> &RhCsr)
{
  if (&RhCsr != this)
    {
      _Patt = RhCsr.Patt();
      _value = RhCsr.value();
    }
  return *this;
};


template<typename PatternType, typename DataType>
void
CSRMatr<PatternType,DataType>::
set_mat(UInt where, DataType loc_val)
{
  _value[where] = loc_val;
  return;
};

template<typename PatternType, typename DataType>
void
CSRMatr<PatternType,DataType>::
set_mat(UInt row, UInt col, DataType loc_val)
{
  pair<UInt,bool> where = *_Patt.locate_index(row,col);
  if (where.second) _value[where.first] = loc_val;
  return;
};

template<typename PatternType, typename DataType>
void
CSRMatr<PatternType,DataType>::
set_mat_inc(UInt row, UInt col, DataType loc_val)
{

  pair<UInt,bool> where = _Patt->locate_index(row,col);
  if (where.second) _value[where.first] += loc_val;

  return;
};


//-------------------------------------------------------------------------------------------------------
// CSR - VALUES WITH CSR Pattern (Usual)
//------------------------------------------------------------------------------------------------------
template<typename DataType>
CSRMatr<CSRPatt,DataType>::CSRMatr():_Patt(0){};


template<typename DataType>
CSRMatr<CSRPatt,DataType>::
CSRMatr(const CSRPatt &ex_pattern)
{
  _Patt = &ex_pattern;
  _value.resize(ex_pattern.nNz());
}
//version for Datatype=Tab2d
CSRMatr<CSRPatt,Tab2d>::
CSRMatr(const CSRPatt &ex_pattern, UInt const nr, UInt const nc);



template<typename DataType>
CSRMatr<CSRPatt,DataType>::
CSRMatr(const CSRPatt &ex_pattern, const vector<DataType> &ex_value)
{
  _value.reserve(ex_pattern.nNz()); // in case of block matrix, there is
  // no default constructor for class KNM
  _Patt = &ex_pattern;
  ASSERT( _Patt->nNz() == ex_value.size(),
	  "Error in CSR Matrix Values LifeV");
  // Warning: if PatternType = CSRPattSymm => _ja.size() != _nnz.===> remember: _ja.size() = (_nnz+_nrows)/2
  _value = ex_value;
}

template<typename DataType>
CSRMatr<CSRPatt,DataType>::
CSRMatr(const CSRMatr<CSRPatt,DataType> &RightHandCSR):
  _value(RightHandCSR.value()),_Patt(RightHandCSR.Patt()) {};

template<typename DataType>
CSRMatr<CSRPatt,DataType>&
CSRMatr<CSRPatt,DataType>::operator= (const CSRMatr<CSRPatt,DataType> &RhCsr)
{
  if (&RhCsr != this)
    {
      _Patt = RhCsr.Patt();
      _value = RhCsr.value();
    }
  return *this;
};


template<typename DataType>
void
CSRMatr<CSRPatt,DataType>::
set_mat(UInt row, UInt col, DataType loc_val)
{
  pair<UInt,bool> where = _Patt->locate_index(row,col);
  if (where.second) _value[where.first] = loc_val;
  return;
};

template<typename DataType>
void
CSRMatr<CSRPatt,DataType>::
set_mat_inc(UInt row, UInt col, DataType loc_val)
{
  pair<UInt,bool> where = _Patt->locate_index(row,col);
  if (where.second) _value[where.first] += loc_val;

  return;
};

// determine the lumped matrix for P1
template<typename DataType>
vector<DataType>
CSRMatr<CSRPatt,DataType>::MassDiagP1() const
{
  UInt nrows = _Patt->nRows();
  UInt ncols = _Patt->nCols();
  ASSERT(ncols == nrows, "The matrix must be square to be lumped");
  vector<DataType> diag(nrows);
  for (UInt nrow = 0; nrow < nrows; ++nrow)
    {
      for (UInt ii = _Patt->ia()[nrow]; ii < _Patt->ia()[nrow+1]; ++ii)
	diag[nrow] += _value[ii];
    };
  return diag;
}


// build the inverse of a diagonal matrix and multiply it by another matrix
void MultInvDiag(const vector<Real> &Diag,
		 const CSRMatr<CSRPatt,Real> &Mat, CSRMatr<CSRPatt,Real> &ans) ;


// transpose-Matrix by vector product
template<typename DataType>
Vector
CSRMatr<CSRPatt,DataType>::
trans_mult(const Vector &v) const
{
  UInt nrows=_Patt->nRows(); // for square matrices...
  ASSERT(nrows==v.size(),"Error in Matrix Vector product");
  Vector ans(nrows);
  ans=0.0;

  for (UInt ir=0+OFFSET;ir<nrows+OFFSET;++ir){
    for (UInt ii=_Patt->ia()[ir]-OFFSET;ii<_Patt->ia()[ir+1]-OFFSET;++ii)
      ans(_Patt->ja()[ii]-OFFSET)+=_value[ii]*v(ir);
  };
  return ans;
}

// version for block matrices
VectorBlock
CSRMatr<CSRPatt,Tab2d>::
trans_mult(const VectorBlock &v);

//
template<typename DataType>
void CSRMatr<CSRPatt,DataType>::
zero_row(UInt const row)
{
  typename vector<DataType>::iterator start= _value.begin()+
    *(_Patt->give_ia().begin()+row-OFFSET);
  typename vector<DataType>::iterator end= _value.begin()+
    *(_Patt->give_ia().begin()+row+1-OFFSET);

  transform(start,end,start,nihil); //nihil is the same used for diagonalize
  //method in MSRMatr class.
}


// Matrix-vector product
template<typename DataType>
Vector  CSRMatr<CSRPatt,DataType>::
operator*(const Vector &v) const
{
  UInt nrows=_Patt->nRows();
  UInt ncols=_Patt->nCols();

  ASSERT(ncols==v.size(),"Error in Matrix Vector product");
  Vector ans(nrows);
  ans=0.0;

  for (UInt ir=0+OFFSET;ir<nrows+OFFSET;++ir){
    for (UInt ii=_Patt->give_ia()[ir]-OFFSET;ii<_Patt->give_ia()[ir+1]-OFFSET;++ii)
      ans(ir)+=_value[ii]*v(_Patt->give_ja()[ii]-OFFSET);
  };
  return ans;
}

// version for block matrices
VectorBlock CSRMatr<CSRPatt,Tab2d>::
operator*(const VectorBlock &v) const;

// Version for C pointer vector. BEWARE: no check on bounds done !
template<typename DataType>
void operMatVec(DataType * const mv,
		const CSRMatr<CSRPatt, DataType> &Mat,
		const DataType *v)
{
  UInt nrows=Mat._Patt->nRows();
  UInt ncols=Mat._Patt->nCols();

  for (UInt ir=0+OFFSET;ir<nrows+OFFSET;++ir){
    mv[ir]= 0.;
    for (UInt ii=Mat._Patt->give_ia()[ir]-OFFSET;ii<Mat._Patt->give_ia()[ir+1]-OFFSET;++ii)
      mv[ir]+=Mat._value[ii]*v[Mat._Patt->give_ja()[ii]-OFFSET];
  };
}

/*! Diagonalization of row r of the system. Done by setting A(r,r) = coeff,
 *  A(r,j) = 0 and A(j,r) = 0 for j!=r, and suitably correcting the right hand
 *  side of the system.
 *  @param r row to diagonalize
 *  @param coeff value to set the diagonal entry A(r,r) to
 *  @param b right hand side vector to be corrected
 *  @param datum value to set the fix the solution entry x(r) at
 */
template<typename DataType>
void CSRMatr<CSRPatt,DataType>::diagonalize(UInt const r, DataType const coeff,
                                            Vector &b, DataType datum) {

    for(UInt j=0; j<_Patt->nCols(); ++j) {
        set_mat(r, j, 0.); // A(r,j) = 0
    }
    for(UInt i=0; i<_Patt->nRows(); ++i) {
        std::pair<UInt,bool> where = _Patt->locate_index(i, r);
        if (where.second) {
            b[i-OFFSET] -= _value[where.first] * datum; // correct rhs
            _value[where.first] = 0.0; // A(j,r) = 0
        }
    }

    set_mat(r, r, coeff); // A(r,r) = coeff

    b[r-OFFSET] = coeff*datum; // correct right hand side for row r
}

/*! Diagonalization of row r of the system. Done by setting A(r,r) = coeff
 *  and A(r,j) = 0 for j!=r without correcting the right hand side
 *  @param r row to diagonalize
 *  @param coeff value to set the diagonal entry A(r,r) to
 */
template<typename DataType>
void CSRMatr<CSRPatt, DataType>::diagonalize_row(UInt const r,
                                                 DataType const coeff) {
    for(UInt j=0; j<_Patt->nCols(); ++j) {
        set_mat(r, j, 0.); // A(r,j) = 0
    }
    set_mat(r, r, coeff); // A(r,r) = coeff
}


// Correction of ShowMe, Alain, 05/02/02.
template<typename DataType>
void
CSRMatr<CSRPatt,DataType>::ShowMe()
{
  UInt i_first,nrows=_Patt->nRows(),ncols=_Patt->nCols(),nnz=_Patt->nNz();
  Container ja = _Patt->ja();

  string pare="[";
  cout << "**************************" << endl;
  cout << "     CSR Matrix           " << endl;
  cout << endl;
  cout << pare;
  for (UInt i_index=0; i_index < nrows; ++i_index)
    {
      cout << pare;
      pare = " [";
      i_first=_Patt->give_ia()[i_index] - OFFSET;

      UInt jj=0;
      for(UInt j=0;j<ncols;j++)
	{
	  if (j==i_index) {cout << " " <<_value[i_first+jj]<<" "; jj++;}
	  else {
	    if (j==ja[i_first+jj]-OFFSET){
	      cout <<" "<< _value[i_first+jj]<<" "; jj++;}
	    else
	      cout << " 0 ";
	  }
	}
      if (i_index==nrows-1)
        cout << " ]] " << endl;
      else
        cout << " ]  " << endl;
    }
  cout << "nnz = " << nnz << ", nrow = " << nrows << ", ncol = " << ncols << endl;
  return;
};

template<typename DataType>
void CSRMatr<CSRPatt,DataType>::
spy(string  const &filename)
{
  // Purpose: Matlab dumping and spy
  string nome=filename, uti=" , ";
  UInt nrows=_Patt->nRows();
  Container ia=_Patt->ia(), ja=_Patt->ja();
  //
  // check on the file name
  //
  int i=filename.find(".");

  if (i<=0) nome=filename+".m";
  else {
    if ((unsigned int)i!=filename.size()-2  || filename[i+1]!='m')
      {cerr << "Wrong file name ";
      nome=filename+".m";}
  }

  ofstream file_out(nome.c_str());
  ASSERT(file_out,"Error: Output Matrix (Values) file cannot be open");


  file_out << "S = [ ";
  for (UInt i=0;i<nrows;++i){
    for (UInt ii=ia[i]-OFFSET;ii<ia[i+1]-OFFSET;++ii)
      file_out << i+1 << uti << ja[ii]+1-OFFSET << uti << _value[ii] << endl; /* */
  }
  file_out << "];" << endl;

  file_out << "I=S(:,1); J=S(:,2); S=S(:,3); A=sparse(I,J,S); spy(A);"<<endl;
};

// the case of block matrices with Tab2d block type.
void CSRMatr<CSRPatt,Tab2d>::
spy(string  const &filename);

//version without using static (I think it is better)
// Modified by A. Gilardi. 03/02.
void colUnify(CSRMatr<CSRPatt,double> &ans, const CSRMatr<CSRPatt,double> &Mat1,
	      const CSRMatr<CSRPatt,double> &Mat2);
// row-concatenation of two blocks of CSRMatr
// Modified by A. Gilardi. 03/02.
void rowUnify(CSRMatr<CSRPatt,double> &ans, const CSRMatr<CSRPatt,double> &Mat1,
	      const CSRMatr<CSRPatt,double> &Mat2);

// set to zero the row trD and the corresponding column=row for D
// passing the known datum to the rhs b.
template<typename VectorType, typename DataType>
void zero_row_col(UInt const row, CSRMatr<CSRPatt,double> &trD,
		  CSRMatr<CSRPatt,double> &D, VectorType &bp,
		  DataType const datum)
{
  // for trD
  trD.zero_row(row);

  // for D
  UInt j;
  //loop on the columns of trD involved
  for ( j= trD._Patt->give_ia()[row-OFFSET];
	j < trD._Patt->give_ia()[row+1-OFFSET]; j++)
    {
      // columns of trD become rows of D
      row_loc= trD._Patt->give_ja()[j];
      bp[row_loc]-= datum * D._value[trD._Patt->give_jaT()[j]];
      D._value[trD._Patt->give_jaT()[j]]= 0.;
    }
}

template<typename DataType>
void CSRMatr<CSRPatt,DataType>::
zeros()
{
  typename vector<DataType>::iterator start= _value.begin();
  typename vector<DataType>::iterator end  = _value.end();
  fill(start, end, 0.0);
}

//---------------------------------------------------------------------
// VBR - VALUES
//---------------------------------------------------------------------

template<typename DataType>
VBRMatr<DataType>::
VBRMatr(const VBRPatt &ex_pattern)
{
  _Patt=&ex_pattern;

  UInt blockSize_r=ex_pattern.rpntr()[1]-ex_pattern.rpntr()[0];
  UInt blockSize_c=blockSize_r; // we work with square blocks !!!

  _value.reserve(ex_pattern.nNz()*blockSize_r*blockSize_c);
}

template<class DataType>
VBRMatr<DataType>::
VBRMatr(const VBRPatt & ex_pattern, const vector<DataType> &ex_value):
  _value(ex_value)
{
  _Patt=&ex_pattern;
  UInt blockSize_r=_Patt->rpntr()[1]-_Patt->rpntr()[0];
  UInt blockSize_c=blockSize_r; // we work with square blocks !!!

  ASSERT( (_Patt->nNz()*blockSize_r*blockSize_c) == ex_value.size(), "Error in VBR Values "); // VBR value has lenghth nnz*blockSize_r*blockSize_c
}

template<class DataType>
VBRMatr<DataType>&
VBRMatr<DataType>::operator=(const VBRMatr<DataType> &RhVbr)
{
  if (&RhVbr != this)
    {
      _value = RhVbr.value();
      _Patt = RhVbr.Patt();
    }
  return *this;
}

template<class DataType>
vector<DataType>
VBRMatr<DataType>::operator*(const vector<DataType> &v) const
{
  UInt blockSize_r=_Patt->rpntr()[1]-_Patt->rpntr()[0];
  //...for square matrices...
  UInt nrows=_Patt->nRows()*blockSize_r;
  UInt ncols=_Patt->nCols()*blockSize_r;
  ASSERT(ncols==v.size(),"Error in Matrix Vector product");
  vector<DataType> ans;
  ans.resize(nrows,0.0);

  for (UInt ir=0+OFFSET;ir<nrows+OFFSET;++ir){
    // column index of non-zero elements and its position respectively
    Container coldata, position;
    UInt nnz_c;
    nnz_c=_Patt->row(ir, coldata, position);
    for (UInt jc=0; jc < nnz_c ; jc++){
      ans[ir]+=_value[position[jc]+OFFSET]*v[coldata[jc]+OFFSET];
    };
  };
  return ans;
}
// version for Vector class
template<class DataType>
Vector
VBRMatr<DataType>::
operator*(const Vector &v) const
{
  UInt blockSize_r=_Patt->rpntr()[1]-_Patt->rpntr()[0];
  //...for square matrices...
  UInt nrows=_Patt->nRows()*blockSize_r;
  UInt ncols=_Patt->nCols()*blockSize_r;
  ASSERT(ncols==v.size(),"Error in Matrix Vector product");
  Vector ans(nrows);
  ans=0.;

  for (UInt ir=0+OFFSET;ir<nrows+OFFSET;++ir){
    // column index of non-zero elements and its position respectively
    Container coldata, position;
    UInt nnz_c;
    nnz_c=_Patt->row(ir, coldata, position);
    for (UInt jc=0; jc < nnz_c ; jc++){
      ans(ir)+=_value[position[jc]+OFFSET]*v(coldata[jc]+OFFSET);
    };
  };
  return ans;
}

//necessary for IML++ library: transpose-matrix by vector product
template<class DataType>
Vector VBRMatr<DataType>::
trans_mult(const Vector &v) const
{
  UInt blockSize_r=_Patt->rpntr()[1]-_Patt->rpntr()[0];
  //...for square matrices...
  UInt nrows=_Patt->nRows()*blockSize_r;
  ASSERT(nrows==v.size(),"Error in Matrix Vector product");
  Vector ans(nrows);
  ans=0.;

  for (UInt ir=0+OFFSET;ir<nrows+OFFSET;++ir){
    // column index of non-zero elements and its position respectively
    Container coldata, position;

    UInt nnz_c;
    nnz_c=_Patt->row(ir, coldata, position);
    for (UInt jc=0; jc < nnz_c ; jc++){
      ans(coldata[jc]+OFFSET)+=_value[position[jc]+OFFSET]*v(ir);
    };
  };
  return ans;
}

template<class DataType>
VBRMatr<DataType>&
VBRMatr<DataType>::operator*=(const DataType num)
{
  UInt blockSize_r=_Patt->rpntr()[1]-_Patt->rpntr()[0];
  UInt blockSize_c=blockSize_r; // we work with square blocks !!!
  UInt stop=_Patt->nNz()*blockSize_r*blockSize_r;
  for (UInt i=0;i<stop;++i)
    _value[i]*=num;

  return *this;
}

template<class DataType>
VBRMatr<DataType>
VBRMatr<DataType>::operator*(const DataType num)
{
  UInt blockSize_r=_Patt->rpntr()[1]-_Patt->rpntr()[0];
  UInt blockSize_c=blockSize_r; // we work with square blocks !!!
  UInt stop=_Patt->nNz()*blockSize_r*blockSize_r;
  VBRMatr<DataType> ans(*this);

  for (UInt i=0;i<stop;++i)
    ans.set_mat(i,num*_value[i]);

  return ans;
}

template<class DataType>
VBRMatr<DataType>&
VBRMatr<DataType>::flop(const DataType num, VBRMatr<DataType>& M,
			VBRMatr<DataType>& A)
{
  /* AIM: Matrix = num*M + A */
  UInt blockSize_r=_Patt->rpntr()[1]-_Patt->rpntr()[0];
  UInt blockSize_c=blockSize_r; // we work with square blocks !!!
  UInt stop=_Patt->nNz()*blockSize_r*blockSize_r;

  for (UInt i=0;i<stop;++i)
    _value[i]=num*M.value()[i]+A.value()[i];

  return *this;
}

template<class DataType>
VBRMatr<DataType>&
VBRMatr<DataType>::flop(const DataType num, VBRMatr<DataType>& M)
{
  /* AIM: Matrix = num*M + Matrix */
  UInt blockSize_r=_Patt->rpntr()[1]-_Patt->rpntr()[0];
  UInt blockSize_c=blockSize_r; // we work with square blocks !!!
  UInt stop=_Patt->nNz()*blockSize_r*blockSize_r;

  for (UInt i=0;i<stop;++i)
    _value[i]+=num*M.value()[i];

  return *this;
}

template<typename DataType>
void VBRMatr<DataType>::
set_mat_inc(UInt row, UInt col, DataType loc_val)
{
  UInt blockSize = _Patt->rpntr()[1]-_Patt->rpntr()[0];
  PatternDefs::Diff_t blocrow = _Patt->rbloc(row+OFFSET);
  PatternDefs::Diff_t bloccol = _Patt->cbloc(col+OFFSET);

  pair<UInt,bool> whereBloc = _Patt->locate_index(blocrow,bloccol);

  UInt ind=_Patt->indx()[whereBloc.first] + _Patt->locr(row+OFFSET) +
    _Patt->locc(col+OFFSET)*(blockSize-OFFSET);

  if (whereBloc.second) _value[ind] += loc_val;

  return;
};

// version for block operation
template<typename DataType>
DataType addition(DataType val1, DataType val2)
{return val1+val2;}

template<typename DataType>
void VBRMatr<DataType>::
set_mat_inc(UInt row, UInt col, vector<DataType> &loc_block)
{
  pair<UInt,bool> whereBloc = _Patt->locate_index(row,col);

  typename vector<DataType>::const_iterator loc_start= loc_block.begin();
  typename vector<DataType>::const_iterator loc_end= loc_block.end();

  UInt ind=_Patt->indx()[whereBloc.first];
  //copy of the block
  typename vector<DataType>::iterator val_start= _value.begin()+ind;
  for (typename vector<DataType>::const_iterator ip= loc_start; ip != loc_end;
       ip++)
    {
      *val_start+=*ip;
      val_start++;
    };
  //I haven't been able to use transform (we would avoid the loop):
  //transform(loc_start, loc_end, val_start, val_start, addition<DataType>);

  return;
};

template<typename DataType>
void
VBRMatr<DataType>::
set_mat(UInt row, UInt col, DataType loc_val)
{
  UInt blockSize=_Patt->rpntr()[1]-_Patt->rpntr()[0];
  PatternDefs::Diff_t blocrow=_Patt->rbloc(row);
  PatternDefs::Diff_t bloccol=_Patt->cbloc(col);

  pair<UInt,bool> whereBloc = _Patt->locate_index(blocrow,bloccol);

  UInt ind=_Patt->indx()[whereBloc.first] + _Patt->locr(row+OFFSET) +
    _Patt->locc(col+OFFSET)*(blockSize-OFFSET);

  if (whereBloc.second) _value[ind] = loc_val;

  return;
};

// version for block operation
template<typename DataType>
void VBRMatr<DataType>::
set_mat(UInt row, UInt col, vector<DataType> &loc_block)
{
  pair<UInt,bool> whereBloc = _Patt->locate_index(row,col);

  typename vector<DataType>::const_iterator loc_start= loc_block.begin();
  typename vector<DataType>::const_iterator loc_end= loc_block.end();

  UInt ind=_Patt->indx()[whereBloc.first];
  //copy of the block
  typename vector<DataType>::iterator val_start= _value.begin()+ind;
  //  _value[ind] += loc_val;
  if (whereBloc.second) copy(loc_start, loc_end, val_start);

  return;
};


template<typename DataType>
void
VBRMatr<DataType>::
set_mat(UInt where, DataType loc_val)
{
  _value[where-OFFSET]=loc_val;

  return;
};

template<typename DataType>
void
VBRMatr<DataType>::
spy(string  const &filename)
{
  // Purpose: Matlab dumping and spy
  string nome=filename, uti=" , ";
  UInt nblocrow=_Patt->nRows(), blocsize=_Patt->rpntr()[1]-_Patt->rpntr()[0];
  Container ia=_Patt->ia(), ja=_Patt->ja(), indx=_Patt->indx();
  //
  // check on the file name
  //
  int i=filename.find(".");

  if (i<=0) nome=filename+".m";
  else {
    if ((unsigned int)i!=filename.size()-2  || filename[i+1]!='m')
      {cerr << "Wrong file name ";
      nome=filename+".m";}
  };

  ofstream file_out(nome.c_str());
  ASSERT(file_out,"Error: Output Matrix (Values) file cannot be open");


  file_out << "S = [ ";
  for (UInt irb=0;irb<nblocrow;++irb){
    for (UInt ic=ia[irb];ic<ia[irb+1];++ic)
      for (UInt i=0;i<blocsize;++i)
	for (UInt j=0;j<blocsize;++j)
	  file_out << irb*blocsize+i-OFFSET+1 << uti <<
	    ja[ic]*blocsize+j-OFFSET+1  << uti <<
	    _value[indx[ic]+i+(j-OFFSET)*blocsize] << endl;
  };
  file_out << "];" << endl;
  file_out << "I=S(:,1); J=S(:,2); S=S(:,3); A=sparse(I,J,S); spy(A);"<<endl;
};

/*
template<typename DataType>
void
VBRMatr<DataType>::
showMe() const
{
  UInt blsize=_Patt->rpntr()[1]-_Patt->rpntr()[0]; // block size
  int i_first;
  string pare="[";
  cout << "**************************" << endl;
  cout << "     VBR Matrix Pattern   " << endl;
  cout << endl;

  cout << pare;
  for (PatternDefs::Diff_t i_index=0; i_index <
	 static_cast<PatternDefs::Diff_t>(_Patt->nRows()); i_index++)
    for (UInt jb=0;jb<blsize;jb++){
      cout << pare;
      pare = " [";
      i_first=_Patt->ia()[i_index]-OFFSET+1; // In _ia[i_index] there is the
                                            // diagonal entry
      UInt jj=0;
      for(PatternDefs::Diff_t j=0;j<static_cast<PatternDefs::Diff_t>(_Patt->nCols());j++)
	{
	  if (j==i_index) for (UInt ib=0;ib<blsize;ib++) cout <<
							   _value[];
	  else {
	    if (j==_Patt->ja()[i_first+jj]-OFFSET){
	      for (UInt ib=0;ib<blsize;ib++) cout << " * ";
	      jj++;
	    };
	    else for (UInt ib=0;ib<blsize;ib++) cout << " 0 ";
	  };
	};
      if (i_index==static_cast<PatternDefs::Diff_t>(_Patt->nRows()-1))
	cout << " ]] " << endl;
      else
	cout << " ]  " << endl;
      cout << endl;
    };
  cout << "**************************" << endl;
  return;
};
*/

//-------------------------------------------------------------------------------------------------------
// MSR - VALUES
//-------------------------------------------------------------------------------------------------------

template<class DataType>
MSRMatr<DataType>::MSRMatr():_Patt(0){};

template<typename DataType>
MSRMatr<DataType>::
MSRMatr(const MSRPatt &ex_pattern)
{
  _Patt = &ex_pattern;
  _value.resize(ex_pattern.nNz()+1);
}

// constructor from CSR Pattern
template<typename DataType>
MSRMatr<DataType>::
MSRMatr(const CSRPatt &ex_pattern)
{
  _Patt = &ex_pattern;
  _value.resize(ex_pattern.nNz()+1);
}

/* template<class DataType> */
/* CSRMatr<DataType>:: */
/* CSRMatr(UInt ex_nnz, UInt ex_nrows, UInt ex_ncols, const vector<UInt> &ex_ia, const vector<UInt> &ex_ja, const vector<DataType> &ex_value): */
/* _value(ex_nnz) */
/*  { */
/*   _Patt = 0; */
/*   Csr_Lv_P Pattern(UInt ex_nnz, UInt ex_nrows, UInt ex_ncols, const vector<UInt> &ex_ia, const vector<UInt> &ex_ja); */
/*   _nnz = ex_nnz; */
/*   _nrows = ex_nrows; */
/*   _ncols = ex_ncols; */
/*   //  Pattern.ShowMe(); */
/*   _Patt = &Pattern; */
/*   assert( ex_nnz == ex_value.size()); */
/*   _value = ex_value;    */
/*  }; */


template<class DataType>
MSRMatr<DataType>::
MSRMatr(const MSRPatt* ex_pattern, const vector<DataType> &ex_value):
  _Patt(ex_pattern),_value(ex_value)
{
  //  _nnz = ex_pattern.nnz();
  //  _nrows = ex_pattern.nrow();
  //  _ncols = ex_pattern.ncol();
  ASSERT( _Patt->nNz() == ex_value.size() - 1, "Error in MSR Values "); // in MSR value has lenghth nnz+1
}

template<class DataType>
MSRMatr<DataType>::
MSRMatr(const MSRMatr<DataType> &RightHandMSR):
  _Patt(RightHandMSR.Patt()),_value(RightHandMSR.value()) {};

template<class DataType>
MSRMatr<DataType>&
MSRMatr<DataType>::operator= (const MSRMatr<DataType> &RhMsr)
{
  if (&RhMsr != this)
    {
      _value = RhMsr.value();
      _Patt = RhMsr.Patt();
    }

  return *this;
};

// conversion CSR->MSR, the pattern must have been already converted
// therefore no check is done on the matrix dimensions !
void CSRmat2MSRmat(MSRMatr<double> &MSRmat,
		   CSRMatr<CSRPatt,double> const &CSRmat);

template<class DataType>
vector<DataType>
MSRMatr<DataType>::operator* (const vector<DataType> &v) const
{
  UInt nrows=_Patt->nRows();
  ASSERT(nrows==v.size(),"Error in Matrix Vector product");
  vector<DataType> ans;
  ans.resize(nrows,0.0);
  for (UInt i=0+OFFSET;i<nrows+OFFSET;++i){
    ans[i]=_value[i]*v[i];
    for (UInt j=_Patt->give_bindx()[i];j < _Patt->give_bindx()[i+1];++j)
      ans[i]+=_value[j]*v[_Patt->give_bindx()[j]];
  }
  return ans;
};

//Matrix-vector product for the class Vector (useful for IML++)
template<class DataType>
Vector
MSRMatr<DataType>::
operator*(const Vector &v) const
{
  UInt nrows=_Patt->nRows();
  ASSERT(nrows==v.size(),"Error in Matrix Vector product");
  Vector ans(nrows);
  ans=0.;

  for (UInt i=0+OFFSET;i<nrows+OFFSET;++i){
    ans(i)=_value[i]*v(i);
    for (UInt j=_Patt->give_bindx()[i];j < _Patt->give_bindx()[i+1];++j)
      ans(i)+=_value[j]*v(_Patt->give_bindx()[j]);
  }
  return ans;
}

// Version for C pointer vector. BEWARE: no check on bounds is done !
template<class DataType>
void operMatVec(DataType * const mv,
		const MSRMatr<DataType> &Mat,
		const DataType *v)
{
  UInt nrows=Mat._Patt->nRows();

  for (UInt i=0+OFFSET;i<nrows+OFFSET;++i){
    mv[i]=Mat._value[i]*v[i];
    for (UInt j=Mat._Patt->give_bindx()[i];j < Mat._Patt->give_bindx()[i+1];++j)
      mv[i]+=Mat._value[j]*v[Mat._Patt->give_bindx()[j]];
  }
}

//necessary for IML++ library: transpose-matrix by vector product
template<class DataType>
Vector
MSRMatr<DataType>::
trans_mult(const Vector &v) const
{
  UInt nrows=_Patt->nRows();
  ASSERT(nrows==v.size(),"Error in Matrix Vector product");
  Vector ans(nrows);
  ans=0.;

  for (UInt i=0+OFFSET;i<nrows+OFFSET;++i){
    ans(i)=_value[i]*v(i);
    for (UInt j=_Patt->give_bindx()[i];j < _Patt->give_bindx()[i+1];++j)
      ans(_Patt->give_bindx()[j])+=_value[j]*v(i);
  }
  return ans;
}

template<class DataType>
MSRMatr<DataType>&
MSRMatr<DataType>::operator*=(const DataType num)
{
  UInt stop=_Patt->nNz()+1;
  for (UInt i=0;i<stop;++i) _value[i]*=num;
  return *this;
};

template<class DataType>
MSRMatr<DataType>
MSRMatr<DataType>::operator*(const DataType num)
{
  UInt stop=_Patt->nNz();
  MSRMatr<DataType> ans(*this);

  for (UInt i=0;i<stop;++i)
    ans.set_mat(i,num*_value[i]);

  return ans;
};

template<class DataType>
MSRMatr<DataType>&
MSRMatr<DataType>::flop(const DataType num, MSRMatr<DataType>& M, MSRMatr<DataType>& A)
{
  /* AIM: Matrix = num*M + A */
  UInt stop=_Patt->nNz();
  //   ASSERT(M.Patt()==A.Patt,"Error in summing matrices");

  for (UInt i=0;i<stop;++i)
    _value[i]=num*M.value()[i]+A.value()[i];

  return *this;
}

template<class DataType>
MSRMatr<DataType>&
MSRMatr<DataType>::flop(const DataType num, MSRMatr<DataType>& M)
{
  /* AIM: Matrix = num*M + Matrix */
  UInt stop=_Patt->nNz();
  //   ASSERT(M.Patt()==A.Patt,"Error in summing matrices");

  for (UInt i=0;i<stop;++i)
    _value[i]+=num*M.value()[i];

  return *this;
}

/*
template<class DataType>
MSRMatr<DataType>&
MSRMatr<DataType>::operator* (const DataType num, MSRMatr<DataType>& M)
  {
    UInt stop=_Patt->nNz();

    for (UInt i=0;i<stop;++i)
      _value[i]=num*M.value()[i];

   return *this;
  };

      */

template<class DataType>
void MSRMatr<DataType>::ShowMe()
{
  Container::iterator found;
  int i_first,i_last;
  vector<UInt> vec_temp = _Patt->bindx();
  UInt _nrows= _Patt->nRows();
  UInt _ncols= _Patt->nCols();
  UInt _nnz  = _Patt->nNz();

  string pare="[";
  cout << "**************************" << endl;
  cout << "     MSR Matrix           " << endl;
  cout << endl;
  cout << pare;
  for (UInt i_index=0; i_index < _nrows; ++i_index)
    {
      cout << pare;
      pare = " [";
      i_first=_Patt->bindx()[i_index];
      i_last =_Patt->bindx()[i_index+1];
      //      cout << i_first << " " << i_last << endl;
      for(int j=0;j<_ncols;++j)
	{
	  if (j==i_index)
	    cout << " " << _value[i_index] << " ";
	  else
	    {
	      found = find(&vec_temp[i_first],&vec_temp[i_last],j); // I am not supposing any given ordering in _ja
	      if (found==&vec_temp[i_last])
		cout << " 0 ";
	      else
		{
		  UInt j_index=found-&vec_temp[0];
		  cout << " " << _value[j_index] << " ";
		}
	    }
        }
      if (i_index==_nrows-1)
        cout << " ]] " << endl;
      else
        cout << " ]  " << endl;
    }
  cout << "nnz = " << _nnz << ", nrow = " << _nrows << ", ncol = " << _ncols << endl;
  return;
};

template<typename DataType>
void
MSRMatr<DataType>::spy(string  const &filename)
{
  // Purpose: Matlab dumping and spy
  string nome=filename, uti=" , ";
  //
  // check on the file name
  //
  int i=filename.find(".");

  if (i<=0) nome=filename+".m";
  else {
    if ((unsigned int)i!=filename.size()-2  || filename[i+1]!='m'){
      cerr << "Wrong file name ";
      nome=filename+".m";}
  }

  ofstream file_out(nome.c_str());

  file_out << "S = [ ";
  for (UInt i=1;i<=_Patt->nRows();++i){
    file_out << i << uti << i << uti << _value[i-1] << endl;
    for (UInt ii=_Patt->give_bindx()[i-1];ii<_Patt->give_bindx()[i];++ii)
      file_out << i << uti << _Patt->give_bindx()[ii]+1-OFFSET <<
	uti << _value[ii] << endl;
  }

    file_out << "];" << endl;

  file_out << "I=S(:,1); J=S(:,2); S=S(:,3); A=sparse(I,J,S); spy(A);"<<endl;
};

template<typename DataType>
void
MSRMatr<DataType>::
set_mat_inc(UInt row, UInt col, DataType loc_val)
{

  pair<UInt,bool> where = _Patt->locate_index(row,col);
  if (where.second)
    _value[where.first] += loc_val;
  else {
    cout << row +1<< "," << col +1<< endl;
    ERROR_MSG("problem in MSR::set_mat_inc");
  }
  return;
};

template<typename DataType>
void
MSRMatr<DataType>::
set_mat(UInt row, UInt col, DataType loc_val)
{

  pair<UInt,bool> where = _Patt->locate_index(row,col);
  if(where.second) _value[where.first] = loc_val;

  return;
};

template<typename DataType>
void
MSRMatr<DataType>::
set_mat(UInt where, DataType loc_val)
{
  _value[where-OFFSET]=loc_val;

  return;
};

// determine the lumped matrix for P1
template<typename DataType>
vector<DataType> MSRMatr<DataType>::MassDiagP1() const
{
  UInt nrows = _Patt->nRows();

  vector<DataType> diag(nrows);

  for (UInt nrow = 0; nrow < nrows; ++nrow)
    {
      diag[nrow] = _value[nrow];
      for (UInt ii = _Patt->bindx()[nrow]; ii < _Patt->bindx()[nrow+1]; ++ii)
	diag[nrow] += _value[ii];
    };

  return diag;
}

// build the inverse of a diagonal matrix and multiply it by another matrix
void MultInvDiag(const vector<Real> &Diag, const MSRMatr<Real> &Mat, MSRMatr<Real> &ans) ;

// give the diagonal of an MSR matrix
template<typename DataType>
vector<DataType> MSRMatr<DataType>::giveDiag() const
{
  UInt nrows = _Patt->nRows();

  vector<DataType> diag(nrows);

  for (UInt nrow = 0; nrow < nrows; ++nrow)
    diag[nrow] = _value[nrow];

  return diag;
}

template<typename DataType>
void
MSRMatr<DataType>::
diagonalize_row(UInt const r, DataType const coeff)
{
  _value[r-OFFSET] = coeff;

  UInt i;
  UInt start= *(_Patt->give_bindx().begin()+r-OFFSET);
  UInt end = *(_Patt->give_bindx().begin()+r+1-OFFSET);
  //  UInt disp = _Patt->nRows()+1;

  for (i=start;i<end;++i) {
    _value[i-OFFSET]=0;
    // Miguel: Why this ?.
    //_value[_Patt->give_ybind()[i-disp]] = 0;
  }

  return;
}

template<typename DataType>
void
MSRMatr<DataType>::
diagonalize(UInt const r, DataType const coeff, vector<DataType> &b, DataType datum)
{
  // AIM: Diagonalization of a row of the system, by setting:
  // A(i,i) = coeff,
  // A(i,j) = 0,  A(j,i) = 0 for j!=i
  // and suitably correcting the right hand side of the system
  _value[r-OFFSET] = coeff;

  UInt istart= *(_Patt->give_bindx().begin()+r-OFFSET);
  UInt iend  = *(_Patt->give_bindx().begin()+r+1-OFFSET);

  typename vector<DataType>::iterator start=_value.begin() + istart;
  typename vector<DataType>::iterator end=_value.begin() + iend;
  UInt disp = _Patt->nRows()+1;
  UInt row,col;

  transform(start,end,start,nihil);

  for (UInt i=istart;i<iend;++i)
    {
      row= _Patt->give_bindx()[i]-OFFSET;
      col= _Patt->give_ybind()[i-disp]-OFFSET;
      b[row] -= _value[col] * datum;
      _value[col]=0.;
    }

  b[r-OFFSET] = coeff*datum;

  //Remark: in processing a list of Dirichlet nodes, there is no need to check a posteriori the right hand side
  return;
}

//version for type Vector
template<typename DataType>
void
MSRMatr<DataType>::
diagonalize(UInt const r, DataType const coeff, Vector &b, DataType datum)
{
  // AIM: Diagonalization of a row of the system, by setting:
  // A(i,i) = coeff,
  // A(i,j) = 0,  A(j,i) = 0 for j!=i
  // and suitably correcting the right hand side of the system
  _value[r-OFFSET] = coeff;

  UInt istart= *(_Patt->give_bindx().begin()+r-OFFSET);
  UInt iend  = *(_Patt->give_bindx().begin()+r+1-OFFSET);

  typename vector<DataType>::iterator start=_value.begin() + istart;
  typename vector<DataType>::iterator end=_value.begin() + iend;

  UInt row,col;

  transform(start,end,start,nihil);


  // Miguel: There is a buh using ybind. Alex, did you fix it?.
  // This code works without ybind.
  //
  for (UInt i=istart;i<iend;++i) {
      row = _Patt->give_bindx()[i]-OFFSET;
      UInt Rstart= *(_Patt->give_bindx().begin()+row-OFFSET);
      UInt Rend  = *(_Patt->give_bindx().begin()+row+1-OFFSET);
      for (UInt j=Rstart;j<Rend;++j) {
	if ( (_Patt->give_bindx()[j]-OFFSET) == r ) {
	  col = j;
	  b[row-OFFSET] -= _value[col] * datum;
	  _value[col]=0.;
	  break;
	}
      }
  }


  b[r-OFFSET] = coeff*datum;

  //Remark: in processing a list of Dirichlet nodes, there is no need to check a posteriori the right hand side
  return;
}


template<typename DataType>
void MSRMatr<DataType>::
zeros()
{
  typename vector<DataType>::iterator start= _value.begin();
  typename vector<DataType>::iterator end  = _value.end();
  fill(start, end, 0.0);
}

double nihil(double val);


//-----------------------------------------------------------------------
// MixedMatr
//-----------------------------------------------------------------------

//Default Constructor
template<UInt BRows, UInt BCols, typename PatternType, typename DataType>
MixedMatr<BRows, BCols, PatternType, DataType>::
MixedMatr()
{
  for (UInt i=0; i<BRows; i++)
    for (UInt j=0; j<BCols; j++)
      _bValues[i][j].resize(0);
}

//Constructor from an existing external pattern.
//version for MSR format: size of values= nnz+1.
template<UInt BRows, UInt BCols, typename PatternType, typename DataType>
MixedMatr<BRows, BCols, PatternType, DataType>::
MixedMatr(const MixedPattern<BRows, BCols, MSRPatt> &ex_pattern)
{
  _Patt = &ex_pattern;

  for (UInt i=0; i<BRows; i++)
    for (UInt j=0; j<BCols; j++)
      _bValues[i][j].resize(ex_pattern.nNz(i,j)+1);
}
//Constructor from an existing external pattern.
//version for CSR format: size of values= nnz.
template<UInt BRows, UInt BCols, typename PatternType, typename DataType>
MixedMatr<BRows, BCols, PatternType, DataType>::
MixedMatr(const MixedPattern<BRows, BCols, CSRPatt> &ex_pattern)
{
  _Patt = &ex_pattern;

  for (UInt i=0; i<BRows; i++)
    for (UInt j=0; j<BCols; j++)
      _bValues[i][j].resize(ex_pattern.nNz(i,j));
}

// Determines the lumped diagonal of P1 mass matrix.
template<UInt BRows, UInt BCols, typename PatternType, typename DataType>
vector<DataType>
MixedMatr<BRows, BCols, PatternType, DataType>::
MassDiagP1() const
{
  UInt nrows = _Patt->nRows();
  UInt nnz = 0, nr_global = 0;
  Container coldata, position;

  vector<DataType> diag;
  diag.resize(nrows,0.0);

  for (UInt ib=0; ib< BRows; ib++){

    nrows = _Patt->nRows(ib,0);

    for (UInt jb=0; jb< BCols; jb++){

      if (_Patt->block_ptr(ib,jb) != 0 ){

	for (UInt nr=0; nr < nrows; nr++){
	  coldata.resize(_Patt->nCols(ib,jb));
	  position.resize(_Patt->nCols(ib,jb));
	  nnz= _Patt->row(ib,jb,nr,coldata.begin(),position.begin());
	  for (UInt jcol=0; jcol < nnz; jcol++)
	    diag[nr_global+nr] +=
	      _bValues[ib][jb][position[jcol]];
	}
      }
    }
    nr_global += nrows;
  }

  return diag;
}

// Inverts the diagonal matrix Diag and mulitply it by the matrix Mat.
template<UInt BR, UInt BC, typename PattType>
void
MultInvDiag(const vector<Real> &Diag,
	    const MixedMatr<BR,BC,PattType,Real> &Mat,
	    MixedMatr<BR,BC,PattType,Real> &ans)
{
  ASSERT(find(Diag.begin(),Diag.end(),0) == Diag.end(), "Diagonal matrix Diag must be invertible");

  ASSERT(Diag.size() == Mat._Patt->nRows(), "Matrix sizes not compatible");

  // Product:
  UInt nrows = 0;
  UInt nnz = 0, nr_global = 0;
  Container coldata, pos;

  for (UInt ib=0; ib< BR; ib++){

    nrows = Mat._Patt->nRows(ib,0);

    for (UInt jb=0; jb< BC; jb++){

      if (Mat._Patt->block_ptr(ib,jb) != 0 ){

	for (UInt nr=0; nr < nrows; nr++){
	  coldata.resize(Mat._Patt->nCols(ib,jb));
	  pos.resize(Mat._Patt->nCols(ib,jb));
	  nnz= Mat._Patt->row(ib,jb,nr,coldata.begin(),pos.begin());

	  for (UInt jcol=0; jcol < nnz; jcol++)
	    ans._bValues[ib][jb][pos[jcol]]=
	      ( Mat._bValues[ib][jb][pos[jcol]]/Diag[nr_global+nr] );
	}
      }
    }
    nr_global += nrows;
  }
}

// Gives the diagonal of a block matrix.
template<UInt BRows, UInt BCols, typename PatternType, typename DataType>
vector<DataType>
MixedMatr<BRows, BCols, PatternType, DataType>::
giveDiag() const
{
  ASSERT_PRE(BRows == BCols, "block matrix must be a square matrix");

  vector<Real> diag;
  diag.resize(_Patt->nRows(),0.0);

  Container coldata, pos;
  UInt nrows= 0, nr_global= 0, nnz= 0;

  for (UInt ib=0; ib < BRows; ib++){
    ASSERT_PRE(_Patt->nRows(ib,ib) == _Patt->nCols(ib,ib),
	       "matrix must have nrows=ncols for diagonal blocks");

    nrows= _Patt->nRows(ib,ib);

    if (_Patt->block_ptr(ib,ib) != 0 )
      for (UInt nr= 0; nr < nrows; nr++){

	coldata.resize(_Patt->nCols(ib,ib));
	pos.resize(_Patt->nCols(ib,ib));

	nnz= _Patt->row(ib,ib,nr,coldata.begin(),pos.begin());

	UInt jcol = 0;

	if (coldata[jcol]-OFFSET == nr)
	  diag[nr_global+nr]= _bValues[ib][ib][pos[jcol]];
	else
	  {
	    while ( coldata[jcol]-OFFSET != nr ) jcol++;
	    diag[nr_global+nr]= _bValues[ib][ib][pos[jcol]];
	  }
      }
    nr_global += nrows;
  }
  return diag;
}


// Warning: the two matrices will point to the same pattern.
template<UInt BRows, UInt BCols, typename PatternType, typename DataType>
MixedMatr<BRows, BCols, PatternType, DataType> &
MixedMatr<BRows, BCols, PatternType, DataType>::
operator=(const MixedMatr<BRows,BCols,PatternType,DataType> &RhMtr)
{
  if (&RhMtr != this)
    {
      UInt ib,jb;
      for (ib= 0; ib < BRows; ib++)
	for (jb= 0; jb < BCols; jb++)
	  _bValues[ib][jb] = RhMtr.bValues(ib,jb);

      _Patt = RhMtr.Patt();
    }

  return *this;

}

// Assigns the matrix to loc_val at the place of index where.
template<UInt BRows, UInt BCols, typename PatternType, typename DataType>
void
MixedMatr<BRows, BCols, PatternType, DataType>::
set_mat(UInt ib, UInt jb, UInt where, DataType loc_val)
{
  _bValues[ib][jb][where-OFFSET] = loc_val;
}

// Assigns the matrix element (row,col) to loc_val.
template<UInt BRows, UInt BCols, typename PatternType, typename DataType>
void
MixedMatr<BRows, BCols, PatternType, DataType>::
set_mat(UInt row, UInt col, DataType loc_val)
{
  UInt m,n,lr,lc;
  extract_pair(_Patt->locateElBlock(row, col),m,n);
  extract_pair(_Patt->localNumber(m, n, row, col), lr, lc);
  set_mat(m, n, lr, lc, loc_val);
}

// Assigns the matrix element (row,col) of block (ib,jb) to loc_val.
template<UInt BRows, UInt BCols, typename PatternType, typename DataType>
void
MixedMatr<BRows, BCols, PatternType, DataType>::
set_mat(UInt ib, UInt jb, UInt row, UInt col, DataType loc_val)
{
  if (_Patt->block_ptr(ib,jb) != 0 )
    {
      pair<UInt,bool> where = _Patt->block_ptr(ib,jb)->locate_index(row,col);
      if (where.second) _bValues[ib][jb][where.first] = loc_val;
    }
}

// Adds loc_val to the matrix element (row,col).
template<UInt BRows, UInt BCols, typename PatternType, typename DataType>
void
MixedMatr<BRows, BCols, PatternType, DataType>::
set_mat_inc(UInt row, UInt col, DataType loc_val)
{
  UInt m,n,lr,lc;
  extract_pair(_Patt->locateElBlock(row, col), m, n);
  extract_pair(_Patt->localNumber(m, n, row, col), lr, lc);
  set_mat_inc(m, n, lr, lc, loc_val);
}

// Adds loc_val to the matrix element (row,col) of block (ib,jb).
template<UInt BRows, UInt BCols, typename PatternType, typename DataType>
void
MixedMatr<BRows, BCols, PatternType, DataType>::
set_mat_inc(UInt ib, UInt jb, UInt row, UInt col, DataType loc_val)
{
  if (_Patt->block_ptr(ib,jb) != 0 )
    {
      pair<UInt,bool> where = _Patt->block_ptr(ib,jb)->locate_index(row,col);
      if (where.second) _bValues[ib][jb][where.first] += loc_val;
    }
}

// Returns the matrix element (i,j) value.
template<UInt BRows, UInt BCols, typename PatternType, typename DataType>
DataType
MixedMatr<BRows, BCols, PatternType, DataType>::
get_value(UInt i, UInt j)
{
  UInt m,n,lr,lc;
  extract_pair(_Patt->locateElBlock(i,j), m, n);
  extract_pair(_Patt->localNumber(m, n, i, j), lr, lc);
  return get_value(m, n, lr, lc);
}
// const qualifyer version
template<UInt BRows, UInt BCols, typename PatternType, typename DataType>
const DataType
MixedMatr<BRows, BCols, PatternType, DataType>::
get_value(UInt i, UInt j) const
{
  UInt m,n,lr,lc;
  extract_pair(_Patt->locateElBlock(i ,j), m, n);
  extract_pair(_Patt->localNumber(m, n, i, j), lr, lc);
  return get_value(m, n, i, j);
}

// Returns the matrix element (i,j) value of block (ib,jb).
template<UInt BRows, UInt BCols, typename PatternType, typename DataType>
DataType
MixedMatr<BRows, BCols, PatternType, DataType>::
get_value(UInt ib, UInt jb, UInt i, UInt j)
{
  if (_Patt->block_ptr(ib,jb) != 0)
    return _bValues[ib][jb][_Patt->block_ptr(ib,jb)->locate_index(i,j).first];
  else
    return 0.0;
}
// const qualifyer version
template<UInt BRows, UInt BCols, typename PatternType, typename DataType>
const DataType
MixedMatr<BRows, BCols, PatternType, DataType>::
get_value(UInt ib, UInt jb, UInt i, UInt j) const
{
  if (_Patt->block_ptr(ib,jb) != 0)
    return _bValues[ib][jb][_Patt->block_ptr(ib,jb)->locate_index(i,j).first];
  else
    return 0.0;
}

// Matrix visualization a la matlab.
template<UInt BRows, UInt BCols, typename PatternType, typename DataType>
void
MixedMatr<BRows, BCols, PatternType, DataType>::
spy(string  const &filename)
{
  // Purpose: Matlab dumping and spy
  string nome=filename, uti=" , ";
  //
  // check on the file name
  //
  UInt i=filename.find(".");

  if (i<=0) nome=filename+".m";
  else {
    if (i!=filename.size()-2  || filename[i+1]!='m'){
      cerr << "Wrong file name ";
      nome=filename+".m";}
  }

  ofstream file_out(nome.c_str());
  UInt nnz, mb, nb;
  Container coldata, pos;
  coldata.resize(_Patt->nCols());
  pos.resize(_Patt->nCols());

  file_out << "S = [ ";
  for (UInt i=0; i<_Patt->nRows(); ++i){
    nnz= _Patt->row(i, coldata.begin(), pos.begin());
    for (UInt j=0; j< nnz; ++j){
      extract_pair(_Patt->locateElBlock(i,coldata[j]-OFFSET),mb,nb);
      file_out << i+1 << uti << coldata[j]+1-OFFSET << uti <<
	_bValues[mb][nb][pos[j]] << endl;
    }
  }
  file_out << "];" << endl;
  file_out << "I=S(:,1); J=S(:,2); S=S(:,3); A=sparse(I,J,S); spy(A);"<<endl;
}

// Assigns matrix diagonal element (r,r) to coeff, other elts
// of row r to zero.
template<UInt BRows, UInt BCols, typename PatternType, typename DataType>
void
MixedMatr<BRows, BCols, PatternType, DataType>::
diagonalize_row( UInt const r, DataType const coeff)
{
  UInt nnz= 0, m=0, n=0;
  Container coldata, pos;
  coldata.resize(_Patt->nCols());
  pos.resize(_Patt->nCols());

  nnz= _Patt->row(r-OFFSET,coldata.begin(),pos.begin());

  UInt jcol=0;
  if (coldata[jcol] == r){
    //diagonal element:
    extract_pair(_Patt->locateElBlock(r-OFFSET ,coldata[jcol]-OFFSET), m, n);
    _bValues[m][n][pos[jcol]]= coeff;
    //other elements:
    for (jcol=1; jcol < nnz; jcol++){
      extract_pair(_Patt->locateElBlock(r-OFFSET ,coldata[jcol]-OFFSET), m, n);
      _bValues[m][n][pos[jcol]]= 0.0;
    }
  }
  else
    for (jcol=0; jcol < nnz; jcol++){
      extract_pair(_Patt->locateElBlock(r-OFFSET ,coldata[jcol]-OFFSET), m, n);
      if (coldata[jcol] == r)
	//diagonal element:
	_bValues[m][n][pos[jcol]]= coeff;
      else
	//other elements
	_bValues[m][n][pos[jcol]]= 0.0;
    }
}

// assign a row to zero. Remark, zero might be defined for any DataType
template<UInt BRows, UInt BCols, typename PatternType, typename DataType>
void
MixedMatr<BRows, BCols, PatternType, DataType>::
zero_row(UInt const row)
{
  diagonalize_row(row, 0.0);
}

// Assigns matrix diagonal element (r,r) to coeff, other elts
// of row r to zero, and vector b element b(r) to coeff*datum.
template<UInt BRows, UInt BCols, typename PatternType, typename DataType>
template< typename VectorType >
void
MixedMatr<BRows, BCols, PatternType, DataType>::
diagonalize( UInt const row, DataType const coeff, VectorType &b,
	     DataType datum)
{
  UInt nnz= 0, m=0, n=0;
  Container coldata, pos;
  coldata.resize(_Patt->nCols());
  pos.resize(_Patt->nCols());

  nnz= _Patt->row(row-OFFSET,coldata.begin(),pos.begin());

  UInt jcol=0;
  if (coldata[jcol] == row){ //case diagfirst
    //diagonal element:
    extract_pair(_Patt->locateElBlock(row-OFFSET ,coldata[jcol]-OFFSET), m, n);
    _bValues[m][n][pos[jcol]]= coeff;
    //other elements:
    for (jcol=1; jcol < nnz; jcol++){
      extract_pair(_Patt->locateElBlock(row-OFFSET ,coldata[jcol]-OFFSET), m, n);
      _bValues[m][n][pos[jcol]]= 0.0;
    }
  }
  else
    for (jcol=0; jcol < nnz; jcol++){
      extract_pair(_Patt->locateElBlock(row-OFFSET ,coldata[jcol]-OFFSET), m, n);
      if (coldata[jcol] == row)
	//diagonal element:
	_bValues[m][n][pos[jcol]]= coeff;
      else
	//other elements
	_bValues[m][n][pos[jcol]]= 0.0;
    }

  //rhs:
  b[r-OFFSET]= coeff*datum;
}

// Matrix-vector product.
template<UInt BRows, UInt BCols, typename PatternType, typename DataType>
vector<DataType>
MixedMatr<BRows, BCols, PatternType, DataType>::
operator*(const vector<DataType> &v) const
{
  ASSERT(_Patt->nCols()==v.size(),"Error in Matrix Vector product");
  vector<DataType> ans;
  ans.resize(_Patt->nRows(),0.0);

  UInt nrows = 0, ncols = 0;
  UInt nnz = 0, nr_global = 0, nc_global;
  Container coldata, pos;

  for (UInt ib=0; ib< BRows; ib++){

    nrows = _Patt->nRows(ib,0);
    nc_global = 0;

    for (UInt jb=0; jb< BCols; jb++){

      ncols = _Patt->nCols(ib,jb);

      if (_Patt->block_ptr(ib,jb) != 0 ){

	for (UInt nr=0; nr < nrows; nr++){
	  coldata.resize(ncols);
	  pos.resize(ncols);
	  nnz= _Patt->row(ib,jb,nr,coldata.begin(),pos.begin());

	  for (UInt jcol=0; jcol < nnz; jcol++)
	    ans[nr_global+nr]+=_bValues[ib][jb][pos[jcol]]*
	      v[nc_global+coldata[jcol]-OFFSET];
	}
      }
      nc_global += ncols;
    }
    nr_global += nrows;
  }
  return ans;
}
// version for type Vector
template<UInt BRows, UInt BCols, typename PatternType, typename DataType>
Vector
MixedMatr<BRows, BCols, PatternType, DataType>::
operator*(const Vector &v) const
{
  ASSERT(_Patt->nCols()==v.size(),"Error in Matrix Vector product");
  Vector ans(_Patt->nRows());
  ans=0.0;

  UInt nrows = 0, ncols = 0;
  UInt nnz = 0, nr_global = 0, nc_global;
  Container coldata, pos;

  for (UInt ib=0; ib< BRows; ib++){

    nrows = _Patt->nRows(ib,0);
    nc_global= 0;

    for (UInt jb=0; jb< BCols; jb++){

      ncols = _Patt->nCols(ib,jb);

      if (_Patt->block_ptr(ib,jb) != 0 ){

	for (UInt nr=0; nr < nrows; nr++){
	  coldata.resize(ncols);
	  pos.resize(ncols);
	  nnz= _Patt->row(ib,jb,nr,coldata.begin(),pos.begin());

	  for (UInt jcol=0; jcol < nnz; jcol++)
	    ans[nr_global+nr]+= _bValues[ib][jb][pos[jcol]]*
	      v[nc_global+coldata[jcol]-OFFSET];
	}
      }
      nc_global += ncols;
    }
    nr_global += nrows;
  }
  return ans;
}
// Version for C pointer vector. BEWARE: no check on bounds is done !
template<typename DataType, UInt BRows, UInt BCols, typename PatternType>
void operMatVec(DataType * const mv,
		const MixedMatr<BRows, BCols, PatternType, DataType> &Mat,
		const DataType *v)
{
  UInt nrows = 0, ncols = 0;
  UInt nnz = 0, nr_global = 0, nc_global;
  Container coldata, pos;

  for (UInt ib=0; ib< BRows; ib++){

    nrows = Mat._Patt->nRows(ib,0);
    nc_global= 0;

    for (UInt nr=0; nr < nrows; nr++){

      //initialize
      mv[nr+ib*nrows]= 0.;

      for (UInt jb=0; jb< BCols; jb++){

	ncols = Mat._Patt->nCols(ib,jb);

	if (Mat._Patt->block_ptr(ib,jb) != 0 ){

	  coldata.resize(ncols);
	  pos.resize(ncols);
	  nnz= Mat._Patt->row(ib,jb,nr,coldata.begin(),pos.begin());

	  for (UInt jcol=0; jcol < nnz; jcol++)
	    mv[nr_global+nr]+= Mat._bValues[ib][jb][pos[jcol]]*
	      v[nc_global+coldata[jcol]-OFFSET];
	}
      }
      nc_global += ncols;
    }
    nr_global += nrows;
  }
}

//necessary for IML++ library: transpose-matrix by vector product.
template<UInt BRows, UInt BCols, typename PatternType, typename DataType>
Vector
MixedMatr<BRows, BCols, PatternType, DataType>::
trans_mult(const Vector &v) const
{
  ASSERT(_Patt->nRows()==v.size(),"Error in Matrix Vector product");
  Vector ans(_Patt->nRows());
  ans=0.0;

  UInt nrows = 0, ncols = 0;
  UInt nnz = 0, nr_global = 0, nc_global;
  Container coldata, pos;

  for (UInt ib=0; ib< BRows; ib++){

    nrows = _Patt->nRows(ib,0);
    nc_global=0;

    for (UInt jb=0; jb< BCols; jb++){

      ncols = _Patt->nCols(ib,jb);

      if (_Patt->block_ptr(ib,jb) != 0 ){

	for (UInt nr=0; nr < nrows; nr++){
	  coldata.resize(ncols);
	  pos.resize(ncols);
	  nnz= _Patt->row(ib,jb,nr,coldata.begin(),pos.begin());

	  for (UInt jcol=0; jcol < nnz; jcol++)
	    ans[nc_global+coldata[jcol]-OFFSET]+=_bValues[ib][jb][pos[jcol]]*
	      v[nr_global+nr];
	}
      }
      nc_global += ncols;
    }
    nr_global += nrows;
  }
  return ans;
}

template<UInt BRows, UInt BCols, typename PatternType, typename DataType>
void MixedMatr<BRows, BCols, PatternType, DataType>::
zeros()
{
  UInt ib,jb;
  for (ib= 0; ib < BRows; ib++)
    for (jb= 0; jb < BCols; jb++)
      if (_Patt->block_ptr(ib,jb) != 0 ){
	typename vector<DataType>::iterator start= _bValues[ib][jb].begin();
	typename vector<DataType>::iterator end  = _bValues[ib][jb].end();
	fill(start, end, 0.0);
      }
}

//-----------------------------------------------------------------------
// MixedMatr SPECIALIZATION FOR MSR
//-----------------------------------------------------------------------

//Default Constructor
template<UInt BRows, UInt BCols>
MixedMatr<BRows, BCols, MSRPatt, double>::
MixedMatr()
{
  for (UInt i=0; i<BRows; i++)
    for (UInt j=0; j<BCols; j++)
      _bValues[i][j].resize(0);
}

//Constructor from an existing external pattern.
//version for MSR format: size of values= nnz+1.
template<UInt BRows, UInt BCols>
MixedMatr<BRows, BCols, MSRPatt, double>::
MixedMatr(const MixedPattern<BRows, BCols, MSRPatt> &ex_pattern)
{
  _Patt = &ex_pattern;

  for (UInt i=0; i<BRows; i++)
    for (UInt j=0; j<BCols; j++)
      _bValues[i][j].resize(ex_pattern.nNz(i,j)+1);
}

// Determines the lumped diagonal of P1 mass matrix.
template<UInt BRows, UInt BCols>
vector<double>
MixedMatr<BRows, BCols, MSRPatt, double>::
MassDiagP1() const
{
  UInt nrows = _Patt->nRows();
  UInt nnz = 0, nr_global = 0;
  Container coldata, position;

  vector<double> diag;
  diag.resize(nrows,0.0);

  for (UInt ib=0; ib< BRows; ib++){

    nrows = _Patt->nRows(ib,0);

    for (UInt jb=0; jb< BCols; jb++){

      if (_Patt->block_ptr(ib,jb) != 0 ){

	for (UInt nr=0; nr < nrows; nr++){
	  coldata.resize(_Patt->nCols(ib,jb));
	  position.resize(_Patt->nCols(ib,jb));
	  nnz= _Patt->row(ib,jb,nr,coldata.begin(),position.begin());
	  for (UInt jcol=0; jcol < nnz; jcol++)
	    diag[nr_global+nr] +=
	      _bValues[ib][jb][position[jcol]];
	}
      }
    }
    nr_global += nrows;
  }

  return diag;
}

// Inverts the diagonal matrix Diag and mulitply it by the matrix Mat.
template<UInt BR, UInt BC>
void
MultInvDiag(const vector<Real> &Diag,
	    const MixedMatr<BR,BC,MSRPatt,Real> &Mat,
	    MixedMatr<BR,BC,MSRPatt,Real> &ans)
{
  ASSERT(find(Diag.begin(),Diag.end(),0) == Diag.end(), "Diagonal matrix Diag must be invertible");

  ASSERT(Diag.size() == Mat._Patt->nRows(), "Matrix sizes not compatible");

  // Product:
  UInt nrows = 0;
  UInt nnz = 0, nr_global = 0;
  Container coldata, pos;

  for (UInt ib=0; ib< BR; ib++){

    nrows = Mat._Patt->nRows(ib,0);

    for (UInt jb=0; jb< BC; jb++){

      if (Mat._Patt->block_ptr(ib,jb) != 0 ){

	for (UInt nr=0; nr < nrows; nr++){
	  coldata.resize(Mat._Patt->nCols(ib,jb));
	  pos.resize(Mat._Patt->nCols(ib,jb));
	  nnz= Mat._Patt->row(ib,jb,nr,coldata.begin(),pos.begin());

	  for (UInt jcol=0; jcol < nnz; jcol++)
	    ans._bValues[ib][jb][pos[jcol]]=
	      ( Mat._bValues[ib][jb][pos[jcol]]/Diag[nr_global+nr] );
	}
      }
    }
    nr_global += nrows;
  }
}

// Gives the diagonal of a block matrix.
template<UInt BRows, UInt BCols>
vector<double>
MixedMatr<BRows, BCols, MSRPatt, double>::
giveDiag() const
{
  ASSERT_PRE(BRows == BCols, "block matrix must be a square matrix");

  vector<Real> diag;
  diag.resize(_Patt->nRows(),0.0);

  Container coldata, pos;
  UInt nrows= 0, nr_global= 0, nnz= 0;

  for (UInt ib=0; ib < BRows; ib++){
    ASSERT_PRE(_Patt->nRows(ib,ib) == _Patt->nCols(ib,ib),
	       "matrix must have nrows=ncols for diagonal blocks");

    nrows= _Patt->nRows(ib,ib);

    if (_Patt->block_ptr(ib,ib) != 0 )
      for (UInt nr= 0; nr < nrows; nr++){

	coldata.resize(_Patt->nCols(ib,ib));
	pos.resize(_Patt->nCols(ib,ib));

	nnz= _Patt->row(ib,ib,nr,coldata.begin(),pos.begin());

	UInt jcol = 0;

	if (coldata[jcol]-OFFSET == nr)
	  diag[nr_global+nr]= _bValues[ib][ib][pos[jcol]];
	else
	  {
	    while ( coldata[jcol]-OFFSET != nr ) jcol++;
	    diag[nr_global+nr]= _bValues[ib][ib][pos[jcol]];
	  }
      }
    nr_global += nrows;
  }
  return diag;
}


// Warning: the two matrices will point to the same pattern.
template<UInt BRows, UInt BCols>
MixedMatr<BRows, BCols, MSRPatt, double> &
MixedMatr<BRows, BCols, MSRPatt, double>::
operator=(const MixedMatr<BRows,BCols,MSRPatt,double> &RhMtr)
{
  if (&RhMtr != this)
    {
      UInt ib,jb;
      for (ib= 0; ib < BRows; ib++)
	for (jb= 0; jb < BCols; jb++)
	  _bValues[ib][jb] = RhMtr.bValues(ib,jb);

      _Patt = RhMtr.Patt();
    }

  return *this;

}

// Assigns the matrix to loc_val at the place of index where.
template<UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, MSRPatt, double>::
set_mat(UInt ib, UInt jb, UInt where, double loc_val)
{
  _bValues[ib][jb][where-OFFSET] = loc_val;
}

// Assigns the matrix element (row,col) to loc_val.
template<UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, MSRPatt, double>::
set_mat(UInt row, UInt col, double loc_val)
{
  UInt m,n,lr,lc;
  extract_pair(_Patt->locateElBlock(row, col),m,n);
  extract_pair(_Patt->localNumber(m, n, row, col), lr, lc);
  set_mat(m, n, lr, lc, loc_val);
}

// Assigns the matrix element (row,col) of block (ib,jb) to loc_val.
template<UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, MSRPatt, double>::
set_mat(UInt ib, UInt jb, UInt row, UInt col, double loc_val)
{
  if (_Patt->block_ptr(ib,jb) != 0 )
    {
      pair<UInt,bool> where = _Patt->block_ptr(ib,jb)->locate_index(row,col);
      if (where.second) _bValues[ib][jb][where.first] = loc_val;
    }
}

// Adds loc_val to the matrix element (row,col).
template<UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, MSRPatt, double>::
set_mat_inc(UInt row, UInt col, double loc_val)
{
  UInt m,n,lr,lc;
  extract_pair(_Patt->locateElBlock(row, col), m, n);
  extract_pair(_Patt->localNumber(m, n, row, col), lr, lc);
  set_mat_inc(m, n, lr, lc, loc_val);
}

// Adds loc_val to the matrix element (row,col) of block (ib,jb).
template<UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, MSRPatt, double>::
set_mat_inc(UInt ib, UInt jb, UInt row, UInt col, double loc_val)
{
  if (_Patt->block_ptr(ib,jb) != 0 )
    {
      pair<UInt,bool> where = _Patt->block_ptr(ib,jb)->locate_index(row,col);
      if (where.second) _bValues[ib][jb][where.first] += loc_val;
    }
}

// Returns the matrix element (i,j) value.
template<UInt BRows, UInt BCols>
double
MixedMatr<BRows, BCols, MSRPatt, double>::
get_value(UInt i, UInt j)
{
  UInt m,n,lr,lc;
  extract_pair(_Patt->locateElBlock(i,j), m, n);
  extract_pair(_Patt->localNumber(m, n, i, j), lr, lc);
  return get_value(m, n, lr, lc);
}
// const qualifyer version
template<UInt BRows, UInt BCols>
const double
MixedMatr<BRows, BCols, MSRPatt, double>::
get_value(UInt i, UInt j) const
{
  UInt m,n,lr,lc;
  extract_pair(_Patt->locateElBlock(i ,j), m, n);
  extract_pair(_Patt->localNumber(m, n, i, j), lr, lc);
  return get_value(m, n, i, j);
}

// Returns the matrix element (i,j) value of block (ib,jb).
template<UInt BRows, UInt BCols>
double
MixedMatr<BRows, BCols, MSRPatt, double>::
get_value(UInt ib, UInt jb, UInt i, UInt j)
{
  if (_Patt->block_ptr(ib,jb) != 0)
    return _bValues[ib][jb][_Patt->block_ptr(ib,jb)->locate_index(i,j).first];
  else
    return 0.0;
}
// const qualifyer version
template<UInt BRows, UInt BCols>
const double
MixedMatr<BRows, BCols, MSRPatt, double>::
get_value(UInt ib, UInt jb, UInt i, UInt j) const
{
  if (_Patt->block_ptr(ib,jb) != 0)
    return _bValues[ib][jb][_Patt->block_ptr(ib,jb)->locate_index(i,j).first];
  else
    return 0.0;
}

// Matrix visualization a la matlab.
template<UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, MSRPatt, double>::
spy(string  const &filename)
{
  // Purpose: Matlab dumping and spy
  string nome=filename, uti=" , ";
  //
  // check on the file name
  //
  UInt i=filename.find(".");

  if (i<=0) nome=filename+".m";
  else {
    if (i!=filename.size()-2  || filename[i+1]!='m'){
      cerr << "Wrong file name ";
      nome=filename+".m";}
  }

  ofstream file_out(nome.c_str());
  UInt nnz, mb, nb;
  Container coldata, pos;
  coldata.resize(_Patt->nCols());
  pos.resize(_Patt->nCols());

  file_out << "S = [ ";
  for (UInt i=0; i<_Patt->nRows(); ++i){
    nnz= _Patt->row(i, coldata.begin(), pos.begin());
    for (UInt j=0; j< nnz; ++j){
      extract_pair(_Patt->locateElBlock(i,coldata[j]-OFFSET),mb,nb);
      file_out << i+1 << uti << coldata[j]+1-OFFSET << uti <<
	_bValues[mb][nb][pos[j]] << endl;
    }
  }
  file_out << "];" << endl;
  file_out << "I=S(:,1); J=S(:,2); S=S(:,3); A=sparse(I,J,S); spy(A);"<<endl;
}

// Assigns matrix diagonal element (r,r) to coeff, other elts
// of row r to zero.
template<UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, MSRPatt, double>::
diagonalize_row( UInt const r, double const coeff)
{
  UInt nnz= 0, m=0, n=0;
  Container coldata, pos;
  coldata.resize(_Patt->nCols());
  pos.resize(_Patt->nCols());

  nnz= _Patt->row(r-OFFSET,coldata.begin(),pos.begin());

  UInt jcol=0;
  if (coldata[jcol] == r){
    //diagonal element:
    extract_pair(_Patt->locateElBlock(r-OFFSET ,coldata[jcol]-OFFSET), m, n);
    _bValues[m][n][pos[jcol]]= coeff;
    //other elements:
    for (jcol=1; jcol < nnz; jcol++){
      extract_pair(_Patt->locateElBlock(r-OFFSET ,coldata[jcol]-OFFSET), m, n);
      _bValues[m][n][pos[jcol]]= 0.0;
    }
  }
  else
    for (jcol=0; jcol < nnz; jcol++){
      extract_pair(_Patt->locateElBlock(r-OFFSET ,coldata[jcol]-OFFSET), m, n);
      if (coldata[jcol] == r)
	//diagonal element:
	_bValues[m][n][pos[jcol]]= coeff;
      else
	//other elements
	_bValues[m][n][pos[jcol]]= 0.0;
    }
}

// assign a row to zero. Remark, zero might be defined for any double
template<UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, MSRPatt, double>::
zero_row(UInt const row)
{
  diagonalize_row(row, 0.0);
}

// Assigns matrix diagonal element (r,r) to coeff, other elts
// of row r to zero, and vector b element b(r) to coeff*datum.
// version of diagonalize for MSR format (performs the symmetric acces
// to the columns)
// Alain. nov. 2002.
template<UInt BRows, UInt BCols>
template< typename VectorType >
void
MixedMatr<BRows, BCols, MSRPatt, double>::
diagonalize( UInt const row, double const coeff, VectorType &b,
	     double datum)
{
  UInt lr,lc, lr_ybind, lc_ybind;
  UInt m, n, disp;
  UInt r_offset, c_offset;
  UInt col, colstart, colend, row_ybind;

  //block row identification
  extract_pair(_Patt->locateElBlock(row,row),m, n);
  //loop on block columns
  for (n= 0; n < BCols; n++)
    if (_Patt->block_ptr(m,n) != 0 ){
      extract_pair(_Patt->blockOffset(m,n),r_offset, c_offset);
      col= c_offset;
      //local row
      extract_pair(_Patt->localNumber(m, n, row, col), lr, lc);
      //diagonal terms
      if (n == m)
	_bValues[m][n][lr]= coeff; //the true diagonal
      else
	_bValues[m][n][lr]= 0.;
      //loop on the columns
      colstart= _Patt->block_ptr(m,n)->give_bindx()[lr];
      colend= _Patt->block_ptr(m,n)->give_bindx()[lr+1];
      disp= _Patt->nRows(m,n) + 1;
      for ( lc= colstart; lc < colend; lc++)
	{
	  _bValues[m][n][lc]= 0;
	  // for the columns update
	  lr_ybind= _Patt->block_ptr(m,n)->give_bindx()[lc];
	  lc_ybind= _Patt->block_ptr(m,n)->give_ybind()[lc-disp];
	  row_ybind= lr_ybind + r_offset;
	  b[row_ybind] -= _bValues[m][n][lc_ybind] * datum;
	  _bValues[m][n][lc_ybind]=0.;
	}
    }
  //rhs:
  b[row]= coeff*datum;
}


// Matrix-vector product.
template<UInt BRows, UInt BCols>
vector<double>
MixedMatr<BRows, BCols, MSRPatt, double>::
operator*(const vector<double> &v) const
{
  ASSERT(_Patt->nCols()==v.size(),"Error in Matrix Vector product");
  vector<double> ans;
  ans.resize(_Patt->nRows(),0.0);

  UInt nrows, istart,iend;
  UInt r_offset, c_offset;
  UInt ib,jb,ir,ii;

  for (ib=0; ib< BRows; ib++){

    nrows = _Patt->nRows(ib,0);

    for (jb=0; jb< BCols; jb++){

      if (_Patt->block_ptr(ib,jb) != 0 ){

	extract_pair(_Patt->blockOffset(ib,jb),r_offset, c_offset);

	for (ir= 0; ir< nrows; ir++){
	  //"diagonal" term
	  ans[ir+r_offset]+= _bValues[ib][jb][ir]*v[ir+c_offset];

	  istart= _Patt->block_ptr(ib,jb)->give_bindx()[ir]-OFFSET;
	  iend  = _Patt->block_ptr(ib,jb)->give_bindx()[ir+1]-OFFSET;
	  for (ii= istart;ii< iend;++ii){
	    ans[ir+r_offset]+= _bValues[ib][jb][ii]*
	      v[c_offset+_Patt->block_ptr(ib,jb)->give_bindx()[ii]-OFFSET];
	  }
	}
      }
    }
  }
  return ans;
}
// version for type Vector
template<UInt BRows, UInt BCols>
Vector
MixedMatr<BRows, BCols, MSRPatt, double>::
operator*(const Vector &v) const
{
  ASSERT(_Patt->nCols()==v.size(),"Error in Matrix Vector product");
  Vector ans(_Patt->nRows());
  ans=0.0;

  UInt nrows, istart,iend;
  UInt r_offset, c_offset;
  UInt ib,jb,ir,ii;

  for (ib=0; ib< BRows; ib++){

    nrows = _Patt->nRows(ib,0);

    for (jb=0; jb< BCols; jb++){

      if (_Patt->block_ptr(ib,jb) != 0 ){

	extract_pair(_Patt->blockOffset(ib,jb),r_offset, c_offset);

	for (ir= 0; ir< nrows; ir++){
	  //"diagonal" term
	  ans[ir+r_offset]+= _bValues[ib][jb][ir]*v[ir+c_offset];
	  //other terms
	  istart= _Patt->block_ptr(ib,jb)->give_bindx()[ir]-OFFSET;
	  iend  = _Patt->block_ptr(ib,jb)->give_bindx()[ir+1]-OFFSET;
	  for (ii= istart;ii< iend;++ii){
	    ans[ir+r_offset]+= _bValues[ib][jb][ii]*
	      v[c_offset+_Patt->block_ptr(ib,jb)->give_bindx()[ii]-OFFSET];
	  }
	}
      }
    }
  }
  return ans;
}
// Version for C pointer vector. BEWARE: no check on bounds is done !
template<UInt BRows, UInt BCols>
void operMatVec(double * const mv,
		const MixedMatr<BRows, BCols, MSRPatt, double> &Mat,
		const double *v)
{
  UInt nrows, istart,iend;
  UInt r_offset, c_offset;
  UInt ib,jb,ir,ii;

  for (ib=0; ib< BRows; ib++){

    nrows = Mat._Patt->nRows(ib,0);

    for (ir= 0; ir< nrows; ir++){

      //initialize
      mv[ir+ib*nrows]= 0.;

      for (jb=0; jb< BCols; jb++){

	if (Mat._Patt->block_ptr(ib,jb) != 0 ){

	  extract_pair(Mat._Patt->blockOffset(ib,jb),r_offset, c_offset);

	  //"diagonal" term
	  mv[ir+r_offset]+= Mat._bValues[ib][jb][ir]*v[ir+c_offset];
	  //other terms
	  istart= Mat._Patt->block_ptr(ib,jb)->give_bindx()[ir]-OFFSET;
	  iend  = Mat._Patt->block_ptr(ib,jb)->give_bindx()[ir+1]-OFFSET;
	  for (ii= istart;ii< iend;++ii){
	    mv[ir+r_offset]+= Mat._bValues[ib][jb][ii]*
	      v[c_offset+Mat._Patt->block_ptr(ib,jb)->give_bindx()[ii]-OFFSET];
	  }
	}
      }
    }
  }
}

//necessary for IML++ library: transpose-matrix by vector product.
template<UInt BRows, UInt BCols>
Vector
MixedMatr<BRows, BCols, MSRPatt, double>::
trans_mult(const Vector &v) const
{
  ASSERT(_Patt->nRows()==v.size(),"Error in Matrix Vector product");
  Vector ans(_Patt->nRows());
  ans=0.0;

  UInt nrows = 0, ncols = 0;
  UInt nnz = 0, nr_global = 0, nc_global;
  Container coldata, pos;

  for (UInt ib=0; ib< BRows; ib++){

    nrows = _Patt->nRows(ib,0);
    nc_global=0;

    for (UInt jb=0; jb< BCols; jb++){

      ncols = _Patt->nCols(ib,jb);

      if (_Patt->block_ptr(ib,jb) != 0 ){

	for (UInt nr=0; nr < nrows; nr++){
	  coldata.resize(ncols);
	  pos.resize(ncols);
	  nnz= _Patt->row(ib,jb,nr,coldata.begin(),pos.begin());

	  for (UInt jcol=0; jcol < nnz; jcol++)
	    ans[nc_global+coldata[jcol]-OFFSET]+=_bValues[ib][jb][pos[jcol]]*
	      v[nr_global+nr];
	}
      }
      nc_global += ncols;
    }
    nr_global += nrows;
  }
  return ans;
}

template<UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, MSRPatt, double>::
zeros()
{
  UInt ib,jb;
  for (ib= 0; ib < BRows; ib++)
    for (jb= 0; jb < BCols; jb++)
      if (_Patt->block_ptr(ib,jb) != 0 ){
	typename vector<double>::iterator start= _bValues[ib][jb].begin();
	typename vector<double>::iterator end  = _bValues[ib][jb].end();
	fill(start, end, 0.0);
      }
}

//-----------------------------------------------------------------------
// MixedMatr SPECIALIZATION FOR CSR
//-----------------------------------------------------------------------

//Default Constructor
template<UInt BRows, UInt BCols>
MixedMatr<BRows, BCols, CSRPatt, double>::
MixedMatr()
{
  for (UInt i=0; i<BRows; i++)
    for (UInt j=0; j<BCols; j++)
      _bValues[i][j].resize(0);
}

//Constructor from an existing external pattern.
//version for CSR format: size of values= nnz.
template<UInt BRows, UInt BCols>
MixedMatr<BRows, BCols, CSRPatt, double>::
MixedMatr(const MixedPattern<BRows, BCols, CSRPatt> &ex_pattern)
{
  _Patt = &ex_pattern;

  for (UInt i=0; i<BRows; i++)
    for (UInt j=0; j<BCols; j++)
      _bValues[i][j].resize(ex_pattern.nNz(i,j));
}

// Determines the lumped diagonal of P1 mass matrix.
template<UInt BRows, UInt BCols>
vector<double>
MixedMatr<BRows, BCols, CSRPatt, double>::
MassDiagP1() const
{
  UInt nrows = _Patt->nRows();
  UInt nnz = 0, nr_global = 0;
  Container coldata, position;

  vector<double> diag;
  diag.resize(nrows,0.0);

  for (UInt ib=0; ib< BRows; ib++){

    nrows = _Patt->nRows(ib,0);

    for (UInt jb=0; jb< BCols; jb++){

      if (_Patt->block_ptr(ib,jb) != 0 ){

	for (UInt nr=0; nr < nrows; nr++){
	  coldata.resize(_Patt->nCols(ib,jb));
	  position.resize(_Patt->nCols(ib,jb));
	  nnz= _Patt->row(ib,jb,nr,coldata.begin(),position.begin());
	  for (UInt jcol=0; jcol < nnz; jcol++)
	    diag[nr_global+nr] +=
	      _bValues[ib][jb][position[jcol]];
	}
      }
    }
    nr_global += nrows;
  }

  return diag;
}

// Inverts the diagonal matrix Diag and mulitply it by the matrix Mat.
template<UInt BR, UInt BC>
void
MultInvDiag(const vector<Real> &Diag,
	    const MixedMatr<BR,BC,CSRPatt,Real> &Mat,
	    MixedMatr<BR,BC,CSRPatt,Real> &ans)
{
  ASSERT(find(Diag.begin(),Diag.end(),0) == Diag.end(), "Diagonal matrix Diag must be invertible");

  ASSERT(Diag.size() == Mat._Patt->nRows(), "Matrices size not compatible");

  // Product:
  UInt nrows = 0;
  UInt nnz = 0, nr_global = 0;
  Container coldata, pos;

  for (UInt ib=0; ib< BR; ib++){

    nrows = Mat._Patt->nRows(ib,0);

    for (UInt jb=0; jb< BC; jb++){

      if (Mat._Patt->block_ptr(ib,jb) != 0 ){

	for (UInt nr=0; nr < nrows; nr++){
	  coldata.resize(Mat._Patt->nCols(ib,jb));
	  pos.resize(Mat._Patt->nCols(ib,jb));
	  nnz= Mat._Patt->row(ib,jb,nr,coldata.begin(),pos.begin());

	  for (UInt jcol=0; jcol < nnz; jcol++)
	    ans._bValues[ib][jb][pos[jcol]]=
	      ( Mat._bValues[ib][jb][pos[jcol]]/Diag[nr_global+nr] );
	}
      }
    }
    nr_global += nrows;
  }
}

// Gives the diagonal of a block matrix.
template<UInt BRows, UInt BCols>
vector<double>
MixedMatr<BRows, BCols, CSRPatt, double>::
giveDiag() const
{
  ASSERT_PRE(BRows == BCols, "block matrix must be a square matrix");

  vector<Real> diag;
  diag.resize(_Patt->nRows(),0.0);

  Container coldata, pos;
  UInt nrows= 0, nr_global= 0, nnz= 0;

  for (UInt ib=0; ib < BRows; ib++){
    ASSERT_PRE(_Patt->nRows(ib,ib) == _Patt->nCols(ib,ib),
	       "matrix must have nrows=ncols for diagonal blocks");

    nrows= _Patt->nRows(ib,ib);

    if (_Patt->block_ptr(ib,ib) != 0 )
      for (UInt nr= 0; nr < nrows; nr++){

	coldata.resize(_Patt->nCols(ib,ib));
	pos.resize(_Patt->nCols(ib,ib));

	nnz= _Patt->row(ib,ib,nr,coldata.begin(),pos.begin());

	UInt jcol = 0;

	if (coldata[jcol]-OFFSET == nr)
	  diag[nr_global+nr]= _bValues[ib][ib][pos[jcol]];
	else
	  {
	    while ( coldata[jcol]-OFFSET != nr ) jcol++;
	    diag[nr_global+nr]= _bValues[ib][ib][pos[jcol]];
	  }
      }
    nr_global += nrows;
  }
  return diag;
}


// Warning: the two matrices will point to the same pattern.
template<UInt BRows, UInt BCols>
MixedMatr<BRows, BCols, CSRPatt, double> &
MixedMatr<BRows, BCols, CSRPatt, double>::
operator=(const MixedMatr<BRows,BCols,CSRPatt,double> &RhMtr)
{
  if (&RhMtr != this)
    {
      UInt ib,jb;
      for (ib= 0; ib < BRows; ib++)
	for (jb= 0; jb < BCols; jb++)
	  _bValues[ib][jb] = RhMtr.bValues(ib,jb);

      _Patt = RhMtr.Patt();
    }

  return *this;

}

// Assigns the matrix to loc_val at the place of index where.
template<UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, CSRPatt, double>::
set_mat(UInt ib, UInt jb, UInt where, double loc_val)
{
  _bValues[ib][jb][where-OFFSET] = loc_val;
}

// Assigns the matrix element (row,col) to loc_val.
template<UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, CSRPatt, double>::
set_mat(UInt row, UInt col, double loc_val)
{
  UInt m,n,lr,lc;
  extract_pair(_Patt->locateElBlock(row, col),m,n);
  extract_pair(_Patt->localNumber(m, n, row, col), lr, lc);
  set_mat(m, n, lr, lc, loc_val);
}

// Assigns the matrix element (row,col) of block (ib,jb) to loc_val.
template<UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, CSRPatt, double>::
set_mat(UInt ib, UInt jb, UInt row, UInt col, double loc_val)
{
  if (_Patt->block_ptr(ib,jb) != 0 )
    {
      pair<UInt,bool> where = _Patt->block_ptr(ib,jb)->locate_index(row,col);
      if (where.second) _bValues[ib][jb][where.first] = loc_val;
    }
}

// Adds loc_val to the matrix element (row,col).
template<UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, CSRPatt, double>::
set_mat_inc(UInt row, UInt col, double loc_val)
{
  UInt m,n,lr,lc;
  extract_pair(_Patt->locateElBlock(row, col), m, n);
  extract_pair(_Patt->localNumber(m, n, row, col), lr, lc);
  set_mat_inc(m, n, lr, lc, loc_val);
}

// Adds loc_val to the matrix element (row,col) of block (ib,jb).
template<UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, CSRPatt, double>::
set_mat_inc(UInt ib, UInt jb, UInt row, UInt col, double loc_val)
{
  if (_Patt->block_ptr(ib,jb) != 0 )
    {
      pair<UInt,bool> where = _Patt->block_ptr(ib,jb)->locate_index(row,col);
      if (where.second) _bValues[ib][jb][where.first] += loc_val;
    }
}

// Returns the matrix element (i,j) value.
template<UInt BRows, UInt BCols>
double
MixedMatr<BRows, BCols, CSRPatt, double>::
get_value(UInt i, UInt j)
{
  UInt m,n,lr,lc;
  extract_pair(_Patt->locateElBlock(i,j), m, n);
  extract_pair(_Patt->localNumber(m, n, i, j), lr, lc);
  return get_value(m, n, lr, lc);
}
// const qualifyer version
template<UInt BRows, UInt BCols>
const double
MixedMatr<BRows, BCols, CSRPatt, double>::
get_value(UInt i, UInt j) const
{
  UInt m,n,lr,lc;
  extract_pair(_Patt->locateElBlock(i ,j), m, n);
  extract_pair(_Patt->localNumber(m, n, i, j), lr, lc);
  return get_value(m, n, i, j);
}

// Returns the matrix element (i,j) value of block (ib,jb).
template<UInt BRows, UInt BCols>
double
MixedMatr<BRows, BCols, CSRPatt, double>::
get_value(UInt ib, UInt jb, UInt i, UInt j)
{
  if (_Patt->block_ptr(ib,jb) != 0)
    return _bValues[ib][jb][_Patt->block_ptr(ib,jb)->locate_index(i,j).first];
  else
    return 0.0;
}
// const qualifyer version
template<UInt BRows, UInt BCols>
const double
MixedMatr<BRows, BCols, CSRPatt, double>::
get_value(UInt ib, UInt jb, UInt i, UInt j) const
{
  if (_Patt->block_ptr(ib,jb) != 0)
    return _bValues[ib][jb][_Patt->block_ptr(ib,jb)->locate_index(i,j).first];
  else
    return 0.0;
}

// Matrix visualization a la matlab.
template<UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, CSRPatt, double>::
spy(string  const &filename)
{
  // Purpose: Matlab dumping and spy
  string nome=filename, uti=" , ";
  //
  // check on the file name
  //
  UInt i=filename.find(".");

  if (i<=0) nome=filename+".m";
  else {
    if (i!=filename.size()-2  || filename[i+1]!='m'){
      cerr << "Wrong file name ";
      nome=filename+".m";}
  }

  ofstream file_out(nome.c_str());
  UInt mb, nb, j;
  UInt r_offset, c_offset, ia_start, ia_end;
  file_out << "S = [ ";
  //loop on block rows
  for (mb=0; mb < BRows; mb++){
    //loop on local rows
    for (i=0; i < _Patt->nRows(mb,0); i++)
      //loop on block columns
      for (nb=0; nb < BCols; nb++)
	if (_Patt->block_ptr(mb,nb) != 0){
	  //row and col offsets
	  extract_pair(_Patt->blockOffset(mb,nb),r_offset, c_offset);
	  //for CSR:
	  ia_start= _Patt->block_ptr(mb,nb)->give_ia()[i]-OFFSET;
	  ia_end= _Patt->block_ptr(mb,nb)->give_ia()[i+1]-OFFSET;
	  for( j= ia_start; j < ia_end; j++){
	    file_out << i+r_offset+1 << uti
		     <<	_Patt->block_ptr(mb,nb)->give_ja()[j]+c_offset+1-OFFSET
		     << uti << _bValues[mb][nb][j] << endl;
	  }
	}
  }

  file_out << "];" << endl;
  file_out << "I=S(:,1); J=S(:,2); S=S(:,3); A=sparse(I,J,S); spy(A);"<<endl;
}

// Assigns matrix diagonal element (r,r) to coeff, other elts
// of row r to zero.
template<UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, CSRPatt, double>::
diagonalize_row( UInt const r, double const coeff)
{
  UInt nnz= 0, m=0, n=0;
  Container coldata, pos;
  coldata.resize(_Patt->nCols());
  pos.resize(_Patt->nCols());

  nnz= _Patt->row(r-OFFSET,coldata.begin(),pos.begin());

  UInt jcol=0;
  if (coldata[jcol] == r){
    //diagonal element:
    extract_pair(_Patt->locateElBlock(r-OFFSET ,coldata[jcol]-OFFSET), m, n);
    _bValues[m][n][pos[jcol]]= coeff;
    //other elements:
    for (jcol=1; jcol < nnz; jcol++){
      extract_pair(_Patt->locateElBlock(r-OFFSET ,coldata[jcol]-OFFSET), m, n);
      _bValues[m][n][pos[jcol]]= 0.0;
    }
  }
  else
    for (jcol=0; jcol < nnz; jcol++){
      extract_pair(_Patt->locateElBlock(r-OFFSET ,coldata[jcol]-OFFSET), m, n);
      if (coldata[jcol] == r)
	//diagonal element:
	_bValues[m][n][pos[jcol]]= coeff;
      else
	//other elements
	_bValues[m][n][pos[jcol]]= 0.0;
    }
}

// assign a row to zero. Remark, zero might be defined for any double
template<UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, CSRPatt, double>::
zero_row(UInt const row)
{
  diagonalize_row(row, 0.0);
}

// Assigns matrix diagonal element (r,r) to coeff, other elts
// of row r to zero, and vector b element b(r) to coeff*datum.
template<UInt BRows, UInt BCols>
template< typename VectorType >
void
MixedMatr<BRows, BCols, CSRPatt, double>::
diagonalize( UInt const row, double const coeff, VectorType &b,
	     double datum)
{
  UInt nnz= 0, m=0, n=0;
  Container coldata, pos;
  coldata.resize(_Patt->nCols());
  pos.resize(_Patt->nCols());

  nnz= _Patt->row(row-OFFSET,coldata.begin(),pos.begin());

  UInt jcol=0;
  if (coldata[jcol] == row){ //case diagfirst
    //diagonal element:
    extract_pair(_Patt->locateElBlock(row-OFFSET ,coldata[jcol]-OFFSET), m, n);
    _bValues[m][n][pos[jcol]]= coeff;
    //other elements:
    for (jcol=1; jcol < nnz; jcol++){
      extract_pair(_Patt->locateElBlock(row-OFFSET ,coldata[jcol]-OFFSET), m, n);
      _bValues[m][n][pos[jcol]]= 0.0;
    }
  }
  else
    for (jcol=0; jcol < nnz; jcol++){
      extract_pair(_Patt->locateElBlock(row-OFFSET ,coldata[jcol]-OFFSET), m, n);
      if (coldata[jcol] == row)
	//diagonal element:
	_bValues[m][n][pos[jcol]]= coeff;
      else
	//other elements
	_bValues[m][n][pos[jcol]]= 0.0;
    }

  //rhs:
  b[r-OFFSET]= coeff*datum;
}

// set to zero the row trD and the corresponding column=row for D
// passing the known datum to the rhs b.
template<UInt BR, UInt BC, typename VectorType, typename DataType>
void zero_row_col(UInt const row, MixedMatr<BR,BC,CSRPatt,Real> &trD,
		  MixedMatr<BC,BR,CSRPatt,Real> &D, VectorType &bp,
		  DataType const datum)
{
  // for trD
  trD.zero_row(row);

  // for D
  UInt lr,lc, row_loc;
  UInt m, n;
  UInt r_offset, c_offset;
  UInt col, start_ia, end_ia;

  //block row identification
  extract_pair(trD._Patt->locateElBlock(row,0),m, n);
  //loop on block columns
  for (n= 0; n < BC; n++)
    if (trD._Patt->block_ptr(m,n) != 0){
      // pattern jaT must have been defined !
      ASSERT(trD._Patt->block_ptr(m,n)->give_jaT().size()>0,"jaT must be built before, cf. buildPattTpatt() function");
      extract_pair(trD._Patt->blockOffset(m,n),r_offset, c_offset);
      col= c_offset;
      //local row
      extract_pair(trD._Patt->localNumber(m, n, row, col), lr, lc);
      //loop on the columns of trD involved
      start_ia= trD._Patt->block_ptr(m,n)->give_ia()[lr];
      end_ia= trD._Patt->block_ptr(m,n)->give_ia()[lr+1];
      for ( lc= start_ia; lc < end_ia; lc++)
	{
	  // columns of trD become rows of D
	  row_loc= trD._Patt->block_ptr(m,n)->give_ja()[lc];
	  bp[row_loc+c_offset]-= datum *
	    D._bValues[n][m][trD._Patt->block_ptr(m,n)->give_jaT()[lc]];
	  D._bValues[n][m][trD._Patt->block_ptr(m,n)->give_jaT()[lc]]= 0.0;
	}
    }
}

// Matrix-vector product.
template<UInt BRows, UInt BCols>
vector<double>
MixedMatr<BRows, BCols, CSRPatt, double>::
operator*(const vector<double> &v) const
{
  ASSERT(_Patt->nCols()==v.size(),"Error in Matrix Vector product");
  vector<double> ans;
  ans.resize(_Patt->nRows(),0.0);

  UInt nrows, iastart,iaend;
  UInt r_offset, c_offset;
  UInt ib,jb,ir,ii;

  for (ib=0; ib< BRows; ib++){

    nrows = _Patt->nRows(ib,0);

    for (jb=0; jb< BCols; jb++){

      if (_Patt->block_ptr(ib,jb) != 0 ){

	extract_pair(_Patt->blockOffset(ib,jb),r_offset, c_offset);

	for (ir= 0; ir< nrows; ir++){
	  iastart= _Patt->block_ptr(ib,jb)->give_ia()[ir]-OFFSET;
	  iaend  = _Patt->block_ptr(ib,jb)->give_ia()[ir+1]-OFFSET;
	  for (ii= iastart;ii< iaend;++ii){
	    ans[ir+r_offset]+= _bValues[ib][jb][ii]*
	      v[c_offset+_Patt->block_ptr(ib,jb)->give_ja()[ii]-OFFSET];
	  }
	}
      }
    }
  }
  return ans;
}
// version for type Vector
template<UInt BRows, UInt BCols>
Vector
MixedMatr<BRows, BCols, CSRPatt, double>::
operator*(const Vector &v) const
{
  ASSERT(_Patt->nCols()==v.size(),"Error in Matrix Vector product");
  Vector ans(_Patt->nRows());
  ans=0.0;

  UInt nrows, iastart,iaend;
  UInt r_offset, c_offset;
  UInt ib,jb,ir,ii;

  for (ib=0; ib< BRows; ib++){

    nrows = _Patt->nRows(ib,0);

    for (jb=0; jb< BCols; jb++){

      if (_Patt->block_ptr(ib,jb) != 0 ){

	extract_pair(_Patt->blockOffset(ib,jb),r_offset, c_offset);

	for (ir= 0; ir< nrows; ir++){
	  iastart= _Patt->block_ptr(ib,jb)->give_ia()[ir]-OFFSET;
	  iaend  = _Patt->block_ptr(ib,jb)->give_ia()[ir+1]-OFFSET;
	  for (ii= iastart;ii< iaend;++ii){
	    ans[ir+r_offset]+= _bValues[ib][jb][ii]*
	      v[c_offset+_Patt->block_ptr(ib,jb)->give_ja()[ii]-OFFSET];
	  }
	}
      }
    }
  }
  return ans;
}
// version C pointer vector. BEWARE: no check on bounds is done !
template<UInt BRows, UInt BCols>
void operMatVec(double * const mv,
		const MixedMatr<BRows, BCols, CSRPatt, double> &Mat,
		const double *v)
{
  UInt nrows, iastart,iaend;
  UInt r_offset, c_offset;
  UInt ib,jb,ir,ii;

  for (ib=0; ib< BRows; ib++){

    nrows = Mat._Patt->nRows(ib,0);

    for (ir= 0; ir< nrows; ir++){

      //initialize
      mv[ir+ib*nrows]= 0.;

      for (jb=0; jb< BCols; jb++){

	if (Mat._Patt->block_ptr(ib,jb) != 0 ){

	  extract_pair(Mat._Patt->blockOffset(ib,jb),r_offset, c_offset);

	  iastart= Mat._Patt->block_ptr(ib,jb)->give_ia()[ir]-OFFSET;
	  iaend  = Mat._Patt->block_ptr(ib,jb)->give_ia()[ir+1]-OFFSET;

	  for (ii= iastart;ii< iaend;++ii){
	    mv[ir+r_offset]+= Mat._bValues[ib][jb][ii]*
	      v[c_offset+Mat._Patt->block_ptr(ib,jb)->give_ja()[ii]-OFFSET];
	  }
	}
      }
    }
  }
}

//necessary for IML++ library: transpose-matrix by vector product.
template<UInt BRows, UInt BCols>
Vector
MixedMatr<BRows, BCols, CSRPatt, double>::
trans_mult(const Vector &v) const
{
  ASSERT(_Patt->nRows()==v.size(),"Error in Matrix Vector product");
  Vector ans(_Patt->nRows());
  ans=0.0;

  UInt nrows = 0, ncols = 0;
  UInt nnz = 0, nr_global = 0, nc_global;
  Container coldata, pos;

  for (UInt ib=0; ib< BRows; ib++){

    nrows = _Patt->nRows(ib,0);
    nc_global=0;

    for (UInt jb=0; jb< BCols; jb++){

      ncols = _Patt->nCols(ib,jb);

      if (_Patt->block_ptr(ib,jb) != 0 ){

	for (UInt nr=0; nr < nrows; nr++){
	  coldata.resize(ncols);
	  pos.resize(ncols);
	  nnz= _Patt->row(ib,jb,nr,coldata.begin(),pos.begin());

	  for (UInt jcol=0; jcol < nnz; jcol++)
	    ans[nc_global+coldata[jcol]-OFFSET]+=_bValues[ib][jb][pos[jcol]]*
	      v[nr_global+nr];
	}
      }
      nc_global += ncols;
    }
    nr_global += nrows;
  }
  return ans;
}

template<UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, CSRPatt, double>::
zeros()
{
  UInt ib,jb;
  for (ib= 0; ib < BRows; ib++)
    for (jb= 0; jb < BCols; jb++)
      if (_Patt->block_ptr(ib,jb) != 0 ){
	typename vector<double>::iterator start= _bValues[ib][jb].begin();
	typename vector<double>::iterator end  = _bValues[ib][jb].end();
	fill(start, end, 0.0);
      }
}

//-----------------------------------------------------------------------
// DiagPreconditioner
//-----------------------------------------------------------------------

//for CSR or MSR normal pattern
DiagPreconditioner<Vector>::
DiagPreconditioner(const CSRMatr<CSRPatt,double> &M);

DiagPreconditioner<Vector>::
DiagPreconditioner(const MSRMatr<double> &M);



//for VBR pattern
DiagPreconditioner<Vector>::
DiagPreconditioner(const VBRMatr<double> &M);


//for CSR block pattern
DiagPreconditioner<VectorBlock>::
DiagPreconditioner(const CSRMatr<CSRPatt,Tab2d> &M);


//solve the diagonal system
Vector
DiagPreconditioner<Vector>::
solve(const Vector &x) const;



VectorBlock
DiagPreconditioner<VectorBlock>::
solve(const VectorBlock &x) const;


//-----------------------------------------------------------------------
// DiagPreconditioner
//-----------------------------------------------------------------------

//for CSR or MSR normal pattern
IDPreconditioner<Vector>::
IDPreconditioner(const CSRMatr<CSRPatt,double> &M);


IDPreconditioner<Vector>::
IDPreconditioner(const MSRMatr<double> &M);


//for VBR pattern
IDPreconditioner<Vector>::
IDPreconditioner(const VBRMatr<double> &M);



//for CSR block pattern
IDPreconditioner<VectorBlock>::
IDPreconditioner(const CSRMatr<CSRPatt,Tab2d> &M);


//solve the diagonal system
Vector
IDPreconditioner<Vector>::
solve(const Vector &x) const;


VectorBlock
IDPreconditioner<VectorBlock>::
solve(const VectorBlock &x) const;
}
#endif
