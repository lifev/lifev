/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2004-10-26

  Copyright (C) 2004 EPFL

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
/**
   \file CSRMatrix.hpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-10-26
 */

namespace LifeV
{
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

template <typename PatternType, typename DataType>
class CSRMatr
{
public:
    //! default constructor. Caution: does not have a pattern.
    CSRMatr();

    // All other constructors MUST be based on an existing pattern

    //! constructor from pattern
    CSRMatr( const PatternType &ex_pattern );

    /*! constructor with given value array
     *  @param ex_pattern the pattern
     *  @param ex_value the value array
     */
    CSRMatr( const PatternType &ex_pattern,
             const std::vector<DataType> &ex_value );

    //! copy constructor
    CSRMatr( const CSRMatr<PatternType, DataType> &RightHandCSR );

    /*! conversion copy constructor.
     *  requires an already converted pattern
     *  @param ex_pattern the converted pattern
     *  @param msrMatr the MSRMatrix to be copied
     */
    CSRMatr( const PatternType& ex_pattern, const MSRMatr<DataType>& msrMatr );

    //! assignment operator for CSR <- MSR.
    //! the two matrices need to have the same dimension and sparsity pattern
    CSRMatr& operator=( const MSRMatr<DataType>& msrMatr );

    //! returns a const * to the pattern
    const PatternType * Patt() const
    {
        return _Patt;
    };

    //! returns the value vector
    std::vector<DataType> & value()
    {
        return _value;
    };

    //! give the value vector in raw format (suitable for C)
    DataType * giveRawCSR_value()
    {
        return & ( _value.front() );
    }

    //! give the value vector in const raw format (suitable for C)
    const DataType * giveRawCSR_value() const
    {
        return & ( _value.front() );
    }

    //! assignment operator.
    //! Warning: the two matrices will point to the same pattern
    CSRMatr& operator= ( const CSRMatr<PatternType, DataType> &RhCsr );

    /*! set single value by position in value vector
     *  @param where position in value vector
     *  @param loc_val new value
     */
    void set_mat( UInt where, DataType loc_val );

    /*! set single value by (row,col)-address
     *  @param row row number
     *  @param col column number
     *  @param loc_val new value
     */
    void set_mat( UInt row, UInt col, DataType loc_val );

    //! add \c loc_val to entry (row,col), useful for assembly
    void set_mat_inc( UInt row, UInt col, DataType loc_val );

    //! return matrix entry (i,j)
    DataType& get_value( UInt i, UInt j ) const
    {
        return _value[ _Patt->locate_index( i, j ).first ];
    };

    //! print matrix to standard output
    void ShowMe();

    //! write matrix in sparse matlab format and spy
    /*! just run the resulting m-file and the matrix is loaded into A
     *  and its sparsity pattern is displayed.
     *  @param filename name of the m-file
     */
    void spy( std::string const &filename );

private:
    std::vector<DataType> _value;
    const PatternType *_Patt; // I want to link the values to a pattern, NOT to change the pattern itself  (which is const)
    //     static const DataType _DefaultValue = 0;
};

////////////////////////////////////////////////////////////////
//
// CSR Format SPECIALIZATION FOR Usual CSR
//
///////////////////////////////////////////////////////////////
//

template <typename DataType>
class CSRMatr<CSRPatt, DataType>
{
public:
    CSRMatr(); //!< default constructor : NULL pattern
    //
    // Note that the constructors MUST be based on an existing pattern
    //
    CSRMatr( const CSRPatt &ex_pattern );
    //! Version for DataType=Tab2d
    CSRMatr( const CSRPatt &ex_pattern, UInt const nr, UInt const nc );
    CSRMatr( const CSRPatt &ex_pattern, const std::vector<DataType> &ex_value );
    CSRMatr( const CSRMatr<CSRPatt, DataType> &RightHandCSR );
    CSRMatr( const CSRPatt &ex_pattern, const MSRMatr<DataType> &msrMatr );
    const CSRPatt * Patt() const
    {
        return _Patt;
    };
    const std::vector<DataType> & value() const
    {
        return _value;
    };
    std::vector<DataType> & value()
    {
        return _value;
    };
    DataType * giveRawCSR_value()
    {
        return & ( _value.front() );
    }
    const DataType * giveRawCSR_value() const
    {
        return & ( _value.front() );
    }

    CSRMatr& operator= ( const CSRMatr<CSRPatt, DataType> &RhCsr );
    // Warning: the two matrices will point to the same pattern

    CSRMatr& operator= ( const MSRMatr<DataType> &msrMatr );


    //! Matrix-vector product
    Vector operator*( const Vector &v ) const;
    //! Version for block matrices
    VectorBlock operator*( const VectorBlock &v ) const;
    //! Version for C pointer vector. BEWARE: no check on bounds done !
    template <typename DataT>
    friend void operMatVec( DataT * const mv,
                            const CSRMatr<CSRPatt, DataT> &Mat,
                            const DataT *v );

    //! necessary for IML++ library: transpose-matrix by vector product
    Vector trans_mult( const Vector &v ) const;
    //! version for block matrices
    VectorBlock trans_mult( const VectorBlock &v );

    //! assign a row to zero. Remark, zero might be defined for any DataType
    void zero_row( UInt const row );

    //! diagonalize a row r and rhs for Dirichlet BC treatment
    void diagonalize( UInt const r, DataType const coeff, Vector &b,
                      DataType datum );
    //! diagonalize the row r of the matrix only (cf diagonalize for taking into
    //! account the rhs of the linear system as well)
    void diagonalize_row( UInt const r, DataType const coeff );

    //! determine the lumped matrix for P1
    std::vector<DataType> MassDiagP1() const;
    //! build the inverse of a diagonal matrix and multiply it by another matrix
    friend void MultInvDiag( const std::vector<Real> &Diag,
                             const CSRMatr<CSRPatt, Real> &Mat,
                             CSRMatr<CSRPatt, Real> &ans );

    void set_mat( UInt row, UInt col, DataType loc_val );
    void set_mat_inc( UInt row, UInt col, DataType loc_val );
    DataType& get_value( UInt i, UInt j )
    {
        return _value[ _Patt->locate_index( i, j ).first ];
    };
    const DataType& get_value( UInt i, UInt j ) const
    {
        return _value[ _Patt->locate_index( i, j ).first ];
    };

    void ShowMe();

    //! write matrix in sparse matlab format and spy
    /*! just run the resulting m-file and the matrix is loaded into A
     *  and its sparsity pattern is displayed.
     *  @param filename name of the m-file
     */
    void spy( std::string const &filename );

    //!column-concatenation of two blocks of CSRMatr
    /*  friend CSRMatr<CSRPatt,double>
        colUnify(const CSRMatr<CSRPatt,double> &Mat1,
        const CSRMatr<CSRPatt,double> &Mat2);*/
    //version without using static : the one to keep (13/12/01)
    friend void
    colUnify( CSRMatr<CSRPatt, double> &ans, const CSRMatr<CSRPatt, double> &Mat1,
              const CSRMatr<CSRPatt, double> &Mat2 );

    //!row-concatenation of two blocks of CSRMatr (without using static)
    friend void
    rowUnify( CSRMatr<CSRPatt, double> &ans, const CSRMatr<CSRPatt, double> &Mat1,
              const CSRMatr<CSRPatt, double> &Mat2 );

    //!set to zero the row trD and the corresponding column=row for D
    template <typename VectorType>
    friend void
    zero_row_col( UInt const row, CSRMatr<CSRPatt, double> &trD,
                  CSRMatr<CSRPatt, double> &D, VectorType &bp,
                  DataType const datum );

    //! set to zero the matrix;
    void zeros();

private:
    std::vector<DataType> _value;
    const CSRPatt *_Patt; // I want to link the values to a pattern, NOT to change the pattern itself  (which is const)
    //  static const DataType _DefaultValue = 0;
};

//-------------------------------------------------------------------------------------------------------
// CSR - VALUES
//------------------------------------------------------------------------------------------------------
template <typename PatternType, typename DataType>
CSRMatr<PatternType, DataType>::CSRMatr()
        :
        _Patt( 0 )
{
    // nothing to do here
}

template <typename PatternType, typename DataType>
CSRMatr<PatternType, DataType>::
CSRMatr( const PatternType &ex_pattern )
{
    _Patt = &ex_pattern;
    _value.reserve( ex_pattern.nNz() ); // no initialization of values !
}


template <typename PatternType, typename DataType>
CSRMatr<PatternType, DataType>::
CSRMatr( const PatternType &ex_pattern, const std::vector<DataType> &ex_value )
{
    _value.reserve( ex_pattern.nNz() ); // in case of block matrix, there is
    // no default constructor for class KNM
    _Patt = &ex_pattern;
    ASSERT( _Patt->nNz() == ex_value.size(),
            "Error in CSR Matrix Values LifeV" );
    // Warning: if PatternType = CSRPattSymm => _ja.size() != _nnz.===> remember: _ja.size() = (_nnz+_nrows)/2
    _value = ex_value;
}

template <typename PatternType, typename DataType>
CSRMatr<PatternType, DataType>::
CSRMatr( const CSRMatr<PatternType, DataType> &RightHandCSR ) :
        _value( RightHandCSR.value() ), _Patt( RightHandCSR.Patt() )
{
    // nothing to do here
}

template <typename PatternType, typename DataType>
CSRMatr<PatternType, DataType>::
CSRMatr( const PatternType &ex_pattern, const MSRMatr<DataType>& msrMatr )
{
    _Patt = &ex_pattern;
    _value.reserve( ex_pattern.nNz() );
    *this = msrMatr;
}

template <typename PatternType, typename DataType>
CSRMatr<PatternType, DataType>&
CSRMatr<PatternType, DataType>::operator= ( const CSRMatr<PatternType, DataType> &RhCsr )
{
    if ( &RhCsr != this )
    {
        _Patt = RhCsr.Patt();
        _value = RhCsr.value();
    }
    return *this;
}

template <typename PatternType, typename DataType>
CSRMatr<PatternType, DataType>&
CSRMatr<PatternType, DataType>::operator= ( const MSRMatr<DataType>& msrMatr )
{
    typename std::vector<DataType>::iterator value = _value.begin();
    const Container& ja = _Patt->ja();
    Container::const_iterator ia = _Patt->ia().begin();
    UInt nrows = _Patt->nRows();
    for ( UInt iRow = 0; iRow < nrows; ++iRow, ++ia )
    {
        for ( UInt i = *ia - OFFSET; i < *( ia + 1 ) - OFFSET; ++i, ++value )
        {
            UInt iCol = ja[ i ] - OFFSET;
            *value = msrMatr.get_value( iRow, iCol );
        }
    }
    return *this;
}

template <typename PatternType, typename DataType>
void
CSRMatr<PatternType, DataType>::
set_mat( UInt where, DataType loc_val )
{
    _value[ where ] = loc_val;
    return ;
}

template <typename PatternType, typename DataType>
void
CSRMatr<PatternType, DataType>::
set_mat( UInt row, UInt col, DataType loc_val )
{
    std::pair<UInt, bool> where = *_Patt.locate_index( row, col );
    if ( where.second )
        _value[ where.first ] = loc_val;
    return ;
}

template <typename PatternType, typename DataType>
void
CSRMatr<PatternType, DataType>::
set_mat_inc( UInt row, UInt col, DataType loc_val )
{

    std::pair<UInt, bool> where = _Patt->locate_index( row, col );
    if ( where.second )
        _value[ where.first ] += loc_val;

    return ;
}


//-------------------------------------------------------------------------------------------------------
// CSR - VALUES WITH CSR Pattern (Usual)
//------------------------------------------------------------------------------------------------------
template <typename DataType>
CSRMatr<CSRPatt, DataType>::CSRMatr()
        :
        _Patt( 0 )
{
    // nothing to do here
}


template <typename DataType>
CSRMatr<CSRPatt, DataType>::
CSRMatr( const CSRPatt &ex_pattern )
{
    _Patt = &ex_pattern;
    _value.resize( ex_pattern.nNz() );
}
//version for Datatype=Tab2d
template <>
CSRMatr<CSRPatt, Tab2d>::CSRMatr( const CSRPatt &ex_pattern,
                                  UInt const nr, UInt const nc );

template <typename DataType>
CSRMatr<CSRPatt, DataType>::CSRMatr( const CSRPatt &ex_pattern,
                                     const std::vector<DataType> &ex_value )
{
    _value.reserve( ex_pattern.nNz() ); // in case of block matrix, there is
    // no default constructor for class KNM
    _Patt = &ex_pattern;
    ASSERT( _Patt->nNz() == ex_value.size(),
            "Error in CSR Matrix Values LifeV" );
    // Warning: if PatternType = CSRPattSymm => _ja.size() != _nnz.===> remember: _ja.size() = (_nnz+_nrows)/2
    _value = ex_value;
}

template <typename DataType>
CSRMatr<CSRPatt, DataType>::
CSRMatr( const CSRMatr<CSRPatt, DataType> &RightHandCSR ) :
        _value( RightHandCSR.value() ), _Patt( RightHandCSR.Patt() )
{
    // nothing to do here
}

template <typename DataType>
CSRMatr<CSRPatt, DataType>::
CSRMatr( const CSRPatt &ex_pattern, const MSRMatr<DataType>& msrMatr )
{
    _Patt = &ex_pattern;
    _value.resize( ex_pattern.nNz() );
    *this = msrMatr;
}

template <typename DataType>
CSRMatr<CSRPatt, DataType>&
CSRMatr<CSRPatt, DataType>::operator= ( const CSRMatr<CSRPatt, DataType> &RhCsr )
{
    if ( &RhCsr != this )
    {
        _Patt = RhCsr.Patt();
        _value = RhCsr.value();
    }
    return *this;
}


template <typename DataType>
CSRMatr<CSRPatt, DataType>&
CSRMatr<CSRPatt, DataType>::operator= ( const MSRMatr<DataType>& msrMatr )
{
    typename std::vector<DataType>::iterator value = _value.begin();
    const Container& ja = _Patt->ja();
    Container::const_iterator ia = _Patt->ia().begin();
    Container::const_iterator endia = _Patt->ia().end();
    UInt nrows = _Patt->nRows();
    for ( UInt iRow = 0; ia != endia && iRow < nrows; ++iRow, ++ia )
    {
        Container::const_iterator nextia = boost::next( ia );
        for ( int i = *ia;
              nextia != endia && i < *nextia;
              ++i, ++value )
        {
            int iCol = ja[ i ];
            *value = msrMatr.get_value( iRow, iCol );
        }
    }
    return *this;
}

template <typename DataType>
void
CSRMatr<CSRPatt, DataType>::
set_mat( UInt row, UInt col, DataType loc_val )
{
    std::pair<UInt, bool> where = _Patt->locate_index( row, col );
    if ( where.second )
        _value[ where.first ] = loc_val;
    return ;
}

template <typename DataType>
void
CSRMatr<CSRPatt, DataType>::
set_mat_inc( UInt row, UInt col, DataType loc_val )
{
    std::pair<UInt, bool> where = _Patt->locate_index( row, col );
    if ( where.second )
        _value[ where.first ] += loc_val;

    return ;
}

// determine the lumped matrix for P1
template <typename DataType>
std::vector<DataType>
CSRMatr<CSRPatt, DataType>::MassDiagP1() const
{
    UInt nrows = _Patt->nRows();
    UInt ncols = _Patt->nCols();
    ASSERT( ncols == nrows, "The matrix must be square to be lumped" );
    std::vector<DataType> diag( nrows );
    for ( UInt nrow = 0; nrow < nrows; ++nrow )
    {
        for ( UInt ii = _Patt->ia() [ nrow ]; ii < _Patt->ia() [ nrow + 1 ]; ++ii )
            diag[ nrow ] += _value[ ii ];
    }
    return diag;
}


// build the inverse of a diagonal matrix and multiply it by another matrix
void MultInvDiag( const std::vector<Real> &Diag,
                  const CSRMatr<CSRPatt, Real> &Mat, CSRMatr<CSRPatt, Real> &ans ) ;


// transpose-Matrix by vector product
template <typename DataType>
Vector
CSRMatr<CSRPatt, DataType>::
trans_mult( const Vector &v ) const
{
    UInt nrows = _Patt->nRows(); // for square matrices...
    ASSERT( nrows == v.size(), "Error in Matrix Vector product" );
    Vector ans( nrows );
    ans = 0.0;

    for ( UInt ir = 0 + OFFSET;ir < nrows + OFFSET;++ir )
    {
        for ( UInt ii = _Patt->ia() [ ir ] - OFFSET;ii < _Patt->ia() [ ir + 1 ] - OFFSET;++ii )
            ans( _Patt->ja() [ ii ] - OFFSET ) += _value[ ii ] * v( ir );
    }
    return ans;
}

// version for block matrices
template <>
VectorBlock
CSRMatr<CSRPatt, Tab2d>::trans_mult( const VectorBlock &v );

//
template <typename DataType>
void CSRMatr<CSRPatt, DataType>::
zero_row( UInt const row )
{
    typename std::vector<DataType>::iterator start = _value.begin() +
            *( _Patt->give_ia().begin() + row - OFFSET );
    typename std::vector<DataType>::iterator end = _value.begin() +
            *( _Patt->give_ia().begin() + row + 1 - OFFSET );

    //nihil is the same used for diagonalize
    //method in MSRMatr class.
    transform( start, end, start, nihil );
}


// Matrix-vector product
template <typename DataType>
Vector CSRMatr<CSRPatt, DataType>::
operator*( const Vector &v ) const
{
    UInt nrows = _Patt->nRows();
    UInt ncols = _Patt->nCols();

    ASSERT( ncols == v.size(), "Error in Matrix Vector product" );
    Vector ans( nrows );
    ans = ZeroVector( nrows );

    for ( UInt ir = 0 + OFFSET;ir < nrows + OFFSET;++ir )
    {
        for ( UInt ii = _Patt->give_ia() [ ir ] - OFFSET;ii < _Patt->give_ia() [ ir + 1 ] - OFFSET;++ii )
            ans( ir ) += _value[ ii ] * v( _Patt->give_ja() [ ii ] - OFFSET );
    }
    return ans;
}

// version for block matrices
template <>
VectorBlock
CSRMatr<CSRPatt, Tab2d>::operator*( const VectorBlock &v ) const;

// Version for C pointer vector. BEWARE: no check on bounds done !
template <typename DataType>
void operMatVec( DataType * const mv,
                 const CSRMatr<CSRPatt, DataType> &Mat,
                 const DataType *v )
{
    UInt nrows = Mat._Patt->nRows();
    UInt ncols = Mat._Patt->nCols();

    for ( UInt ir = 0 + OFFSET;ir < nrows + OFFSET;++ir )
    {
        mv[ ir ] = 0.;
        for ( UInt ii = Mat._Patt->give_ia() [ ir ] - OFFSET;ii < Mat._Patt->give_ia() [ ir + 1 ] - OFFSET;++ii )
            mv[ ir ] += Mat._value[ ii ] * v[ Mat._Patt->give_ja() [ ii ] - OFFSET ];
    }
}

/*! Diagonalization of row r of the system. Done by setting A(r,r) = coeff,
 *  A(r,j) = 0 and A(j,r) = 0 for j!=r, and suitably correcting the right hand
 *  side of the system.
 *  @param r row to diagonalize
 *  @param coeff value to set the diagonal entry A(r,r) to
 *  @param b right hand side vector to be corrected
 *  @param datum value to set the fix the solution entry x(r) at
 */
template <typename DataType>
void CSRMatr<CSRPatt, DataType>::diagonalize( UInt const r, DataType const coeff,
        Vector &b, DataType datum )
{

    for ( UInt j = 0; j < _Patt->nCols(); ++j )
    {
        set_mat( r, j, 0. ); // A(r,j) = 0
    }
    for ( UInt i = 0; i < _Patt->nRows(); ++i )
    {
        std::pair<UInt, bool> where = _Patt->locate_index( i, r );
        if ( where.second )
        {
            b[ i - OFFSET ] -= _value[ where.first ] * datum; // correct rhs
            _value[ where.first ] = 0.0; // A(j,r) = 0
        }
    }

    set_mat( r, r, coeff ); // A(r,r) = coeff

    b[ r - OFFSET ] = coeff * datum; // correct right hand side for row r
}

/*! Diagonalization of row r of the system. Done by setting A(r,r) = coeff
 *  and A(r,j) = 0 for j!=r without correcting the right hand side
 *  @param r row to diagonalize
 *  @param coeff value to set the diagonal entry A(r,r) to
 */
template <typename DataType>
void CSRMatr<CSRPatt, DataType>::diagonalize_row( UInt const r,
        DataType const coeff )
{
    for ( UInt j = 0; j < _Patt->nCols(); ++j )
    {
        set_mat( r, j, 0. ); // A(r,j) = 0
    }
    set_mat( r, r, coeff ); // A(r,r) = coeff
}


// Correction of ShowMe, Alain, 05/02/02.
template <typename DataType>
void
CSRMatr<CSRPatt, DataType>::ShowMe()
{
    UInt i_first, nrows = _Patt->nRows(), ncols = _Patt->nCols(), nnz = _Patt->nNz();
    Container ja = _Patt->ja();

    std::string pare = "[";
    std::cout << "**************************" << std::endl;
    std::cout << "     CSR Matrix           " << std::endl;
    std::cout << std::endl;
    std::cout << pare;
    for ( UInt i_index = 0; i_index < nrows; ++i_index )
    {
        std::cout << pare;
        pare = " [";
        i_first = _Patt->give_ia() [ i_index ] - OFFSET;

        UInt jj = 0;
        for ( UInt j = 0;j < ncols;j++ )
        {
            if ( j == i_index )
            {
                std::cout << " " << _value[ i_first + jj ] << " ";
                jj++;
            }
            else
            {
                if ( j == ja[ i_first + jj ] - OFFSET )
                {
                    std::cout << " " << _value[ i_first + jj ] << " ";
                    jj++;
                }
                else
                    std::cout << " 0 ";
            }
        }
        if ( i_index == nrows - 1 )
            std::cout << " ]] " << std::endl;
        else
            std::cout << " ]  " << std::endl;
    }
    std::cout << "nnz = " << nnz << ", nrow = " << nrows << ", ncol = " << ncols << std::endl;
    return ;
}

template <typename DataType>
void CSRMatr<CSRPatt, DataType>::
spy( std::string const &filename )
{
    // Purpose: Matlab dumping and spy
    std::string nome = filename, uti = " , ";
    UInt nrows = _Patt->nRows();
    Container ia = _Patt->ia(), ja = _Patt->ja();
    //
    // check on the file name
    //
    int i = filename.find( "." );

    if ( i <= 0 )
        nome = filename + ".m";
    else
    {
        if ( ( unsigned int ) i != filename.size() - 2 || filename[ i + 1 ] != 'm' )
        {
            std::cerr << "Wrong file name ";
            nome = filename + ".m";
        }
    }

    std::ofstream file_out( nome.c_str() );
    ASSERT( file_out, "Error: Output Matrix (Values) file cannot be open" );


    file_out << "S = [ ";
    for ( UInt i = 0;i < nrows;++i )
    {
        for ( UInt ii = ia[ i ] - OFFSET;ii < ia[ i + 1 ] - OFFSET;++ii )
            file_out << i + 1 << uti << ja[ ii ] + 1 - OFFSET << uti << _value[ ii ] << std::endl; /* */
    }
    file_out << "];" << std::endl;

    file_out << "I=S(:,1); J=S(:,2); S=S(:,3); A=sparse(I,J,S); spy(A);" << std::endl;
}

// the case of block matrices with Tab2d block type.
template <>
void
CSRMatr<CSRPatt, Tab2d>::spy( std::string const &filename );

//version without using static (I think it is better)
// Modified by A. Gilardi. 03/02.
void colUnify( CSRMatr<CSRPatt, double> &ans, const CSRMatr<CSRPatt, double> &Mat1,
               const CSRMatr<CSRPatt, double> &Mat2 );
// row-concatenation of two blocks of CSRMatr
// Modified by A. Gilardi. 03/02.
void rowUnify( CSRMatr<CSRPatt, double> &ans, const CSRMatr<CSRPatt, double> &Mat1,
               const CSRMatr<CSRPatt, double> &Mat2 );

// set to zero the row trD and the corresponding column=row for D
// passing the known datum to the rhs b.
template <typename VectorType, typename DataType>
void zero_row_col( UInt const row, CSRMatr<CSRPatt, double> &trD,
                   CSRMatr<CSRPatt, double> &D, VectorType &bp,
                   DataType const datum )
{
    // for trD
    trD.zero_row( row );

    // for D
    UInt j;
    //loop on the columns of trD involved
    for ( j = trD._Patt->give_ia() [ row - OFFSET ];
            j < trD._Patt->give_ia() [ row + 1 - OFFSET ]; j++ )
    {
        // columns of trD become rows of D
        UInt row_loc = trD._Patt->give_ja() [ j ];
        bp[ row_loc ] -= datum * D._value[ trD._Patt->give_jaT() [ j ] ];
        D._value[ trD._Patt->give_jaT() [ j ] ] = 0.;
    }
}

template <typename DataType>
void CSRMatr<CSRPatt, DataType>::
zeros()
{
    typename std::vector<DataType>::iterator start = _value.begin();
    typename std::vector<DataType>::iterator end = _value.end();
    fill( start, end, 0.0 );
}

}
