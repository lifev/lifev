/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2004-10-27

  Copyright (C) 2004 EPFL, INRIA, Politecnico di Milano

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
   \file MixedMatrix.hpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-10-27
 */

namespace LifeV
{
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
template <UInt BRows, UInt BCols, typename PatternType, typename DataType>
class MixedMatr
{
public:
    //! Default constructor.
    MixedMatr();
    //! Constructor from an existing MixedPattern pattern.
    //! Version for MSR format.
    MixedMatr( const MixedPattern<BRows, BCols, MSRPatt> &ex_pattern );
    //! Version for CSR format.
    MixedMatr( const MixedPattern<BRows, BCols, CSRPatt> &ex_pattern );
    //! gives the pattern pointer.
    const MixedPattern<BRows, BCols, PatternType> * Patt() const
    {
        return _Patt;
    }
    //! points to the values vector array.
    std::vector<DataType> * bValues() const
    {
        return _bValues;
    }
    //! The container of block (i,j).
    //! BEWARE: no check is done on the array size.
    const std::vector<DataType> & bValues( UInt i, UInt j ) const
    {
        return _bValues[ i ][ j ];
    }
    std::vector<DataType> & bValues( UInt i, UInt j )
    {
        return _bValues[ i ][ j ];
    }
    //! give the block (i,j) value vector in Raw format (suitable for C)
    DataType * giveRaw_value( UInt i, UInt j )
    {
        return & ( _bValues[ i ][ j ].front() );
    }
    //! give the block (i,j) value vector in Raw format (suitable for C)
    DataType const * giveRaw_value( UInt i, UInt j ) const
    {
        return & ( _bValues[ i ][ j ].front() );
    }

    //! Determines the lumped diagonal of P1 mass matrix.
    std::vector<DataType> MassDiagP1() const;
    //! Inverts the diagonal matrix Diag and mulitply it by the matrix Mat.
    template <UInt BR, UInt BC, typename PattType>
    friend void MultInvDiag( const std::vector<Real> &Diag,
                             const MixedMatr<BR, BC, PattType, Real> &Mat,
                             MixedMatr<BR, BC, PattType, Real> &ans );
    //! Gives the diagonal of a block matrix.
    std::vector<DataType> giveDiag() const;

    //! Warning: the two matrices will point to the same pattern.
    MixedMatr & operator= ( const MixedMatr<BRows, BCols, PatternType, DataType>
                            &RhMtr );
    //! Assigns the matrix to loc_val at the place of index where in the
    //! block (ib,jb).
    void set_mat( UInt ib, UInt jb, UInt where, DataType loc_val );
    //! Assigns the matrix element (row,col) to loc_val.
    void set_mat( UInt row, UInt col, DataType loc_val );
    //! Assigns the matrix element (row,col) of block (ib,jb) to loc_val.
    void set_mat( UInt ib, UInt jb, UInt row, UInt col, DataType loc_val );
    //! Adds loc_val to the matrix element (row,col).
    void set_mat_inc( UInt row, UInt col, DataType loc_val );
    //! Adds loc_val to the matrix element (row,col) of block (ib,jb).
    void set_mat_inc( UInt ib, UInt jb, UInt row, UInt col, DataType loc_val );
    //! Returns the matrix element (i,j) value.
    DataType& get_value( UInt i, UInt j );
    const DataType& get_value( UInt i, UInt j ) const;
    //! Returns the matrix element (i,j) value of block (ib,jb).
    DataType& get_value( UInt ib, UInt jb, UInt i, UInt j );
    const DataType& get_value( UInt ib, UInt jb, UInt i, UInt j ) const;
    //! Shows the matrix (only the pattern here).
    void ShowMe()
    {
        return _Patt->showMe( true );
    }
    //! write matrix in sparse matlab format and spy
    /*! just run the resulting m-file and the matrix is loaded into A
     *  and its sparsity pattern is displayed.
     *  @param filename name of the m-file
     */
    void spy( std::string const &filename );
    //! Assigns matrix diagonal element (r,r) to coeff, other elts
    //! of row r to zero.
    void diagonalize_row( UInt const r, DataType const coeff );
    //! assign a row to zero. Remark, zero might be defined for any DataType
    void zero_row( UInt const row );

    //! set entries (rVec(i),rVec(i)) to coeff and rest of row r(i) to zero
    template < typename VectorType >
    void diagonalize( std::vector<UInt> const rVec, DataType const coeff, VectorType &b,
                      std::vector<DataType> datum );
    //! Assigns matrix diagonal element (r,r) to coeff, other elts
    //! of row r to zero, and vector b element b(r) to coeff*datum.
    template < typename VectorType >
    void diagonalize( UInt const row, DataType const coeff, VectorType &b,
                      DataType datum );

    //! Matrix-vector product.
    std::vector<DataType> operator*( const std::vector<DataType> &v ) const;
    //! Version for type Vector.
    Vector operator*( const Vector &v ) const;
    //! Version for C pointer vector. BEWARE: no check on bounds done !
    template <typename DataT, UInt BR, UInt BC, typename PatternT>
    friend void operMatVec( DataT * const mv,
                            const MixedMatr<BR, BC, PatternT, DataT> &Mat,
                            const DataT *v );
    //! necessary for IML++ library: transpose-matrix by vector product.
    Vector trans_mult( const Vector &v ) const;

    //! set to zero the matrix;
    void zeros();

private:
    //! A pointer to the pattern.
    const MixedPattern<BRows, BCols, PatternType> * _Patt;
    //! An array of values vector: there is one vector of values for
    //! each block.
    std::vector<DataType> _bValues[ BRows ][ BCols ];

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
template <UInt BRows, UInt BCols>
class MixedMatr<BRows, BCols, MSRPatt, double>
{
public:
    //! Default constructor.
    MixedMatr();
    //! Constructor from an existing MixedPattern pattern.
    //! Version for MSR format.
    MixedMatr( const MixedPattern<BRows, BCols, MSRPatt> &ex_pattern );
    //! Version for CSR format.
    MixedMatr( const MixedPattern<BRows, BCols, CSRPatt> &ex_pattern );
    //! gives the pattern pointer.
    const MixedPattern<BRows, BCols, MSRPatt> * Patt() const
    {
        return _Patt;
    }
    //! points to the values vector array.
    std::vector<double> * bValues() const
    {
        return _bValues;
    }
    //! The container of block (i,j).
    //! BEWARE: no check is done on the array size.
    const std::vector<double> & bValues( UInt i, UInt j ) const
    {
        return _bValues[ i ][ j ];
    }
    std::vector<double> & bValues( UInt i, UInt j )
    {
        return _bValues[ i ][ j ];
    }
    //! give the block (i,j) value vector in Raw format (suitable for C)
    double * giveRaw_value( UInt i, UInt j )
    {
        return & ( _bValues[ i ][ j ].front() );
    }
    //! give the block (i,j) value vector in Raw format (suitable for C)
    double const * giveRaw_value( UInt i, UInt j ) const
    {
        return & ( _bValues[ i ][ j ].front() );
    }

    //! Determines the lumped diagonal of P1 mass matrix.
    std::vector<double> MassDiagP1() const;
    //! Inverts the diagonal matrix Diag and mulitply it by the matrix Mat.
    template <UInt BR, UInt BC>
    friend void MultInvDiag( const std::vector<Real> &Diag,
                             const MixedMatr<BR, BC, MSRPatt, Real> &Mat,
                             MixedMatr<BR, BC, MSRPatt, Real> &ans );
    //! Gives the diagonal of a block matrix.
    std::vector<double> giveDiag() const;

    //! Warning: the two matrices will point to the same pattern.
    MixedMatr & operator= ( const MixedMatr<BRows, BCols, MSRPatt, double>
                            &RhMtr );
    //! Assigns the matrix to loc_val at the place of index where in the
    //! block (ib,jb).
    void set_mat( UInt ib, UInt jb, UInt where, double loc_val );
    //! Assigns the matrix element (row,col) to loc_val.
    void set_mat( UInt row, UInt col, double loc_val );
    //! Assigns the matrix element (row,col) of block (ib,jb) to loc_val.
    void set_mat( UInt ib, UInt jb, UInt row, UInt col, double loc_val );
    //! Adds loc_val to the matrix element (row,col).
    void set_mat_inc( UInt row, UInt col, double loc_val );
    //! Adds loc_val to the matrix element (row,col) of block (ib,jb).
    void set_mat_inc( UInt ib, UInt jb, UInt row, UInt col, double loc_val );
    //! Returns the matrix element (i,j) value.
    double& get_value( UInt i, UInt j );
    const double get_value( UInt i, UInt j ) const;
    //! Returns the matrix element (i,j) value of block (ib,jb).
    double& get_value( UInt ib, UInt jb, UInt i, UInt j );
    const double get_value( UInt ib, UInt jb, UInt i, UInt j ) const;
    //! Shows the matrix (only the pattern here).
    void ShowMe()
    {
        return _Patt->showMe( true );
    }
    //! write matrix in sparse matlab format and spy
    /*! just run the resulting m-file and the matrix is loaded into A
     *  and its sparsity pattern is displayed.
     *  @param filename name of the m-file
     */
    void spy( std::string const &filename );
    //! Assigns matrix diagonal element (r,r) to coeff, other elts
    //! of row r to zero.
    void diagonalize_row( UInt const r, double const coeff );
    //! assign a row to zero. Remark, zero might be defined for any double
    void zero_row( UInt const row );

    //! set entries (rVec(i),rVec(i)) to coeff and rest of row r(i) to zero
    template < typename VectorType >
    void diagonalize( std::vector<UInt> const rVec, double const coeff, VectorType &b,
		      std::vector<double> datum );

    //! Assigns matrix diagonal element (r,r) to coeff, other elts
    //! of row r to zero, and vector b element b(r) to coeff*datum.
    template < typename VectorType >
    void diagonalize( UInt const row, double const coeff, VectorType &b,
                      double datum );

    //! Matrix-vector product.
    std::vector<double> operator*( const std::vector<double> &v ) const;
    //! Version for type Vector.
    Vector operator*( const Vector &v ) const;
    //! Version for C pointer vector. BEWARE: no check on bounds done !
    template <UInt BR, UInt BC>
    friend void operMatVec( double * const mv,
                            const MixedMatr<BR, BC, MSRPatt, double> &Mat,
                            const double *v );
    //! necessary for IML++ library: transpose-matrix by vector product.
    Vector trans_mult( const Vector &v ) const;
    //! set to zero the matrix;
    void zeros();

private:
    //! A pointer to the pattern.
    const MixedPattern<BRows, BCols, MSRPatt> * _Patt;
    //! An array of values vector: there is one vector of values for
    //! each block.
    std::vector<double> _bValues[ BRows ][ BCols ];

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
template <UInt BRows, UInt BCols>
class MixedMatr<BRows, BCols, CSRPatt, double>
{
public:
    //! Default constructor.
    MixedMatr();
    //! Constructor from an existing MixedPattern pattern.
    //! Version for CSR format.
    MixedMatr( const MixedPattern<BRows, BCols, CSRPatt> &ex_pattern );
    //! gives the pattern pointer.
    const MixedPattern<BRows, BCols, CSRPatt> * Patt() const
    {
        return _Patt;
    }
    //! points to the values vector array.
    std::vector<double> * bValues() const
    {
        return _bValues;
    }
    //! The container of block (i,j).
    //! BEWARE: no check is done on the array size.
    const std::vector<double> & bValues( UInt i, UInt j ) const
    {
        return _bValues[ i ][ j ];
    }
    std::vector<double> & bValues( UInt i, UInt j )
    {
        return _bValues[ i ][ j ];
    }
    //! give the block (i,j) value vector in Raw format (suitable for C)
    double * giveRaw_value( UInt i, UInt j )
    {
        return & ( _bValues[ i ][ j ].front() );
    }
    //! give the block (i,j) value vector in Raw format (suitable for C)
    double const * giveRaw_value( UInt i, UInt j ) const
    {
        return & ( _bValues[ i ][ j ].front() );
    }

    //! Determines the lumped diagonal of P1 mass matrix.
    std::vector<double> MassDiagP1() const;
    //! Inverts the diagonal matrix Diag and mulitply it by the matrix Mat.
    template <UInt BR, UInt BC>
    friend void MultInvDiag( const std::vector<Real> &Diag,
                             const MixedMatr<BR, BC, CSRPatt, Real> &Mat,
                             MixedMatr<BR, BC, CSRPatt, Real> &ans );

    //! Gives the diagonal of a block matrix.
    std::vector<double> giveDiag() const;

    //! Warning: the two matrices will point to the same pattern.
    MixedMatr & operator= ( const MixedMatr<BRows, BCols, CSRPatt, double>
                            &RhMtr );
    //! Assigns the matrix to loc_val at the place of index where in the
    //! block (ib,jb).
    void set_mat( UInt ib, UInt jb, UInt where, double loc_val );
    //! Assigns the matrix element (row,col) to loc_val.
    void set_mat( UInt row, UInt col, double loc_val );
    //! Assigns the matrix element (row,col) of block (ib,jb) to loc_val.
    void set_mat( UInt ib, UInt jb, UInt row, UInt col, double loc_val );
    //! Adds loc_val to the matrix element (row,col).
    void set_mat_inc( UInt row, UInt col, double loc_val );
    //! Adds loc_val to the matrix element (row,col) of block (ib,jb).
    void set_mat_inc( UInt ib, UInt jb, UInt row, UInt col, double loc_val );
    //! Returns the matrix element (i,j) value.
    double& get_value( UInt i, UInt j );
    const double get_value( UInt i, UInt j ) const;
    //! Returns the matrix element (i,j) value of block (ib,jb).
    double& get_value( UInt ib, UInt jb, UInt i, UInt j );
    const double get_value( UInt ib, UInt jb, UInt i, UInt j ) const;
    //! Shows the matrix (only the pattern here).
    void ShowMe()
    {
        return _Patt->showMe( true );
    }
    //! write matrix in sparse matlab format and spy
    /*! just run the resulting m-file and the matrix is loaded into A
     *  and its sparsity pattern is displayed.
     *  @param filename name of the m-file
     */
    void spy( std::string const &filename );
    //! Assigns matrix diagonal element (r,r) to coeff, other elts
    //! of row r to zero.
    void diagonalize_row( UInt const r, double const coeff );
    //! assign a row to zero. Remark, zero might be defined for any double
    void zero_row( UInt const row );

    //! set entries (rVec(i),rVec(i)) to coeff and rest of row r(i) to zero
    template < typename VectorType >
    void diagonalize( std::vector<UInt> const rVec, double const coeff, VectorType &b,
                      std::vector<double> datum );

    //! Assigns matrix diagonal element (r,r) to coeff, other elts
    //! of row r to zero, and vector b element b(r) to coeff*datum.
    template < typename VectorType >
    void diagonalize( UInt const row, double const coeff, VectorType &b,
                      double datum );

    //!set to zero the row trD and the corresponding column=row for D
    template <UInt BR, UInt BC, typename VectorType, typename DataType>
    friend void
    zero_row_col( UInt const row, MixedMatr<BR, BC, CSRPatt, Real> &trD,
                  MixedMatr<BC, BR, CSRPatt, Real> &D, VectorType &bp,
                  DataType const datum );

    //! Matrix-vector product.
    std::vector<double> operator*( const std::vector<double> &v ) const;
    //! Version for type Vector.
    Vector operator*( const Vector &v ) const;
    //! Version for C pointer vector. BEWARE: no check on bounds done !
    template <UInt BR, UInt BC>
    friend void operMatVec( double * const mv,
                            const MixedMatr<BR, BC, CSRPatt, double> &Mat,
                            const double *v );
    //! necessary for IML++ library: transpose-matrix by vector product.
    Vector trans_mult( const Vector &v ) const;
    //! set to zero the matrix;
    void zeros();

private:
    //! A pointer to the pattern.
    const MixedPattern<BRows, BCols, CSRPatt> * _Patt;
    //! An array of values vector: there is one vector of values for
    //! each block.
    std::vector<double> _bValues[ BRows ][ BCols ];

};
// MixedMatr
//-----------------------------------------------------------------------

//Default Constructor
template <UInt BRows, UInt BCols, typename PatternType, typename DataType>
MixedMatr<BRows, BCols, PatternType, DataType>::
MixedMatr()
{
    for ( UInt i = 0; i < BRows; i++ )
        for ( UInt j = 0; j < BCols; j++ )
            _bValues[ i ][ j ].resize( 0 );
}

//Constructor from an existing external pattern.
//version for MSR format: size of values= nnz+1.
template <UInt BRows, UInt BCols, typename PatternType, typename DataType>
MixedMatr<BRows, BCols, PatternType, DataType>::
MixedMatr( const MixedPattern<BRows, BCols, MSRPatt> &ex_pattern )
{
    _Patt = &ex_pattern;

    for ( UInt i = 0; i < BRows; i++ )
        for ( UInt j = 0; j < BCols; j++ )
            _bValues[ i ][ j ].resize( ex_pattern.nNz( i, j ) + 1 );
}
//Constructor from an existing external pattern.
//version for CSR format: size of values= nnz.
template <UInt BRows, UInt BCols, typename PatternType, typename DataType>
MixedMatr<BRows, BCols, PatternType, DataType>::
MixedMatr( const MixedPattern<BRows, BCols, CSRPatt> &ex_pattern )
{
    _Patt = &ex_pattern;

    for ( UInt i = 0; i < BRows; i++ )
        for ( UInt j = 0; j < BCols; j++ )
            _bValues[ i ][ j ].resize( ex_pattern.nNz( i, j ) );
}

// Determines the lumped diagonal of P1 mass matrix.
template <UInt BRows, UInt BCols, typename PatternType, typename DataType>
std::vector<DataType>
MixedMatr<BRows, BCols, PatternType, DataType>::
MassDiagP1() const
{
    UInt nrows = _Patt->nRows();
    UInt nnz = 0, nr_global = 0;
    Container coldata, position;

    std::vector<DataType> diag;
    diag.resize( nrows, 0.0 );

    for ( UInt ib = 0; ib < BRows; ib++ )
    {

        nrows = _Patt->nRows( ib, 0 );

        for ( UInt jb = 0; jb < BCols; jb++ )
        {

            if ( _Patt->block_ptr( ib, jb ) != 0 )
            {

                for ( UInt nr = 0; nr < nrows; nr++ )
                {
                    coldata.resize( _Patt->nCols( ib, jb ) );
                    position.resize( _Patt->nCols( ib, jb ) );
                    nnz = _Patt->row( ib, jb, nr, coldata.begin(), position.begin() );
                    for ( UInt jcol = 0; jcol < nnz; jcol++ )
                        diag[ nr_global + nr ] +=
                            _bValues[ ib ][ jb ][ position[ jcol ] ];
                }
            }
        }
        nr_global += nrows;
    }

    return diag;
}

// Inverts the diagonal matrix Diag and mulitply it by the matrix Mat.
template <UInt BR, UInt BC, typename PattType>
void
MultInvDiag( const std::vector<Real> &Diag,
             const MixedMatr<BR, BC, PattType, Real> &Mat,
             MixedMatr<BR, BC, PattType, Real> &ans )
{
    ASSERT( find( Diag.begin(), Diag.end(), 0 ) == Diag.end(), "Diagonal matrix Diag must be invertible" );

    ASSERT( Diag.size() == Mat._Patt->nRows(), "Matrix sizes not compatible" );

    // Product:
    UInt nrows = 0;
    UInt nnz = 0, nr_global = 0;
    Container coldata, pos;

    for ( UInt ib = 0; ib < BR; ib++ )
    {

        nrows = Mat._Patt->nRows( ib, 0 );

        for ( UInt jb = 0; jb < BC; jb++ )
        {

            if ( Mat._Patt->block_ptr( ib, jb ) != 0 )
            {

                for ( UInt nr = 0; nr < nrows; nr++ )
                {
                    coldata.resize( Mat._Patt->nCols( ib, jb ) );
                    pos.resize( Mat._Patt->nCols( ib, jb ) );
                    nnz = Mat._Patt->row( ib, jb, nr, coldata.begin(), pos.begin() );

                    for ( UInt jcol = 0; jcol < nnz; jcol++ )
                        ans._bValues[ ib ][ jb ][ pos[ jcol ] ] =
                            ( Mat._bValues[ ib ][ jb ][ pos[ jcol ] ] / Diag[ nr_global + nr ] );
                }
            }
        }
        nr_global += nrows;
    }
}

// Gives the diagonal of a block matrix.
template <UInt BRows, UInt BCols, typename PatternType, typename DataType>
std::vector<DataType>
MixedMatr<BRows, BCols, PatternType, DataType>::
giveDiag() const
{
    ASSERT_PRE( BRows == BCols, "block matrix must be a square matrix" );

    std::vector<Real> diag;
    diag.resize( _Patt->nRows(), 0.0 );

    Container coldata, pos;
    UInt nrows = 0, nr_global = 0, nnz = 0;

    for ( UInt ib = 0; ib < BRows; ib++ )
    {
        ASSERT_PRE( _Patt->nRows( ib, ib ) == _Patt->nCols( ib, ib ),
                    "matrix must have nrows=ncols for diagonal blocks" );

        nrows = _Patt->nRows( ib, ib );

        if ( _Patt->block_ptr( ib, ib ) != 0 )
            for ( UInt nr = 0; nr < nrows; nr++ )
            {

                coldata.resize( _Patt->nCols( ib, ib ) );
                pos.resize( _Patt->nCols( ib, ib ) );

                nnz = _Patt->row( ib, ib, nr, coldata.begin(), pos.begin() );

                UInt jcol = 0;

                if ( coldata[ jcol ] - OFFSET == nr )
                    diag[ nr_global + nr ] = _bValues[ ib ][ ib ][ pos[ jcol ] ];
                else
                {
                    while ( coldata[ jcol ] - OFFSET != nr )
                        jcol++;
                    diag[ nr_global + nr ] = _bValues[ ib ][ ib ][ pos[ jcol ] ];
                }
            }
        nr_global += nrows;
    }
    return diag;
}


// Warning: the two matrices will point to the same pattern.
template <UInt BRows, UInt BCols, typename PatternType, typename DataType>
MixedMatr<BRows, BCols, PatternType, DataType> &
MixedMatr<BRows, BCols, PatternType, DataType>::
operator=( const MixedMatr<BRows, BCols, PatternType, DataType> &RhMtr )
{
    if ( &RhMtr != this )
    {
        UInt ib, jb;
        for ( ib = 0; ib < BRows; ib++ )
            for ( jb = 0; jb < BCols; jb++ )
                _bValues[ ib ][ jb ] = RhMtr.bValues( ib, jb );

        _Patt = RhMtr.Patt();
    }

    return *this;

}

// Assigns the matrix to loc_val at the place of index where.
template <UInt BRows, UInt BCols, typename PatternType, typename DataType>
void
MixedMatr<BRows, BCols, PatternType, DataType>::
set_mat( UInt ib, UInt jb, UInt where, DataType loc_val )
{
    _bValues[ ib ][ jb ][ where - OFFSET ] = loc_val;
}

// Assigns the matrix element (row,col) to loc_val.
template <UInt BRows, UInt BCols, typename PatternType, typename DataType>
void
MixedMatr<BRows, BCols, PatternType, DataType>::
set_mat( UInt row, UInt col, DataType loc_val )
{
    UInt m, n, lr, lc;
    extract_pair( _Patt->locateElBlock( row, col ), m, n );
    extract_pair( _Patt->localNumber( m, n, row, col ), lr, lc );
    set_mat( m, n, lr, lc, loc_val );
}

// Assigns the matrix element (row,col) of block (ib,jb) to loc_val.
template <UInt BRows, UInt BCols, typename PatternType, typename DataType>
void
MixedMatr<BRows, BCols, PatternType, DataType>::
set_mat( UInt ib, UInt jb, UInt row, UInt col, DataType loc_val )
{
    if ( _Patt->block_ptr( ib, jb ) != 0 )
    {
        std::pair<UInt, bool> where = _Patt->block_ptr( ib, jb ) ->locate_index( row, col );
        if ( where.second )
            _bValues[ ib ][ jb ][ where.first ] = loc_val;
    }
}

// Adds loc_val to the matrix element (row,col).
template <UInt BRows, UInt BCols, typename PatternType, typename DataType>
void
MixedMatr<BRows, BCols, PatternType, DataType>::
set_mat_inc( UInt row, UInt col, DataType loc_val )
{
    UInt m, n, lr, lc;
    extract_pair( _Patt->locateElBlock( row, col ), m, n );
    extract_pair( _Patt->localNumber( m, n, row, col ), lr, lc );
    set_mat_inc( m, n, lr, lc, loc_val );
}

// Adds loc_val to the matrix element (row,col) of block (ib,jb).
template <UInt BRows, UInt BCols, typename PatternType, typename DataType>
void
MixedMatr<BRows, BCols, PatternType, DataType>::
set_mat_inc( UInt ib, UInt jb, UInt row, UInt col, DataType loc_val )
{
    if ( _Patt->block_ptr( ib, jb ) != 0 )
    {
        std::pair<UInt, bool> where = _Patt->block_ptr( ib, jb ) ->locate_index( row, col );
        if ( where.second )
            _bValues[ ib ][ jb ][ where.first ] += loc_val;
    }
}

// Returns the matrix element (i,j) value.
template <UInt BRows, UInt BCols, typename PatternType, typename DataType>
DataType&
MixedMatr<BRows, BCols, PatternType, DataType>::
get_value( UInt i, UInt j )
{
    UInt m, n, lr, lc;
    extract_pair( _Patt->locateElBlock( i, j ), m, n );
    extract_pair( _Patt->localNumber( m, n, i, j ), lr, lc );
    return get_value( m, n, lr, lc );
}
// const qualifyer version
template <UInt BRows, UInt BCols, typename PatternType, typename DataType>
const DataType&
MixedMatr<BRows, BCols, PatternType, DataType>::
get_value( UInt i, UInt j ) const
{
    UInt m, n, lr, lc;
    extract_pair( _Patt->locateElBlock( i , j ), m, n );
    extract_pair( _Patt->localNumber( m, n, i, j ), lr, lc );
    return get_value( m, n, i, j );
}

// Returns the matrix element (i,j) value of block (ib,jb).
template <UInt BRows, UInt BCols, typename PatternType, typename DataType>
DataType&
MixedMatr<BRows, BCols, PatternType, DataType>::
get_value( UInt ib, UInt jb, UInt i, UInt j )
{
    if ( _Patt->block_ptr( ib, jb ) != 0 )
        return _bValues[ ib ][ jb ][ _Patt->block_ptr( ib, jb ) ->locate_index( i, j ).first ];
    else
        return 0.0;
}
// const qualifyer version
template <UInt BRows, UInt BCols, typename PatternType, typename DataType>
const DataType&
MixedMatr<BRows, BCols, PatternType, DataType>::
get_value( UInt ib, UInt jb, UInt i, UInt j ) const
{
    if ( _Patt->block_ptr( ib, jb ) != 0 )
        return _bValues[ ib ][ jb ][ _Patt->block_ptr( ib, jb ) ->locate_index( i, j ).first ];
    else
        return 0.0;
}

// Matrix visualization a la matlab.
template <UInt BRows, UInt BCols, typename PatternType, typename DataType>
void
MixedMatr<BRows, BCols, PatternType, DataType>::
spy( std::string const &filename )
{
    // Purpose: Matlab dumping and spy
    std::string nome = filename, uti = " , ";
    //
    // check on the file name
    //
    UInt i = filename.find( "." );

    if ( i <= 0 )
        nome = filename + ".m";
    else
    {
        if ( i != filename.size() - 2 || filename[ i + 1 ] != 'm' )
        {
            std::cerr << "Wrong file name ";
            nome = filename + ".m";
        }
    }

    std::ofstream file_out( nome.c_str() );
    UInt nnz, mb, nb;
    Container coldata, pos;
    coldata.resize( _Patt->nCols() );
    pos.resize( _Patt->nCols() );

    file_out << "S = [ ";
    for ( UInt i = 0; i < _Patt->nRows(); ++i )
    {
        nnz = _Patt->row( i, coldata.begin(), pos.begin() );
        for ( UInt j = 0; j < nnz; ++j )
        {
            extract_pair( _Patt->locateElBlock( i, coldata[ j ] - OFFSET ), mb, nb );
            file_out << i + 1 << uti << coldata[ j ] + 1 - OFFSET << uti <<
            _bValues[ mb ][ nb ][ pos[ j ] ] << std::endl;
        }
    }
    file_out << "];" << std::endl;
    file_out << "I=S(:,1); J=S(:,2); S=S(:,3); A=sparse(I,J,S); spy(A);" << std::endl;
}

// Assigns matrix diagonal element (r,r) to coeff, other elts
// of row r to zero.
template <UInt BRows, UInt BCols, typename PatternType, typename DataType>
void
MixedMatr<BRows, BCols, PatternType, DataType>::
diagonalize_row( UInt const r, DataType const coeff )
{
    UInt nnz = 0, m = 0, n = 0;
    Container coldata, pos;
    coldata.resize( _Patt->nCols() );
    pos.resize( _Patt->nCols() );

    nnz = _Patt->row( r - OFFSET, coldata.begin(), pos.begin() );

    UInt jcol = 0;
    if ( coldata[ jcol ] == r )
    {
        //diagonal element:
        extract_pair( _Patt->locateElBlock( r - OFFSET , coldata[ jcol ] - OFFSET ), m, n );
        _bValues[ m ][ n ][ pos[ jcol ] ] = coeff;
        //other elements:
        for ( jcol = 1; jcol < nnz; jcol++ )
        {
            extract_pair( _Patt->locateElBlock( r - OFFSET , coldata[ jcol ] - OFFSET ), m, n );
            _bValues[ m ][ n ][ pos[ jcol ] ] = 0.0;
        }
    }
    else
        for ( jcol = 0; jcol < nnz; jcol++ )
        {
            extract_pair( _Patt->locateElBlock( r - OFFSET , coldata[ jcol ] - OFFSET ), m, n );
            if ( coldata[ jcol ] == r )
                //diagonal element:
                _bValues[ m ][ n ][ pos[ jcol ] ] = coeff;
            else
                //other elements
                _bValues[ m ][ n ][ pos[ jcol ] ] = 0.0;
        }
}

// assign a row to zero. Remark, zero might be defined for any DataType
template <UInt BRows, UInt BCols, typename PatternType, typename DataType>
void
MixedMatr<BRows, BCols, PatternType, DataType>::
zero_row( UInt const row )
{
    diagonalize_row( row, 0.0 );
}

template <UInt BRows, UInt BCols, typename PatternType, typename DataType>
template < typename VectorType >
void
MixedMatr<BRows, BCols, PatternType, DataType>::
diagonalize( std::vector<UInt> const rVec,
				DataType const coeff,
				VectorType &b,
				std::vector<DataType> datumVec )
{
     UInt sizeVec(rVec.size());
     if ( sizeVec != datumVec.size()) 
     { //! vectors must be of the same size
         ERROR_MSG( "diagonalize: vectors must be of the same size\n" );
     }
      
     for (UInt i=0; i < sizeVec; i++)
       diagonalize( rVec[i], coeff, b, datumVec[i]);
      
}

// Assigns matrix diagonal element (r,r) to coeff, other elts
// of row r to zero, and vector b element b(r) to coeff*datum.
template <UInt BRows, UInt BCols, typename PatternType, typename DataType>
template < typename VectorType >
void
MixedMatr<BRows, BCols, PatternType, DataType>::
diagonalize( UInt const row, DataType const coeff, VectorType &b,
             DataType datum )
{
    UInt nnz = 0, m = 0, n = 0;
    Container coldata, pos;
    coldata.resize( _Patt->nCols() );
    pos.resize( _Patt->nCols() );

    nnz = _Patt->row( row - OFFSET, coldata.begin(), pos.begin() );

    UInt jcol = 0;
    if ( coldata[ jcol ] == row )
    { //case diagfirst
        //diagonal element:
        extract_pair( _Patt->locateElBlock( row - OFFSET , coldata[ jcol ] - OFFSET ), m, n );
        _bValues[ m ][ n ][ pos[ jcol ] ] = coeff;
        //other elements:
        for ( jcol = 1; jcol < nnz; jcol++ )
        {
            extract_pair( _Patt->locateElBlock( row - OFFSET , coldata[ jcol ] - OFFSET ), m, n );
            _bValues[ m ][ n ][ pos[ jcol ] ] = 0.0;
        }
    }
    else
        for ( jcol = 0; jcol < nnz; jcol++ )
        {
            extract_pair( _Patt->locateElBlock( row - OFFSET , coldata[ jcol ] - OFFSET ), m, n );
            if ( coldata[ jcol ] == row )
                //diagonal element:
                _bValues[ m ][ n ][ pos[ jcol ] ] = coeff;
            else
                //other elements
                _bValues[ m ][ n ][ pos[ jcol ] ] = 0.0;
        }

    //rhs:
    b[ row - OFFSET ] = coeff * datum;
}

// Matrix-vector product.
template <UInt BRows, UInt BCols, typename PatternType, typename DataType>
std::vector<DataType>
MixedMatr<BRows, BCols, PatternType, DataType>::operator*( const std::vector<DataType> &v ) const
{
    ASSERT( _Patt->nCols() == v.size(), "Error in Matrix Vector product" );
    std::vector<DataType> ans;
    ans.resize( _Patt->nRows(), 0.0 );

    UInt nrows = 0, ncols = 0;
    UInt nnz = 0, nr_global = 0, nc_global;
    Container coldata, pos;

    for ( UInt ib = 0; ib < BRows; ib++ )
    {

        nrows = _Patt->nRows( ib, 0 );
        nc_global = 0;

        for ( UInt jb = 0; jb < BCols; jb++ )
        {

            ncols = _Patt->nCols( ib, jb );

            if ( _Patt->block_ptr( ib, jb ) != 0 )
            {

                for ( UInt nr = 0; nr < nrows; nr++ )
                {
                    coldata.resize( ncols );
                    pos.resize( ncols );
                    nnz = _Patt->row( ib, jb, nr, coldata.begin(), pos.begin() );

                    for ( UInt jcol = 0; jcol < nnz; jcol++ )
                        ans[ nr_global + nr ] += _bValues[ ib ][ jb ][ pos[ jcol ] ] *
                                                 v[ nc_global + coldata[ jcol ] - OFFSET ];
                }
            }
            nc_global += ncols;
        }
        nr_global += nrows;
    }
    return ans;
}
// version for type Vector
template <UInt BRows, UInt BCols, typename PatternType, typename DataType>
Vector
MixedMatr<BRows, BCols, PatternType, DataType>::
operator*( const Vector &v ) const
{
    ASSERT( _Patt->nCols() == v.size(), "Error in Matrix Vector product" );
    Vector ans( _Patt->nRows() );
    ans = 0.0;

    UInt nrows = 0, ncols = 0;
    UInt nnz = 0, nr_global = 0, nc_global;
    Container coldata, pos;

    for ( UInt ib = 0; ib < BRows; ib++ )
    {

        nrows = _Patt->nRows( ib, 0 );
        nc_global = 0;

        for ( UInt jb = 0; jb < BCols; jb++ )
        {

            ncols = _Patt->nCols( ib, jb );

            if ( _Patt->block_ptr( ib, jb ) != 0 )
            {

                for ( UInt nr = 0; nr < nrows; nr++ )
                {
                    coldata.resize( ncols );
                    pos.resize( ncols );
                    nnz = _Patt->row( ib, jb, nr, coldata.begin(), pos.begin() );

                    for ( UInt jcol = 0; jcol < nnz; jcol++ )
                        ans[ nr_global + nr ] += _bValues[ ib ][ jb ][ pos[ jcol ] ] *
                                                 v[ nc_global + coldata[ jcol ] - OFFSET ];
                }
            }
            nc_global += ncols;
        }
        nr_global += nrows;
    }
    return ans;
}
// Version for C pointer vector. BEWARE: no check on bounds is done !
template <typename DataType, UInt BRows, UInt BCols, typename PatternType>
void operMatVec( DataType * const mv,
                 const MixedMatr<BRows, BCols, PatternType, DataType> &Mat,
                 const DataType *v )
{
    UInt nrows = 0, ncols = 0;
    UInt nnz = 0, nr_global = 0, nc_global;
    Container coldata, pos;

    for ( UInt ib = 0; ib < BRows; ib++ )
    {

        nrows = Mat._Patt->nRows( ib, 0 );
        nc_global = 0;

        for ( UInt nr = 0; nr < nrows; nr++ )
        {

            //initialize
            mv[ nr + ib * nrows ] = 0.;

            for ( UInt jb = 0; jb < BCols; jb++ )
            {

                ncols = Mat._Patt->nCols( ib, jb );

                if ( Mat._Patt->block_ptr( ib, jb ) != 0 )
                {

                    coldata.resize( ncols );
                    pos.resize( ncols );
                    nnz = Mat._Patt->row( ib, jb, nr, coldata.begin(), pos.begin() );

                    for ( UInt jcol = 0; jcol < nnz; jcol++ )
                        mv[ nr_global + nr ] += Mat._bValues[ ib ][ jb ][ pos[ jcol ] ] *
                                                v[ nc_global + coldata[ jcol ] - OFFSET ];
                }
            }
            nc_global += ncols;
        }
        nr_global += nrows;
    }
}

//necessary for IML++ library: transpose-matrix by vector product.
template <UInt BRows, UInt BCols, typename PatternType, typename DataType>
Vector
MixedMatr<BRows, BCols, PatternType, DataType>::
trans_mult( const Vector &v ) const
{
    ASSERT( _Patt->nRows() == v.size(), "Error in Matrix Vector product" );
    Vector ans( _Patt->nRows() );
    ans = 0.0;

    UInt nrows = 0, ncols = 0;
    UInt nnz = 0, nr_global = 0, nc_global;
    Container coldata, pos;

    for ( UInt ib = 0; ib < BRows; ib++ )
    {

        nrows = _Patt->nRows( ib, 0 );
        nc_global = 0;

        for ( UInt jb = 0; jb < BCols; jb++ )
        {

            ncols = _Patt->nCols( ib, jb );

            if ( _Patt->block_ptr( ib, jb ) != 0 )
            {

                for ( UInt nr = 0; nr < nrows; nr++ )
                {
                    coldata.resize( ncols );
                    pos.resize( ncols );
                    nnz = _Patt->row( ib, jb, nr, coldata.begin(), pos.begin() );

                    for ( UInt jcol = 0; jcol < nnz; jcol++ )
                        ans[ nc_global + coldata[ jcol ] - OFFSET ] += _bValues[ ib ][ jb ][ pos[ jcol ] ] *
                                v[ nr_global + nr ];
                }
            }
            nc_global += ncols;
        }
        nr_global += nrows;
    }
    return ans;
}

template <UInt BRows, UInt BCols, typename PatternType, typename DataType>
void MixedMatr<BRows, BCols, PatternType, DataType>::
zeros()
{
    UInt ib, jb;
    for ( ib = 0; ib < BRows; ib++ )
        for ( jb = 0; jb < BCols; jb++ )
            if ( _Patt->block_ptr( ib, jb ) != 0 )
            {
                typename std::vector<DataType>::iterator start = _bValues[ ib ][ jb ].begin();
                typename std::vector<DataType>::iterator end = _bValues[ ib ][ jb ].end();
                fill( start, end, 0.0 );
            }
}

//-----------------------------------------------------------------------
// MixedMatr SPECIALIZATION FOR MSR
//-----------------------------------------------------------------------

//Default Constructor
template <UInt BRows, UInt BCols>
MixedMatr<BRows, BCols, MSRPatt, double>::
MixedMatr()
{
    for ( UInt i = 0; i < BRows; i++ )
        for ( UInt j = 0; j < BCols; j++ )
            _bValues[ i ][ j ].resize( 0 );
}

//Constructor from an existing external pattern.
//version for MSR format: size of values= nnz+1.
template <UInt BRows, UInt BCols>
MixedMatr<BRows, BCols, MSRPatt, double>::
MixedMatr( const MixedPattern<BRows, BCols, MSRPatt> &ex_pattern )
{
    _Patt = &ex_pattern;

    for ( UInt i = 0; i < BRows; i++ )
        for ( UInt j = 0; j < BCols; j++ )
            _bValues[ i ][ j ].resize( ex_pattern.nNz( i, j ) + 1 );
}

// Determines the lumped diagonal of P1 mass matrix.
template <UInt BRows, UInt BCols>
std::vector<double>
MixedMatr<BRows, BCols, MSRPatt, double>::
MassDiagP1() const
{
    UInt nrows = _Patt->nRows();
    UInt nnz = 0, nr_global = 0;
    Container coldata, position;

    std::vector<double> diag;
    diag.resize( nrows, 0.0 );

    for ( UInt ib = 0; ib < BRows; ib++ )
    {

        nrows = _Patt->nRows( ib, 0 );

        for ( UInt jb = 0; jb < BCols; jb++ )
        {

            if ( _Patt->block_ptr( ib, jb ) != 0 )
            {

                for ( UInt nr = 0; nr < nrows; nr++ )
                {
                    coldata.resize( _Patt->nCols( ib, jb ) );
                    position.resize( _Patt->nCols( ib, jb ) );
                    nnz = _Patt->row( ib, jb, nr, coldata.begin(), position.begin() );
                    for ( UInt jcol = 0; jcol < nnz; jcol++ )
                        diag[ nr_global + nr ] +=
                            _bValues[ ib ][ jb ][ position[ jcol ] ];
                }
            }
        }
        nr_global += nrows;
    }

    return diag;
}

// Inverts the diagonal matrix Diag and mulitply it by the matrix Mat.
template <UInt BR, UInt BC>
void
MultInvDiag( const std::vector<Real> &Diag,
             const MixedMatr<BR, BC, MSRPatt, Real> &Mat,
             MixedMatr<BR, BC, MSRPatt, Real> &ans )
{
    ASSERT( find( Diag.begin(), Diag.end(), 0 ) == Diag.end(), "Diagonal matrix Diag must be invertible" );

    ASSERT( Diag.size() == Mat._Patt->nRows(), "Matrix sizes not compatible" );

    // Product:
    UInt nrows = 0;
    UInt nnz = 0, nr_global = 0;
    Container coldata, pos;

    for ( UInt ib = 0; ib < BR; ib++ )
    {

        nrows = Mat._Patt->nRows( ib, 0 );

        for ( UInt jb = 0; jb < BC; jb++ )
        {

            if ( Mat._Patt->block_ptr( ib, jb ) != 0 )
            {

                for ( UInt nr = 0; nr < nrows; nr++ )
                {
                    coldata.resize( Mat._Patt->nCols( ib, jb ) );
                    pos.resize( Mat._Patt->nCols( ib, jb ) );
                    nnz = Mat._Patt->row( ib, jb, nr, coldata.begin(), pos.begin() );

                    for ( UInt jcol = 0; jcol < nnz; jcol++ )
                        ans._bValues[ ib ][ jb ][ pos[ jcol ] ] =
                            ( Mat._bValues[ ib ][ jb ][ pos[ jcol ] ] / Diag[ nr_global + nr ] );
                }
            }
        }
        nr_global += nrows;
    }
}

// Gives the diagonal of a block matrix.
template <UInt BRows, UInt BCols>
std::vector<double>
MixedMatr<BRows, BCols, MSRPatt, double>::
giveDiag() const
{
    ASSERT_PRE( BRows == BCols, "block matrix must be a square matrix" );

    std::vector<Real> diag;
    diag.resize( _Patt->nRows(), 0.0 );

    Container coldata, pos;
    UInt nrows = 0, nr_global = 0, nnz = 0;

    for ( UInt ib = 0; ib < BRows; ib++ )
    {
        ASSERT_PRE( _Patt->nRows( ib, ib ) == _Patt->nCols( ib, ib ),
                    "matrix must have nrows=ncols for diagonal blocks" );

        nrows = _Patt->nRows( ib, ib );

        if ( _Patt->block_ptr( ib, ib ) != 0 )
            for ( UInt nr = 0; nr < nrows; nr++ )
            {

                coldata.resize( _Patt->nCols( ib, ib ) );
                pos.resize( _Patt->nCols( ib, ib ) );

                nnz = _Patt->row( ib, ib, nr, coldata.begin(), pos.begin() );

                UInt jcol = 0;

                if ( coldata[ jcol ] - OFFSET == nr )
                    diag[ nr_global + nr ] = _bValues[ ib ][ ib ][ pos[ jcol ] ];
                else
                {
                    while ( coldata[ jcol ] - OFFSET != nr )
                        jcol++;
                    diag[ nr_global + nr ] = _bValues[ ib ][ ib ][ pos[ jcol ] ];
                }
            }
        nr_global += nrows;
    }
    return diag;
}


// Warning: the two matrices will point to the same pattern.
template <UInt BRows, UInt BCols>
MixedMatr<BRows, BCols, MSRPatt, double> &
MixedMatr<BRows, BCols, MSRPatt, double>::
operator=( const MixedMatr<BRows, BCols, MSRPatt, double> &RhMtr )
{
    if ( &RhMtr != this )
    {
        UInt ib, jb;
        for ( ib = 0; ib < BRows; ib++ )
            for ( jb = 0; jb < BCols; jb++ )
                _bValues[ ib ][ jb ] = RhMtr.bValues( ib, jb );

        _Patt = RhMtr.Patt();
    }

    return *this;

}

// Assigns the matrix to loc_val at the place of index where.
template <UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, MSRPatt, double>::
set_mat( UInt ib, UInt jb, UInt where, double loc_val )
{
    _bValues[ ib ][ jb ][ where - OFFSET ] = loc_val;
}

// Assigns the matrix element (row,col) to loc_val.
template <UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, MSRPatt, double>::
set_mat( UInt row, UInt col, double loc_val )
{
    UInt m, n, lr, lc;
    extract_pair( _Patt->locateElBlock( row, col ), m, n );
    extract_pair( _Patt->localNumber( m, n, row, col ), lr, lc );
    set_mat( m, n, lr, lc, loc_val );
}

// Assigns the matrix element (row,col) of block (ib,jb) to loc_val.
template <UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, MSRPatt, double>::
set_mat( UInt ib, UInt jb, UInt row, UInt col, double loc_val )
{
    if ( _Patt->block_ptr( ib, jb ) != 0 )
    {
        std::pair<UInt, bool> where = _Patt->block_ptr( ib, jb ) ->locate_index( row, col );
        if ( where.second )
            _bValues[ ib ][ jb ][ where.first ] = loc_val;
    }
}

// Adds loc_val to the matrix element (row,col).
template <UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, MSRPatt, double>::
set_mat_inc( UInt row, UInt col, double loc_val )
{
    UInt m, n, lr, lc;
    extract_pair( _Patt->locateElBlock( row, col ), m, n );
    extract_pair( _Patt->localNumber( m, n, row, col ), lr, lc );
    set_mat_inc( m, n, lr, lc, loc_val );
}

// Adds loc_val to the matrix element (row,col) of block (ib,jb).
template <UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, MSRPatt, double>::
set_mat_inc( UInt ib, UInt jb, UInt row, UInt col, double loc_val )
{
    if ( _Patt->block_ptr( ib, jb ) != 0 )
    {
        std::pair<UInt, bool> where = _Patt->block_ptr( ib, jb ) ->locate_index( row, col );
        if ( where.second )
            _bValues[ ib ][ jb ][ where.first ] += loc_val;
    }
}

// Returns the matrix element (i,j) value.
template <UInt BRows, UInt BCols>
double&
MixedMatr<BRows, BCols, MSRPatt, double>::
get_value( UInt i, UInt j )
{
    UInt m, n, lr, lc;
    extract_pair( _Patt->locateElBlock( i, j ), m, n );
    extract_pair( _Patt->localNumber( m, n, i, j ), lr, lc );
    return get_value( m, n, lr, lc );
}
// const qualifyer version
template <UInt BRows, UInt BCols>
const double
MixedMatr<BRows, BCols, MSRPatt, double>::
get_value( UInt i, UInt j ) const
{
    UInt m, n, lr, lc;
    extract_pair( _Patt->locateElBlock( i , j ), m, n );
    extract_pair( _Patt->localNumber( m, n, i, j ), lr, lc );
    return get_value( m, n, i, j );
}

// Returns the matrix element (i,j) value of block (ib,jb).
template <UInt BRows, UInt BCols>
double&
MixedMatr<BRows, BCols, MSRPatt, double>::
get_value( UInt ib, UInt jb, UInt i, UInt j )
{
    if ( _Patt->block_ptr( ib, jb ) != 0 )
        return _bValues[ ib ][ jb ][ _Patt->block_ptr( ib, jb ) ->locate_index( i, j ).first ];
    else
        return 0.0;
}
// const qualifyer version
template <UInt BRows, UInt BCols>
const double
MixedMatr<BRows, BCols, MSRPatt, double>::
get_value( UInt ib, UInt jb, UInt i, UInt j ) const
{
    if ( _Patt->block_ptr( ib, jb ) != 0 )
        return _bValues[ ib ][ jb ][ _Patt->block_ptr( ib, jb ) ->locate_index( i, j ).first ];
    else
        return 0.0;
}

// Matrix visualization a la matlab.
template <UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, MSRPatt, double>::
spy( std::string const &filename )
{
    // Purpose: Matlab dumping and spy
    std::string nome = filename, uti = " , ";
    //
    // check on the file name
    //
    UInt i = filename.find( "." );

    if ( i <= 0 )
        nome = filename + ".m";
    else
    {
        if ( i != filename.size() - 2 || filename[ i + 1 ] != 'm' )
        {
            std::cerr << "Wrong file name ";
            nome = filename + ".m";
        }
    }

    std::ofstream file_out( nome.c_str() );
    UInt nnz, mb, nb;
    Container coldata, pos;
    coldata.resize( _Patt->nCols() );
    pos.resize( _Patt->nCols() );

    file_out << "S = [ ";
    for ( UInt i = 0; i < _Patt->nRows(); ++i )
    {
        nnz = _Patt->row( i, coldata.begin(), pos.begin() );
        for ( UInt j = 0; j < nnz; ++j )
        {
            extract_pair( _Patt->locateElBlock( i, coldata[ j ] - OFFSET ), mb, nb );
            file_out << i + 1 << uti << coldata[ j ] + 1 - OFFSET << uti <<
            _bValues[ mb ][ nb ][ pos[ j ] ] << std::endl;
        }
    }
    file_out << "];" << std::endl;
    file_out << "I=S(:,1); J=S(:,2); S=S(:,3); A=sparse(I,J,S); spy(A);" << std::endl;
}

// Assigns matrix diagonal element (r,r) to coeff, other elts
// of row r to zero.
template <UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, MSRPatt, double>::
diagonalize_row( UInt const r, double const coeff )
{
    UInt nnz = 0, m = 0, n = 0;
    Container coldata, pos;
    coldata.resize( _Patt->nCols() );
    pos.resize( _Patt->nCols() );

    nnz = _Patt->row( r - OFFSET, coldata.begin(), pos.begin() );

    UInt jcol = 0;
    if ( coldata[ jcol ] == r )
    {
        //diagonal element:
        extract_pair( _Patt->locateElBlock( r - OFFSET , coldata[ jcol ] - OFFSET ), m, n );
        _bValues[ m ][ n ][ pos[ jcol ] ] = coeff;
        //other elements:
        for ( jcol = 1; jcol < nnz; jcol++ )
        {
            extract_pair( _Patt->locateElBlock( r - OFFSET , coldata[ jcol ] - OFFSET ), m, n );
            _bValues[ m ][ n ][ pos[ jcol ] ] = 0.0;
        }
    }
    else
        for ( jcol = 0; jcol < nnz; jcol++ )
        {
            extract_pair( _Patt->locateElBlock( r - OFFSET , coldata[ jcol ] - OFFSET ), m, n );
            if ( coldata[ jcol ] == r )
                //diagonal element:
                _bValues[ m ][ n ][ pos[ jcol ] ] = coeff;
            else
                //other elements
                _bValues[ m ][ n ][ pos[ jcol ] ] = 0.0;
        }
}

// assign a row to zero. Remark, zero might be defined for any double
template <UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, MSRPatt, double>::
zero_row( UInt const row )
{
    diagonalize_row( row, 0.0 );
}

template <UInt BRows, UInt BCols>
template < typename VectorType >
void
MixedMatr<BRows, BCols, MSRPatt, double>::
diagonalize( std::vector<UInt> const rVec,
				double const coeff,
				VectorType &b,
				std::vector<double> datumVec )
{
     UInt sizeVec(rVec.size());
     if ( sizeVec != datumVec.size()) 
     { //! vectors must be of the same size
         ERROR_MSG( "diagonalize: vectors must be of the same size\n" );
     }
      

     for (UInt i=0; i < sizeVec; i++)
       diagonalize( rVec[i], coeff, b, datumVec[i]);
      
}

// Assigns matrix diagonal element (r,r) to coeff, other elts
// of row r to zero, and vector b element b(r) to coeff*datum.
// version of diagonalize for MSR format (performs the symmetric acces
// to the columns)
// Alain. nov. 2002.
template <UInt BRows, UInt BCols>
template < typename VectorType >
void
MixedMatr<BRows, BCols, MSRPatt, double>::
diagonalize( UInt const row, double const coeff, VectorType &b,
             double datum )
{
    UInt lr, lc, lr_ybind, lc_ybind;
    UInt m, n, disp;
    UInt r_offset, c_offset;
    UInt col, colstart, colend, row_ybind;

    //block row identification
    extract_pair( _Patt->locateElBlock( row, row ), m, n );
    //loop on block columns
    for ( n = 0; n < BCols; n++ )
        if ( _Patt->block_ptr( m, n ) != 0 )
        {
            extract_pair( _Patt->blockOffset( m, n ), r_offset, c_offset );
            col = c_offset;
            //local row
            extract_pair( _Patt->localNumber( m, n, row, col ), lr, lc );
            //diagonal terms
            if ( n == m )
                _bValues[ m ][ n ][ lr ] = coeff; //the true diagonal
            else
                _bValues[ m ][ n ][ lr ] = 0.;
            //loop on the columns
            colstart = _Patt->block_ptr( m, n ) ->give_bindx() [ lr ];
            colend = _Patt->block_ptr( m, n ) ->give_bindx() [ lr + 1 ];
            disp = _Patt->nRows( m, n ) + 1;
            for ( lc = colstart; lc < colend; lc++ )
            {
                _bValues[ m ][ n ][ lc ] = 0;
                // for the columns update
                lr_ybind = _Patt->block_ptr( m, n ) ->give_bindx() [ lc ];
                lc_ybind = _Patt->block_ptr( m, n ) ->give_ybind() [ lc - disp ];
                row_ybind = lr_ybind + r_offset;
                b[ row_ybind ] -= _bValues[ m ][ n ][ lc_ybind ] * datum;
                _bValues[ m ][ n ][ lc_ybind ] = 0.;
            }
        }
    //rhs:
    b[ row ] = coeff * datum;
}


// Matrix-vector product.
template <UInt BRows, UInt BCols>
std::vector<double>
MixedMatr<BRows, BCols, MSRPatt, double>::
operator*( const std::vector<double> &v ) const
{
    ASSERT( _Patt->nCols() == v.size(), "Error in Matrix Vector product" );
    std::vector<double> ans;
    ans.resize( _Patt->nRows(), 0.0 );

    UInt nrows, istart, iend;
    UInt r_offset, c_offset;
    UInt ib, jb, ir, ii;

    for ( ib = 0; ib < BRows; ib++ )
    {

        nrows = _Patt->nRows( ib, 0 );

        for ( jb = 0; jb < BCols; jb++ )
        {

            if ( _Patt->block_ptr( ib, jb ) != 0 )
            {

                extract_pair( _Patt->blockOffset( ib, jb ), r_offset, c_offset );

                for ( ir = 0; ir < nrows; ir++ )
                {
                    //"diagonal" term
                    ans[ ir + r_offset ] += _bValues[ ib ][ jb ][ ir ] * v[ ir + c_offset ];

                    istart = _Patt->block_ptr( ib, jb ) ->give_bindx() [ ir ] - OFFSET;
                    iend = _Patt->block_ptr( ib, jb ) ->give_bindx() [ ir + 1 ] - OFFSET;
                    for ( ii = istart;ii < iend;++ii )
                    {
                        ans[ ir + r_offset ] += _bValues[ ib ][ jb ][ ii ] *
                                                v[ c_offset + _Patt->block_ptr( ib, jb ) ->give_bindx() [ ii ] - OFFSET ];
                    }
                }
            }
        }
    }
    return ans;
}
// version for type Vector
template <UInt BRows, UInt BCols>
Vector
MixedMatr<BRows, BCols, MSRPatt, double>::
operator*( const Vector &v ) const
{
    ASSERT( _Patt->nCols() == v.size(), "Error in Matrix Vector product" );
    Vector ans( ZeroVector(_Patt->nRows() ) );

    UInt nrows, istart, iend;
    UInt r_offset, c_offset;
    UInt ib, jb, ir, ii;

    for ( ib = 0; ib < BRows; ib++ )
    {

        nrows = _Patt->nRows( ib, 0 );

        for ( jb = 0; jb < BCols; jb++ )
        {

            if ( _Patt->block_ptr( ib, jb ) != 0 )
            {

                extract_pair( _Patt->blockOffset( ib, jb ), r_offset, c_offset );

                for ( ir = 0; ir < nrows; ir++ )
                {
                    //"diagonal" term
                    ans[ ir + r_offset ] += _bValues[ ib ][ jb ][ ir ] * v[ ir + c_offset ];
                    //other terms
                    istart = _Patt->block_ptr( ib, jb ) ->give_bindx() [ ir ] - OFFSET;
                    iend = _Patt->block_ptr( ib, jb ) ->give_bindx() [ ir + 1 ] - OFFSET;
                    for ( ii = istart;ii < iend;++ii )
                    {
                        ans[ ir + r_offset ] += _bValues[ ib ][ jb ][ ii ] *
                                                v[ c_offset + _Patt->block_ptr( ib, jb ) ->give_bindx() [ ii ] - OFFSET ];
                    }
                }
            }
        }
    }
    return ans;
}
// Version for C pointer vector. BEWARE: no check on bounds is done !
template <UInt BRows, UInt BCols>
void operMatVec( double * const mv,
                 const MixedMatr<BRows, BCols, MSRPatt, double> &Mat,
                 const double *v )
{
    UInt nrows, istart, iend;
    UInt r_offset, c_offset;
    UInt ib, jb, ir, ii;

    for ( ib = 0; ib < BRows; ib++ )
    {

        nrows = Mat._Patt->nRows( ib, 0 );

        for ( ir = 0; ir < nrows; ir++ )
        {

            //initialize
            mv[ ir + ib * nrows ] = 0.;

            for ( jb = 0; jb < BCols; jb++ )
            {

                if ( Mat._Patt->block_ptr( ib, jb ) != 0 )
                {

                    extract_pair( Mat._Patt->blockOffset( ib, jb ), r_offset, c_offset );

                    //"diagonal" term
                    mv[ ir + r_offset ] += Mat._bValues[ ib ][ jb ][ ir ] * v[ ir + c_offset ];
                    //other terms
                    istart = Mat._Patt->block_ptr( ib, jb ) ->give_bindx() [ ir ] - OFFSET;
                    iend = Mat._Patt->block_ptr( ib, jb ) ->give_bindx() [ ir + 1 ] - OFFSET;
                    for ( ii = istart;ii < iend;++ii )
                    {
                        mv[ ir + r_offset ] += Mat._bValues[ ib ][ jb ][ ii ] *
                                               v[ c_offset + Mat._Patt->block_ptr( ib, jb ) ->give_bindx() [ ii ] - OFFSET ];
                    }
                }
            }
        }
    }
}

//necessary for IML++ library: transpose-matrix by vector product.
template <UInt BRows, UInt BCols>
Vector
MixedMatr<BRows, BCols, MSRPatt, double>::
trans_mult( const Vector &v ) const
{
    ASSERT( _Patt->nRows() == v.size(), "Error in Matrix Vector product" );
    Vector ans( _Patt->nRows() );
    ans = 0.0;

    UInt nrows = 0, ncols = 0;
    UInt nnz = 0, nr_global = 0, nc_global;
    Container coldata, pos;

    for ( UInt ib = 0; ib < BRows; ib++ )
    {

        nrows = _Patt->nRows( ib, 0 );
        nc_global = 0;

        for ( UInt jb = 0; jb < BCols; jb++ )
        {

            ncols = _Patt->nCols( ib, jb );

            if ( _Patt->block_ptr( ib, jb ) != 0 )
            {

                for ( UInt nr = 0; nr < nrows; nr++ )
                {
                    coldata.resize( ncols );
                    pos.resize( ncols );
                    nnz = _Patt->row( ib, jb, nr, coldata.begin(), pos.begin() );

                    for ( UInt jcol = 0; jcol < nnz; jcol++ )
                        ans[ nc_global + coldata[ jcol ] - OFFSET ] += _bValues[ ib ][ jb ][ pos[ jcol ] ] *
                                v[ nr_global + nr ];
                }
            }
            nc_global += ncols;
        }
        nr_global += nrows;
    }
    return ans;
}

template <UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, MSRPatt, double>::
zeros()
{
    UInt ib, jb;
    for ( ib = 0; ib < BRows; ib++ )
        for ( jb = 0; jb < BCols; jb++ )
            if ( _Patt->block_ptr( ib, jb ) != 0 )
            {
                typename std::vector<double>::iterator start = _bValues[ ib ][ jb ].begin();
                typename std::vector<double>::iterator end = _bValues[ ib ][ jb ].end();
                fill( start, end, 0.0 );
            }
}

//-----------------------------------------------------------------------
// MixedMatr SPECIALIZATION FOR CSR
//-----------------------------------------------------------------------

//Default Constructor
template <UInt BRows, UInt BCols>
MixedMatr<BRows, BCols, CSRPatt, double>::
MixedMatr()
{
    for ( UInt i = 0; i < BRows; i++ )
        for ( UInt j = 0; j < BCols; j++ )
            _bValues[ i ][ j ].resize( 0 );
}

//Constructor from an existing external pattern.
//version for CSR format: size of values= nnz.
template <UInt BRows, UInt BCols>
MixedMatr<BRows, BCols, CSRPatt, double>::
MixedMatr( const MixedPattern<BRows, BCols, CSRPatt> &ex_pattern )
{
    _Patt = &ex_pattern;

    for ( UInt i = 0; i < BRows; i++ )
        for ( UInt j = 0; j < BCols; j++ )
            _bValues[ i ][ j ].resize( ex_pattern.nNz( i, j ) );
}

// Determines the lumped diagonal of P1 mass matrix.
template <UInt BRows, UInt BCols>
std::vector<double>
MixedMatr<BRows, BCols, CSRPatt, double>::
MassDiagP1() const
{
    UInt nrows = _Patt->nRows();
    UInt nnz = 0, nr_global = 0;
    Container coldata, position;

    std::vector<double> diag;
    diag.resize( nrows, 0.0 );

    for ( UInt ib = 0; ib < BRows; ib++ )
    {

        nrows = _Patt->nRows( ib, 0 );

        for ( UInt jb = 0; jb < BCols; jb++ )
        {

            if ( _Patt->block_ptr( ib, jb ) != 0 )
            {

                for ( UInt nr = 0; nr < nrows; nr++ )
                {
                    coldata.resize( _Patt->nCols( ib, jb ) );
                    position.resize( _Patt->nCols( ib, jb ) );
                    nnz = _Patt->row( ib, jb, nr, coldata.begin(), position.begin() );
                    for ( UInt jcol = 0; jcol < nnz; jcol++ )
                        diag[ nr_global + nr ] +=
                            _bValues[ ib ][ jb ][ position[ jcol ] ];
                }
            }
        }
        nr_global += nrows;
    }

    return diag;
}

// Inverts the diagonal matrix Diag and mulitply it by the matrix Mat.
template <UInt BR, UInt BC>
void
MultInvDiag( const std::vector<Real> &Diag,
             const MixedMatr<BR, BC, CSRPatt, Real> &Mat,
             MixedMatr<BR, BC, CSRPatt, Real> &ans )
{
    ASSERT( find( Diag.begin(), Diag.end(), 0 ) == Diag.end(), "Diagonal matrix Diag must be invertible" );

    ASSERT( Diag.size() == Mat._Patt->nRows(), "Matrices size not compatible" );

    // Product:
    UInt nrows = 0;
    UInt nnz = 0, nr_global = 0;
    Container coldata, pos;

    for ( UInt ib = 0; ib < BR; ib++ )
    {

        nrows = Mat._Patt->nRows( ib, 0 );

        for ( UInt jb = 0; jb < BC; jb++ )
        {

            if ( Mat._Patt->block_ptr( ib, jb ) != 0 )
            {

                for ( UInt nr = 0; nr < nrows; nr++ )
                {
                    coldata.resize( Mat._Patt->nCols( ib, jb ) );
                    pos.resize( Mat._Patt->nCols( ib, jb ) );
                    nnz = Mat._Patt->row( ib, jb, nr, coldata.begin(), pos.begin() );

                    for ( UInt jcol = 0; jcol < nnz; jcol++ )
                        ans._bValues[ ib ][ jb ][ pos[ jcol ] ] =
                            ( Mat._bValues[ ib ][ jb ][ pos[ jcol ] ] / Diag[ nr_global + nr ] );
                }
            }
        }
        nr_global += nrows;
    }
}

// Gives the diagonal of a block matrix.
template <UInt BRows, UInt BCols>
std::vector<double>
MixedMatr<BRows, BCols, CSRPatt, double>::
giveDiag() const
{
    ASSERT_PRE( BRows == BCols, "block matrix must be a square matrix" );

    std::vector<Real> diag;
    diag.resize( _Patt->nRows(), 0.0 );

    Container coldata, pos;
    UInt nrows = 0, nr_global = 0, nnz = 0;

    for ( UInt ib = 0; ib < BRows; ib++ )
    {
        ASSERT_PRE( _Patt->nRows( ib, ib ) == _Patt->nCols( ib, ib ),
                    "matrix must have nrows=ncols for diagonal blocks" );

        nrows = _Patt->nRows( ib, ib );

        if ( _Patt->block_ptr( ib, ib ) != 0 )
            for ( UInt nr = 0; nr < nrows; nr++ )
            {

                coldata.resize( _Patt->nCols( ib, ib ) );
                pos.resize( _Patt->nCols( ib, ib ) );

                nnz = _Patt->row( ib, ib, nr, coldata.begin(), pos.begin() );

                UInt jcol = 0;

                if ( coldata[ jcol ] - OFFSET == nr )
                    diag[ nr_global + nr ] = _bValues[ ib ][ ib ][ pos[ jcol ] ];
                else
                {
                    while ( coldata[ jcol ] - OFFSET != nr )
                        jcol++;
                    diag[ nr_global + nr ] = _bValues[ ib ][ ib ][ pos[ jcol ] ];
                }
            }
        nr_global += nrows;
    }
    return diag;
}


// Warning: the two matrices will point to the same pattern.
template <UInt BRows, UInt BCols>
MixedMatr<BRows, BCols, CSRPatt, double> &
MixedMatr<BRows, BCols, CSRPatt, double>::
operator=( const MixedMatr<BRows, BCols, CSRPatt, double> &RhMtr )
{
    if ( &RhMtr != this )
    {
        UInt ib, jb;
        for ( ib = 0; ib < BRows; ib++ )
            for ( jb = 0; jb < BCols; jb++ )
                _bValues[ ib ][ jb ] = RhMtr.bValues( ib, jb );

        _Patt = RhMtr.Patt();
    }

    return *this;

}

// Assigns the matrix to loc_val at the place of index where.
template <UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, CSRPatt, double>::
set_mat( UInt ib, UInt jb, UInt where, double loc_val )
{
    _bValues[ ib ][ jb ][ where - OFFSET ] = loc_val;
}

// Assigns the matrix element (row,col) to loc_val.
template <UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, CSRPatt, double>::
set_mat( UInt row, UInt col, double loc_val )
{
    UInt m, n, lr, lc;
    extract_pair( _Patt->locateElBlock( row, col ), m, n );
    extract_pair( _Patt->localNumber( m, n, row, col ), lr, lc );
    set_mat( m, n, lr, lc, loc_val );
}

// Assigns the matrix element (row,col) of block (ib,jb) to loc_val.
template <UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, CSRPatt, double>::
set_mat( UInt ib, UInt jb, UInt row, UInt col, double loc_val )
{
    if ( _Patt->block_ptr( ib, jb ) != 0 )
    {
        std::pair<UInt, bool> where = _Patt->block_ptr( ib, jb ) ->locate_index( row, col );
        if ( where.second )
            _bValues[ ib ][ jb ][ where.first ] = loc_val;
    }
}

// Adds loc_val to the matrix element (row,col).
template <UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, CSRPatt, double>::
set_mat_inc( UInt row, UInt col, double loc_val )
{
    UInt m, n, lr, lc;
    extract_pair( _Patt->locateElBlock( row, col ), m, n );
    extract_pair( _Patt->localNumber( m, n, row, col ), lr, lc );
    set_mat_inc( m, n, lr, lc, loc_val );
}

// Adds loc_val to the matrix element (row,col) of block (ib,jb).
template <UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, CSRPatt, double>::
set_mat_inc( UInt ib, UInt jb, UInt row, UInt col, double loc_val )
{
    if ( _Patt->block_ptr( ib, jb ) != 0 )
    {
        std::pair<UInt, bool> where = _Patt->block_ptr( ib, jb ) ->locate_index( row, col );
        if ( where.second )
            _bValues[ ib ][ jb ][ where.first ] += loc_val;
    }
}

// Returns the matrix element (i,j) value.
template <UInt BRows, UInt BCols>
double&
MixedMatr<BRows, BCols, CSRPatt, double>::
get_value( UInt i, UInt j )
{
    UInt m, n, lr, lc;
    extract_pair( _Patt->locateElBlock( i, j ), m, n );
    extract_pair( _Patt->localNumber( m, n, i, j ), lr, lc );
    return get_value( m, n, lr, lc );
}
// const qualifyer version
template <UInt BRows, UInt BCols>
const double
MixedMatr<BRows, BCols, CSRPatt, double>::
get_value( UInt i, UInt j ) const
{
    UInt m, n, lr, lc;
    extract_pair( _Patt->locateElBlock( i , j ), m, n );
    extract_pair( _Patt->localNumber( m, n, i, j ), lr, lc );
    return get_value( m, n, i, j );
}

// Returns the matrix element (i,j) value of block (ib,jb).
template <UInt BRows, UInt BCols>
double&
MixedMatr<BRows, BCols, CSRPatt, double>::
get_value( UInt ib, UInt jb, UInt i, UInt j )
{
    if ( _Patt->block_ptr( ib, jb ) != 0 )
        return _bValues[ ib ][ jb ][ _Patt->block_ptr( ib, jb ) ->locate_index( i, j ).first ];
    else
        return 0.0;
}
// const qualifyer version
template <UInt BRows, UInt BCols>
const double
MixedMatr<BRows, BCols, CSRPatt, double>::
get_value( UInt ib, UInt jb, UInt i, UInt j ) const
{
    if ( _Patt->block_ptr( ib, jb ) != 0 )
        return _bValues[ ib ][ jb ][ _Patt->block_ptr( ib, jb ) ->locate_index( i, j ).first ];
    else
        return 0.0;
}

// Matrix visualization a la matlab.
template <UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, CSRPatt, double>::
spy( std::string const &filename )
{
    // Purpose: Matlab dumping and spy
    std::string nome = filename, uti = " , ";
    //
    // check on the file name
    //
    UInt i = filename.find( "." );

    if ( i <= 0 )
        nome = filename + ".m";
    else
    {
        if ( i != filename.size() - 2 || filename[ i + 1 ] != 'm' )
        {
            std::cerr << "Wrong file name ";
            nome = filename + ".m";
        }
    }

    std::ofstream file_out( nome.c_str() );
    UInt mb, nb, j;
    UInt r_offset, c_offset, ia_start, ia_end;
    file_out << "S = [ ";
    //loop on block rows
    for ( mb = 0; mb < BRows; mb++ )
    {
        //loop on local rows
        for ( i = 0; i < _Patt->nRows( mb, 0 ); i++ )
            //loop on block columns
            for ( nb = 0; nb < BCols; nb++ )
                if ( _Patt->block_ptr( mb, nb ) != 0 )
                {
                    //row and col offsets
                    extract_pair( _Patt->blockOffset( mb, nb ), r_offset, c_offset );
                    //for CSR:
                    ia_start = _Patt->block_ptr( mb, nb ) ->give_ia() [ i ] - OFFSET;
                    ia_end = _Patt->block_ptr( mb, nb ) ->give_ia() [ i + 1 ] - OFFSET;
                    for ( j = ia_start; j < ia_end; j++ )
                    {
                        file_out << i + r_offset + 1 << uti
                        << _Patt->block_ptr( mb, nb ) ->give_ja() [ j ] + c_offset + 1 - OFFSET
                        << uti << _bValues[ mb ][ nb ][ j ] << std::endl;
                    }
                }
    }

    file_out << "];" << std::endl;
    file_out << "I=S(:,1); J=S(:,2); S=S(:,3); A=sparse(I,J,S); spy(A);" << std::endl;
}

// Assigns matrix diagonal element (r,r) to coeff, other elts
// of row r to zero.
template <UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, CSRPatt, double>::
diagonalize_row( UInt const r, double const coeff )
{
    UInt nnz = 0, m = 0, n = 0;
    Container coldata, pos;
    coldata.resize( _Patt->nCols() );
    pos.resize( _Patt->nCols() );

    nnz = _Patt->row( r - OFFSET, coldata.begin(), pos.begin() );

    UInt jcol = 0;
    if ( coldata[ jcol ] == r )
    {
        //diagonal element:
        extract_pair( _Patt->locateElBlock( r - OFFSET , coldata[ jcol ] - OFFSET ), m, n );
        _bValues[ m ][ n ][ pos[ jcol ] ] = coeff;
        //other elements:
        for ( jcol = 1; jcol < nnz; jcol++ )
        {
            extract_pair( _Patt->locateElBlock( r - OFFSET , coldata[ jcol ] - OFFSET ), m, n );
            _bValues[ m ][ n ][ pos[ jcol ] ] = 0.0;
        }
    }
    else
        for ( jcol = 0; jcol < nnz; jcol++ )
        {
            extract_pair( _Patt->locateElBlock( r - OFFSET , coldata[ jcol ] - OFFSET ), m, n );
            if ( coldata[ jcol ] == r )
                //diagonal element:
                _bValues[ m ][ n ][ pos[ jcol ] ] = coeff;
            else
                //other elements
                _bValues[ m ][ n ][ pos[ jcol ] ] = 0.0;
        }
}

// assign a row to zero. Remark, zero might be defined for any double
template <UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, CSRPatt, double>::
zero_row( UInt const row )
{
    diagonalize_row( row, 0.0 );
}

template <UInt BRows, UInt BCols>
template < typename VectorType >
void
MixedMatr<BRows, BCols, CSRPatt, double>::
diagonalize( std::vector<UInt> const rVec,
	     double const coeff,
	     VectorType &b,
	     std::vector<double> datumVec )
{
     UInt sizeVec(rVec.size());
     if ( sizeVec != datumVec.size()) 
     { //! vectors must be of the same size
         ERROR_MSG( "diagonalize: vectors must be of the same size\n" );
     }
      
     for (UInt i=0; i < sizeVec; i++)
       diagonalize( rVec[i], coeff, b, datumVec[i]);
      
}

// Assigns matrix diagonal element (r,r) to coeff, other elts
// of row r to zero, and vector b element b(r) to coeff*datum.
template <UInt BRows, UInt BCols>
template < typename VectorType >
void
MixedMatr<BRows, BCols, CSRPatt, double>::
diagonalize( UInt const row, double const coeff, VectorType &b,
             double datum )
{
    UInt nnz = 0, m = 0, n = 0;
    Container coldata, pos;
    coldata.resize( _Patt->nCols() );
    pos.resize( _Patt->nCols() );

    nnz = _Patt->row( row - OFFSET, coldata.begin(), pos.begin() );

    UInt jcol = 0;
    if ( coldata[ jcol ] == row )
    { //case diagfirst
        //diagonal element:
        extract_pair( _Patt->locateElBlock( row - OFFSET , coldata[ jcol ] - OFFSET ), m, n );
        _bValues[ m ][ n ][ pos[ jcol ] ] = coeff;
        //other elements:
        for ( jcol = 1; jcol < nnz; jcol++ )
        {
            extract_pair( _Patt->locateElBlock( row - OFFSET , coldata[ jcol ] - OFFSET ), m, n );
            _bValues[ m ][ n ][ pos[ jcol ] ] = 0.0;
        }
    }
    else
        for ( jcol = 0; jcol < nnz; jcol++ )
        {
            extract_pair( _Patt->locateElBlock( row - OFFSET , coldata[ jcol ] - OFFSET ), m, n );
            if ( coldata[ jcol ] == row )
                //diagonal element:
                _bValues[ m ][ n ][ pos[ jcol ] ] = coeff;
            else
                //other elements
                _bValues[ m ][ n ][ pos[ jcol ] ] = 0.0;
        }

    //rhs:
    b[ row - OFFSET ] = coeff * datum;
}

// set to zero the row trD and the corresponding column=row for D
// passing the known datum to the rhs b.
template <UInt BR, UInt BC, typename VectorType, typename DataType>
void zero_row_col( UInt const row, MixedMatr<BR, BC, CSRPatt, Real> &trD,
                   MixedMatr<BC, BR, CSRPatt, Real> &D, VectorType &bp,
                   DataType const datum )
{
    // for trD
    trD.zero_row( row );

    // for D
    UInt lr, lc, row_loc;
    UInt m, n;
    UInt r_offset, c_offset;
    UInt col, start_ia, end_ia;

    //block row identification
    extract_pair( trD._Patt->locateElBlock( row, 0 ), m, n );
    //loop on block columns
    for ( n = 0; n < BC; n++ )
        if ( trD._Patt->block_ptr( m, n ) != 0 )
        {
            // pattern jaT must have been defined !
            ASSERT( trD._Patt->block_ptr( m, n ) ->give_jaT().size() > 0, "jaT must be built before, cf. buildPattTpatt() function" );
            extract_pair( trD._Patt->blockOffset( m, n ), r_offset, c_offset );
            col = c_offset;
            //local row
            extract_pair( trD._Patt->localNumber( m, n, row, col ), lr, lc );
            //loop on the columns of trD involved
            start_ia = trD._Patt->block_ptr( m, n ) ->give_ia() [ lr ];
            end_ia = trD._Patt->block_ptr( m, n ) ->give_ia() [ lr + 1 ];
            for ( lc = start_ia; lc < end_ia; lc++ )
            {
                // columns of trD become rows of D
                row_loc = trD._Patt->block_ptr( m, n ) ->give_ja() [ lc ];
                bp[ row_loc + c_offset ] -= datum *
                                            D._bValues[ n ][ m ][ trD._Patt->block_ptr( m, n ) ->give_jaT() [ lc ] ];
                D._bValues[ n ][ m ][ trD._Patt->block_ptr( m, n ) ->give_jaT() [ lc ] ] = 0.0;
            }
        }
}

// Matrix-vector product.
template <UInt BRows, UInt BCols>
std::vector<double>
MixedMatr<BRows, BCols, CSRPatt, double>::
operator*( const std::vector<double> &v ) const
{
    ASSERT( _Patt->nCols() == v.size(), "Error in Matrix Vector product" );
    std::vector<double> ans;
    ans.resize( _Patt->nRows(), 0.0 );

    UInt nrows, iastart, iaend;
    UInt r_offset, c_offset;
    UInt ib, jb, ir, ii;

    for ( ib = 0; ib < BRows; ib++ )
    {

        nrows = _Patt->nRows( ib, 0 );

        for ( jb = 0; jb < BCols; jb++ )
        {

            if ( _Patt->block_ptr( ib, jb ) != 0 )
            {

                extract_pair( _Patt->blockOffset( ib, jb ), r_offset, c_offset );

                for ( ir = 0; ir < nrows; ir++ )
                {
                    iastart = _Patt->block_ptr( ib, jb ) ->give_ia() [ ir ] - OFFSET;
                    iaend = _Patt->block_ptr( ib, jb ) ->give_ia() [ ir + 1 ] - OFFSET;
                    for ( ii = iastart;ii < iaend;++ii )
                    {
                        ans[ ir + r_offset ] += _bValues[ ib ][ jb ][ ii ] *
                                                v[ c_offset + _Patt->block_ptr( ib, jb ) ->give_ja() [ ii ] - OFFSET ];
                    }
                }
            }
        }
    }
    return ans;
}
// version for type Vector
template <UInt BRows, UInt BCols>
Vector
MixedMatr<BRows, BCols, CSRPatt, double>::
operator*( const Vector &v ) const
{
    ASSERT( _Patt->nCols() == v.size(), "Error in Matrix Vector product" );
    Vector ans( ZeroVector( _Patt->nRows() ) );

    UInt nrows, iastart, iaend;
    UInt r_offset, c_offset;
    UInt ib, jb, ir, ii;

    for ( ib = 0; ib < BRows; ib++ )
    {

        nrows = _Patt->nRows( ib, 0 );

        for ( jb = 0; jb < BCols; jb++ )
        {

            if ( _Patt->block_ptr( ib, jb ) != 0 )
            {

                extract_pair( _Patt->blockOffset( ib, jb ), r_offset, c_offset );

                for ( ir = 0; ir < nrows; ir++ )
                {
                    iastart = _Patt->block_ptr( ib, jb ) ->give_ia() [ ir ] - OFFSET;
                    iaend = _Patt->block_ptr( ib, jb ) ->give_ia() [ ir + 1 ] - OFFSET;
                    for ( ii = iastart;ii < iaend;++ii )
                    {
                        ans[ ir + r_offset ] += _bValues[ ib ][ jb ][ ii ] *
                                                v[ c_offset + _Patt->block_ptr( ib, jb ) ->give_ja() [ ii ] - OFFSET ];
                    }
                }
            }
        }
    }
    return ans;
}
// version C pointer vector. BEWARE: no check on bounds is done !
template <UInt BRows, UInt BCols>
void operMatVec( double * const mv,
                 const MixedMatr<BRows, BCols, CSRPatt, double> &Mat,
                 const double *v )
{
    UInt nrows, iastart, iaend;
    UInt r_offset, c_offset;
    UInt ib, jb, ir, ii;

    for ( ib = 0; ib < BRows; ib++ )
    {

        nrows = Mat._Patt->nRows( ib, 0 );

        for ( ir = 0; ir < nrows; ir++ )
        {

            //initialize
            mv[ ir + ib * nrows ] = 0.;

            for ( jb = 0; jb < BCols; jb++ )
            {

                if ( Mat._Patt->block_ptr( ib, jb ) != 0 )
                {

                    extract_pair( Mat._Patt->blockOffset( ib, jb ), r_offset, c_offset );

                    iastart = Mat._Patt->block_ptr( ib, jb ) ->give_ia() [ ir ] - OFFSET;
                    iaend = Mat._Patt->block_ptr( ib, jb ) ->give_ia() [ ir + 1 ] - OFFSET;

                    for ( ii = iastart;ii < iaend;++ii )
                    {
                        mv[ ir + r_offset ] += Mat._bValues[ ib ][ jb ][ ii ] *
                                               v[ c_offset + Mat._Patt->block_ptr( ib, jb ) ->give_ja() [ ii ] - OFFSET ];
                    }
                }
            }
        }
    }
}

//necessary for IML++ library: transpose-matrix by vector product.
template <UInt BRows, UInt BCols>
Vector
MixedMatr<BRows, BCols, CSRPatt, double>::
trans_mult( const Vector &v ) const
{
    ASSERT( _Patt->nRows() == v.size(), "Error in Matrix Vector product" );
    Vector ans( _Patt->nRows() );
    ans = 0.0;

    UInt nrows = 0, ncols = 0;
    UInt nnz = 0, nr_global = 0, nc_global;
    Container coldata, pos;

    for ( UInt ib = 0; ib < BRows; ib++ )
    {

        nrows = _Patt->nRows( ib, 0 );
        nc_global = 0;

        for ( UInt jb = 0; jb < BCols; jb++ )
        {

            ncols = _Patt->nCols( ib, jb );

            if ( _Patt->block_ptr( ib, jb ) != 0 )
            {

                for ( UInt nr = 0; nr < nrows; nr++ )
                {
                    coldata.resize( ncols );
                    pos.resize( ncols );
                    nnz = _Patt->row( ib, jb, nr, coldata.begin(), pos.begin() );

                    for ( UInt jcol = 0; jcol < nnz; jcol++ )
                        ans[ nc_global + coldata[ jcol ] - OFFSET ] += _bValues[ ib ][ jb ][ pos[ jcol ] ] *
                                v[ nr_global + nr ];
                }
            }
            nc_global += ncols;
        }
        nr_global += nrows;
    }
    return ans;
}

template <UInt BRows, UInt BCols>
void
MixedMatr<BRows, BCols, CSRPatt, double>::
zeros()
{
    UInt ib, jb;
    for ( ib = 0; ib < BRows; ib++ )
        for ( jb = 0; jb < BCols; jb++ )
            if ( _Patt->block_ptr( ib, jb ) != 0 )
            {
                typename std::vector<double>::iterator start = _bValues[ ib ][ jb ].begin();
                typename std::vector<double>::iterator end = _bValues[ ib ][ jb ].end();
                fill( start, end, 0.0 );
            }
}

}
