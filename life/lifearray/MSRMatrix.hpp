/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2004-10-26

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
   \file MSRMatrix.hpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-10-26
 */

#ifndef _MSRMATRIX_HPP_
#define _MSRMATRIX_HPP_

namespace LifeV
{

////////////////////////////////////////////////////////////////
//
// MSR Format
//
///////////////////////////////////////////////////////////////

template <typename DataType>
class MSRMatr
{
public:
    //! default constructor. Caution: does not have a pattern.
    MSRMatr();

    // All other constructors MUST be based on an existing pattern

    //! constructor from MSR pattern
    MSRMatr( const MSRPatt &ex_pattern );

    /*! constructor from CSR pattern
     *  @author Alain Gauthier
     */
    MSRMatr( const CSRPatt &ex_pattern );

    /*! constructor with given value array
     *  @param ex_pattern the MSR pattern
     *  @param ex_value the value array
     */
    MSRMatr( const MSRPatt* ex_pattern, const std::vector<DataType> &ex_value );

    //! copy constructor
    MSRMatr( const MSRMatr<DataType> &RightHandMSR );

    /*! conversion copy constructor.
     *  requires an already converted pattern
     *  @param ex_pattern the converted pattern
     *  @param RightHandCSR the CSRMatrix to be copied
     */
    MSRMatr( const MSRPatt &ex_pattern,
             const CSRMatr<CSRPatt, DataType> &RightHandCSR );

    //! returns a const * to the pattern
    const MSRPatt * Patt() const
    {
        return _Patt;
    };

    //! returns the value vector read only
    const std::vector<DataType> & value() const
    {
        return _value;
    };

    //! returns the value vector writable
    std::vector<DataType> & value()
    {
        return _value;
    };

    //! give the value vector writable in raw format (suitable for C)
    DataType * giveRaw_value()
    {
        return & ( _value.front() );
    }

    //! give the value vector read only in raw format (suitable for C)
    DataType const * giveRaw_value() const
    {
        return & ( _value.front() );
    }

    //! determine the lumped matrix for P1
    std::vector<DataType> MassDiagP1() const;

    //! build the inverse of a diagonal matrix and multiply it by another matrix
    friend void MultInvDiag( const std::vector<Real> &Diag,
                             const MSRMatr<Real> &Mat, MSRMatr<Real> &ans );

    //! give the diagonal of an MSR matrix
    std::vector<DataType> giveDiag() const;

    //! assignment operator.
    //! Warning: the two matrices will point to the same pattern
    MSRMatr & operator= ( const MSRMatr<DataType> &RhMsr );

    //! assignment operator for MSR <- CSR.
    //! the two matrices need to have the same dimension and sparsity pattern
    MSRMatr & operator= ( const CSRMatr<CSRPatt, DataType> &RhCsr );

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
    DataType& get_value( UInt i, UInt j )
    {
        return _value[ _Patt->locate_index( i, j ).first ];
    };

    //! return matrix entry (i,j)
    const DataType& get_value( UInt i, UInt j ) const
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

    //! set entry (r,r) to coeff and rest of row r to zero
    void diagonalize_row ( UInt const r, DataType const coeff );

    /*! apply constraint on row r
     *  @param r row number
     *  @param coeff value to set entry (r,r) at
     *  @param b right hand side vector of the system to be adapted accordingly
     *  @param datum value to constrain entry r of the solution at
     */
    void diagonalize ( UInt const r, DataType const coeff,
                       std::vector<DataType> &b, DataType datum );
    /*! apply constraint on row r
     *  @param r row number
     *  @param coeff value to set entry (r,r) at
     *  @param b right hand side Vector of the system to be adapted accordingly
     *  @param datum value to constrain entry r of the solution at
     */
    void diagonalize( UInt const r, DataType const coeff, Vector &b,
                      DataType datum );

    //! matrix vector product
    std::vector<DataType> operator* ( const std::vector<DataType> &v ) const;

    //! Matrix-vector product for the class Vector (useful for IML++)
    Vector operator*( const Vector &v ) const;

    //! Version for C pointer vector. BEWARE: no check on bounds done !
    template <typename DataT>
    friend void operMatVec( DataT * const mv,
                            const MSRMatr<DataT> &Mat,
                            const DataT *v );

    //! transpose-matrix by vector product. necessary for IML++ library
    Vector trans_mult( const Vector &v ) const;

    //! scaling operator
    MSRMatr operator*( const DataType num );

    //! in place scaling operator
    MSRMatr & operator*=( const DataType num );

    /*! set this matrix to numM*M + numA*A of two matrices scaling plus addition
     *  @return numM*M + numA*A
     */
    MSRMatr & flop( const DataType numM, MSRMatr<DataType> & M,
                    const DataType numA,MSRMatr<DataType> & A );

    /*! set this matrix to num*M + A of two matrices scaling plus addition
     *  @return num*M + A
     */
    MSRMatr & flop( const DataType num, MSRMatr<DataType> & M,
                    MSRMatr<DataType> & A );

    /*! add num*M to this matrix
     *  @return new value of this matrix
     */
    MSRMatr & flop( const DataType num, MSRMatr<DataType> & M );

    //    MSRMatr & operator* (const DataType num,const MSRMatr<DataType> &M); //Real*Matrix product

    //! set the matrix to zero
    void zeros();

private:
    std::vector<DataType> _value;
    const MSRPatt *_Patt; // I want to link the values to a pattern, NOT to change the pattern itself  (which is const)
    //     static const DataType _DefaultValue = 0;
}; // class MSRMatr

//-------------------------------------------------------------------------------------------------------
// MSR - VALUES
//-------------------------------------------------------------------------------------------------------

template <class DataType>
MSRMatr<DataType>::MSRMatr()
        :
        _Patt( 0 )
{
    // nothing to do here
}

template <typename DataType>
MSRMatr<DataType>::
MSRMatr( const MSRPatt &ex_pattern )
{
    _Patt = &ex_pattern;
    _value.resize( ex_pattern.nNz() + 1 );
}

// constructor from CSR Pattern
template <typename DataType>
MSRMatr<DataType>::
MSRMatr( const CSRPatt &ex_pattern )
{
    _Patt = &ex_pattern;
    _value.resize( ex_pattern.nNz() + 1 );
}

/* template<class DataType> */
/* CSRMatr<DataType>:: */
/* CSRMatr(UInt ex_nnz, UInt ex_nrows, UInt ex_ncols, const std::vector<UInt> &ex_ia, const std::vector<UInt> &ex_ja, const std::vector<DataType> &ex_value): */
/* _value(ex_nnz) */
/*  { */
/*   _Patt = 0; */
/*   Csr_Lv_P Pattern(UInt ex_nnz, UInt ex_nrows, UInt ex_ncols, const std::vector<UInt> &ex_ia, const std::vector<UInt> &ex_ja); */
/*   _nnz = ex_nnz; */
/*   _nrows = ex_nrows; */
/*   _ncols = ex_ncols; */
/*   //  Pattern.ShowMe(); */
/*   _Patt = &Pattern; */
/*   assert( ex_nnz == ex_value.size()); */
/*   _value = ex_value;    */
/*  }; */


template <class DataType>
MSRMatr<DataType>::
MSRMatr( const MSRPatt* ex_pattern, const std::vector<DataType> &ex_value ) :
        _Patt( ex_pattern ), _value( ex_value )
{
    //  _nnz = ex_pattern.nnz();
    //  _nrows = ex_pattern.nrow();
    //  _ncols = ex_pattern.ncol();
    ASSERT( _Patt->nNz() == ex_value.size() - 1, "Error in MSR Values " ); // in MSR value has lenghth nnz+1
}

template <class DataType>
MSRMatr<DataType>::
MSRMatr( const MSRMatr<DataType> &RightHandMSR ) :
        _value( RightHandMSR.value() ), _Patt( RightHandMSR.Patt() )
{
    // nothing to do here
}

template <class DataType>
MSRMatr<DataType>&
MSRMatr<DataType>::operator= ( const MSRMatr<DataType> &RhMsr )
{
    if ( &RhMsr != this )
    {
        _value = RhMsr.value();
        _Patt = RhMsr.Patt();
    }

    return *this;
}

//! convert CSR matrix to MSR matrix. provided for backward compatibility only
template <class DataType>
void CSRmat2MSRmat( MSRMatr<DataType> &MSRmat,
                    CSRMatr<CSRPatt, DataType> const &CSRmat )
{
    MSRmat = CSRmat;
}

template <typename DataType>
MSRMatr<DataType>::MSRMatr( const MSRPatt &ex_pattern,
                            const CSRMatr<CSRPatt, DataType> &RightHandCSR )
{
    _Patt = &ex_pattern;
    _value.resize( ex_pattern.nNz() + 1 );
    *this = RightHandCSR;
}

template <typename DataType>
MSRMatr<DataType>& MSRMatr<DataType>::
operator= ( const CSRMatr<CSRPatt, DataType> &CSRmat )
{
    Container::const_iterator bindx = Patt() ->give_bindx().begin();
    Container::const_iterator ia = CSRmat.Patt() ->give_ia().begin();
    Container::const_iterator ja = CSRmat.Patt() ->give_ja().begin();
    std::vector<double>::const_iterator valueCSR = CSRmat.value().begin();

    // copy of CSR_value into MSR_value
    UInt nrows = _Patt->nRows();

    for ( UInt i = 0; i < nrows; ++i )
    {
        UInt j_msr = *( bindx + i - OFFSET );
        UInt ifirst_csr = *( ia + i - OFFSET );
        UInt ilast_csr = *( ia + i + 1 - OFFSET );

        for ( UInt i_csr = ifirst_csr; i_csr < ilast_csr; i_csr++ )
        {
            UInt j_csr = *( ja + i_csr );
            if ( i == j_csr )
            {
                _value[ i ] = *( valueCSR + i_csr );
            }
            else
            {
                _value[ j_msr++ ] = *( valueCSR + i_csr );
            }
        }
    }
    return *this;
}

template <class DataType>
std::vector<DataType>
MSRMatr<DataType>::operator* ( const std::vector<DataType> &v ) const
{
    UInt nrows = _Patt->nRows();
    ASSERT( nrows == v.size(), "Error in Matrix Vector product" );
    std::vector<DataType> ans;
    ans.resize( nrows, 0.0 );
    for ( UInt i = 0 + OFFSET;i < nrows + OFFSET;++i )
    {
        ans[ i ] = _value[ i ] * v[ i ];
        for ( UInt j = _Patt->give_bindx() [ i ];j < _Patt->give_bindx() [ i + 1 ];++j )
            ans[ i ] += _value[ j ] * v[ _Patt->give_bindx() [ j ] ];
    }
    return ans;
}

//Matrix-vector product for the class Vector (useful for IML++)
template <class DataType>
Vector
MSRMatr<DataType>::
operator*( const Vector &v ) const
{
    UInt nrows = _Patt->nRows();
    ASSERT( nrows == v.size(), "Error in Matrix Vector product" );
    Vector ans( nrows );
    ans = ZeroVector( nrows );

    for ( UInt i = 0 + OFFSET;i < nrows + OFFSET;++i )
    {
        ans( i ) = _value[ i ] * v( i );
        for ( UInt j = _Patt->give_bindx() [ i ];j < _Patt->give_bindx() [ i + 1 ];++j )
            ans( i ) += _value[ j ] * v( _Patt->give_bindx() [ j ] );
    }
    return ans;
}

// Version for C pointer vector. BEWARE: no check on bounds is done !
template <class DataType>
void operMatVec( DataType * const mv,
                 const MSRMatr<DataType> &Mat,
                 const DataType *v )
{
    UInt nrows = Mat._Patt->nRows();

    for ( UInt i = 0 + OFFSET;i < nrows + OFFSET;++i )
    {
        mv[ i ] = Mat._value[ i ] * v[ i ];
        for ( UInt j = Mat._Patt->give_bindx() [ i ];j < Mat._Patt->give_bindx() [ i + 1 ];++j )
            mv[ i ] += Mat._value[ j ] * v[ Mat._Patt->give_bindx() [ j ] ];
    }
}

//necessary for IML++ library: transpose-matrix by vector product
template <class DataType>
Vector
MSRMatr<DataType>::
trans_mult( const Vector &v ) const
{
    UInt nrows = _Patt->nRows();
    ASSERT( nrows == v.size(), "Error in Matrix Vector product" );
    Vector ans( nrows );
    ans = 0.;

    for ( UInt i = 0 + OFFSET;i < nrows + OFFSET;++i )
    {
        ans( i ) = _value[ i ] * v( i );
        for ( UInt j = _Patt->give_bindx() [ i ];j < _Patt->give_bindx() [ i + 1 ];++j )
            ans( _Patt->give_bindx() [ j ] ) += _value[ j ] * v( i );
    }
    return ans;
}

template <class DataType>
MSRMatr<DataType>&
MSRMatr<DataType>::operator*=( const DataType num )
{
    UInt stop = _Patt->nNz() + 1;
    for ( UInt i = 0;i < stop;++i )
        _value[ i ] *= num;
    return *this;
}

template <class DataType>
MSRMatr<DataType>
MSRMatr<DataType>::operator*( const DataType num )
{
    UInt stop = _Patt->nNz()+1;
    MSRMatr<DataType> ans( *this );

    for ( UInt i = 0;i < stop;++i )
        ans.set_mat( i, num * _value[ i ] );

    return ans;
}
template <class DataType>
MSRMatr<DataType>&
MSRMatr<DataType>::flop( const DataType numM, MSRMatr<DataType> & M,
                    const DataType numA,MSRMatr<DataType> & A )
{
    UInt stop = _Patt->nNz()+1;
    //   ASSERT(M.Patt()==A.Patt,"Error in summing matrices");

    for ( UInt i = 0;i < stop;++i )
        _value[ i ] = numM * M.value() [ i ] + numA*A.value() [ i ];

    return *this;
}
template <class DataType>
MSRMatr<DataType>&
MSRMatr<DataType>::flop( const DataType num, MSRMatr<DataType>& M, MSRMatr<DataType>& A )
{
    /* AIM: Matrix = num*M + A */
    UInt stop = _Patt->nNz()+1;
    //   ASSERT(M.Patt()==A.Patt,"Error in summing matrices");

    for ( UInt i = 0;i < stop;++i )
        _value[ i ] = num * M.value() [ i ] + A.value() [ i ];

    return *this;
}

template <class DataType>
MSRMatr<DataType>&
MSRMatr<DataType>::flop( const DataType num, MSRMatr<DataType>& M )
{
    /* AIM: Matrix = num*M + Matrix */
    UInt stop = _Patt->nNz()+1;
    //   ASSERT(M.Patt()==A.Patt,"Error in summing matrices");

    for ( UInt i = 0;i < stop;++i )
        _value[ i ] += num * M.value() [ i ];

    return *this;
}

/*
template<class DataType>
MSRMatr<DataType>&
MSRMatr<DataType>::operator* (const DataType num, MSRMatr<DataType>& M)
  {
    UInt stop=_Patt->nNz()+1;

    for (UInt i=0;i<stop;++i)
      _value[i]=num*M.value()[i];

   return *this;
  }

      */

template <class DataType>
void MSRMatr<DataType>::ShowMe()
{
    Container::iterator found;
    int i_first, i_last;
    std::vector<UInt> vec_temp = _Patt->bindx();
    UInt _nrows = _Patt->nRows();
    UInt _ncols = _Patt->nCols();
    UInt _nnz = _Patt->nNz();

    std::string pare = "[";
    std::cout << "**************************" << std::endl;
    std::cout << "     MSR Matrix           " << std::endl;
    std::cout << std::endl;
    std::cout << pare;
    for ( UInt i_index = 0; i_index < _nrows; ++i_index )
    {
        std::cout << pare;
        pare = " [";
        i_first = _Patt->bindx() [ i_index ];
        i_last = _Patt->bindx() [ i_index + 1 ];
        //      std::cout << i_first << " " << i_last << std::endl;
        for ( int j = 0;j < _ncols;++j )
        {
            if ( j == i_index )
                std::cout << " " << _value[ i_index ] << " ";
            else
            {
                // I am not supposing any given ordering in _ja
                found = std::find( &vec_temp[ i_first ], &vec_temp[ i_last ], j );
                if ( found == ( Container::iterator ) & vec_temp[ i_last ] )
                    std::cout << " 0 ";
                else
                {
                    UInt j_index = found - ( &vec_temp.front() );
                    std::cout << " " << _value[ j_index ] << " ";
                }
            }
        }
        if ( i_index == _nrows - 1 )
            std::cout << " ]] " << std::endl;
        else
            std::cout << " ]  " << std::endl;
    }
    std::cout << "nnz = " << _nnz << ", nrow = " << _nrows << ", ncol = " << _ncols << std::endl;
    return ;
}

template <typename DataType>
void
MSRMatr<DataType>::spy( std::string const &filename )
{
    // Purpose: Matlab dumping and spy
    std::string nome = filename, uti = " , ";
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

    file_out << "S = [ ";
    for ( UInt i = 1;i <= _Patt->nRows();++i )
    {
        file_out << i << uti << i << uti << _value[ i - 1 ] << std::endl;
        for ( UInt ii = _Patt->give_bindx() [ i - 1 ];ii < _Patt->give_bindx() [ i ];++ii )
            file_out << i << uti << _Patt->give_bindx() [ ii ] + 1 - OFFSET <<
            uti << _value[ ii ] << std::endl;
    }

    file_out << "];" << std::endl;

    file_out << "I=S(:,1); J=S(:,2); S=S(:,3); A=sparse(I,J,S); spy(A);" << std::endl;
}

template <typename DataType>
void
MSRMatr<DataType>::set_mat_inc( UInt row, UInt col, DataType loc_val )
{

    std::pair<UInt, bool> where = _Patt->locate_index( row, col );
    if ( where.second )
        _value[ where.first ] += loc_val;
    else
    {
        std::cout << row + 1 << "," << col + 1 << std::endl;
        ERROR_MSG( "problem in MSR::set_mat_inc" );
    }
    return ;
}

template <typename DataType>
void
MSRMatr<DataType>::set_mat( UInt row, UInt col, DataType loc_val )
{

    std::pair<UInt, bool> where = _Patt->locate_index( row, col );
    if ( where.second )
        _value[ where.first ] = loc_val;

    return ;
}

template <typename DataType>
void
MSRMatr<DataType>::set_mat( UInt where, DataType loc_val )
{
    _value[ where - OFFSET ] = loc_val;

    return ;
}

// determine the lumped matrix for P1
template <typename DataType>
std::vector<DataType>
MSRMatr<DataType>::MassDiagP1() const
{
    UInt nrows = _Patt->nRows();

    std::vector<DataType> diag( nrows );

    for ( UInt nrow = 0; nrow < nrows; ++nrow )
    {
        diag[ nrow ] = _value[ nrow ];
        for ( UInt ii = _Patt->bindx() [ nrow ]; ii < _Patt->bindx() [ nrow + 1 ]; ++ii )
            diag[ nrow ] += _value[ ii ];
    }

    return diag;
}

// build the inverse of a diagonal matrix and multiply it by another matrix
void MultInvDiag( const std::vector<Real> &Diag,
                  const MSRMatr<Real> &Mat,
                  MSRMatr<Real> &ans ) ;

// give the diagonal of an MSR matrix
template <typename DataType>
std::vector<DataType>
MSRMatr<DataType>::giveDiag() const
{
    UInt nrows = _Patt->nRows();

    std::vector<DataType> diag( nrows );

    for ( UInt nrow = 0; nrow < nrows; ++nrow )
        diag[ nrow ] = _value[ nrow ];

    return diag;
}

template <typename DataType>
void
MSRMatr<DataType>::diagonalize_row( UInt const r, DataType const coeff )
{
    _value[ r - OFFSET ] = coeff;

    UInt i;
    UInt start = *( _Patt->give_bindx().begin() + r - OFFSET );
    UInt end = *( _Patt->give_bindx().begin() + r + 1 - OFFSET );
    //  UInt disp = _Patt->nRows()+1;

    for ( i = start;i < end;++i )
    {
        _value[ i - OFFSET ] = 0;
        // Miguel: Why this ?.
        //_value[_Patt->give_ybind()[i-disp]] = 0;
    }

    return ;
}

template <typename DataType>
void
MSRMatr<DataType>::diagonalize( UInt const r, DataType const coeff, std::vector<DataType> &b,
                                DataType datum )
{
    // AIM: Diagonalization of a row of the system, by setting:
    // A(i,i) = coeff,
    // A(i,j) = 0,  A(j,i) = 0 for j!=i
    // and suitably correcting the right hand side of the system
    _value[ r - OFFSET ] = coeff;

    UInt istart = *( _Patt->give_bindx().begin() + r - OFFSET );
    UInt iend = *( _Patt->give_bindx().begin() + r + 1 - OFFSET );

    typename std::vector<DataType>::iterator start = _value.begin() + istart;
    typename std::vector<DataType>::iterator end = _value.begin() + iend;
    UInt disp = _Patt->nRows() + 1;
    UInt row, col;

    transform( start, end, start, zero );

    for ( UInt i = istart;i < iend;++i )
    {
        row = _Patt->give_bindx() [ i ] - OFFSET;
        col = _Patt->give_ybind() [ i - disp ] - OFFSET;
        b[ row ] -= _value[ col ] * datum;
        _value[ col ] = 0.;
    }

    b[ r - OFFSET ] = coeff * datum;

    //Remark: in processing a list of Dirichlet nodes, there is no need to check a posteriori the right hand side
    return ;
}

//version for type Vector
template <typename DataType>
void
MSRMatr<DataType>::diagonalize( UInt const r, DataType const coeff, Vector &b, DataType datum )
{
    // AIM: Diagonalization of a row of the system, by setting:
    // A(i,i) = coeff,
    // A(i,j) = 0,  A(j,i) = 0 for j!=i
    // and suitably correcting the right hand side of the system
    _value[ r - OFFSET ] = coeff;

    UInt istart = *( _Patt->give_bindx().begin() + r - OFFSET );
    UInt iend = *( _Patt->give_bindx().begin() + r + 1 - OFFSET );

    typename std::vector<DataType>::iterator start = _value.begin() + istart;
    typename std::vector<DataType>::iterator end = _value.begin() + iend;

    UInt row, col;

    transform( start, end, start, zero );


    // Miguel: There is a buh using ybind. Alex, did you fix it?.
    // This code works without ybind.
    //
    for ( UInt i = istart;i < iend;++i )
    {
        row = _Patt->give_bindx() [ i ] - OFFSET;
        UInt Rstart = *( _Patt->give_bindx().begin() + row - OFFSET );
        UInt Rend = *( _Patt->give_bindx().begin() + row + 1 - OFFSET );
        for ( UInt j = Rstart;j < Rend;++j )
        {
            if ( ( _Patt->give_bindx() [ j ] - OFFSET ) == r )
            {
                col = j;
                b[ row - OFFSET ] -= _value[ col ] * datum;
                _value[ col ] = 0.;
                break;
            }
        }
    }


    b[ r - OFFSET ] = coeff * datum;

    //Remark: in processing a list of Dirichlet nodes, there is no need to check a posteriori the right hand side
    return ;
}


template <typename DataType>
void
MSRMatr<DataType>::zeros()
{
    typename std::vector<DataType>::iterator start = _value.begin();
    typename std::vector<DataType>::iterator end = _value.end();
    fill( start, end, 0.0 );
}

} // namespace LifeV

#endif /* _MSRMATRIX_HPP_ */
