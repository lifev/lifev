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
   \file VBRMatrix.hpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-10-26
 */

namespace LifeV
{
////////////////////////////////////////////////////////////////
//
// VBR Format based on CSR pattern for SQUARE blocks !!!
//
///////////////////////////////////////////////////////////////

template <typename DataType>
class VBRMatr
{
public:
    VBRMatr() : _Patt( 0 )
    {}
    ; // default constructor : NULL pattern
    //
    // Note that the constructors MUST be based on an existing pattern
    //
    VBRMatr( const VBRPatt &ex_pattern );
    VBRMatr( const VBRPatt& ex_pattern, const std::vector<DataType> &ex_value );
    VBRMatr( const VBRMatr<DataType> &RightHandVBR ) :
            _Patt( RightHandVBR.Patt() ), _value( RightHandVBR.value() )
    {}

    const VBRPatt *Patt() const
    {
        return _Patt;
    }
    std::vector<DataType> value() const
    {
        return _value;
    }
    DataType* giveRaw_value()
    {
        return & ( _value.front() );
    } // give the
    // value vector in Raw format (suitable for C)

    VBRMatr& operator=( const VBRMatr<DataType> &RhVbr ); // Warning: the
    // two matrices will point to the same pattern

    std::vector<DataType> operator*( const std::vector<DataType> &v ) const; //Matrix-vector product
    //Matrix vector product for IML++ (uses the Vector class)
    Vector operator*( const Vector &v ) const;
    //necessary for IML++ library: transpose-matrix by vector product
    Vector trans_mult( const Vector &v ) const;

    VBRMatr operator*( const DataType num );
    VBRMatr & operator*=( const DataType num ); //Real*Matrix product
    VBRMatr & flop( const DataType num, VBRMatr<DataType> & M, VBRMatr<DataType> & A ); //Matrix=num*M+A
    VBRMatr & flop( const DataType num, VBRMatr<DataType> & M ); //Matrix=num*M+Matrix
    void set_mat( UInt where, DataType loc_val );
    void set_mat( UInt row, UInt col, DataType loc_val );
    void set_mat_inc( UInt row, UInt col, DataType loc_val );
    // version for block operation
    void set_mat( UInt row, UInt col, std::vector<DataType> &loc_block );
    void set_mat_inc( UInt row, UInt col, std::vector<DataType> &loc_block );

    DataType& get_value( UInt i, UInt j )
    {
        return _value[ _Patt->locate_index( i, j ).first ];
    };
    const DataType& get_value( UInt i, UInt j ) const
    {
        return _value[ _Patt->locate_index( i, j ).first ];
    };
    //  void ShowMe() const;

    //! write matrix in sparse matlab format and spy
    /*! just run the resulting m-file and the matrix is loaded into A
     *  and its sparsity pattern is displayed.
     *  @param filename name of the m-file
     */
    void spy( std::string const &filename );
    //  void diagonalize_row ( UInt const r, DataType const coeff);
    //  void diagonalize ( UInt const r, DataType const coeff, std::vector<DataType> &b, DataType datum);

private:
    std::vector<DataType> _value;
    const VBRPatt *_Patt; // I want to link the values to a pattern, NOT
    // changing the pattern itself  (which is const)
    //  static const DataType _DefaultValue = 0;
};
//---------------------------------------------------------------------
// VBR - VALUES
//---------------------------------------------------------------------

template <typename DataType>
VBRMatr<DataType>::
VBRMatr( const VBRPatt &ex_pattern )
{
    _Patt = &ex_pattern;

    UInt blockSize_r = ex_pattern.rpntr() [ 1 ] - ex_pattern.rpntr() [ 0 ];
    UInt blockSize_c = blockSize_r; // we work with square blocks !!!

    _value.reserve( ex_pattern.nNz() * blockSize_r * blockSize_c );
}

template <class DataType>
VBRMatr<DataType>::
VBRMatr( const VBRPatt & ex_pattern, const std::vector<DataType> &ex_value ) :
        _value( ex_value )
{
    _Patt = &ex_pattern;
    UInt blockSize_r = _Patt->rpntr() [ 1 ] - _Patt->rpntr() [ 0 ];
    UInt blockSize_c = blockSize_r; // we work with square blocks !!!

    ASSERT( ( _Patt->nNz() * blockSize_r * blockSize_c ) == ex_value.size(), "Error in VBR Values " ); // VBR value has lenghth nnz*blockSize_r*blockSize_c
}

template <class DataType>
VBRMatr<DataType>&
VBRMatr<DataType>::operator=( const VBRMatr<DataType> &RhVbr )
{
    if ( &RhVbr != this )
    {
        _value = RhVbr.value();
        _Patt = RhVbr.Patt();
    }
    return *this;
}

template <class DataType>
std::vector<DataType>
VBRMatr<DataType>::operator*( const std::vector<DataType> &v ) const
{
    UInt blockSize_r = _Patt->rpntr() [ 1 ] - _Patt->rpntr() [ 0 ];
    //...for square matrices...
    UInt nrows = _Patt->nRows() * blockSize_r;
    UInt ncols = _Patt->nCols() * blockSize_r;
    ASSERT( ncols == v.size(), "Error in Matrix Vector product" );
    std::vector<DataType> ans;
    ans.resize( nrows, 0.0 );

    for ( UInt ir = 0 + OFFSET;ir < nrows + OFFSET;++ir )
    {
        // column index of non-zero elements and its position respectively
        Container coldata, position;
        UInt nnz_c;
        nnz_c = _Patt->row( ir, coldata, position );
        for ( UInt jc = 0; jc < nnz_c ; jc++ )
        {
            ans[ ir ] += _value[ position[ jc ] + OFFSET ] * v[ coldata[ jc ] + OFFSET ];
        }
    }
    return ans;
}
// version for Vector class
template <class DataType>
Vector
VBRMatr<DataType>::
operator*( const Vector &v ) const
{
    UInt blockSize_r = _Patt->rpntr() [ 1 ] - _Patt->rpntr() [ 0 ];
    //...for square matrices...
    UInt nrows = _Patt->nRows() * blockSize_r;
    UInt ncols = _Patt->nCols() * blockSize_r;
    ASSERT( ncols == v.size(), "Error in Matrix Vector product" );
    Vector ans( nrows );
    ans = 0.;

    for ( UInt ir = 0 + OFFSET;ir < nrows + OFFSET;++ir )
    {
        // column index of non-zero elements and its position respectively
        Container coldata, position;
        UInt nnz_c;
        nnz_c = _Patt->row( ir, coldata, position );
        for ( UInt jc = 0; jc < nnz_c ; jc++ )
        {
            ans( ir ) += _value[ position[ jc ] + OFFSET ] * v( coldata[ jc ] + OFFSET );
        }
    }
    return ans;
}

//necessary for IML++ library: transpose-matrix by vector product
template <class DataType>
Vector VBRMatr<DataType>::
trans_mult( const Vector &v ) const
{
    UInt blockSize_r = _Patt->rpntr() [ 1 ] - _Patt->rpntr() [ 0 ];
    //...for square matrices...
    UInt nrows = _Patt->nRows() * blockSize_r;
    ASSERT( nrows == v.size(), "Error in Matrix Vector product" );
    Vector ans( nrows );
    ans = 0.;

    for ( UInt ir = 0 + OFFSET;ir < nrows + OFFSET;++ir )
    {
        // column index of non-zero elements and its position respectively
        Container coldata, position;

        UInt nnz_c;
        nnz_c = _Patt->row( ir, coldata, position );
        for ( UInt jc = 0; jc < nnz_c ; jc++ )
        {
            ans( coldata[ jc ] + OFFSET ) += _value[ position[ jc ] + OFFSET ] * v( ir );
        }
    }
    return ans;
}

template <class DataType>
VBRMatr<DataType>&
VBRMatr<DataType>::operator*=( const DataType num )
{
    UInt blockSize_r = _Patt->rpntr() [ 1 ] - _Patt->rpntr() [ 0 ];
    UInt blockSize_c = blockSize_r; // we work with square blocks !!!
    UInt stop = _Patt->nNz() * blockSize_r * blockSize_r;
    for ( UInt i = 0;i < stop;++i )
        _value[ i ] *= num;

    return *this;
}

template <class DataType>
VBRMatr<DataType>
VBRMatr<DataType>::operator*( const DataType num )
{
    UInt blockSize_r = _Patt->rpntr() [ 1 ] - _Patt->rpntr() [ 0 ];
    UInt blockSize_c = blockSize_r; // we work with square blocks !!!
    UInt stop = _Patt->nNz() * blockSize_r * blockSize_r;
    VBRMatr<DataType> ans( *this );

    for ( UInt i = 0;i < stop;++i )
        ans.set_mat( i, num * _value[ i ] );

    return ans;
}

template <class DataType>
VBRMatr<DataType>&
VBRMatr<DataType>::flop( const DataType num, VBRMatr<DataType>& M,
                         VBRMatr<DataType>& A )
{
    /* AIM: Matrix = num*M + A */
    UInt blockSize_r = _Patt->rpntr() [ 1 ] - _Patt->rpntr() [ 0 ];
    UInt blockSize_c = blockSize_r; // we work with square blocks !!!
    UInt stop = _Patt->nNz() * blockSize_r * blockSize_r;

    for ( UInt i = 0;i < stop;++i )
        _value[ i ] = num * M.value() [ i ] + A.value() [ i ];

    return *this;
}

template <class DataType>
VBRMatr<DataType>&
VBRMatr<DataType>::flop( const DataType num, VBRMatr<DataType>& M )
{
    /* AIM: Matrix = num*M + Matrix */
    UInt blockSize_r = _Patt->rpntr() [ 1 ] - _Patt->rpntr() [ 0 ];
    UInt blockSize_c = blockSize_r; // we work with square blocks !!!
    UInt stop = _Patt->nNz() * blockSize_r * blockSize_r;

    for ( UInt i = 0;i < stop;++i )
        _value[ i ] += num * M.value() [ i ];

    return *this;
}

template <typename DataType>
void VBRMatr<DataType>::
set_mat_inc( UInt row, UInt col, DataType loc_val )
{
    UInt blockSize = _Patt->rpntr() [ 1 ] - _Patt->rpntr() [ 0 ];
    PatternDefs::Diff_t blocrow = _Patt->rbloc( row + OFFSET );
    PatternDefs::Diff_t bloccol = _Patt->cbloc( col + OFFSET );

    std::pair<UInt, bool> whereBloc = _Patt->locate_index( blocrow, bloccol );

    UInt ind = _Patt->indx() [ whereBloc.first ] + _Patt->locr( row + OFFSET ) +
               _Patt->locc( col + OFFSET ) * ( blockSize - OFFSET );

    if ( whereBloc.second )
        _value[ ind ] += loc_val;

    return ;
}

// version for block operation
template <typename DataType>
DataType addition( DataType val1, DataType val2 )
{
    return val1 + val2;
}

template <typename DataType>
void VBRMatr<DataType>::
set_mat_inc( UInt row, UInt col, std::vector<DataType> &loc_block )
{
    std::pair<UInt, bool> whereBloc = _Patt->locate_index( row, col );

    typename std::vector<DataType>::const_iterator loc_start = loc_block.begin();
    typename std::vector<DataType>::const_iterator loc_end = loc_block.end();

    UInt ind = _Patt->indx() [ whereBloc.first ];
    //copy of the block
    typename std::vector<DataType>::iterator val_start = _value.begin() + ind;
    for ( typename std::vector<DataType>::const_iterator ip = loc_start;
            ip != loc_end; ip++ )
    {
        *val_start += *ip;
        val_start++;
    }
    //I haven't been able to use transform (we would avoid the loop):
    //transform(loc_start, loc_end, val_start, val_start, addition<DataType>);

    return ;
}

template <typename DataType>
void
VBRMatr<DataType>::
set_mat( UInt row, UInt col, DataType loc_val )
{
    UInt blockSize = _Patt->rpntr() [ 1 ] - _Patt->rpntr() [ 0 ];
    PatternDefs::Diff_t blocrow = _Patt->rbloc( row );
    PatternDefs::Diff_t bloccol = _Patt->cbloc( col );

    std::pair<UInt, bool> whereBloc = _Patt->locate_index( blocrow, bloccol );

    UInt ind = _Patt->indx() [ whereBloc.first ] + _Patt->locr( row + OFFSET ) +
               _Patt->locc( col + OFFSET ) * ( blockSize - OFFSET );

    if ( whereBloc.second )
        _value[ ind ] = loc_val;

    return ;
}

// version for block operation
template <typename DataType>
void VBRMatr<DataType>::
set_mat( UInt row, UInt col, std::vector<DataType> &loc_block )
{
    std::pair<UInt, bool> whereBloc = _Patt->locate_index( row, col );

    typename std::vector<DataType>::const_iterator loc_start = loc_block.begin();
    typename std::vector<DataType>::const_iterator loc_end = loc_block.end();

    UInt ind = _Patt->indx() [ whereBloc.first ];
    //copy of the block
    typename std::vector<DataType>::iterator val_start = _value.begin() + ind;
    //  _value[ind] += loc_val;
    if ( whereBloc.second )
        copy( loc_start, loc_end, val_start );

    return ;
}


template <typename DataType>
void
VBRMatr<DataType>::
set_mat( UInt where, DataType loc_val )
{
    _value[ where - OFFSET ] = loc_val;

    return ;
}

template <typename DataType>
void
VBRMatr<DataType>::
spy( std::string const &filename )
{
    // Purpose: Matlab dumping and spy
    std::string nome = filename, uti = " , ";
    UInt nblocrow = _Patt->nRows(), blocsize = _Patt->rpntr() [ 1 ] - _Patt->rpntr() [ 0 ];
    Container ia = _Patt->ia(), ja = _Patt->ja(), indx = _Patt->indx();
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
    for ( UInt irb = 0;irb < nblocrow;++irb )
    {
        for ( UInt ic = ia[ irb ];ic < ia[ irb + 1 ];++ic )
            for ( UInt i = 0;i < blocsize;++i )
                for ( UInt j = 0;j < blocsize;++j )
                    file_out << irb*blocsize + i - OFFSET + 1 << uti <<
                    ja[ ic ] * blocsize + j - OFFSET + 1 << uti <<
                    _value[ indx[ ic ] + i + ( j - OFFSET ) * blocsize ] << std::endl;
    }
    file_out << "];" << std::endl;
    file_out << "I=S(:,1); J=S(:,2); S=S(:,3); A=sparse(I,J,S); spy(A);" << std::endl;
}

/*
template<typename DataType>
void
VBRMatr<DataType>::
showMe() const
{
  UInt blsize=_Patt->rpntr()[1]-_Patt->rpntr()[0]; // block size
  int i_first;
  std::string pare="[";
  std::cout << "**************************" << std::endl;
  std::cout << "     VBR Matrix Pattern   " << std::endl;
  std::cout << std::endl;

  std::cout << pare;
  for (PatternDefs::Diff_t i_index=0; i_index <
  static_cast<PatternDefs::Diff_t>(_Patt->nRows()); i_index++)
    for (UInt jb=0;jb<blsize;jb++){
      std::cout << pare;
      pare = " [";
      i_first=_Patt->ia()[i_index]-OFFSET+1; // In _ia[i_index] there is the
                                            // diagonal entry
      UInt jj=0;
      for(PatternDefs::Diff_t j=0;j<static_cast<PatternDefs::Diff_t>(_Patt->nCols());j++)
 {
   if (j==i_index) for (UInt ib=0;ib<blsize;ib++) std::cout <<
          _value[];
   else {
     if (j==_Patt->ja()[i_first+jj]-OFFSET){
       for (UInt ib=0;ib<blsize;ib++) std::cout << " * ";
       jj++;
     }
     else for (UInt ib=0;ib<blsize;ib++) std::cout << " 0 ";
   }
 }
      if (i_index==static_cast<PatternDefs::Diff_t>(_Patt->nRows()-1))
 std::cout << " ]] " << std::endl;
      else
 std::cout << " ]  " << std::endl;
      std::cout << std::endl;
    }
  std::cout << "**************************" << std::endl;
  return;
}
*/

}
