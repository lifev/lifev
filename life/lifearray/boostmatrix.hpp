/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2005-02-22

  Copyright (C) 2005 EPFL

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
   \file boostmatrix.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2005-02-22
 */

#ifndef _BOOSTMATRIX_HPP_
#define _BOOSTMATRIX_HPP_

#include <fstream>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/lambda/lambda.hpp>

#include <life/lifefem/pattern.hpp>
#include <life/lifearray/tab.hpp>

namespace LifeV
{

/**
   \class BoostMatrix

   Wrapper class which allows for using boost uBlas matrices in existing LifeV
   routines. WARNING: Quite a few parts are slow at the moment.
*/
template<typename storage_scheme>
class BoostMatrix
    : public boost::numeric::ublas::compressed_matrix<double, storage_scheme>
{
public:
    //! Constructor from a pattern object

    template<typename PatternType>
    BoostMatrix( const PatternType& pattern )
        : boost::numeric::ublas::compressed_matrix<double, storage_scheme>
    ( pattern.nRows(), pattern.nCols() ) {
        UInt __nnz = pattern.nNz();

        // Save current non-zero entries of the matrix in the vector val

        boost::numeric::ublas::unbounded_array<double> __val( __nnz );

        std::for_each( __val.begin(), __val.end(), boost::lambda::_1 = 0.0 );
        std::copy( this->value_data().begin(), this->value_data().end(), __val.begin() );

        // Make room for non-zero elements

        this->reserve( __nnz, false );

        // Fill matrix from pattern

        UInt  __max_nnz_per_line = 0;
        UInt __nnz_entry = 0;

        this->filled1() = 1;
        this->filled2() = 0;

        for( UInt __row = 0; __row < this->size1(); __row++ ) {
            this->index1_data()[__row] = __nnz_entry;
            ++( this->filled1() );

            UInt __nnz_line = 0;

            // Find number of non-zero entries in __row and fill the matrix with
            // old values

            for( UInt __col = 0; __col < this->size2(); __col++) {
                if( pattern.isThere( __row, __col ) ) {
                    this->index2_data()[__nnz_entry] = __col;
                    ++( this->filled2() );
                    this->value_data()[__nnz_entry] = __val[__nnz_entry];
                    ++__nnz_entry;
                    ++__nnz_line;
                }
                __max_nnz_per_line = std::max( __nnz_line, __max_nnz_per_line );
            }
        }
        this->index1_data()[this->size1()] = this->non_zeros();
    }

    //! equivalent of += for assembly. WARNING: slow
    void set_mat_inc( typename BoostMatrix::size_type i,
                      typename BoostMatrix::size_type j, double inc )
        {
            operator()( i, j ) += inc;
        }

    //! matrix vector multiplication
    Vector operator*( const Vector& b ) const
        {
            return prod(*this, b);
        }

    //! set to zero. WARNING: erases the pattern also
    void zeros()
        {
            clear();
        }

    /** Diagonalization of row r of the system. Done by setting A(r,r) = coeff,
     *  A(r,j) = 0 for j!=r.
     *  @param iRow row to diagonalize
     *  @param coeff value to set the diagonal entry A(r,r) to
     *  @param b right hand side vector to be corrected
     *  @param datum value to set the fix the solution entry x(r) at
     */
    void diagonalize_row( typename BoostMatrix::size_type iRow, double coeff )
        {
            for ( typename BoostMatrix::size_type iCol=0; iCol<iRow; ++iCol )
            {
                erase( iRow, iCol );
            }
            operator()( iRow, iRow ) = coeff;
            for ( typename BoostMatrix::size_type iCol=iRow+1;
                  iCol<size2(); ++iCol )
            {
                erase( iRow, iCol );
            }
        }

    /** Set all elements in a row to zero
     *  @param iRow row to diagonalize
     */
    void zero_row( typename BoostMatrix::size_type iRow) {
        for (typename BoostMatrix::size_type iCol = 0; iCol < size2(); ++iCol) {
            this->erase(iRow, iCol);
        }
    }

    /** Diagonalization of row r of the system. Done by setting A(r,r) = coeff,
     *  A(r,j) = 0 and A(j,r) = 0 for j!=r, and suitably correcting the right
     *  hand side of the system.
     *  @param r row to diagonalize
     *  @param coeff value to set the diagonal entry A(r,r) to
     *  @param b right hand side vector to be corrected
     *  @param datum value to set the fix the solution entry x(r) at
     */
    void diagonalize( typename BoostMatrix::size_type r,
                      double coeff, Vector &b, double datum )
        {
            for ( UInt j = 0; j < size2(); ++j )
            {
                erase( r, j ); // A(r,j) = 0
            }
            for ( UInt i = 0; i < size1(); ++i )
            {
                b[ i ] -= operator()( i, r ) * datum;
                erase( i, r ); // A(i,r) = 0
            }

            operator()( r, r ) = coeff; // A(r,r) = coeff
            b[ r ] = coeff * datum; // correct right hand side for row r
        }

    //! save matrix to file. Should be matlab format, but is not yet.
    //! Matlab format can be achieved though.
    void spy( std::string const &filename )
        {
            std::ofstream file_out( filename.c_str() );
            ASSERT( file_out,
                    "ERROR in BoostMatrix::spy(): File " << filename <<
                    " cannot be opened for writing.");
            file_out << *this << std::endl;
        }
}; // class BoostMatrix

} // namespace LifeV

#endif /* _BOOSTMATRIX_HPP_ */
