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

        //UInt  __max_nnz_per_line = 0;
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
                //__max_nnz_per_line = std::max( __nnz_line, __max_nnz_per_line );
            }
        }
        this->index1_data()[this->size1()] = this->non_zeros();
    }

    //! equivalent of += for assembly.
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
     *  @param r row to diagonalize
     *  @param coeff value to set the diagonal entry A(r,r) to
     *  @param b right hand side vector to be corrected
     *  @param datum value to set the fix the solution entry x(r) at
     */
    void diagonalize_row( typename BoostMatrix::size_type r, double coeff )
        {
            zero_row( r );
            operator()( r, r ) = coeff;
        }

    /** Set all elements in a row to zero
     *  @param iRow row to diagonalize
     */
    void zero_row( typename BoostMatrix::size_type iRow )
        {
            typename BoostMatrix::iterator1 row = begin1();
            for ( UInt i=0; i<r; ++i, ++row );
            std::for_each( row.begin(), row.end(), boost::lambda::_1 = 0.0 );
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
            zero_row( r );
            typename BoostMatrix::iterator2 col = begin2();
            for ( UInt i=0; i<r; ++i, ++col );
            for ( typename BoostMatrix::iterator1 iRow = col.begin();
                  iRow != col.end(); ++iRow )
            {
                b[ iRow.index1() ] -= *iRow * datum;
                *iRow = 0;
            }
            operator()( r, r ) = coeff; // A(r,r) = coeff
            b[ r ] = coeff * datum; // correct right hand side for row r
        }

    //! Dump matrix to file in Matlab format and spy
    void spy( std::string const &filename )
        {
            std::string name = filename;
            std::string separator = " , ";

            // check on the file name
            int i = filename.find( "." );

            if ( i <= 0 )
                name = filename + ".m";
            else
            {
                if ( ( unsigned int ) i != filename.size() - 2 ||
                     filename[ i + 1 ] != 'm' )
                {
                    std::cerr << "Wrong file name ";
                    name = filename + ".m";
                }
            }

            std::ofstream file_out( name.c_str() );
            ASSERT( file_out, "[BoostMatrix::spy] ERROR: File " << filename <<
                    " cannot be opened for writing.");

            file_out << "S = [ ";
            for ( typename BoostMatrix::iterator1 i1=begin1();
                  i1!=end1(); ++i1 )
            {
                for ( typename BoostMatrix::iterator2 i2=i1.begin();
                      i2!=i1.end(); ++i2 )
                    file_out << i2.index1() + 1 << separator
                             << i2.index2() + 1 << separator
                             << *i2  << std::endl;
            }
            file_out << "];" << std::endl;
            file_out << "I=S(:,1); J=S(:,2); S=S(:,3);" << std::endl;
            file_out << "A=sparse(I,J,S); spy(A);" << std::endl;
        }

}; // class BoostMatrix

} // namespace LifeV

#endif /* _BOOSTMATRIX_HPP_ */
