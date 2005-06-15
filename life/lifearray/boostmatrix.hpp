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
   \author Daniele A. Di Pietro <dipietro@unibg.it>
   \date 2005-02-22
*/

#ifndef _BOOSTMATRIX_HPP_
#define _BOOSTMATRIX_HPP_

#include <fstream>

#include <boost/numeric/ublas/matrix_sparse.hpp>
/* #include <boost/numeric/ublas/banded.hpp> */
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
        : public boost::numeric::ublas::compressed_matrix<double, storage_scheme, 0, boost::numeric::ublas::unbounded_array<int> >
    {
    public:
        typedef boost::numeric::ublas::compressed_matrix<double, storage_scheme, 0, boost::numeric::ublas::unbounded_array<int> > super;
        //! empty constructor
        BoostMatrix( typename BoostMatrix::size_type size1,
                     typename BoostMatrix::size_type size2 )
            : super( size1, size2 ) { }

        //! Copy constructor

        template<typename __E>
        BoostMatrix( const boost::numeric::ublas::matrix_expression<__E>& M)
            : super(M) { }

        //! Constructor from a pattern object

        BoostMatrix( const BasePattern& pattern )
            : super( pattern.nRows(), pattern.nCols() ) {
            UInt __nnz = pattern.nNz();

            // Save current non-zero entries of the matrix in vector __val

            boost::numeric::ublas::unbounded_array<double> __val( __nnz );

            std::for_each( __val.begin(), __val.end(), boost::lambda::_1 = 0.0 );
            std::copy( this->value_data().begin(), this->value_data().end(), __val.begin() );

            // Make room for non-zero elements

            this->reserve( __nnz, false );

            // Fill matrix from pattern

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
                }
            }
            this->index1_data()[this->size1()] = this->non_zeros();
        }

        template <UInt BROWS, UInt BCOLS, typename PATTERN>
        BoostMatrix( const MixedPattern<BROWS, BCOLS, PATTERN>& pattern );

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

        //! set to zero without erasing pattern
        void zeros();

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
            for ( UInt i=0; i<iRow; ++i, ++row );
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

        //! Operator =
        template<typename __storage_scheme, typename __E>
        BoostMatrix<__storage_scheme> operator=(const boost::numeric::ublas::matrix_expression<__E>& M) {
            return super::operator=(M);
        }

        //! Dump matrix to file in Matlab format and spy
        void spy( std::string const &filename ) const
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
            for ( typename BoostMatrix::const_iterator1 i1=begin1();
                  i1!=end1(); ++i1 )
                {
                    for ( typename BoostMatrix::const_iterator2 i2=i1.begin();
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

    template<>
    template<UInt BROWS, UInt BCOLS, typename PATTERN>
    BoostMatrix<boost::numeric::ublas::row_major>::BoostMatrix( const MixedPattern<BROWS, BCOLS, PATTERN>& pattern )
        : super( pattern.nRows(), pattern.nCols() )
    {
        UInt __nnz = pattern.nNz();

        // Save current non-zero entries of the matrix in vector __val

        boost::numeric::ublas::unbounded_array<double> __val( __nnz );

        std::for_each( __val.begin(), __val.end(), boost::lambda::_1 = 0.0 );
        std::copy( this->value_data().begin(), this->value_data().end(), __val.begin() );

        // Make room for non-zero elements

        this->reserve( __nnz, false );

        // Fill matrix from pattern

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

    template<>
    template<UInt BROWS, UInt BCOLS, typename PATTERN>
    BoostMatrix<boost::numeric::ublas::column_major>::BoostMatrix( const MixedPattern<BROWS, BCOLS, PATTERN>& pattern )
        : super( pattern.nRows(), pattern.nCols() )
    {
        UInt __nnz = pattern.nNz();

        // Save current non-zero entries of the matrix in vector __val

        boost::numeric::ublas::unbounded_array<double> __val( __nnz );

        std::for_each( __val.begin(), __val.end(), boost::lambda::_1 = 0.0 );
        std::copy( this->value_data().begin(), this->value_data().end(), __val.begin() );

        // Make room for non-zero elements

        this->reserve( __nnz, false );

        // Fill matrix from pattern

        UInt __nnz_entry = 0;

        this->filled1() = 1;
        this->filled2() = 0;

        for( UInt __col = 0; __col < this->size2(); __col++ ) {
            this->index1_data()[__col] = __nnz_entry;
            ++( this->filled1() );

            UInt __nnz_line = 0;

            // Find number of non-zero entries in __row and fill the matrix with
            // old values

            for( UInt __row = 0; __row < this->size1(); __row++) {
                if( pattern.isThere( __row, __col ) ) {
                    this->index2_data()[__nnz_entry] = __row;
                    ++( this->filled2() );
                    this->value_data()[__nnz_entry] = __val[__nnz_entry];
                    ++__nnz_entry;
                    ++__nnz_line;
                }
            }
        }
        this->index1_data()[this->size2()] = this->non_zeros();
    }

    template<>
    void BoostMatrix<boost::numeric::ublas::row_major>::zeros()
    {
        for ( BoostMatrix::iterator1 i1=begin1();
              i1!=end1(); ++i1 )
            {
                for ( BoostMatrix::iterator2 i2=i1.begin();
                      i2!=i1.end(); ++i2 )
                    {
                        *i2 = 0;
                    }
            }
    }

    template <>
    void BoostMatrix<boost::numeric::ublas::column_major>::zeros()
    {
        for ( BoostMatrix::iterator2 i2=begin2();
              i2!=end2(); ++i2 )
            {
                for ( BoostMatrix::iterator1 i1=i2.begin();
                      i1!=i2.end(); ++i1 )
                    {
                        *i1 = 0;
                    }
            }
    }

    class DiagonalBoostMatrix : public boost::numeric::ublas::compressed_matrix<double, boost::numeric::ublas::row_major, 0, boost::numeric::ublas::unbounded_array<int> >
    {
    public:
        typedef boost::numeric::ublas::compressed_matrix<double, boost::numeric::ublas::row_major, 0, boost::numeric::ublas::unbounded_array<int> > super;
        DiagonalBoostMatrix( size_type n )
            : super( n, n, n )
        {
            for( size_type i=0; i<this->size1(); ++i )
                {
                    this->push_back( i, i, 0 );
                }
        }

        //! Copy constructor

        template<typename __E>
        DiagonalBoostMatrix( const boost::numeric::ublas::matrix_expression<__E>& M)
            : super(M) { }

        void invert()
        {
            for ( size_type i=0; i<this->size1(); ++i )
                {
                    this->value_data()[ i ] = 1. / this->value_data()[ i ];
                }
        }

        //! Operator =
        template<typename __E>
        DiagonalBoostMatrix operator=(const boost::numeric::ublas::matrix_expression<__E>& M) {
            return super::operator=(M);
        }

        template<typename matrix_type>
        void lumpRowSum( const matrix_type& m )
        {
            for ( typename matrix_type::const_iterator1 i1=m.begin1();
                  i1!=m.end1(); ++i1 )
                {
                    double rowSum = 0;
                    for ( typename matrix_type::const_iterator2 i2=i1.begin();
                          i2!=i1.end(); ++i2)
                        {
                            rowSum += *i2;
                        }
                    this->value_data()[ i1.index1() ] = rowSum;
                }
        }

        template<typename matrix_type>
        void lumpDiagonalScaling( const matrix_type& m )
        {
            double totalMass = 0;
            double diagonalMass = 0;
            for ( typename matrix_type::const_iterator1 i1=m.begin1();
                  i1!=m.end1(); ++i1 )
                {
                    for ( typename matrix_type::const_iterator2 i2=i1.begin();
                          i2!=i1.end(); ++i2)
                        {
                            totalMass += *i2;
                        }
                    diagonalMass += m( i1.index1(), i1.index1() );
                }
            double factor = totalMass / diagonalMass;
            for ( size_type i=0; i<size1(); ++i )
                {
                    this->operator()( i, i ) = m( i, i ) * factor;
                }
        }

        //! Dump matrix to file in Matlab format and spy
        void spy( std::string const &filename ) const
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
            ASSERT( file_out, "[DiagonalBoostMatrix::spy] ERROR: File " << filename <<
                    " cannot be opened for writing.");

            file_out << "S = [ ";

            for ( DiagonalBoostMatrix::const_iterator1 i1=begin1();
                  i1!=end1(); ++i1 )
                {
                    for ( DiagonalBoostMatrix::const_iterator2 i2=i1.begin();
                          i2!=i1.end(); ++i2 )
                        file_out << i2.index1() + 1 << separator
                                 << i2.index2() + 1 << separator
                                 << *i2  << std::endl;
                }

            file_out << "];" << std::endl;
            file_out << "I=S(:,1); J=S(:,2); S=S(:,3);" << std::endl;
            file_out << "A=sparse(I,J,S); spy(A);" << std::endl;
        }

    }; // class DiagonalBoostMatrix

    /** efficient (Schur) product of sparse, diagonal, and sparse matrix D*H*G
        @param D first matrix
        @param H second matrix, diagonal
        @param G third matrix
        @param S schur product matrix, empty sparse matrix at input
    */
    void schurProduct( const BoostMatrix<boost::numeric::ublas::row_major>& D,
                       const DiagonalBoostMatrix& H,
                       const BoostMatrix<boost::numeric::ublas::column_major>& G,
                       boost::numeric::ublas::compressed_matrix<double, boost::numeric::ublas::row_major>& S )
    {
        using namespace boost::numeric::ublas;
        ASSERT( D.size2() == H.size1() && H.size2() == G.size1() &&
                D.size1() == S.size1() && G.size2() == S.size2(),
                "[schurProduct] ERROR: cannot multiply " <<
                D.size1() << "x" << D.size2() << " matrix by " <<
                H.size1() << "x" << H.size2() << " matrix by " <<
                G.size1() << "x" << G.size2() << " matrix into " <<
                S.size1() << "x" << S.size2() << " matrix." );

        // Chain multiplication
        S = prod(D,  prod<compressed_matrix<double, column_major> >(H, G) );
        /*
          for ( BoostMatrix<row_major>::const_iterator1 iD1=D.begin1(); iD1!=D.end1(); ++iD1 ) {
          for ( BoostMatrix<column_major>::const_iterator2 iG2=G.begin2(); iG2!=G.end2(); ++iG2 ) {
          double sij = 0;
          BoostMatrix<column_major>::const_iterator1 iG1=iG2.begin();
          for ( BoostMatrix<row_major>::const_iterator2 iD2 = iD1.begin(); iD2!=iD1.end(); ++iD2 ) {
          while ( iG1.index1() < iD2.index2() ) ++iG1;
          if ( iD2.index2() == iG1.index1() )
          sij += (*iD2) * H.value_data()[iD2.index2()] * (*iG1);

          }
          S.insert( iD1.index1(), iG2.index2(), sij );
          }
          }
        */
    }
} // namespace LifeV

#endif /* _BOOSTMATRIX_HPP_ */
