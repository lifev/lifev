/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2005-02-03

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
   \file blockMatrix.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2005-02-03
 */

#ifndef BLOCKMATRIX_HPP
#define BLOCKMATRIX_HPP

// forward declarations for aztec
// struct AZ_MATRIX_STRUCT;
// typedef struct AZ_MATRIX_STRUCT AZ_MATRIX;
// void* AZ_get_matvec_data( AZ_MATRIX *aMat );

namespace LifeV
{

//! matrix consisting of 4 blocks
class BlockMatrix
{
public:
    //! constructor
    BlockMatrix( const MSRMatr<double>* matrix11,
                 const CSRMatr<CSRPatt, double>* matrix12,
                 const CSRMatr<CSRPatt, double>* matrix21,
                 const MSRMatr<double>* matrix22 = 0);
    //! matrix vector product
    Vector operator*( const Vector& v ) const;
    //! matrix vector product for aztec
    void matrixVectorProductAztec( double* p, double* ap ) const;
    //! (total) number of rows
    UInt nRows() const { return M_nRows1 + M_nRows2; }
    //! (total) number of columns
    UInt nCols() const { return M_nCols1 + M_nCols2; }
private:
    //! submatrix 1,1
    const MSRMatr<double>* M_matrix11;

    //! submatrix 1,2
    const CSRMatr<CSRPatt, double>* M_matrix12;

    //! submatrix 2,1
    const CSRMatr<CSRPatt, double>* M_matrix21;

    //! submatrix 2,2
    const MSRMatr<double>* M_matrix22;

    //! number of rows of submatrices 1,1 and 1,2
    const UInt M_nRows1;

    //! number of rows of submatrices 2,1 and 2,2
    const UInt M_nRows2;

    //! number of columns of submatrices 1,1 and 2,1
    const UInt M_nCols1;

    //! nubmer of columns of submatrices 1,2 and 2,2
    const UInt M_nCols2;
}; // class BlockMatrix

//! matrix vector product function with interface for Aztec
//! usage: BlockMatrix a(a11, a12, a21, a22);
//!        SolverAztec.setMatrixFree(a.nRows(), &a, &blockMatrixVectorProduct);
void blockMatrixVectorProduct( double* p, double* ap,
                               AZ_MATRIX* aMat, int proc_config[] );

// IMPLEMENTATIONS

BlockMatrix::BlockMatrix( const MSRMatr<double>* matrix11,
                          const CSRMatr<CSRPatt, double>* matrix12,
                          const CSRMatr<CSRPatt, double>* matrix21,
                          const MSRMatr<double>* matrix22 )
    : M_matrix11( matrix11 ), M_matrix12( matrix12 ),
      M_matrix21( matrix21 ), M_matrix22( matrix22 ),
      M_nRows1( matrix11->Patt()->nRows() ),
      M_nRows2( matrix21->Patt()->nRows() ),
      M_nCols1( matrix11->Patt()->nCols() ),
      M_nCols2( matrix12->Patt()->nCols() ) { }


Vector BlockMatrix::operator*( const Vector& v ) const
{
    //const boost::numeric::ublas::vector_range<const Vector>
    //    v1( v, boost::numeric::ublas::range ( 0, M_nCols1 ) );
    //const boost::numeric::ublas::vector_range<const Vector>
    //    v2( v, boost::numeric::ublas::range ( M_nCols1, M_nCols1 + M_nCols2 ) );
    Vector v1( M_nCols1 );
    for ( UInt i=0; i<M_nCols1; ++i)
        v1[i] = v[i];
    Vector v2( M_nCols2 );
    for ( UInt i=0; i<M_nCols2; ++i)
        v2[i] = v[i+M_nCols1];
    Vector mv( M_nRows1+M_nRows2 );
    boost::numeric::ublas::vector_range<Vector>
        mv1( mv, boost::numeric::ublas::range ( 0, M_nRows1 ) );
    boost::numeric::ublas::vector_range<Vector>
        mv2( mv, boost::numeric::ublas::range ( M_nRows1,
                                                M_nRows1 + M_nRows2 ) );
    mv1 = (*M_matrix11) * v1 + (*M_matrix12) * v2;
    if ( M_matrix22 )
        mv2 = (*M_matrix21) * v1 + (*M_matrix22) * v2;
    else
        mv2 = (*M_matrix21) * v1;
    return mv;
}

void BlockMatrix::matrixVectorProductAztec( double* p, double* ap) const
{
    Vector v1( M_nCols1 );
    for ( UInt i=0; i<M_nCols1; ++i)
        v1[i] = p[i];

    Vector v2( M_nCols2 );
    for ( UInt i=0; i<M_nCols2; ++i)
        v2[i] = p[i+M_nCols1];

    Vector av1( M_nRows1 );
    av1 = (*M_matrix11) * v1 + (*M_matrix12) * v2;
    for ( UInt i=0; i<M_nRows1; ++i )
        ap[i] = av1[i];

    Vector av2( M_nRows2 );
    if ( M_matrix22 )
        av2 = (*M_matrix21) * v1 + (*M_matrix22) * v2;
    else
        av2 = (*M_matrix21) * v1;
    for ( UInt i=0; i<M_nRows2; ++i )
        ap[i+M_nRows1] = av2[i];
}

void blockMatrixVectorProduct( double* p, double* ap,
                               AZ_MATRIX* aMat, int proc_config[] )
{
    BlockMatrix& a = *(static_cast<BlockMatrix*>(AZ_get_matvec_data( aMat )));
    a.matrixVectorProductAztec( p, ap );
}

} // namespace LifeV

#endif
