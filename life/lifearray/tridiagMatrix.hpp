/* -*- mode: c++ -*-
  This program is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
/*!
  \file tridiagMatrix.hpp
  \author Vincent Martin
  \date 07/2004
  \version 1.0

  \brief This file contains a class for tridiagonal matrices.
*/

#ifndef _TRIDIAGMATRIX_H_
#define _TRIDIAGMATRIX_H_

#include <cmath>

#include <life/lifearray/RNM.hpp>


namespace LifeV
{
#define const_R   R


/*!
  \class TriDiagMatrix

   A class for tridiagonal matrices.

*/
template < class R >
class TriDiagMatrix
{
protected:
    UInt _M_order; //!< order of the matrix

    KN<R> _M_diag;    //!< diagonal
    KN<R> _M_updiag;  //!< upper diagonal (size: _M_order-1)
    KN<R> _M_lowdiag; //!< lower diagonal (size: _M_order-1)


public:

    //! Constructor
    TriDiagMatrix( UInt n_mat );

    //! order of the matrix
    int OrderMatrix() const
    {
        return _M_order;
    }

    //! set to zero the matrix (all 3 diagonals)
    void zero();

    //! set the value loc_val at (row , col) in the matrix (may be slow!)
    void set_mat( UInt row, UInt col, const_R loc_val );

    //! add the value loc_val at (row , col) in the matrix (may be slow!)
    void set_mat_inc( UInt row , UInt col , const_R loc_val );

    //! output
    void showMe( std::ostream& c = std::cout, UInt verbose = 0 );

    //! return the three diagonals (const)
    const KN<R>& Diag() const
    {
        return _M_diag;
    }
    const KN<R>& UpDiag() const
    {
        return _M_updiag;
    }
    const KN<R>& LowDiag() const
    {
        return _M_lowdiag;
    }

    //! return the three diagonals
    KN<R>& Diag()
    {
        return _M_diag;
    }
    KN<R>& UpDiag()
    {
        return _M_updiag;
    }
    KN<R>& LowDiag()
    {
        return _M_lowdiag;
    }

    const TriDiagMatrix& operator= ( const TriDiagMatrix<const_R> & mattrid );

    /*!
      simple axpy product.

     @@ It SHOULD use the blas2 routine : (but it does NOT!)
     SUBROUTINE DGBMV ( TRANS, M, N, KL, KU, ALPHA, A, LDA, X, INCX,
                        BETA, Y, INCY )
    *  DGBMV  performs one of the matrix-vector operations
    *
    *     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
    *
    *  where alpha and beta are scalars, x and y are vectors and A is an
    *  m by n band matrix, with kl sub-diagonals and ku super-diagonals.

     One should first transform the three vectors into an adapted array for dgbmv...
    */
    void Axpy( const_R alpha, const ScalUnknown<Vector>& x ,
               const_R beta, ScalUnknown<Vector>& y );
};


//----------------------------------------------------------------------
//! IMPLEMENTATION for TriDiagMatrix
//----------------------------------------------------------------------
//! constructor
template < class R >
TriDiagMatrix<R>::TriDiagMatrix( UInt n_mat ) :
        _M_order( n_mat ),
        _M_diag( n_mat ),
        _M_updiag( n_mat - 1 ),
        _M_lowdiag( n_mat - 1 )
{
    ASSERT_PRE( n_mat > 1,
                "Tridiag matrix class does not accept 1x1 matrices (nor 0x0)." );
}

template < class R >
const TriDiagMatrix<R>&
TriDiagMatrix<R>::operator= ( const TriDiagMatrix<const_R> & mattrid )
{
    _M_order = mattrid.OrderMatrix();
    _M_diag = mattrid.Diag();
    _M_updiag = mattrid.UpDiag();
    _M_lowdiag = mattrid.LowDiag();
    return *this;
}

//! set to zero the matrix (all 3 diagonals)
template < class R >
void TriDiagMatrix<R>::zero()
{
    _M_diag = 0.;
    _M_updiag = 0.;
    _M_lowdiag = 0.;

}

//! set the value loc_val at (row , col) in the matrix (may be slow!)
template < class R >
void TriDiagMatrix<R>::set_mat( UInt row, UInt col, const_R loc_val )
{
    switch ( row - col )
    {
    case 0:
        _M_diag( row ) = loc_val;
        break;
    case 1:
        _M_lowdiag( col ) = loc_val;
        break;
    case - 1:
        _M_updiag( row ) = loc_val;
        break;
    default:
        ERROR_MSG( "Not a tridiagonal matrix. Check the indices." );
    }
}

//! add the value loc_val at (row , col) in the matrix (may be slow!)
template < class R >
void TriDiagMatrix<R>::set_mat_inc( UInt row, UInt col, const_R loc_val )
{
    switch ( (int)row - (int)col )
    {
    case 0:
        _M_diag( row ) += loc_val;
        break;
    case 1:
        _M_lowdiag( col ) += loc_val;
        break;
    case -1:
        _M_updiag( row ) += loc_val;
        break;
    default:
        ERROR_MSG( "Not a tridiagonal matrix. Check the indices." );
    }
}

/*!
  simple axpy product.

  @@ It SHOULD use the blas2 routine : (but it does NOT!)
  SUBROUTINE DGBMV ( TRANS, M, N, KL, KU, ALPHA, A, LDA, X, INCX,
                      BETA, Y, INCY )
*  DGBMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n band matrix, with kl sub-diagonals and ku super-diagonals.

One should first transform the three vectors into an adapted array for dgbmv...
*/
template < class R >
void TriDiagMatrix<R>::Axpy( const_R alpha, const ScalUnknown<Vector>& x ,
                             const_R beta, ScalUnknown<Vector>& y )
{
    //! wrong if _M_order = 1 (but the constructor requires >1)
    y( 0 ) = alpha * ( _M_diag( 0 ) * x( 0 ) +
                       _M_updiag( 0 ) * x( 1 ) )
             + beta * y( 0 );

    UInt last = _M_order - 1;
    y( last ) = alpha * ( _M_lowdiag( last - 1 ) * x( last - 1 ) +
                          _M_diag( last ) * x( last ) )
                + beta * y( last );

    for ( UInt ii = 1; ii < _M_order - 1; ii++ )
    {
        y( ii ) = alpha * ( _M_lowdiag( ii - 1 ) * x( ii - 1 ) +
                            _M_diag( ii ) * x( ii ) +
                            _M_updiag( ii ) * x( ii + 1 ) )
                  + beta * y( ii );
    }
}

//! output
template < class R >
void TriDiagMatrix<R>::showMe( std::ostream& c, UInt verbose )
{
    c << "\n--- TriDiagonal Matrix of order " << OrderMatrix() << "\n";
    if ( verbose > 1 )
    {
        c << "diagonal:\n";
        c << _M_diag << std::endl;
    }
    if ( verbose > 1 )
    {
        c << "upper diagonal:\n";
        c << _M_updiag << std::endl;
    }
    if ( verbose > 1 )
    {
        c << "lower diagonal:\n";
        c << _M_lowdiag << std::endl;
    }
    c << "\n--- End of TriDiagonal Matrix " << std::endl;
}
}
#endif
