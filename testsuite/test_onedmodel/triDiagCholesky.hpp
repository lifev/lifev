/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Vincent Martin <vincent.martin@mate.polimi.it>
       Date: 2004-11-09

  Copyright (C) 2004 Politecnico di Milano

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
   \file triDiagCholesky.hpp
   \author Vincent Martin <vincent.martin@mate.polimi.it>
   \date 2004-11-09
 */


#ifndef __TRIDIAGCHOLESKY_H_
#define __TRIDIAGCHOLESKY_H_

#include <life/lifecore/life.hpp>
#include <cmath>

namespace LifeV
{
#define const_R   R

/*!
  \class TriDiagCholesky

  This class contains a basic Cholesky algorithm

*/
template < class R , class TriDiagMat, class TriDiagVect >
class TriDiagCholesky
{
private:
    bool _M_isFactorized;
public:
    TriDiagCholesky( UInt __unused ):
        _M_isFactorized( false ) {};

    void Factor( TriDiagMat & __mat );

    void Solve( TriDiagMat const& __mat,
                        TriDiagVect& __x );

};


/*!
  Cholesky factorization for a symmetric tridiagonal
  matrix A [n x n] defined through:
  @  A_diag(n):       main diagonal
  @  A_diag_sup(n-1): upper diagonale

  -  CholeskyFactor:  Cholesky factorization
*/
template < class R , class TriDiagMat, class TriDiagVect >
void TriDiagCholesky<R, TriDiagMat, TriDiagVect>::
Factor( TriDiagMat & __mat )
{
    ASSERT_PRE( !_M_isFactorized, "Cholesky factorization already performed!");

    KN<R> & m_diag = __mat.Diag();
    KN<R> & m_updiag = __mat.UpDiag();

    m_diag( 0 ) = std::sqrt( m_diag( 0 ) );
    for( Int ii = 1 ; ii < __mat.OrderMatrix() ; ii ++ ) {
        R val_ud = m_updiag( ii - 1) / m_diag( ii - 1 );
        m_updiag( ii - 1 ) = val_ud;
        m_diag( ii ) = std::sqrt( m_diag(ii) - val_ud * val_ud );
    }
    _M_isFactorized = true;
}

/*!
  Cholesky solve (AFTER FACTORIZATION) for a symmetric tridiagonal
  matrix A [n x n] defined through:
  @  A_diag(n):       main diagonal
  @  A_diag_sup(n-1): upper diagonale

  \param TriDiagMat const& __mat : factorized symmetric
         tridiagonal matrix (only the upper part is used)
  \param TriDiagVect& __x : input : right hand side
                            output: solution
     solve __x = __mat^{-1} * __x

  -  backsub_tri: alg. of back substitution for triang. sup matrices
    solve the linear system Ax=b with A symmetric lower triangular
    with the main diagonal and the upper diag (symmetry!)
    (algorithm of forward-substitution)
  -  forsub_tri: alg. of forward substitution for triang. inf matrices
    solve the linear system Ax=b with A symmetric upper triangular
    with the main diagonal and the upper diag.
    (algorithm of backward-substitution)
*/
template < class R , class TriDiagMat, class TriDiagVect >
void TriDiagCholesky<R, TriDiagMat, TriDiagVect>::
Solve( TriDiagMat const& __mat,
               TriDiagVect& __x )
{
    ASSERT_PRE( _M_isFactorized, "Cholesky factorization must be performed before!");

    KN<R> const& m_diag = __mat.Diag();
    KN<R> const& m_updiag = __mat.UpDiag();

    int dim = __mat.OrderMatrix();

    //! forward-substitution
    __x( 0 ) = __x( 0 ) / m_diag( 0 );
    for( Int ii = 1 ; ii < dim ; ii ++ ) {
        __x( ii ) = ( __x( ii ) - m_updiag( ii - 1 ) * __x( ii - 1 ) )
            / m_diag( ii );
    }
    //! backward-substitution
    __x( dim - 1 ) = __x( dim - 1 ) / m_diag( dim - 1 );
    for( Int ii = dim - 2 ; ii >= 0 ; ii -- ) {
        __x( ii ) = ( __x( ii ) - m_updiag( ii ) * __x( ii + 1 ) )
            / m_diag( ii );
    }
}

}

#endif
