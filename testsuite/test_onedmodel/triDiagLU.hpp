/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Vincent Martin <vincent.martin@mate.polimi.it>
       Date: 2004-11-10

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
   \file triDiagLU.hpp
   \author Vincent Martin <vincent.martin@mate.polimi.it>
   \date 2004-11-10
 */

#ifndef __TRIDIAGLU_H_
#define __TRIDIAGLU_H_


#include <lifeV.hpp>
#include <tab.hpp>
#include <clapack.h>


namespace LifeV
{
#define const_R   R

/*!
  \class TriDiagLU

  This class contains an interface for lapack LU solver

*/
template < class R , class TriDiagMat, class TriDiagVect >
class TriDiagLU
{
private:
    //! order of the system
    UInt _M_order;

    //! vectors used by lapack for factorization:
    Tab1d _M_massupdiag2; //!< second upper diagonal (used by lapack) (size _M_order-2)
    Tab1dInt  _M_massipiv;   //!< indices of pivot in the lapack LU (size _M_order)

    bool _M_isFactorized;
public:
    //! constructor
    TriDiagLU( UInt __size ):
        _M_order( __size ),
        _M_massupdiag2( _M_order - 2 ),
        _M_massipiv( _M_order ),
        _M_isFactorized( false )
    {};

    void Factor( TriDiagMat & __mat );

    void Solve( TriDiagMat const& __mat, TriDiagVect& __x );

};


/*!
  LAPACK LU factorization for a tridiagonal
  matrix A [n x n] defined through:
  @  A_diag(n):       main diagonal
  @  A_diag_up(n-1): upper diagonale
  @  A_diag_low(n-1): lower diagonale

  \param TriDiagMat & __mat : tridiagonal matrix

  LU factorize with lapack _M_factorMassMatrix
  SUBROUTINE DGTTRF( N, DL, D, DU, DU2, IPIV, INFO )
*/
template < class R , class TriDiagMat, class TriDiagVect >
void TriDiagLU<R, TriDiagMat, TriDiagVect>::
Factor( TriDiagMat & __mat )
{
    int INFO = 0;
    int OrderMat = __mat.OrderMatrix();

    ASSERT_PRE( !_M_isFactorized, "Lapack factorization already performed!");

    //! solve with lapack (for tridiagonal matrices)
    dgttrf_( &OrderMat,
             __mat.LowDiag(), __mat.Diag(), __mat.UpDiag(),
             _M_massupdiag2.data().begin(), _M_massipiv.data().begin(),
             &INFO);
    ASSERT_PRE(!INFO,"Lapack factorization of tridiagonal matrix not achieved.");

    _M_isFactorized = true;
}

/*! lapack LU solve AFTER FACTORIZATION of __mat
  SUBROUTINE DGTTRS( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB, INFO )

  *  TRANS   (input) CHARACTER
  *          Specifies the form of the system of equations.
  *          = 'N':  A * X = B  (No transpose)
  *          = 'T':  A'* X = B  (Transpose)
  *          = 'C':  A'* X = B  (Conjugate transpose = Transpose)

  tridiagonal matrix A [n x n] defined through:
  @  A_diag(n):       main diagonal
  @  A_diag_sup(n-1): upper diagonal
  @  A_diag_low(n-1): lower diagonal

  \param TriDiagMat const& __mat : factorized tridiagonal matrix
  \param TriDiagVect& __x : input : right hand side
                            output: solution
     solve __x = __mat^{-1} * __x

*/
template < class R , class TriDiagMat, class TriDiagVect >
void TriDiagLU<R, TriDiagMat, TriDiagVect>::
Solve( TriDiagMat const& __mat, TriDiagVect& __x )
{
    int INFO = 0;
    int NBRHS = 1;//  nb columns of the vec := 1.
    int OrderMat = __mat.OrderMatrix();

    ASSERT_PRE( _M_isFactorized, "Lapack factorization must be performed before!");

    ASSERT_PRE( OrderMat == (int) __x.size() ,
                "The right-hand side must have the same dimensions as the tridiag matrix.");

    //! solve with lapack (for tridiagonal matrices)
    dgttrs_( "N", &OrderMat, &NBRHS,
             __mat.LowDiag(), __mat.Diag(), __mat.UpDiag(),
             _M_massupdiag2.data().begin(), _M_massipiv.data().begin(),
             __x.data().begin(), &OrderMat, &INFO);
    ASSERT_PRE(!INFO,"Lapack solve of tridiagonal matrix not achieved.");

}

}

#endif

