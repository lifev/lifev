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


#include <lifeV.hpp>
#include <cmath>
#include <tab.hpp>

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

public:
    //! constructor
    TriDiagLU( UInt __size ):
        _M_order( __size ),
        _M_massupdiag2( _M_order - 2 ),
        _M_massipiv( _M_order )
    {};

    void LUFactor( TriDiagMat const& __mat );

    void LUSolve( TriDiagMat const& __mat, TriDiagVect& __x );

};


/*!
  LAPACK LU factorization for a tridiagonal
  matrix A [n x n] defined through:
  @  A_diag(n):       main diagonal
  @  A_diag_up(n-1): upper diagonale
  @  A_diag_low(n-1): lower diagonale

  -  LUFactor:  LU factorization
*/
template < class R , class TriDiagMat, class TriDiagVect >
void TriDiagLU<R, TriDiagMat, TriDiagVect>::

/*!
  LAPACK LU solve (AFTER FACTORIZATION) for a tridiagonal
  matrix A [n x n] defined through:
  @  A_diag(n):       main diagonal
  @  A_diag_sup(n-1): upper diagonale
  @  A_diag_low(n-1): lower diagonale

  \param TriDiagMat const& __mat : factorized symmetric
         tridiagonal matrix (only the upper part is used)
  \param TriDiagVect& __x : input : right hand side
                            output: solution
     solve __x = __mat^{-1} * __x

*/
template < class R , class TriDiagMat, class TriDiagVect >
void TriDiagLU<R, TriDiagMat, TriDiagVect>::

}

#endif

