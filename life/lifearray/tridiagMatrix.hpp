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

#include "clapack.h"
#include "RNM.hpp"


#define const_R   R


/*! 
  \class TriDiagMatrix

   A class for tridiagonal matrices.

*/
template< class R >
class TriDiagMatrix
{
private:
  UInt _M_order; //!< order of the matrix

  KN<R> _M_diag;    //!< diagonal 
  KN<R> _M_updiag;  //!< upper diagonal (size: _M_order-1)
  KN<R> _M_lowdiag; //!< lower diagonal (size: _M_order-1)


public:
  
  //! Constructor 
  TriDiagMatrix( UInt n_mat );

  //! order of the matrix
  int OrderMatrix() const {return _M_order;}

  //! set to zero the matrix (all 3 diagonals)
  void zero();

  //! output
  void showMe(std::ostream& c = std::cout, UInt verbose = 0);
  
  //! return the three diagonals
  const KN<R>& Diag()    const { return _M_diag; }
  const KN<R>& UpDiag()  const { return _M_updiag; }
  const KN<R>& LowDiag() const { return _M_lowdiag; }

  const TriDiagMatrix& operator= (const TriDiagMatrix<const_R> & mattrid);

};


//----------------------------------------------------------------------
//! IMPLEMENTATION 
//----------------------------------------------------------------------
//! constructor
template< class R >
TriDiagMatrix<R>::TriDiagMatrix( UInt n_mat ):
  _M_order(n_mat), 
  _M_diag(n_mat),
  _M_updiag(n_mat - 1),
  _M_lowdiag(n_mat - 1)
{}

template< class R >
const TriDiagMatrix<R>& 
TriDiagMatrix<R>::operator= (const TriDiagMatrix<const_R> & mattrid)
{
  _M_order   = mattrid.OrderMatrix();
  _M_diag    = mattrid.Diag();
  _M_updiag  = mattrid.UpDiag();
  _M_lowdiag = mattrid.LowDiag();
  return *this;
}

//! set to zero the matrix (all 3 diagonals)
template< class R > 
void TriDiagMatrix<R>::zero()
{
  _M_diag    = 0.;
  _M_updiag  = 0.;
  _M_lowdiag = 0.;

}
//! set to zero the matrix (all 3 diagonals)
template< class R > 
void TriDiagMatrix<R>::showMe(std::ostream& c, UInt verbose)
{
  c << "\n--- TriDiagonal Matrix of order " << OrderMatrix() << "\n";
  if ( verbose > 1) {
    c << "diagonal:\n";
    c << _M_diag << std::endl;
  }
  if ( verbose > 1) {
    c << "upper diagonal:\n";
    c << _M_updiag << std::endl;
  }
  if ( verbose > 1) {
    c << "lower diagonal:\n";
    c << _M_lowdiag << std::endl;
  }
  c << "\n--- End of TriDiagonal Matrix " << std::endl;
}

#endif
