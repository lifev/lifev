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
  \file oneDModelSolver.hpp
  \author Vincent Martin
  \date 07/2004
  \version 1.0

  \brief This file contains a solver class for the 1D model. 
*/

#ifndef _ONEDMODELSOLVER_H_
#define _ONEDMODELSOLVER_H_

#include <string>

#include "clapack.h"

#include "oneDModelHandler.hpp"
#include "elemMat.hpp"
#include "elemVec.hpp"
#include "elemOper.hpp"
#include "RNM.hpp"

#include "values.hpp"
#include "assemb.hpp"
#include "chrono.hpp"

#include "tridiagMatrix.hpp"

/*! 
  \class OneDModelSolver

   This class contains a solver class for the 1D model.

*/
class OneDModelSolver:
  public OneDModelHandler 
{
 
public:
  
  //! Constructor 
  /*!
    \param data_file GetPot data file
  */
  OneDModelSolver(const GetPot& data_file);

  //! return the solution at current time step
  const ScalUnknown<Vector>& U_thistime() const { return _M_U_thistime;}

  //! return the solution at next time step
  const ScalUnknown<Vector>& U_nexttime() const{ return _M_U_nexttime;}

  //! Update the right  hand side  for time advancing 
  void timeAdvance();

  //! Update convective term, bc treatment and solve the linearized ns system
  void iterate();

  //! plotting 
  void gplot();

private:

  //! coefficient in front of the corresponding _M_elmat*
  double _M_coeffMass;  
  double _M_coeffStiff;
  double _M_coeffGrad;
  double _M_coeffDiv;

 
  ElemMat _M_elmatMass;  //!< element mass matrix
  ElemMat _M_elmatStiff; //!< element stiffness matrix
  ElemMat _M_elmatGrad;  //!< element gradient matrix
  ElemMat _M_elmatDiv;   //!< element divergence matrix

  ElemVec _M_elvec; // Elementary right hand side

  //! Unknown at present time step
  ScalUnknown<Vector> _M_U_thistime; //! not used??
  //! Unknown at next time step
  ScalUnknown<Vector> _M_U_nexttime;
  //! Flux (at present time step: F_h(U_h^n) )
  ScalUnknown<Vector> _M_FluxU;
  //! Right hand side 
  ScalUnknown<Vector> _M_rhs;
      
  //! tridiagonal mass matrix
  TriDiagMatrix<double> _M_massMatrix;
   //! tridiagonal stiffness matrix
  TriDiagMatrix<double> _M_stiffMatrix;
  //! tridiagonal gradient matrix
  TriDiagMatrix<double> _M_gradMatrix;
  //! tridiagonal divergence matrix
  TriDiagMatrix<double> _M_divMatrix;

  //! lapack LU factorized tridiagonal mass matrix 
  TriDiagMatrix<double> _M_factorMassMatrix;
  //! vectors used by lapack for factorization:
  KN<double> _M_massupdiag2; //!< second upper diagonal (used by lapack) (size _M_order-2)
  KN<int>    _M_massipiv;   //!< indices of pivot in the lapack LU (size _M_order)

  //! Update the element matrices with the current element (from 1)
  void _updateElemMatrices( const UInt& iedge );

  /*! modify the matrix to take into account 
    the Dirichlet boundary conditions 
    (works for P1Seg and canonic numbering!)
  */
  void _updateBCDirichletMatrix( TriDiagMatrix<double>& mat );

  /*! modify the vector to take into account 
    the Dirichlet boundary conditions 
    (works for P1Seg and canonic numbering!)
 
     \param vec : the rhs vector   
     \param val_left  : Dirichlet value inserted to the left
     \param val_right : Dirichlet value inserted to the right
 */
  void _updateBCDirichletVector( ScalUnknown<Vector>& vec, 
				 const double& val_left, 
				 const double& val_right );

  //! update the flux from the current unknown: FfluxU = F_h(U_h^n)
  void _updateFlux();

  //! lapack LU factorization for tridiag matrices 
  //! (modifies _M_factorMassMatrix, _M_massupdiag2 and _M_massipiv.)
  void _factorizeMassMatrix();
  //! lapack LU solve for tridiag matrices (AFTER factorization)
  void _solveMassMatrix( ScalUnknown<Vector>& vec );

  //! direct lapack LU solve for tridiag matrices 
  //! quite useless as you can call it only once! 
  //! REMOVE it some day... (just an example)
  void _directsolveMassMatrix( ScalUnknown<Vector>& vec );

};


#endif
