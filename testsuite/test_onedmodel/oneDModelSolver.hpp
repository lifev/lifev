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

  //! Update the right  hand side  for time advancing 
  /*! 
    \param source volumic source  
    \param time present time
  */
  void timeAdvance(const Function source, const double& time);

  //! Update convective term, bc treatment and solve the linearized ns system
  void iterate(const double& time, PhysVectUnknown<Vector> & u);

 private:

  //! coefficient in the corresponding _M_elmat*
  double _M_coeffMass;  
  double _M_coeffStiff;
  double _M_coeffGrad;
  double _M_coeffDiv;

 
  ElemMat _M_elmatMass;  //!< element mass matrix
  ElemMat _M_elmatStiff; //!< element stiffness matrix
  ElemMat _M_elmatGrad;  //!< element gradient matrix
  ElemMat _M_elmatDiv;   //!< element divergence matrix

  ElemVec _M_elvec; // Elementary right hand side

  //! Right  hand  side 
  ScalUnknown<Vector> _M_rhs;
  
  //! tridiagonal mass matrix
  TriDiagMatrix<double> _M_massMatrix;

  //! Update the element matrices with the current element (from 1)
  void _updateElemMatrices( const UInt& iedge );

};


#endif
