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
  \file oneDModelHandler.hpp
  \author Vincent Martin
  \date 07/2004 
  \version 1.0

  \brief This file contains a class for the OneD model solver.  

*/

#ifndef _ONEDMODELHANDLER_H_
#define _ONEDMODELHANDLER_H_


#include <cmath>
#include "lifeV.hpp"
#include "dataOneDModel.hpp"
#include "basicOneDMesh.hpp"
#include "dataAztec.hpp"
#include "geoMap.hpp"
#include "currentFE.hpp"
#include "refFE.hpp"
#include "dof.hpp"
#include "bcCond.hpp"
#include "medit_wrtrs.hpp"



/*! 
  \class OneDModelHandler
  
*/

class OneDModelHandler:
  public DataOneDModel,
  public DataAztec
{ 
public:

  typedef Real (*Function)(const Real&, const Real&, const Real&, const Real&, const ID&);

  //! Constructor   
  /*!
    \param data_file GetPot data file
    \param refFE_c reference FE for the concentration
    \param Qr_c volumic quadrature rule for the concentration
    \param bdQr_c surface quadrature rule for the concentration
    \param BCh_c boundary conditions for the concentration
  */
  OneDModelHandler(const GetPot& data_file);
    
  //! Do nothing destructor
  ~OneDModelHandler() {};

   //! Sets initial condition for the concentration 
  //! (incremental approach): the initial time is t0, the time step dt
  void initialize(const Function& c0, Real t0, Real dt); 

 //! Sets initial condition for the concentration from file
  void initialize(const std::string & vname);

  //! Update the right  hand side  for time advancing   
  /*! 
    \param source volumic source  
    \param time present time
  */
  //  virtual void timeAdvance(const Function source, const Real& time) = 0; 

  //! Update convective term, bc treatment and solve the linearized cdr system
  //  virtual void iterate(const Real& time, PhysVectUnknown<Vector> & u) = 0;

  //! Returns the Area vector
  ScalUnknown<Vector>& AreaUnkn();
  
  //! Returns the dof 
  //  const Dof& theDof() const;

   //! Output 
  void showMeHandler(std::ostream& c=std::cout, UInt verbose=0);

 protected:
  
  const UInt _M_nbCoor; //!< := 1 (ONE dimensional model...)

  BasicOneDMesh _M_mesh;

  const GeoMap&   _M_geoMap; 
  const QuadRule& _M_qr;   //!< Quadrature rule for segment elementary computations

  const RefFE& _M_refFE;   //!< Reference Finite Element

  //  Dof  _M_dof;       //!< the degree of freedom  (useless without the 1D mesh)
  UInt _M_dimDof;   //!< number of dof  := nb of vertices

  CurrentFE _M_fe;     //!< current finite element

  int _M_nb_BC;                //!< number of boundary conditions
  BC_Handler _M_BCh;            //!< boundary conditions handler

  //!!! Why  BC_Handler& _BCh_c;??


  //! Area unknown
  ScalUnknown<Vector> _M_AreaUnkn;
  //! Flux unknown
  ScalUnknown<Vector> _M_FluxUnkn;

  //! boundary conditions functions
  BCFunction_Base _M_bc_fct1;  //!< low  X
  BCFunction_Base _M_bc_fct2;  //!< high X
   

};


#endif
