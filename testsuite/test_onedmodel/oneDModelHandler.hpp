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

  \brief This file contains a class for the OneD model handler.

*/

#ifndef _ONEDMODELHANDLER_H_
#define _ONEDMODELHANDLER_H_


#include <cmath>
#include <utility>

#include <life/lifecore/life.hpp>
#include <dataOneDModel.hpp>
#include <basicOneDMesh.hpp>
#include <life/lifefem/geoMap.hpp>
#include <life/lifefem/currentFE.hpp>
#include <life/lifefem/refFE.hpp>
#include <dofOneD.hpp>

#include <gracePlot.hpp>

namespace LifeV
{
/*!
  \class OneDModelHandler

*/

class OneDModelHandler:
  public DataOneDModel
{
public:

  typedef pair< Real, Real > Vec2D;


  typedef Real (*Function)(const Real&, const Real&, const Real&, const Real&, const ID&);

  //! Constructor
  /*!
    \param data_file GetPot data file
    \param refFE reference FE
    \param Qr volumic quadrature rule
    \param bdQr surface quadrature rule
    \param BCh boundary conditions
  */
  OneDModelHandler(const GetPot& data_file);

  //! Do nothing destructor
  ~OneDModelHandler() {};

  //! Update the right  hand side  for time advancing
  /*!
    \param source volumic source
    \param time present time
  */
  // virtual void timeAdvance(const Function source, const Real& time) = 0;

  //! Update convective term, bc treatment and solve the linearized cdr system
  //  virtual void iterate(const Real& time, PhysVectUnknown<Vector> & u) = 0;

   //! Output
  void showMeHandler(std::ostream& c=std::cout, UInt verbose=0);

protected:

  const UInt _M_nbCoor; //!< := 1 (ONE dimensional model...)

  BasicOneDMesh _M_mesh;

  const GeoMap&   _M_geoMap;
  const QuadRule& _M_qr;   //!< Quadrature rule for segment elementary computations

  const RefFE& _M_refFE;   //!< Reference Finite Element

  DofOneD _M_dof1D; //!< the simplified degrees of freedom
  const UInt _M_dimDof;   //!< number of dof  := nb of vertices

  CurrentFE _M_fe;     //!< current finite element


  //! Dirichlet boundary value at left and right boundaries (NO BCHandler)
  Vec2D _M_bcDirLeft; //! first -> U1, second ->U2
  Vec2D _M_bcDirRight;

  GracePlot _M_GracePlot; //!< for plotting

};
}

#endif
