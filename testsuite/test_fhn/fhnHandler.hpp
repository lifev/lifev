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
#ifndef _FHNHANDLER_H_
#define _FHNHANDLER_H_
#include "dataFhN.hpp"
#include <life/lifealg/dataAztec.hpp>
#include <life/lifefem/dataTransient.hpp>
#include <life/lifefilters/readMesh3D.hpp>
#include <life/lifefem/refFE.hpp>
#include <life/lifefem/refHdivFE.hpp>
#include <life/lifefem/refHybridFE.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/bcCond.hpp>
#include <life/lifefem/dof.hpp>
#include <life/lifearray/pattern.hpp>
#include <life/lifefem/values.hpp>
#include <life/lifefem/assemb.hpp>
#include <life/lifefilters/vtk_wrtrs.hpp>
#include "user_fct.hpp"

namespace LifeV
{
/*!
  \brief A simple Fitzhugh-Nagumo solver
  \file fhnSolver.hpp
  \author J.-F. Gerbeau
  \date 09/2004

  \par The equations:

   du/dt - div( fhn_diff \grad u) = fhn_f0 u(1-u)(u-fhn_alpha) - v
   dv/dt = fhn_eps (fhn_beta u - fhn_gamma v)

   where fhn_diff,fhn_f0,fhn_eps,fhn_alpha,fhn_beta and fhn_gamma are parameters
   defined in DataFhN
   
   \par Numerical method:

   *** Time:
   
   First order time scheme 
   All the linear terms are implicit, and
   u(1-u)(u-alpha) is explicit (in such a way, the matrix
   does not change during the simulation)

   *** Space:

   P1 tetra
   This can easily be changed. But be careful: we often assume
   that points=nodes (to be improved).

*/
class FhNHandler:
  public DataFhN,
  public DataAztec,
  public DataTransient
{
public:
  const UInt nbCoor; //!< = 3 in 3D, 2 in 2D
  //
  const GeoMap& geoMap;
  const QuadRule& qr;
  const RefFE& refFE;
  CurrentFE fe;
  //
  const GeoMap& geoMapBd;
  const QuadRule& qrBd;
  const RefFE& refBdFE;
  CurrentBdFE feBd;
  //
  double time;
  Dof dof; 
  UInt dimdof; //! number of dof
  RegionMesh3D<LinearTetra> mesh; // the mesh
  int nb_bc;                //!< number of boundary conditions
  BCHandler bc;            //!< boundary conditions handler
  //! boundary conditions functions
  BCFunctionBase bc_fct;
  BCFunctionMixte bc_fct_rob;  //!< a mixte (or Robin) bc function

public:
  FhNHandler(const GetPot& data_file);
};
}
#endif

