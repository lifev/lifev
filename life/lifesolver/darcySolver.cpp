/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2004-10-25

  Copyright (C) 2004 EPFL

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
   \file darcySolver.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-10-25
 */

#include <darcySolver.hpp>


namespace LifeV
{
typedef DarcySolverBase* darcy_solver;

darcy_solver
createDarcyTetra()
{
    return darcy_solver( new DarcySolver<RegionMesh3D<LinearTetra> >( GetPot( "data" ),
                                                                      feTetraRT0,
                                                                      feTetraP0,
                                                                      feTetraRT0Hyb,
                                                                      feTetraRT0VdotNHyb,
                                                                      feTetraP1,
                                                                      quadRuleTetra15pt, quadRuleTria4pt ) );
}

static const bool __rt =  FactoryDarcy::instance().registerProduct( "darcy_tetra", &createDarcyTetra );

darcy_solver
createDarcyHexa()
{
    return darcy_solver( new DarcySolver<RegionMesh3D<LinearHexa> >( GetPot( "data" ),
                                                                     feHexaRT0,
                                                                     feHexaQ0,
                                                                     feHexaRT0Hyb,
                                                                     feHexaRT0VdotNHyb,
                                                                     feHexaQ1,
                                                                     quadRuleHexa8pt, quadRuleQuad4pt) );
}

static const bool __rh =  FactoryDarcy::instance().registerProduct( "darcy_hexa", &createDarcyHexa );

}
