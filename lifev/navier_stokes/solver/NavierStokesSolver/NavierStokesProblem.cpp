//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
    @file NavierStokesProblem abstract class
    @brief This class contains all the informations necessary to generate a Navier-Stokes problem

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 23-10-2011
 */

#include <lifev/navier_stokes/solver/NavierStokesSolver/NavierStokesProblem.hpp>

namespace LifeV
{

NavierStokesProblem::NavierStokesProblem()
    : M_refinement ( 0 ), M_resourcesPath ( "" ), M_viscosity ( 1.0 ), M_density ( 1.0 )
{

}

NavierStokesProblem::~NavierStokesProblem()
{

}

bool
NavierStokesProblem::hasExactSolution() const
{
    return false;
}

NavierStokesProblem::function_Type
NavierStokesProblem::xexact()
{
    return 0;
}

NavierStokesProblem::function_Type
NavierStokesProblem::uexact()
{
    return 0;
}

NavierStokesProblem::function_Type
NavierStokesProblem::uderexact()
{
    return 0;
}

NavierStokesProblem::function_Type
NavierStokesProblem::pexact()
{
    return 0;
}

void
NavierStokesProblem::setMesh ( const UInt& refinement,
                               const std::string& resourcesPath )
{
    M_refinement    = refinement;
    M_resourcesPath = resourcesPath;
}

void
NavierStokesProblem::setViscosity ( const Real& viscosity )
{
    M_viscosity = viscosity;
}

void
NavierStokesProblem::setDensity ( const Real& density )
{
    M_density = density;
}

Real
NavierStokesProblem::viscosity() const
{
    return M_viscosity;
}

Real
NavierStokesProblem::density() const
{
    return M_density;
}

} // namespace LifeV
