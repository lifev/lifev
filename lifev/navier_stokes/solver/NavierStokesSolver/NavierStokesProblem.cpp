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

template< class mesh_Type >
Real
NavierStokesProblem< mesh_Type >::nullFunction( const Real&, const Real&, const Real&,
                                   const Real&, const ID& )
{
	return 0.0;
}

template< class mesh_Type >
NavierStokesProblem< mesh_Type >::NavierStokesProblem()
    : M_refinement ( 0 ), M_resourcesPath ( "" ), M_viscosity ( 1.0 ), M_density ( 1.0 )
{

}

template< class mesh_Type >
NavierStokesProblem< mesh_Type >::~NavierStokesProblem()
{

}

template< class mesh_Type >
bool
NavierStokesProblem< mesh_Type >::hasExactSolution() const
{
    return false;
}

template< class mesh_Type >
typename NavierStokesProblem< mesh_Type >::function_Type
NavierStokesProblem< mesh_Type >::xexact()
{
    return 0;
}

template< class mesh_Type >
typename NavierStokesProblem< mesh_Type >::function_Type
NavierStokesProblem< mesh_Type >::uexact()
{
    return 0;
}

template< class mesh_Type >
typename NavierStokesProblem< mesh_Type >::function_Type
NavierStokesProblem< mesh_Type >::uderexact()
{
    return 0;
}

template< class mesh_Type >
typename NavierStokesProblem< mesh_Type >::function_Type
NavierStokesProblem< mesh_Type >::pexact()
{
    return 0;
}

template< class mesh_Type >
void
NavierStokesProblem< mesh_Type >::setMesh ( const UInt& refinement,
                               const std::string& resourcesPath )
{
    M_refinement    = refinement;
    M_resourcesPath = resourcesPath;
}

template< class mesh_Type >
void
NavierStokesProblem< mesh_Type >::setViscosity ( const Real& viscosity )
{
    M_viscosity = viscosity;
}

template< class mesh_Type >
void
NavierStokesProblem< mesh_Type >::setDensity ( const Real& density )
{
    M_density = density;
}

template< class mesh_Type >
Real
NavierStokesProblem< mesh_Type >::viscosity() const
{
    return M_viscosity;
}

template< class mesh_Type >
Real
NavierStokesProblem< mesh_Type >::density() const
{
    return M_density;
}

template< class mesh_Type >
typename NavierStokesProblem< mesh_Type >::function_Type
NavierStokesProblem< mesh_Type >::force()
{
    return nullFunction;
}

} // namespace LifeV
