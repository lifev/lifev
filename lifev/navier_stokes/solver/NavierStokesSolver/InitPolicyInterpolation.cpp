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
    @file InitPolicyInterpolation class
    @brief This class is a strategy to initialize a Navier-Stokes problem

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 11-12-2012
 */

#include <lifev/navier_stokes/solver/NavierStokesSolver/InitPolicyInterpolation.hpp>

#include <string>
#include <lifev/core/util/LifeChrono.hpp>

namespace LifeV
{

void
InitPolicyInterpolation::
setupInit ( Teuchos::ParameterList& /*list*/ )
{

}

void
InitPolicyInterpolation::
initSimulation ( bcContainerPtr_Type /*bchandler*/,
                 vectorPtr_Type solution )
{
    ASSERT ( problem()->hasExactSolution(), "Error: You cannot use the interpolation method if the problem has not an analytical solution." );

    Real currentTime = initialTime() - timestep() * bdf()->order();
    UInt pressureOffset = uFESpace()->fieldDim() * uFESpace()->dof().numTotalDof();

    vectorPtr_Type velocity;
    velocity.reset ( new vector_Type ( uFESpace()->map(), Unique ) );

    vectorPtr_Type pressure;
    pressure.reset ( new vector_Type ( pFESpace()->map(), Unique ) );

    *solution = 0.0;
    uFESpace()->interpolate ( problem()->uexact(), *velocity, currentTime );
    pFESpace()->interpolate ( problem()->pexact(), *pressure, currentTime );
    solution->add ( *velocity );
    solution->add ( *pressure, pressureOffset );
    bdf()->setInitialCondition ( *solution );

    currentTime += timestep();
    for ( ; currentTime <=  initialTime() + timestep() / 2.0; currentTime += timestep() )
    {
        *solution = 0.0;
        *velocity = 0.0;
        *pressure = 0.0;

        uFESpace()->interpolate ( problem()->uexact(), *velocity, currentTime );
        pFESpace()->interpolate ( problem()->pexact(), *pressure, currentTime );
        solution->add ( *velocity );
        solution->add ( *pressure, pressureOffset );

        // Updating bdf
        bdf()->shiftRight ( *solution );
    }
}

} // namespace LifeV
