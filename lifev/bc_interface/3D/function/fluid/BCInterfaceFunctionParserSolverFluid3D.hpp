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
 *  @file
 *  @brief File containing the BCInterfaceFunctionParserSolver class
 *
 *  @date 24-08-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterfaceFunctionParserSolverFluid3D_H
#define BCInterfaceFunctionParserSolverFluid3D_H 1

// Oseen includes
#include <lifev/navier_stokes/solver/OseenSolverShapeDerivative.hpp>

// BCInterface includes
#include <lifev/bc_interface/function/BCInterfaceFunctionParserSolver.hpp>

namespace LifeV
{

// ===================================================
// Methods
// ===================================================
template< >
void
BCInterfaceFunctionParserSolver< OseenSolver< RegionMesh< LinearTetra > > >::updatePhysicalSolverVariables();

template< >
void
BCInterfaceFunctionParserSolver< OseenSolverShapeDerivative< RegionMesh< LinearTetra > > >::updatePhysicalSolverVariables();

// ===================================================
// Protected Methods
// ===================================================
template< >
void
BCInterfaceFunctionParserSolver< OseenSolver< RegionMesh< LinearTetra > > >::createAccessList ( const BCInterfaceData& data );

template< >
void
BCInterfaceFunctionParserSolver< OseenSolverShapeDerivative< RegionMesh< LinearTetra > > >::createAccessList ( const BCInterfaceData& data );

} // Namespace LifeV

#endif /* BCInterfaceFunctionParserSolverFluid3D_H */
