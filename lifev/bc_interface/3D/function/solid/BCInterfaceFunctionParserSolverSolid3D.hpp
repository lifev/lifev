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

#ifndef BCInterfaceFunctionParserSolverSolid3D_H
#define BCInterfaceFunctionParserSolverSolid3D_H 1

// Structural solver includes
#include <lifev/structure/solver/StructuralOperator.hpp>

// BCInterface includes
#include <lifev/bc_interface/3D/function/solid/BCInterfaceFunctionParserSolid3D.hpp>
#include <lifev/bc_interface/core/function/BCInterfaceFunctionParserSolver.hpp>

namespace LifeV
{

// ===================================================
// Methods
// ===================================================
template< >
void
BCInterfaceFunctionParserSolver< BCHandler, StructuralOperator<RegionMesh <LinearTetra> > >::updatePhysicalSolverVariables();

// ===================================================
// Protected Methods
// ===================================================
template< >
void
BCInterfaceFunctionParserSolver< BCHandler, StructuralOperator<RegionMesh <LinearTetra> > >::createAccessList ( const std::shared_ptr< BCInterfaceData >& data );

} // Namespace LifeV

#endif /* BCInterfaceFunctionParserSolverSolid3D_H */
