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
    @file
    @brief File containing a base class for the source function of the 1D hyperbolic problem.

    @date 15-04-2010
    @author Cristiano Malossi <cristiano.malossi@epfl.ch>

    @contributor Simone Rossi <simone.rossi@epfl.ch>

    @mantainer  Cristiano Malossi <cristiano.malossi@epfl.ch>
 */
#include "OneDimensionalModel_Source.hpp"

namespace LifeV
{

std::map< std::string, OneDimensionalModel_SourceTypes > OneDimensionalModel_SourceMap;

// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalModel_Source::OneDimensionalModel_Source() :
    M_physics   ()
{}

OneDimensionalModel_Source::OneDimensionalModel_Source( const physicsPtr_Type Physics ) :
    M_physics   ( Physics )
{}

// ===================================================
// Set Methods
// ===================================================
void
OneDimensionalModel_Source::SetPhysics( const physicsPtr_Type& Physics )
{
    M_physics = Physics;
}

// ===================================================
// Get Methods
// ===================================================
OneDimensionalModel_Source::physicsPtr_Type
OneDimensionalModel_Source::Physics() const
{
    return M_physics;
}

}
