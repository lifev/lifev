//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief File containing a base class for the source function of the 1D hyperbolic problem.
 *
 *  @version 1.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 15-04-2010
 */

#include "OneDimensionalModel_Source.hpp"

namespace LifeV {

std::map< std::string, OneDimensionalModel_SourceTypes > OneDimensionalModel_SourceMap;

// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalModel_Source::OneDimensionalModel_Source() :
    M_Physics   ()
{}

OneDimensionalModel_Source::OneDimensionalModel_Source( const Physics_PtrType Physics ) :
    M_Physics   ( Physics )
{}

// ===================================================
// Set Methods
// ===================================================
void
OneDimensionalModel_Source::SetPhysics( const Physics_PtrType& Physics )
{
    M_Physics = Physics;
}

// ===================================================
// Get Methods
// ===================================================
OneDimensionalModel_Source::Physics_PtrType
OneDimensionalModel_Source::Physics() const
{
    return M_Physics;
}

}
