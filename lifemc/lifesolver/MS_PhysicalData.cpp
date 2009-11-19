//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

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
 *  @brief MultiScale Physical Data
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 09-09-2009
 */

#include <lifemc/lifesolver/MS_PhysicalData.hpp>

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
MS_PhysicalData::MS_PhysicalData() :
    M_fluidDensity  (),
    M_fluidViscosity()
{
}

MS_PhysicalData::MS_PhysicalData( const MS_PhysicalData& PhysicalData ) :
    M_fluidDensity  ( PhysicalData.M_fluidDensity ),
    M_fluidViscosity( PhysicalData.M_fluidViscosity )
{
}

// ===================================================
// Operators
// ===================================================
MS_PhysicalData&
MS_PhysicalData::operator=( const MS_PhysicalData& PhysicalData )
{
    if ( this != &PhysicalData )
    {
        M_fluidDensity   = PhysicalData.M_fluidDensity;
        M_fluidViscosity = PhysicalData.M_fluidViscosity;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
void
MS_PhysicalData::ReadData( const GetPot& dataFile )
{
    M_fluidDensity   = dataFile( "Physics/fluidDensity", 0. );
    M_fluidViscosity = dataFile( "Physics/fluidViscosity", 0. );
}

void
MS_PhysicalData::ShowMe()
{
    std::cout << "Fluid density     = " << M_fluidDensity << std::endl
              << "Fluid viscosity   = " << M_fluidViscosity << std::endl << std::endl;
}

// ===================================================
// Get Methods
// ===================================================
const Real&
MS_PhysicalData::GetFluidDensity() const
{
    return M_fluidDensity;
}

const Real&
MS_PhysicalData::GetFluidViscosity() const
{
    return M_fluidViscosity;
}

} // Namespace LifeV
