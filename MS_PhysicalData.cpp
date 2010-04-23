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
    M_DataTime      (),
    M_FluidDensity  (),
    M_FluidViscosity()
{
}

MS_PhysicalData::MS_PhysicalData( const MS_PhysicalData& PhysicalData ) :
    M_DataTime      ( PhysicalData.M_DataTime ),
    M_FluidDensity  ( PhysicalData.M_FluidDensity ),
    M_FluidViscosity( PhysicalData.M_FluidViscosity )
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
        M_DataTime       = PhysicalData.M_DataTime;
        M_FluidDensity   = PhysicalData.M_FluidDensity;
        M_FluidViscosity = PhysicalData.M_FluidViscosity;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
void
MS_PhysicalData::ReadData( const GetPot& dataFile )
{
    M_DataTime.reset( new Time_Type( dataFile, "Solver/time_discretization" ) );
    M_FluidDensity   = dataFile( "Physics/fluidDensity", 0. );
    M_FluidViscosity = dataFile( "Physics/fluidViscosity", 0. );
}

void
MS_PhysicalData::ShowMe()
{
    std::cout << "Fluid density       = " << M_FluidDensity << std::endl
              << "Fluid viscosity     = " << M_FluidViscosity << std::endl << std::endl;

    std::cout << "Initial time        = " << M_DataTime->getInitialTime() << std::endl
              << "End time            = " << M_DataTime->getEndTime() << std::endl
              << "TimeStep            = " << M_DataTime->getTimeStep() << std::endl << std::endl;
}

// ===================================================
// Get Methods
// ===================================================
MS_PhysicalData::Time_ptrType
MS_PhysicalData::GetDataTime() const
{
    return M_DataTime;
}

const Real&
MS_PhysicalData::GetFluidDensity() const
{
    return M_FluidDensity;
}

const Real&
MS_PhysicalData::GetFluidViscosity() const
{
    return M_FluidViscosity;
}

} // Namespace LifeV
