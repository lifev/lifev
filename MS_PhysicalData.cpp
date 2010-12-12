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
 *  @brief File containing the MultiScale Physical Data
 *
 *  @date 09-09-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/MS_PhysicalData.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
MS_PhysicalData::MS_PhysicalData() :
        M_dataTime                      (),
        M_fluidDensity                  (),
        M_fluidViscosity                (),
        M_fluidReferencePressure        (),
        M_structureDensity              (),
        M_structurePoissonCoefficient   (),
        //M_structureThickness            (),
        M_structureYoungModulus         ()
{
}

MS_PhysicalData::MS_PhysicalData( const MS_PhysicalData& data ) :
        M_dataTime                      ( data.M_dataTime ),
        M_fluidDensity                  ( data.M_fluidDensity ),
        M_fluidViscosity                ( data.M_fluidViscosity ),
        M_fluidReferencePressure        ( data.M_fluidReferencePressure ),
        M_structureDensity              ( data.M_structureDensity ),
        M_structurePoissonCoefficient   ( data.M_structurePoissonCoefficient ),
        //M_structureThickness            ( data.M_structureThickness ),
        M_structureYoungModulus         ( data.M_structureYoungModulus )
{
}

// ===================================================
// Operators
// ===================================================
MS_PhysicalData&
MS_PhysicalData::operator=( const MS_PhysicalData& data )
{
    if ( this != &data )
    {
        M_dataTime                      = data.M_dataTime;
        M_fluidDensity                  = data.M_fluidDensity;
        M_fluidViscosity                = data.M_fluidViscosity;
        M_fluidReferencePressure        = data.M_fluidReferencePressure;
        M_structureDensity              = data.M_structureDensity;
        M_structurePoissonCoefficient   = data.M_structurePoissonCoefficient;
        //M_structureThickness            = data.M_structureThickness;
        M_structureYoungModulus         = data.M_structureYoungModulus;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
void
MS_PhysicalData::readData( const GetPot& dataFile )
{
    M_dataTime.reset( new time_Type( dataFile, "Solver/time_discretization" ) );
    M_fluidDensity                  = dataFile( "Physics/FluidDensity", 0. );
    M_fluidViscosity                = dataFile( "Physics/FluidViscosity", 0. );
    M_fluidReferencePressure        = dataFile( "Physics/FluidReferencePressure", 0. );
    M_structureDensity              = dataFile( "Physics/StructureDensity", 0. );
    M_structurePoissonCoefficient   = dataFile( "Physics/StructurePoissonCoefficient", 0. );
    //M_structureThickness            = dataFile( "Physics/StructureThickness", 0. );
    M_structureYoungModulus         = dataFile( "Physics/StructureYoungModulus", 0. );
}

void
MS_PhysicalData::showMe()
{
    std::cout << "Fluid density                 = " << M_fluidDensity << std::endl
              << "Fluid viscosity               = " << M_fluidViscosity << std::endl
              << "Fluid reference pressure      = " << M_fluidReferencePressure << std::endl << std::endl;

    std::cout << "Structure density coefficient = " << M_structureDensity << std::endl
              << "Structure Poisson coefficient = " << M_structurePoissonCoefficient << std::endl
              //<< "Structure Thickness           = " << M_structureThickness << std::endl
              << "Structure Young modulus       = " << M_structureYoungModulus << std::endl << std::endl;

    std::cout << "Initial time                  = " << M_dataTime->getInitialTime() << std::endl
              << "End time                      = " << M_dataTime->getEndTime() << std::endl
              << "TimeStep                      = " << M_dataTime->getTimeStep() << std::endl << std::endl;
}

// ===================================================
// Get Methods
// ===================================================
MS_PhysicalData::timePtr_Type
MS_PhysicalData::dataTime() const
{
    return M_dataTime;
}

const Real&
MS_PhysicalData::fluidDensity() const
{
    return M_fluidDensity;
}

const Real&
MS_PhysicalData::fluidViscosity() const
{
    return M_fluidViscosity;
}

const Real&
MS_PhysicalData::fluidReferencePressure() const
{
    return M_fluidReferencePressure;
}

const Real&
MS_PhysicalData::structureDensity() const
{
    return M_structureDensity;
}

const Real&
MS_PhysicalData::structurePoissonCoefficient() const
{
    return M_structurePoissonCoefficient;
}

//const Real&
//MS_PhysicalData::GetStructureThickness() const
//{
//    return M_structureThickness;
//}

const Real&
MS_PhysicalData::structureYoungModulus() const
{
    return M_structureYoungModulus;
}

} // Namespace LifeV
