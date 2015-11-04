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
 *  @brief File containing the Multiscale Global Physical Data
 *
 *  @date 09-09-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/multiscale/framework/MultiscaleGlobalData.hpp>

namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleGlobalData::MultiscaleGlobalData() :
    M_timeData                      (),
    M_fluidDensity                  (),
    M_fluidViscosity                (),
    M_fluidVenousPressure           (),
    M_solidExternalPressure         (),
    M_solidDensity                  (),
    M_solidPoissonCoefficient       (),
    //M_solidThickness                (),
    M_solidYoungModulus             (),
    M_scalingFactorResistance       (),
    M_scalingFactorCompliance       ()
{
}

MultiscaleGlobalData::MultiscaleGlobalData ( const MultiscaleGlobalData& data ) :
    M_timeData                      ( data.M_timeData ),
    M_fluidDensity                  ( data.M_fluidDensity ),
    M_fluidViscosity                ( data.M_fluidViscosity ),
    M_fluidVenousPressure           ( data.M_fluidVenousPressure ),
    M_solidExternalPressure         ( data.M_solidExternalPressure ),
    M_solidDensity                  ( data.M_solidDensity ),
    M_solidPoissonCoefficient       ( data.M_solidPoissonCoefficient ),
    //M_solidThickness                ( data.M_solidThickness ),
    M_solidYoungModulus             ( data.M_solidYoungModulus ),
    M_scalingFactorResistance       ( data.M_scalingFactorResistance ),
    M_scalingFactorCompliance       ( data.M_scalingFactorCompliance )
{
}

// ===================================================
// Operators
// ===================================================
MultiscaleGlobalData&
MultiscaleGlobalData::operator= ( const MultiscaleGlobalData& data )
{
    if ( this != &data )
    {
        M_timeData                      = data.M_timeData;
        M_fluidDensity                  = data.M_fluidDensity;
        M_fluidViscosity                = data.M_fluidViscosity;
        M_fluidVenousPressure           = data.M_fluidVenousPressure;
        M_solidExternalPressure         = data.M_solidExternalPressure;
        M_solidDensity                  = data.M_solidDensity;
        M_solidPoissonCoefficient       = data.M_solidPoissonCoefficient;
        //M_solidThickness                = data.M_solidThickness;
        M_solidYoungModulus             = data.M_solidYoungModulus;
        M_scalingFactorResistance       = data.M_scalingFactorResistance;
        M_scalingFactorCompliance       = data.M_scalingFactorCompliance;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
void
MultiscaleGlobalData::readData ( const GetPot& dataFile )
{
    M_timeData.reset ( new time_Type ( dataFile, "Solver/time_discretization" ) );
    M_fluidDensity                  = dataFile ( "Physics/FluidDensity", 0. );
    M_fluidViscosity                = dataFile ( "Physics/FluidViscosity", 0. );
    M_fluidVenousPressure           = dataFile ( "Physics/FluidVenousPressure", 0. );
    M_solidExternalPressure         = dataFile ( "Physics/SolidExternalPressure", 0. );
    M_solidDensity                  = dataFile ( "Physics/SolidDensity", 0. );
    M_solidPoissonCoefficient       = dataFile ( "Physics/SolidPoissonCoefficient", 0. );
    //M_solidThickness                = dataFile( "Physics/SolidThickness", 0. );
    M_solidYoungModulus             = dataFile ( "Physics/SolidYoungModulus", 0. );
    M_scalingFactorResistance       = dataFile ( "Physics/ScalingFactorResistance", 1. );
    M_scalingFactorCompliance       = dataFile ( "Physics/ScalingFactorCompliance", 1. );
}

void
MultiscaleGlobalData::showMe()
{
    std::cout << "Fluid density                 = " << M_fluidDensity << std::endl
              << "Fluid viscosity               = " << M_fluidViscosity << std::endl
              << "Fluid venous pressure         = " << M_fluidVenousPressure << std::endl << std::endl;

    std::cout << "Solid reference pressure      = " << M_solidExternalPressure << std::endl
              << "Solid density coefficient     = " << M_solidDensity << std::endl
              << "Solid Poisson coefficient     = " << M_solidPoissonCoefficient << std::endl
              //<< "Solid Thickness               = " << M_solidThickness << std::endl
              << "Solid Young modulus           = " << M_solidYoungModulus << std::endl << std::endl;

    std::cout << "Resistance scaling factor     = " << M_scalingFactorResistance << std::endl
              << "Compliance scaling factor     = " << M_scalingFactorCompliance << std::endl << std::endl;

    std::cout << "Initial time                  = " << M_timeData->initialTime() << std::endl
              << "End time                      = " << M_timeData->endTime() << std::endl
              << "TimeStep                      = " << M_timeData->timeStep() << std::endl << std::endl;
}

} // Namespace Multiscale
} // Namespace LifeV
