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

// BCInterface includes
#include <lifev/bc_interface/1D/function/BCInterfaceFunctionParserSolver1D.hpp>

namespace LifeV
{

// ===================================================
// Methods
// ===================================================
template< >
void
BCInterfaceFunctionParserSolver< OneDFSIBCHandler, OneDFSISolver >::updatePhysicalSolverVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5023 ) << "BCInterfaceFunctionSolver<FSI>::updatePhysicalSolverVariables" << "\n";
#endif

    OneDFSI::bcSide_Type side = ( M_boundaryID == 0 ) ? OneDFSI::left : OneDFSI::right;
    // Create/Update variables
    for ( std::set< physicalSolverList >::iterator j = M_list.begin(); j != M_list.end(); ++j )
        switch ( *j )
        {
                // f_ -> FLUID
            case f_timeStep:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              f_timeStep(): " << M_physicalSolver->physics()->data()->dataTime()->timeStep() << "\n";
#endif
                setVariable ( "f_timeStep", M_physicalSolver->physics()->data()->dataTime()->timeStep() );

                break;

            case f_area:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              f_area(" << static_cast<Real> (side) << "): " << M_physicalSolver->boundaryValue ( *M_solution, OneDFSI::A, side ) << "\n";
#endif
                setVariable ( "f_area", M_physicalSolver->boundaryValue ( *M_solution, OneDFSI::A, side ) );

                break;

            case f_density:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              f_density: " << M_physicalSolver->physics()->data()->densityRho() << "\n";
#endif
                setVariable ( "f_density", M_physicalSolver->physics()->data()->densityRho() );

                break;

            case f_flux:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              f_flux(" << static_cast<Real> (side) << "): " << M_physicalSolver->boundaryValue ( *M_solution, OneDFSI::Q, side ) << "\n";
#endif

                setVariable ( "f_flux", M_physicalSolver->boundaryValue ( *M_solution, OneDFSI::Q, side ) );

                break;

            case f_pressure:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              f_pressure(" << static_cast<Real> (side) << "): " << M_physicalSolver->boundaryValue ( *M_solution, OneDFSI::P, side ) << "\n";
#endif

                setVariable ( "f_pressure", M_physicalSolver->boundaryValue ( *M_solution, OneDFSI::P, side ) );

                break;

            case f_viscosity:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              f_viscosity: " << M_physicalSolver->physics()->data()->viscosity() << "\n";
#endif
                setVariable ( "f_viscosity", M_physicalSolver->physics()->data()->viscosity() );

                break;

            case f_venousPressure:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              f_venousPressure: " << M_physicalSolver->physics()->data()->venousPressure() << "\n";
#endif
                setVariable ( "f_venousPressure", M_physicalSolver->physics()->data()->venousPressure() );

                break;

                // s_ -> SOLID
            case s_density:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              s_density: " << M_physicalSolver->physics()->data()->densityWall() << "\n";
#endif

                setVariable ( "s_density", M_physicalSolver->physics()->data()->densityWall() );

                break;

            case s_poisson:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              s_poisson: " << M_physicalSolver->physics()->data()->poisson() << "\n";
#endif

                setVariable ( "s_poisson", M_physicalSolver->physics()->data()->poisson() );

                break;

            case s_thickness:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              s_thickness: " << M_physicalSolver->physics()->data()->thickness ( M_physicalSolver->boundaryDOF ( side ) ) << "\n";
#endif

                setVariable ( "s_thickness", M_physicalSolver->physics()->data()->thickness ( M_physicalSolver->boundaryDOF ( side ) ) );

                break;

            case s_young:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              s_young: " << M_physicalSolver->physics()->data()->young() << "\n";
#endif

                setVariable ( "s_young", M_physicalSolver->physics()->data()->young() );

                break;

            case s_externalPressure:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              s_externalPressure: " << M_physicalSolver->physics()->data()->externalPressure() << "\n";
#endif

                setVariable ( "s_externalPressure", M_physicalSolver->physics()->data()->externalPressure() );

                break;

            default:
                switchErrorMessage ( "OneDFSIModel_Solver" );

                break;
        }
}



// ===================================================
// Protected Methods
// ===================================================
template< >
void
BCInterfaceFunctionParserSolver< OneDFSIBCHandler, OneDFSISolver >::createAccessList ( const std::shared_ptr< BCInterfaceData >& data )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5023 ) << "BCInterfaceFunctionSolver<OneDimensionaSolver>::createAccessList( data )" << "\n";
#endif

    std::map< std::string, physicalSolverList > mapList;

    createFluidMap ( mapList );
    createSolidMap ( mapList );
    createList ( mapList, data );

    if ( M_physicalSolver.get() )
    {
        updatePhysicalSolverVariables();
    }
}

} // Namespace LifeV
