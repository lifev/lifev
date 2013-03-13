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
#include <lifev/bc_interface/3D/function/fluid/BCInterfaceFunctionParserSolverFluid3D.hpp>

namespace LifeV
{

// ===================================================
// Methods
// ===================================================
template< >
void
BCInterfaceFunctionParserSolver< BCHandler, OseenSolver< RegionMesh< LinearTetra > > >::updatePhysicalSolverVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5023 ) << "BCInterfaceFunctionSolver<BCHandler, OseenSolver>::updatePhysicalSolverVariables" << "\n";
#endif

    // Create/Update variables
    for ( std::set< physicalSolverList >::iterator j = M_list.begin(); j != M_list.end(); ++j )
        switch ( *j )
        {
                // f_ -> FLUID
            case f_timeStep:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              f_timeStep(): " << M_physicalSolver->data()->dataTime()->timeStep() << "\n";
#endif
                setVariable ( "f_timeStep", M_physicalSolver->data()->dataTime()->timeStep() );

                break;

            case f_area:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              f_area(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->area ( M_flag ) << "\n";
#endif
                setVariable ( "f_area", M_physicalSolver->area ( M_flag ) );

                break;

            case f_density:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              f_density: " << M_physicalSolver->density() << "\n";
#endif
                setVariable ( "f_density", M_physicalSolver->density() );

                break;

            case f_flux:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              f_flux(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->flux ( M_flag ) << "\n";
#endif

                setVariable ( "f_flux", M_physicalSolver->flux ( M_flag ) );

                break;

            case f_pressure:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              f_pressure(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->pressure ( M_flag ) << "\n";
#endif

                setVariable ( "f_pressure", M_physicalSolver->pressure ( M_flag ) );

                break;

            case f_viscosity:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              f_viscosity: " << M_physicalSolver->viscosity() << "\n";
#endif
                setVariable ( "f_viscosity", M_physicalSolver->viscosity() );

                break;

            default:

                switchErrorMessage ( "OSEEN" );

                break;
        }
}

template< >
void
BCInterfaceFunctionParserSolver< BCHandler, OseenSolverShapeDerivative< RegionMesh< LinearTetra > > >::updatePhysicalSolverVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5023 ) << "BCInterfaceFunctionSolver<BCHandler, OseenSolverShapeDerivative>::updatePhysicalSolverVariables" << "\n";
#endif

    // Create/Update variables
    for ( std::set< physicalSolverList >::iterator j = M_list.begin(); j != M_list.end(); ++j )
        switch ( *j )
        {
                // f_ -> FLUID
            case f_timeStep:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              f_timeStep(): " << M_physicalSolver->data()->dataTime()->timeStep() << "\n";
#endif
                setVariable ( "f_timeStep", M_physicalSolver->data()->dataTime()->timeStep() );

                break;

            case f_area:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              f_area(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->area ( M_flag ) << "\n";
#endif
                setVariable ( "f_area", M_physicalSolver->area ( M_flag ) );

                break;

            case f_density:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              f_density(): " << M_physicalSolver->density() << "\n";
#endif
                setVariable ( "f_density", M_physicalSolver->density() );

                break;

            case f_flux:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              f_flux(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->flux ( M_flag ) << "\n";
#endif

                setVariable ( "f_flux", M_physicalSolver->flux ( M_flag ) );

                break;

            case f_pressure:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              f_pressure(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->pressure ( M_flag ) << "\n";
#endif

                setVariable ( "f_pressure", M_physicalSolver->pressure ( M_flag ) );

                break;

            case f_viscosity:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              f_viscosity(): " << M_physicalSolver->viscosity() << "\n";
#endif
                setVariable ( "f_viscosity", M_physicalSolver->viscosity() );

                break;

            default:

                switchErrorMessage ( "OSEENSHAPEDERIVATIVE" );

                break;
        }
}

// ===================================================
// Protected Methods
// ===================================================
template< >
void
BCInterfaceFunctionParserSolver< BCHandler, OseenSolver< RegionMesh< LinearTetra > > >::createAccessList ( const BCInterfaceData& data )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5023 ) << "BCInterfaceFunctionSolver<BCHandler, OseenSolver>::createAccessList( data )" << "\n";
#endif

    std::map< std::string, physicalSolverList > mapList;

    createFluidMap ( mapList );
    createList ( mapList, data );

    if ( M_physicalSolver.get() )
    {
        updatePhysicalSolverVariables();
    }
}

template< >
void
BCInterfaceFunctionParserSolver< BCHandler, OseenSolverShapeDerivative< RegionMesh< LinearTetra > > >::createAccessList ( const BCInterfaceData& data )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5023 ) << "BCInterfaceFunctionSolver<BCHandler, OseenSolverShapeDerivative>::createAccessList( data )" << "\n";
#endif

    std::map< std::string, physicalSolverList > mapList;

    createFluidMap ( mapList );
    createList ( mapList, data );

    if ( M_physicalSolver.get() )
    {
        updatePhysicalSolverVariables();
    }
}

} // Namespace LifeV
