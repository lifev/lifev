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
#include <lifev/bc_interface/function/BCInterfaceFunctionParserSolverFSI3D.hpp>

namespace LifeV
{

// ===================================================
// Methods
// ===================================================
template< >
void
BCInterfaceFunctionParserSolver< FSIOperator >::updatePhysicalSolverVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5023 ) << "BCInterfaceFunctionSolver<FSIOperator>::updatePhysicalSolverVariables" << "\n";
#endif

    // Create/Update variables
    for ( std::set< physicalSolverList >::iterator j = M_list.begin(); j != M_list.end(); ++j )
        switch ( *j )
        {
                // f_ -> FLUID
            case f_timeStep:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              f_timeStep(): " << M_physicalSolver->data().dataFluid()->dataTime()->timeStep() << "\n";
#endif
                setVariable ( "f_timeStep", M_physicalSolver->data().dataFluid()->dataTime()->timeStep() );

                break;

            case f_area:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              f_area(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->fluid().area ( M_flag ) << "\n";
#endif
                setVariable ( "f_area", M_physicalSolver->fluid().area ( M_flag ) );

                break;

            case f_density:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              f_density: " << M_physicalSolver->fluid().density() << "\n";
#endif
                setVariable ( "f_density", M_physicalSolver->fluid().density() );

                break;

            case f_flux:

                if ( M_physicalSolver->isFluid() )
                {
#ifdef HAVE_LIFEV_DEBUG
                    debugStream ( 5023 ) << "!!! Warning: fluid not initialized yet, setting flux = 0 in BCInterface !!!\n";

#endif
                    setVariable ( "f_flux", 0.0 );
                }
                else
                {
#ifdef HAVE_LIFEV_DEBUG
                    debugStream ( 5023 ) << "                                              f_flux(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->fluid().flux ( M_flag, *M_physicalSolver->fluid().solution() ) << "\n";
#endif
                    setVariable ( "f_flux", M_physicalSolver->fluid().flux ( M_flag, *M_physicalSolver->fluid().solution() ) );
                }

                break;

            case f_pressure:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              f_pressure(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->fluid().pressure ( M_flag, *M_physicalSolver->fluid().solution() ) << "\n";
#endif

                setVariable ( "f_pressure", M_physicalSolver->fluid().pressure ( M_flag, *M_physicalSolver->fluid().solution() ) );

                break;

            case f_viscosity:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              f_viscosity: " << M_physicalSolver->fluid().viscosity() << "\n";
#endif
                setVariable ( "f_viscosity", M_physicalSolver->fluid().viscosity() );

                break;

                // s_ -> SOLID
            case s_density:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              s_density: " << M_physicalSolver->solid().rho() << "\n";
#endif

                setVariable ( "s_density", M_physicalSolver->solid().rho() );

                break;

            case s_poisson:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              s_poisson: " << M_physicalSolver->solid().poisson() << "\n";
#endif

                setVariable ( "s_poisson", M_physicalSolver->solid().poisson (1) );

                break;

            case s_thickness:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              s_thickness: " << M_physicalSolver->solid().thickness() << "\n";
#endif

                setVariable ( "s_thickness", M_physicalSolver->solid().thickness() );

                break;

            case s_young:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              s_young: " << M_physicalSolver->solid().young() << "\n";
#endif

                setVariable ( "s_young", M_physicalSolver->solid().young (1) );

                break;

            case s_externalPressure:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              s_externalPressure: " << M_physicalSolver->solid().data()->externalPressure() << "\n";
#endif

                setVariable ( "s_externalPressure", M_physicalSolver->solid().data()->externalPressure() );

                break;

            default:

                switchErrorMessage ( "FSIOperator" );

                break;
        }
}



// ===================================================
// Protected Methods
// ===================================================
template< >
void
BCInterfaceFunctionParserSolver< FSIOperator >::createAccessList ( const BCInterfaceData& data )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5023 ) << "BCInterfaceFunctionSolver<FSIOperator>::createAccessList( data )" << "\n";
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
