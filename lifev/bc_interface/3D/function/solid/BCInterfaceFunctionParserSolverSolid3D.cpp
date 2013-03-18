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
#include <lifev/bc_interface/3D/function/solid/BCInterfaceFunctionParserSolverSolid3D.hpp>

namespace LifeV
{

// ===================================================
// Methods
// ===================================================
template< >
void
BCInterfaceFunctionParserSolver< BCHandler, StructuralOperator<RegionMesh <LinearTetra> > >::updatePhysicalSolverVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5023 ) << "BCInterfaceFunctionSolver<BCHandler, StructuralOperator>::updatePhysicalSolverVariables" << "\n";
#endif

    // Create/Update variables
    for ( std::set< physicalSolverList >::iterator j = M_list.begin(); j != M_list.end(); ++j )
        switch ( *j )
        {
                // s_ -> SOLID
            case s_density:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              s_density: " << M_physicalSolver->data()->rho() << "\n";
#endif

                setVariable ( "s_density", M_physicalSolver->data()->rho() );

                break;

            case s_poisson:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              s_poisson: " << M_physicalSolver->data()->poisson (1) << "\n";
#endif

                setVariable ( "s_poisson", M_physicalSolver->data()->poisson (1) );

                break;

            case s_thickness:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              s_thickness: " << M_physicalSolver->data()->thickness() << "\n";
#endif

                setVariable ( "s_thickness", M_physicalSolver->data()->thickness() );

                break;

            case s_young:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              s_young: " << M_physicalSolver->data()->young (1) << "\n";
#endif

                setVariable ( "s_young", M_physicalSolver->data()->young (1) );

                break;

            case s_externalPressure:

#ifdef HAVE_LIFEV_DEBUG
                debugStream ( 5023 ) << "                                              s_externalPressure: " << M_physicalSolver->data()->externalPressure() << "\n";
#endif

                setVariable ( "s_externalPressure", M_physicalSolver->data()->externalPressure() );

                break;

            default:

                switchErrorMessage ( "StructuralOperator" );

                break;
        }
}



// ===================================================
// Protected Methods
// ===================================================
template< >
void
BCInterfaceFunctionParserSolver< BCHandler, StructuralOperator<RegionMesh <LinearTetra> > >::createAccessList ( const boost::shared_ptr< BCInterfaceData >& data )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5023 ) << "BCInterfaceFunctionSolver<BCHandler, StructuralOperator>::createAccessList( data )" << "\n";
#endif

    std::map< std::string, physicalSolverList > mapList;

    createSolidMap ( mapList );
    createList ( mapList, data );

    if ( M_physicalSolver.get() )
    {
        updatePhysicalSolverVariables();
    }
}

} // Namespace LifeV
