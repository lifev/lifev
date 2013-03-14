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
#include <lifev/bc_interface/3D/function/fsi/BCInterfaceFunctionParserFSI3D.hpp>

namespace LifeV
{

// ===================================================
// Methods
// ===================================================
template< >
void
BCInterfaceFunctionParser< BCHandler, FSIOperator >::assignFunction ( bcBase_Type& base )
{
    base.setFunction ( functionSelectorTimeSpaceID() );
}

// ===================================================
// Set Methods
// ===================================================
template< >
void
BCInterfaceFunctionParser< BCHandler, FSIOperator >::setData ( const boost::shared_ptr< BCInterfaceData >& data )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5022 ) << "BCInterfaceFunction::setData" << "\n";
#endif

    setupParser ( data );

    boost::shared_ptr< BCInterfaceData3D > castedData = boost::dynamic_pointer_cast< BCInterfaceData3D > ( data );

    if ( castedData != 0 )
    {

        /*
         * MODE          COMPONENT     FUNCTION      |      COMV.SIZE     ARGUMENTS     INTERFACEFUNCTION
         * ------------------------------------------|---------------------------------------------------
         *                                           |
         * COMPONENT     2             x*y*z         |      1             1             function
         * FULL          3             x*y*z         |      1             1             function
         * FULL          1             x*y*z         |      1             1             function
         * FULL          3             (y*z,x*z,x*y) |      1             3             functionID
         * FULL          2             (x,y)         |      1             2             functionID
         * COMPONENT     '1 3'         (x,y)         |      2             2             functionID
         */

    #ifdef HAVE_LIFEV_DEBUG
        debugStream ( 5021 ) << "BCInterfaceFunction::setData                arguments: " << M_parser->countSubstring ( "," ) << "\n";
    #endif

        // Note: the map ID is used only for 3D handler.
        if ( M_parser->countSubstring ( "," ) )
        {
            //Create the ID map
            if ( castedData->componentsVector().size() > 1 ) // Component
                for ( ID i ( 0 ); i < static_cast< ID > ( castedData->componentsVector().size() ); ++i )
                {
                    M_mapID[castedData->componentsVector() [i]] = i + 1;
                }
            else
                // if ( castedData->componentsVector().front() == arguments )  Full
                for ( ID i ( 0 ); i < castedData->componentsVector().front(); ++i )
                {
                    M_mapID[i + 1] = i + 1;
                }
        }
    }
    else
        std::cerr << "!!! ERROR: BCInterface wrong data cast !!!" << std::endl;
}

} // Namespace LifeV
