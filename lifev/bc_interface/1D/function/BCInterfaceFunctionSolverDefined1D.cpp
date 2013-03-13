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
 *  @brief File containing the BCInterfaceFunctionSolverDefined class and specializations
 *
 *  @date 23-04-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/bc_interface/1D/function/BCInterfaceFunctionSolverDefined1D.hpp>

namespace LifeV
{

// ===================================================
// Constructors
// ===================================================
BCInterfaceFunctionSolverDefined< OneDFSIBCHandler, OneDFSISolver >::BCInterfaceFunctionSolverDefined() :
    M_defaultFunction (),
    M_function        ()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5025 ) << "BCInterfaceFunctionSolverDefined::BCInterfaceFunctionSolverDefined()" << "\n";
#endif

}

// ===================================================
// Methods
// ===================================================
void
BCInterfaceFunctionSolverDefined< OneDFSIBCHandler, OneDFSISolver >::assignFunction ( OneDFSIFunction& base )
{
    switch ( M_defaultFunction )
    {
        case Riemann:

            base.setFunction ( boost::bind ( &OneDFSIFunctionSolverDefinedRiemann::operator(),
                                             dynamic_cast<OneDFSIFunctionSolverDefinedRiemann*> ( & ( *M_function ) ), _1, _2 ) );

            break;

        case Compatibility:

            base.setFunction ( boost::bind ( &OneDFSIFunctionSolverDefinedCompatibility::operator(),
                                             dynamic_cast<OneDFSIFunctionSolverDefinedCompatibility*> ( & ( *M_function ) ), _1, _2 ) );

            break;

        case Absorbing:

            base.setFunction ( boost::bind ( &OneDFSIFunctionSolverDefinedAbsorbing::operator(),
                                             dynamic_cast<OneDFSIFunctionSolverDefinedAbsorbing*> ( & ( *M_function ) ), _1, _2 ) );

            break;

        case Resistance:

            base.setFunction ( boost::bind ( &OneDFSIFunctionSolverDefinedResistance::operator(),
                                             dynamic_cast<OneDFSIFunctionSolverDefinedResistance*> ( & ( *M_function ) ), _1, _2 ) );

            break;
    }
}

// ===================================================
// Set Methods
// ===================================================
void
BCInterfaceFunctionSolverDefined< OneDFSIBCHandler, OneDFSISolver >::setData ( const BCInterfaceData1D& data )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5025 ) << "BCInterfaceFunctionSolverDefined::setData( data )" << "\n";
#endif

    //Set mapFunction
    std::map< std::string, solverDefinedFunctions > mapFunction;
    mapFunction["Riemann"]       = Riemann;
    mapFunction["Compatibility"] = Compatibility;
    mapFunction["Absorbing"]     = Absorbing;
    mapFunction["Resistance"]    = Resistance;

    M_defaultFunction = mapFunction[data.baseString()];

    switch ( M_defaultFunction )
    {
        case Riemann:

            M_function.reset ( new OneDFSIFunctionSolverDefinedRiemann ( data.side(), data.quantity() ) );

            break;

        case Compatibility:

            M_function.reset ( new OneDFSIFunctionSolverDefinedCompatibility ( data.side(), data.quantity() ) );

            break;

        case Absorbing:

            M_function.reset ( new OneDFSIFunctionSolverDefinedAbsorbing ( data.side(), data.quantity() ) );

            break;

        case Resistance:

            M_function.reset ( new OneDFSIFunctionSolverDefinedResistance ( data.side(), data.quantity(), data.resistance() [0] ) );

            break;
    }
}

} // Namespace LifeV
