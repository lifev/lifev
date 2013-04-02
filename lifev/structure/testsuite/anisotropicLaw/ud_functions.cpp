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
 *  @brief File containing the boundary conditions for the Monolithic Test
 *
 *  @date 2009-04-09
 *  @author Paolo Tricerri <paolo.tricerri@epfl.ch>
 *
 *  @maintainer Paolo Tricerri <paolo.tricerri@epfl.ch>
 *
 *  Contains the functions to be assigned as boundary conditions, in the file boundaryConditions.hpp . The functions
 *  can depend on time and space, while they can take in input an ID specifying one of the three principal axis
 *  if the functions to assign is vectorial and the boundary condition is of type \c Full \c.
 */

#include "ud_functions.hpp"

namespace LifeV
{

Real f (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    switch (i)
    {
        case 0:
            return 0.0;
            break;
        case 1:
            return 0.0;
            break;
        case 2:
            return 0.0;
            break;
        default:
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }
}


Real fzero_scalar (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}

Real InternalPressure (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -1e+5;
    //return -260000*sin(80*3.141592*t);
}

// Initial displacement and velocity
Real d0 (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{

    switch (i)
    {
    case 0:
        return  0.088002 * ( x + 0.5 );
        break;
    case 1:
        return - ( 0.02068 * 2.0 ) * ( y );
        break;
    case 2:
        return - ( 0.02068 * 2.0 ) * ( z );
        break;
    default:
        ERROR_MSG ("This entry is not allowed: ud_functions.hpp");
        return 0.;
        break;
    }

}

Real w0 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{

    switch (i)
    {
        case 0:
            return 0.0;
            break;
        case 1:
            return 0.0;
            break;
        case 2:
            return 0.0;
            break;
        default:
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }
}

Real a0 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{

    switch (i)
    {
        case 0:
            return 0.0;
            break;
        case 1:
            return 0.0;
            break;
        case 2:
            return 0.0;
            break;
        default:
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }
}


//----------------------------------------------Boundary Conditions--------------
Real bcZero (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
    return  0.;
}

Real bcNonZero (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
    return  100000.;
}

Real bcNonZeroSecondOrderExponential (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
    return  19180.;
}


//----------------------------------------------Fibers Directions--------------
Real Family1 ( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{
    switch (i)
    {
        case 0:
            return 1.0;
            break;
        case 1:
            return 0.0;
            break;
        case 2:
            return 0.0;
            break;
        default:
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }
}

Real Family2 ( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{

    switch (i)
    {
        case 0:
            return 0.0;
            break;
        case 1:
            return 1.0;
            break;
        case 2:
            return 0.0;
            break;
        default:
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }
}

Real Family3 ( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{

    switch (i)
    {
        case 0:
            return 0.0;
            break;
        case 1:
            return 0.0;
            break;
        case 2:
            return 0.0;
            break;
        default:
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }
}


Real Family4 ( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{

    switch (i)
    {
        case 0:
            return 0.0;
            break;
        case 1:
            return 0.0;
            break;
        case 2:
            return 0.0;
            break;
        default:
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }
}


Real Family5 ( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{

    switch (i)
    {
        case 0:
            return 0.0;
            break;
        case 1:
            return 0.0;
            break;
        case 2:
            return 0.0;
            break;
        default:
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }
}

Real Family6 ( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{

    switch (i)
    {
        case 0:
            return 0.0;
            break;
        case 1:
            return 0.0;
            break;
        case 2:
            return 0.0;
            break;
        default:
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }
}

// Method for the definition of the fibers
fibersDirectionList::fibersDirectionList() :
    M_mapNameDefinition( )
{}

fibersDirectionList::~fibersDirectionList()
{}

void fibersDirectionList::setupFiberDefinitions( const UInt nbFamilies )
{
    // At the moment the creation of the table of fiber functions is done
    // manually. There should be a way to make it automatically. Btw, only
    // the first nbFamilies that are set in the data file are taken into account

    ASSERT( nbFamilies < 6, "At the moment, a maximum number = 6 of families can be used! If you want more \n modifiy the file ud_functions.hpp in the application folder." );

    // Creation of the database of functions
    fiberFunctionPtr_Type pointerToFunction( new fiberFunction_Type( Family1 ) );
    M_mapNameDefinition.insert( std::pair<std::string, fiberFunctionPtr_Type>
                                  ( "Family1", pointerToFunction ) );

    pointerToFunction.reset( new fiberFunction_Type( Family2 ) );
    M_mapNameDefinition.insert( std::pair<std::string, fiberFunctionPtr_Type>
                                  ( "Family2", pointerToFunction ) );

    pointerToFunction.reset( new fiberFunction_Type( Family3 ) );
    M_mapNameDefinition.insert( std::pair<std::string, fiberFunctionPtr_Type>
                                  ( "Family3", pointerToFunction ) );

    pointerToFunction.reset( new fiberFunction_Type( Family4 ) );
    M_mapNameDefinition.insert( std::pair<std::string, fiberFunctionPtr_Type>
                                  ( "Family4", pointerToFunction ) );

    pointerToFunction.reset( new fiberFunction_Type( Family5 ) );
    M_mapNameDefinition.insert( std::pair<std::string, fiberFunctionPtr_Type>
                                  ( "Family5", pointerToFunction ) );

    pointerToFunction.reset( new fiberFunction_Type( Family6 ) );
    M_mapNameDefinition.insert( std::pair<std::string, fiberFunctionPtr_Type>
                                  ( "Family6", pointerToFunction ) );


}

fibersDirectionList::fiberFunctionPtr_Type fibersDirectionList::fiberDefinition( const std::string nameFamily )
{

    mapNameDefinitionFiberFunction_Type::const_iterator IT;

    IT = M_mapNameDefinition.find ( nameFamily );

    if ( IT != M_mapNameDefinition.end() )
    {
        return IT->second;
    }
    else
    {
        std::cout << " Wrong identification of the fiber function! " << std::endl;
        fiberFunctionPtr_Type pointerToFunction( new fiberFunction_Type() );

        return pointerToFunction;
    }
}

}
