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
 *  @author Davide Forti <davide.forti@epfl.ch>
 *
 *  Contains the functions to be assigned as boundary conditions, in the file boundaryConditions.hpp . The functions
 *  can depend on time and space, while they can take in input an ID specifying one of the three principal axis
 *  if the functions to assign is vectorial and the boundary condition is of type \c Full \c.
 */

#ifndef BCNS_HPP
#define BCNS_HPP

// LifeV includes
#include <lifev/core/LifeV.hpp>
#include <lifev/core/fem/BCHandler.hpp>

#include "ud_functions.hpp"

namespace LifeV
{

typedef boost::shared_ptr<BCHandler> bcPtr_Type;

bcPtr_Type BCh_fluid ()
{
    BCFunctionBase zero (fZero);

    bcPtr_Type bc (new BCHandler );

    bc->addBC( "Inflow",         3, Essential,      Full,        zero, 3 );
    bc->addBC( "Walls",        200, Essential, 	    Full, 	     zero, 3 );
    bc->addBC( "Outflow2",       2, Natural, 	  Normal,    	 zero    );
    bc->addBC( "Outflow4",       4, Essential,      Full,        zero, 3 );
    bc->addBC( "Outflow5",       5, Essential,      Full,        zero, 3 );
    bc->addBC( "Outflow6",       6, Essential,      Full,        zero, 3 );
    bc->addBC( "Outflow7",       7, Essential,      Full,        zero, 3 );
    bc->addBC( "Outflow8",       8, Essential,      Full,        zero, 3 );
    bc->addBC( "Outflow9",       9, Essential,      Full,        zero, 3 );

    return bc;
}

bcPtr_Type BCh_preprocessing ()
{
    BCFunctionBase zero (fZero);

    bcPtr_Type bc (new BCHandler );

    bc->addBC( "Walls", 200, Essential, Full, zero, 1 );

    return bc;
}

}

#endif
