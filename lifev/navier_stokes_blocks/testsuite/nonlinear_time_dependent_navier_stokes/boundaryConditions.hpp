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

#define OUTLET     5
#define INLET_UP   3
#define INLET_DOWN 4
#define WALL       210
#define INTERFACE  200
#define CONSTRAINT 6

namespace LifeV
{

typedef boost::shared_ptr<BCHandler> bcPtr_Type;

bcPtr_Type BCh_fluid ()
{
    BCFunctionBase zero_function (fZero);
    BCFunctionBase inflow_function_up (inflow_up);
    BCFunctionBase inflow_function_down (inflow_down);

    bcPtr_Type bc (new BCHandler );

    bc->addBC ("InletUp",  	 INLET_UP,    Essential,    Full,   inflow_function_up,     3);
    bc->addBC ("InletDown",  INLET_DOWN,  Essential, 	Full,   inflow_function_down,   3);
    bc->addBC ("Wall",       WALL,        Essential,    Full,   zero_function,          3);
    bc->addBC ("Interface",  INTERFACE,   Essential,    Full,   zero_function,          3);

    return bc;
}

bcPtr_Type BCh_PCD ()
{
    BCFunctionBase zero_function (fZero);

    bcPtr_Type bc (new BCHandler );

    bc->addBC ("Outflow", OUTLET, Essential, Full, zero_function, 3);

    return bc;
}

}

#endif
