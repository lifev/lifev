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

#ifndef BCFSITUBE_HPP
#define BCFSITUBE_HPP

// LifeV includes
#include <lifev/core/LifeV.hpp>
#include <lifev/core/fem/BCHandler.hpp>

#include "ud_functions.hpp"

#define OUTLET 3
#define INLET  2
#define WALL   1
#define INOUTEDGE 20
#define OUTERWALL 10
#define FLUIDINTERFACE 100

namespace LifeV
{

typedef boost::shared_ptr<BCHandler> bcPtr_Type;

bcPtr_Type BCh_fluid ()
{
    BCFunctionBase zero_function (fZero);
    BCFunctionBase inflow_function (inflow);

    bcPtr_Type bc (new BCHandler );

    //bc->addBC ("Inflow",  	 INLET,      Essential,      Full,   inflow_function, 3);
    // bc->addBC ("Inflow",  	 INLET,      Natural,        Normal, inflow_function, 3);
    bc->addBC ("INOUTEDGE",  INOUTEDGE,  EssentialEdges, Full,   zero_function,   3);
    bc->addBC ("Outflow",    OUTLET,     Natural,        Normal, zero_function);

    return bc;
}

bcPtr_Type BCh_fluid_residual ()
{
    BCFunctionBase zero_function (fZero);
    BCFunctionBase inflow_function (inflow);
    BCFunctionBase pressure_wave (pressure);

    bcPtr_Type bc (new BCHandler );

    bc->addBC ("Inflow",  	 INLET,      Natural,        Normal, pressure_wave);
    //bc->addBC ("Inflow",  	 INLET,      Essential,      Full,   zero_function, 3);
    bc->addBC ("INOUTEDGE",  INOUTEDGE,  EssentialEdges, Full,   zero_function, 3);
    bc->addBC ("Outflow",    OUTLET,     Natural,        Normal, zero_function);

    return bc;
}

bcPtr_Type BCh_structure ()
{
    BCFunctionBase zero_function (fZero);
    BCFunctionBase inflow_function (inflow);

    bcPtr_Type bc (new BCHandler );

    bc->addBC ("Inflow",     INLET,  	 Essential,      Full, zero_function, 3);
    bc->addBC ("Outflow",    OUTLET, 	 Essential,      Full, zero_function, 3);

    return bc;
}

bcPtr_Type BCh_ale ()
{
    BCFunctionBase zero_function (fZero);

    bcPtr_Type bc_ale (new BCHandler );

    bc_ale->addBC ("Inflow",  INLET,     Essential, Full, zero_function, 3);
    bc_ale->addBC ("Outflow", OUTLET,    Essential, Full, zero_function, 3);
    //bc_ale->addBC ("Interface",    1,    Essential, Full, zero_function, 3);
    bc_ale->addBC ("Gamma",  FLUIDINTERFACE, Essential, Full, zero_function,   3);

    return bc_ale;
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
