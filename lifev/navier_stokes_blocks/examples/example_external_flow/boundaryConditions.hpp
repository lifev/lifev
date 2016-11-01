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
	BCFunctionBase uZero( zeroFunction );
	BCFunctionBase uInflow( inflowFunction );

	std::vector<LifeV::ID> xComp(1), yComp(1), zComp(1);
	xComp[0] = 0;
	yComp[0] = 1;
	zComp[0] = 2;

    bcPtr_Type bcH (new BCHandler );

    bcH->addBC( "Inflow",       1, Essential, Full,      	uInflow, 3 );
    bcH->addBC( "Outflow",      3, Natural,   Full,      	uZero,   3 );
    bcH->addBC( "TopBottom",    2, Essential, Component, 	uZero,   yComp  );
    bcH->addBC( "LeftRight",    5, Essential, Component, 	uZero,   zComp  );
    bcH->addBC ("Cylinder" ,    4, Essential, Full,    	    uZero,	 3 );

    return bcH;
}

bcPtr_Type BCh_drag ()
{
	bcPtr_Type bcH (new BCHandler );

	BCFunctionBase uOneX( oneFunctionX );

	bcH->addBC( "CylinderDrag",   4, Essential, Full, uOneX, 3 );

	return bcH;
}

bcPtr_Type BCh_lift ()
{
	bcPtr_Type bcH (new BCHandler );

	BCFunctionBase uOneY( oneFunctionY );

	bcH->addBC( "CylinderLift",   4, Essential, Full, uOneY, 3 );

	return bcH;
}

}

#endif
