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
 *  @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
 *
 *  @contributor Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @maintainer Paolo Crosetto <crosetto@iacspc70.epfl.ch>
 *
 *  Contains the functions to be assigned as boundary conditions, in the file boundaryConditions.hpp . The functions
 *  can depend on time and space, while they can take in input an ID specifying one of the three principal axis
 *  if the functions to assign is vectorial and the boundary condition is of type \c Full \c.
 */

#ifndef BC_HPP
#define BC_HPP

// LifeV includes
#include <lifev/core/LifeV.hpp>
#include <lifev/core/fem/BCHandler.hpp>

#include <lifev/fsi/solver/FSIMonolithicGE.hpp>
#include <lifev/fsi/solver/FSIMonolithicGI.hpp>

#include "flowConditions.hpp"
#include "resistance.hpp"
#include "ud_functions.hpp"


//Regular mesh
#define INLET 2
#define INLETRING 20
#define OUTLET 3
#define OUTLETRING 30
#define FLUIDINTERFACE 100

#define SOLIDINTERFACE 100
#define INLETWALL 2
#define INLETWALL_INTRING 20

#define OUTLETWALL 3
#define OUTLETWALL_INTRING 30

#define OUTERWALL 1000



namespace LifeV
{

typedef FSIOperator::fluid_Type fluid;
typedef FSIOperator::solid_Type solid;

FSIOperator::fluidBchandlerPtr_Type BCh_harmonicExtension (FSIOperator& /*_oper*/)
{

    // Boundary condition for the mesh
    debugStream ( 10000 ) << "Boundary condition for the harmonic extension\n";

    BCFunctionBase bcf (fZero);

    FSISolver::fluidBchandlerPtr_Type BCh_he (new FSIOperator::fluidBchandler_Type );

    BCh_he->addBC ("in", INLET, Essential, Full, bcf,   3);
    BCh_he->addBC ("in", OUTLET, Essential, Full, bcf,   3);

    // Rings of the fluid domain
    // BCh_he->addBC ("inRing", INLETRING,  EssentialVertices, Full, bcf,   3);
    // BCh_he->addBC ("inRing", OUTLETRING, EssentialVertices, Full, bcf,   3);

    return BCh_he;
}


FSIOperator::fluidBchandlerPtr_Type BCh_monolithicFlux (bool /*isOpen=true*/)
{
    FSIOperator::fluidBchandlerPtr_Type BCh_fluid ( new FSIOperator::fluidBchandler_Type );

    BCFunctionBase flowAneurysm (fluxFunctionAneurysm);
    BCFunctionBase bcf      (fZero);
    //uncomment  to use fluxes

    //BCh_fluid->addBC("InFlow" , INLET,  Flux, Normal, flowAneurysm);
    //   if(!isOpen)
    //       BCh_fluid->addBC("InFlow" , INLET,  Flux,   Normal, bcf);

    //uncomment  to use fluxes
    //BCh_fluid->addBC("InFlow" , INLET,  Flux, Normal, flowAneurysm);

    return BCh_fluid;
}

FSIOperator::fluidBchandlerPtr_Type BCh_monolithicFluid (FSIOperator& _oper, bool const& /*isOpen=true*/)
{
    // Boundary conditions for the fluid velocity
    debugStream ( 10000 ) << "Boundary condition for the fluid\n";

    if (! _oper.isFluid() )
    {
        return FSIOperator::fluidBchandlerPtr_Type();
    }

    FSIOperator::fluidBchandlerPtr_Type BCh_fluid ( new FSIOperator::fluidBchandler_Type );

    BCFunctionBase bcf      (fZero);
    BCFunctionBase in_flow  (uInterpolated);
    //    BCFunctionBase out_flow (fZero);

    BCFunctionBase out_press3 (ResistanceBCs::outPressure0);

    BCFunctionBase InletVect (aneurismFluxInVectorial);
    //BCFunctionBase bcfw0 (w0);

    //Inlets
    BCh_fluid->addBC ("InFlow" , INLET,  Essential, Full, InletVect, 3);

    //Outlets

    //Absorbing BC seemed not to work
    //Absorbing BC on outlet 2and3 caused instabilities
    BCh_fluid->addBC ("out3", OUTLET, Natural,  Normal, out_press3);
    //BCh_fluid->addBC("out3", OUTLET, Natural,  Normal, bcf);

    return BCh_fluid;
}

FSIOperator::solidBchandlerPtr_Type BCh_monolithicSolid (FSIOperator& _oper)
{

    if (! _oper.isSolid() )
    {
        return FSIOperator::solidBchandlerPtr_Type();
    }

    // Boundary conditions for the solid displacement
    debugStream ( 10000 ) << "Boundary condition for the solid\n";
    FSIOperator::solidBchandlerPtr_Type BCh_solid ( new FSIOperator::solidBchandler_Type );

    BCFunctionBase bcf (fZero);

    //Inlets & Outlets
    BCh_solid->addBC ("BORDERS",   INLETWALL, Essential, Full, bcf,  3);
    //BCh_solid->addBC ("BORDERS-RIN",   INLETWALL_INTRING, EssentialVertices, Full, bcf,  3);
    BCh_solid->addBC ("BORDERS",   OUTLETWALL, Essential, Full, bcf,  3);
    //BCh_solid->addBC ("BORDERS-rin",   OUTLETWALL_INTRING, EssentialVertices, Full, bcf,  3);

    //Robin BC
    BCFunctionBase hyd (fZero);
    BCFunctionBase young (E);
    BCFunctionBase externalPressure (outerWallPressure);
    //robin condition on the outer wall
    _oper.setRobinOuterWall (externalPressure, young);

    //BCh_solid->addBC ("OuterWall", OUTERWALL, Robin, Normal, _oper.bcfRobinOuterWall() );
    //First try: Homogeneous Neumann
    //BCh_solid->addBC ("OuterWall", OUTERWALL, Natural, Normal, bcf);

    //Constant pressure  Neumann
    BCh_solid->addBC ("OuterWall", OUTERWALL, Natural, Normal, externalPressure);

    return BCh_solid;
}

}

#endif
