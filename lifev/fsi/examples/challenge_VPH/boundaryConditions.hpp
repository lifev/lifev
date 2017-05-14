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


#include "ud_functions.hpp"
// LifeV includes
#include "lifev/core/LifeV.hpp"
#include "lifev/core/fem/BCHandler.hpp"

// Mathcard includes
#include "lifev/fsi/solver/FSIMonolithicGE.hpp"
#include "lifev/fsi/solver/FSIMonolithicGI.hpp"

#define OUTLET 3
#define INLET 2
#define FLUIDINTERFACE 200
#define OUTERWALL 210
#define SOLIDINTERFACE 200
//#define RING2 22
#define RING 2
#define RING3 3
#define RING4 4
#define RING5 5
#define RING6 6
#define INOUTEDGE 20

namespace LifeV
{

typedef FSIOperator::fluid_Type fluid;
typedef FSIOperator::solid_Type solid;

FSIOperator::fluidBchandlerPtr_Type BCh_harmonicExtension (FSIOperator& _oper)
{

    // Boundary condition for the mesh
    debugStream ( 10000 ) << "Boundary condition for the harmonic extension\n";

    BCFunctionBase bcf (fZero);

    FSISolver::fluidBchandlerPtr_Type BCh_he (new FSIOperator::fluidBchandler_Type );


    //BCh_he->addBC("Base",  INLET,     Essential, Full, bcf,   3);

    if (_oper.data().method() == "monolithicGE")
    {
        debugStream (10000) << "FSIMonolithic GCE harmonic extension\n";
        FSIMonolithicGE* MOper = dynamic_cast<FSIMonolithicGE*> (&_oper);
        MOper->setStructureDispToHarmonicExtension (_oper.lambdaFluidRepeated() );
        BCh_he->addBC ("Interface", SOLIDINTERFACE, Essential, Full,
                       *MOper->bcvStructureDispToHarmonicExtension(), 3);
    }

    return BCh_he;
}


FSIOperator::fluidBchandlerPtr_Type BCh_monolithicFlux()
{
    // Boundary conditions for the fluid velocity

    FSIOperator::fluidBchandlerPtr_Type BCh_fluid ( new FSIOperator::fluidBchandler_Type );

    // BCFunctionBase flow_in (aortaFluxIn);
    // BCFunctionBase flow_3 (linearFlux3_);
    // BCFunctionBase flow_4 (linearFlux4);
    // BCFunctionBase flow_5 (linearFlux5);
    // BCFunctionBase flow_6 (linearFlux6_);
    // BCFunctionBase flow_7 (linearFlux7);
    // BCFunctionBase flow_8 (linearFlux8);
    // BCFunctionBase flow_9 (linearFlux9);

    BCFunctionBase flow_in (abdominalAorta);

    //     BCFunctionBase flow_jean (aortaFluxJean);

    //uncomment  to use fluxes
    BCh_fluid->addBC ("InFlow" , INLET,  Flux, /*Full*/Normal, flow_in);
    // BCh_fluid->addBC("OutFlow" , OUTLET,  Flux/*Essential*/, Normal, flow_3);

    // BCh_fluid->addBC("Flow4" , 4,  Flux/*Essential*/, Normal, flow_4);
    // BCh_fluid->addBC("Flow7" , 7,  Flux/*Essential*/, Normal, flow_7);
    // BCh_fluid->addBC("Flow6" , 6,  Flux/*Essential*/, Normal, flow_6);
    // BCh_fluid->addBC("Flow5" , 5,  Flux/*Essential*/, Normal, flow_5);
    // BCh_fluid->addBC("Flow8" , 8,  Flux/*Essential*/, Normal, flow_8);
    // BCh_fluid->addBC("Flow9" , 9,  Flux/*Essential*/, Normal, flow_9);

    return BCh_fluid;
}

FSIOperator::fluidBchandlerPtr_Type BCh_monolithicFluid (FSIOperator& _oper)
{
    // Boundary conditions for the fluid velocity
    debugStream ( 10000 ) << "Boundary condition for the fluid\n";

    if (! _oper.isFluid() )
    {
        return FSIOperator::fluidBchandlerPtr_Type();
    }

    FSIOperator::fluidBchandlerPtr_Type BCh_fluid ( new FSIOperator::fluidBchandler_Type );

    BCFunctionBase bcf      (fZero);
    //BCFunctionBase pressure_out(FlowConditions::outPressure0);

    BCFunctionBase in_flow  (aortaPhisPress);

    BCFunctionBase out_flow (fZero);


    BCh_fluid->addBC ("InFlow" , INLET,  Natural/*Essential*/, Normal, in_flow/*, 3*/);
    return BCh_fluid;
}


FSIOperator::solidBchandlerPtr_Type BCh_monolithicSolid (FSIOperator& _oper)
{
    // Boundary conditions for the solid displacement

    if (! _oper.isSolid() )
    {
        return FSIOperator::solidBchandlerPtr_Type();
    }

    FSIOperator::solidBchandlerPtr_Type BCh_solid ( new FSIOperator::solidBchandler_Type );

    BCFunctionBase bcf (fZero);


    BCFunctionBase young (E);

    //robin condition on the outer wall
    _oper.setRobinOuterWall (bcf, young);
    BCh_solid->addBC ("OuterWall", OUTERWALL, Robin, Full, _oper.bcfRobinOuterWall(), 3);

    //BCh_solid->addBC("Top",  OUTLET, Essential, Full, bcf, 3);
    //BCh_solid->addBC("Base",  INLET, Essential, Full, bcf, 3);

    return BCh_solid;
}


}

#endif
