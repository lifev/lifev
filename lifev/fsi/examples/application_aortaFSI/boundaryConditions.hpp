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

//#ifndef BL
//Fluid Mesh
#define OUTLET 5
#define OUTLETRING 50
#define INLET 6
#define INLETRING 60

#define RCCA 8 //Right common carotid artery
#define RCCARING 80
#define RCSUBCA 2 //Right subclavial artery
#define RCSUBCARING 20
#define RCSUBSUBCA 9 //Right sub-subclavial artery
#define RCSUBSUBCARING 90
#define LSUBCA 4 //Left subclavial artery
#define LSSUBCARING 40
#define LSUBSUBCA 3 //Left sub-subclavial artery
#define LSSUBSUBCARING 30
#define LCCA 7 //Left common carotid artery
#define LCCARING 70
#define FLUIDINTERFACE 200

//Solid Mesh
#define OUTLETWALL 5
#define OUTLETWALLRINGIN 50
#define OUTLETWALLRINGOUT 51
#define INLETWALL 6
#define INLETWALLRINGIN 60
#define INLETWALLRINGOUT 61

#define RCCAWALL 8 //Right common carotid artery
#define RCCAWALLRINGIN 80
#define RCCAWALLRINGOUT 81
#define RCSUBCAWALL 2 //Right subclavial artery
#define RCSUBCAWALLRINGIN 20
#define RCSUBCAWALLRINGOUT 21
#define RCSUBSUBCAWALL 9 //Right sub-subclavial artery
#define RCSUBSUBCAWALLRINGIN 90
#define RCSUBSUBCAWALLRINGOUT 91
#define LSUBCAWALL 4 //Left subclavial artery
#define LSSUBCAWALLRINGIN 40
#define LSSUBCAWALLRINGOUT 41
#define LSUBSUBCAWALL 3 //Left sub-subclavial artery
#define LSSUBSUBCAWALLRINGIN 30
#define LSSUBSUBCAWALLRINGOUT 31
#define LCCAWALL 7 //Left common carotid artery
#define LCCAWALLRINGIN 70
#define LCCAWALLRINGOUT 70
#define SOLIDINTERFACE 200
#define OUTERWALL 210

//#endif

/*
//#define RING2 22
#define RING 2
//thoracic aorta,
#define RING3 3
//21, L. Brachia, bhanch 3_2
#define RING4 4
//first branchstd::placeholders::_1,
#define RING5 5
//branch 1_2 smallest
#define RING6 6
//R. Brachia, branch 1_3
#define RING7 7
// 15, LCCA, branch 2
#define RING8 8
// 20 LVA branch 3_1
#define RING9 9
#define INOUTEDGE 20
#else
#define OUTLET 5
#define INLET 6
#define FLUIDINTERFACE 200
#define OUTERWALL 201
#define SOLIDINTERFACE 200
//#define RING2 22
#define RING4 4
#define RING5 8
#define RING6 9
#define RING7 2
#define RING8 7
#define RING6 4
//#define INOUTEDGE 20
*/

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


    BCh_he->addBC ("Base",  INLET,     Essential, Full, bcf,   3);

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

    BCFunctionBase flow_in (aortaFluxIn);
    BCFunctionBase flow_3 (linearFlux3_);
    BCFunctionBase flow_4 (linearFlux4);
    BCFunctionBase flow_5 (linearFlux5);
    BCFunctionBase flow_6 (linearFlux6_);
    BCFunctionBase flow_7 (linearFlux7);
    BCFunctionBase flow_8 (linearFlux8);
    BCFunctionBase flow_9 (linearFlux9);

    //     BCFunctionBase flow_jean (aortaFluxJean);

    //uncomment  to use fluxes
    BCh_fluid->addBC ("InFlow" , INLET,  Flux, /*Full/**/Normal, flow_in);
    BCh_fluid->addBC ("OutFlow" , OUTLET,  Flux/*Essential*/, Normal, flow_3);

    //BCh_fluid->addBC("Flow4" , 4,  Flux/*Essential*/, Normal, flow_4);
    BCh_fluid->addBC ("Flow4" , LSUBCA,  Flux/*Essential*/, Normal, flow_4);
    BCh_fluid->addBC ("Flow7" , RCSUBCA,  Flux/*Essential*/, Normal, flow_7);
    BCh_fluid->addBC ("Flow6" , RCSUBSUBCA,  Flux/*Essential*/, Normal, flow_6);
    BCh_fluid->addBC ("Flow5" , RCCA,  Flux/*Essential*/, Normal, flow_5);
    BCh_fluid->addBC ("Flow8" , LCCA,  Flux/*Essential*/, Normal, flow_8);
    BCh_fluid->addBC ("Flow9" , LSUBSUBCA,  Flux/*Essential*/, Normal, flow_9);

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


    //BCh_fluid->addBC("InFlow" , INLET,  Natural/*Essential*/, Normal, in_flow/*, 3*/);
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
    BCh_solid->addBC ("OuterWall", OUTERWALL, Robin, Normal, _oper.bcfRobinOuterWall() );

    //BCh_solid->addBC("Top",  OUTLET, Essential, Full, bcf, 3);
    BCh_solid->addBC ("Base",  INLET, Essential, Full, bcf, 3);

    return BCh_solid;
}


}

#endif
