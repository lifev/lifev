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
#include "life/lifecore/life.hpp"
#include "life/lifefem/bcHandler.hpp"
#include "life/lifefem/bcFunction.hpp"

// Mathcard includes
#include "lifemc/lifesolver/FSIMonolithicGE.hpp"
#include "lifemc/lifesolver/FSIMonolithicGI.hpp"

#include "flowConditions.hpp"
//#include "lumpedHeart.hpp"
#include "ud_functions.hpp"

#define OUTLET 3
#define INLET 2
#define FLUIDINTERFACE 1
#define SOLIDINTERFACE 1
#define OUTERWALL 10
#define RING  2
#define RING2 3
#define INOUTEDGE 20
#define INEDGE 30

namespace LifeV
{

typedef FSI::fluid_Type fluid;
typedef FSI::solid_Type solid;

FSI::fluidBchandlerPtr_Type BCh_harmonicExtension(FSI &_oper)
{

    // Boundary condition for the mesh
    Debug( 10000 ) << "Boundary condition for the harmonic extension\n";

    BCFunctionBase bcf(fZero);

    FSISolver::fluidBchandlerPtr_Type BCh_he(new FSI::fluidBchandler_Type );

    BCh_he->addBC("Edges", INOUTEDGE, Essential, Full, bcf,   3);
    BCh_he->addBC("Edges", INEDGE, Essential, Full, bcf,   3);
    BCh_he->addBC("Base",  INLET,     Essential, Full, bcf,   3);

    if (_oper.data().method() == "monolithicGE")
    {
        Debug(10000) << "FSIMonolithic GCE harmonic extension\n";
        FSIMonolithicGE *MOper = dynamic_cast<FSIMonolithicGE *>(&_oper);
        MOper->setStructureDispToHarmonicExtension(_oper.lambdaFluidRepeated());
        BCh_he->addBC("Interface", SOLIDINTERFACE, Essential, Full,
                      *MOper->bcvStructureDispToHarmonicExtension(), 3);
    }
    else if (_oper.data().method() == "monolithicGI")
    {

    }

    return BCh_he;
}


FSI::fluidBchandlerPtr_Type BCh_monolithicFlux(bool /*isOpen=true*/)
{
    FSI::fluidBchandlerPtr_Type BCh_fluid( new FSI::fluidBchandler_Type );

    BCFunctionBase flow_3 (fluxFunction);
    BCFunctionBase bcf      (fZero);
    //uncomment  to use fluxes

    //  BCh_fluid->addBC("InFlow" , INLET,  Flux, Normal, flow_3);
//   if(!isOpen)
//       BCh_fluid->addBC("InFlow" , INLET,  Flux,   Normal, bcf);

    //uncomment  to use fluxes
    BCh_fluid->addBC("InFlow" , INLET,  Flux, Normal, flow_3);

    return BCh_fluid;
}

FSI::fluidBchandlerPtr_Type BCh_monolithicFluid(FSI &_oper, bool const & /*isOpen=true*/)
{
    // Boundary conditions for the fluid velocity
    Debug( 10000 ) << "Boundary condition for the fluid\n";

    if (! _oper.isFluid() )
        return FSI::fluidBchandlerPtr_Type();

    FSI::fluidBchandlerPtr_Type BCh_fluid( new FSI::fluidBchandler_Type );

    BCFunctionBase bcf      (fZero);
    BCFunctionBase in_flow  (/*uInterpolated*/u2normal/*aortaPhisPress*/);
    //    BCFunctionBase out_flow (fZero);
    //BCFunctionBase in_flow  (LumpedHeart::outPressure);

    BCFunctionBase out_press (FlowConditions::outPressure0);


//     if(isOpen)
    //BCh_fluid->addBC("InFlow" , INLET,  Natural,   Normal, in_flow);

    BCh_fluid->addBC("OutFlow", OUTLET,  Natural,  Normal, out_press);
    return BCh_fluid;
}

FSI::solidBchandlerPtr_Type BCh_monolithicRobin(FSI &_oper)
{

    FSI::solidBchandlerPtr_Type BCh_solid( new FSI::solidBchandler_Type );
    BCFunctionBase hyd(fZero);
    BCFunctionBase young (E);

    //robin condition on the outer wall
    _oper.setRobinOuterWall(hyd, young);
    BCh_solid->addBC("OuterWall", OUTERWALL, Robin, Normal, _oper.bcfRobinOuterWall());

    return BCh_solid;
}

FSI::solidBchandlerPtr_Type BCh_monolithicSolid(FSI &_oper)
{

    if (! _oper.isSolid() )
        return FSI::solidBchandlerPtr_Type();

    // Boundary conditions for the solid displacement
    Debug( 10000 ) << "Boundary condition for the solid\n";
    FSI::solidBchandlerPtr_Type BCh_solid( new FSI::solidBchandler_Type );

    BCFunctionBase bcf(fZero);
    //BCFunctionBase young (E);

    BCh_solid->addBC("Top",   RING, Essential, Full, bcf,  3);
    BCh_solid->addBC("Base",  RING2, Essential, Full, bcf,  3);

    BCh_solid->addBC("OuterWall", OUTERWALL, Natural, Full, bcf,  3);
    BCh_solid->addBC("Edges", INOUTEDGE, Essential, Full, bcf,  3);
    BCh_solid->addBC("Edges", INEDGE, Essential, Full, bcf,  3);

    return BCh_solid;
}

}

#endif
