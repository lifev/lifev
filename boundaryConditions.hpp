
/**
   \file boundaryConditions.hpp
   \author Paolo Crosetto <paolo.crosetto@epfl.ch>
   \date 2009-04-09
*/

/**

contains the functions to be assigned as boundary conditions, in the file boundaryConditions.hpp . The functions
can depend on time and space, while they can take in input an ID specifying one of the three principal axis
if the functions to
assign is vectorial and the boundary condition is of type \c Full \c.
 */

#ifndef BC_HPP
#define BC_HPP

#include "life/lifecore/life.hpp"
#include "flowConditions.hpp"
//#include "lumpedHeart.hpp"
#include "ud_functions.hpp"
#include "life/lifefem/bcHandler.hpp"
#include "life/lifefem/bcFunction.hpp"

#include "lifemc/lifesolver/Monolithic.hpp"
#include "lifemc/lifesolver/fullMonolithic.hpp"



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

typedef FSIOperator::fluid_raw_type fluid;
typedef FSIOperator::solid_raw_type solid;

FSIOperator::fluid_bchandler_type BCh_harmonicExtension(FSIOperator &_oper)
{

    // Boundary condition for the mesh
    Debug( 10000 ) << "Boundary condition for the harmonic extension\n";

    BCFunctionBase bcf(fZero);

    FSISolver::fluid_bchandler_type BCh_he(new FSIOperator::fluid_bchandler_raw_type );

    BCh_he->addBC("Edges", INOUTEDGE, Essential, Full, bcf,   3);
    BCh_he->addBC("Edges", INEDGE, Essential, Full, bcf,   3);
    BCh_he->addBC("Base",  INLET,     Essential, Full, bcf,   3);

    if (_oper.data().method() == "monolithic")
    {
        Debug(10000) << "Monolithic GCE harmonic extension\n";
        Monolithic *MOper = dynamic_cast<Monolithic *>(&_oper);
        MOper->setStructureDispToHarmonicExtension(_oper.lambdaFluidRepeated());
        BCh_he->addBC("Interface", SOLIDINTERFACE, Essential, Full,
                      *MOper->bcvStructureDispToHarmonicExtension(), 3);
    }
    else if (_oper.data().method() == "fullMonolithic")
    {

        BCh_he->addBC("Interface", SOLIDINTERFACE, Essential, Full,
                      bcf, 3);
    }

    return BCh_he;
}


FSIOperator::fluid_bchandler_type BCh_monolithicFlux(bool isOpen=true)
{
    FSIOperator::fluid_bchandler_type BCh_fluid( new FSIOperator::fluid_bchandler_raw_type );

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

FSIOperator::fluid_bchandler_type BCh_monolithicFluid(FSIOperator &_oper, bool const & isOpen=true)
{
    // Boundary conditions for the fluid velocity
    Debug( 10000 ) << "Boundary condition for the fluid\n";

    if (! _oper.isFluid() )
        return FSIOperator::fluid_bchandler_type();

    FSIOperator::fluid_bchandler_type BCh_fluid( new FSIOperator::fluid_bchandler_raw_type );

    BCFunctionBase bcf      (fZero);
    BCFunctionBase in_flow  (/*uInterpolated*/u2/*aortaPhisPress*/);
    //    BCFunctionBase out_flow (fZero);
    //BCFunctionBase in_flow  (LumpedHeart::outPressure);

    BCFunctionBase out_press (FlowConditions::outPressure0);


//     if(isOpen)
    BCh_fluid->addBC("InFlow" , INLET,  Natural,   Normal, in_flow);

    BCh_fluid->addBC("OutFlow", OUTLET,  Natural,  Normal, out_press);
    return BCh_fluid;
}

FSIOperator::solid_bchandler_type BCh_monolithicRobin(FSIOperator &_oper)
{

    FSIOperator::solid_bchandler_type BCh_solid( new FSIOperator::solid_bchandler_raw_type );
    BCFunctionBase hyd(fZero);
    BCFunctionBase young (E);

    //robin condition on the outer wall
    _oper.setMixteOuterWall(hyd, young);
    BCh_solid->addBC("OuterWall", OUTERWALL, Mixte, Normal, _oper.bcfMixteOuterWall());

    return BCh_solid;
}

FSIOperator::solid_bchandler_type BCh_monolithicSolid(FSIOperator &_oper)
{

    if (! _oper.isSolid() )
        return FSIOperator::solid_bchandler_type();

    // Boundary conditions for the solid displacement
    Debug( 10000 ) << "Boundary condition for the solid\n";
    FSIOperator::solid_bchandler_type BCh_solid( new FSIOperator::solid_bchandler_raw_type );

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
