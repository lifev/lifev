/* -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-

  <short description here>

  Gilles Fourestey gilles.fourestey@epfl.ch

*/
/** \file boundaryConditions.hpp
*/

#ifndef BC_HPP
#define BC_HPP

#include "life/lifecore/life.hpp"
#include "ud_functions.hpp"
#include "life/lifefem/bcHandler.hpp"
#include "life/lifefem/bcFunction.hpp"
//#include "life/lifesolver/reducedLinFluid.hpp"

// #include "life/lifesolver/steklovPoincareBase.hpp"
#include "lifemc/lifesolver/Monolithic.hpp"
#include "lifemc/lifesolver/fullMonolithic.hpp"
#include "lifemc/lifesolver/fixedPointBase.hpp"
//#define BLUFF_CYLINDER
//#define ANEURISM
//#define ANEURISM100170
//#define AORTA
#ifdef ANEURISM        // ifdef ANEURISM

#ifdef ANEURISM100170       // ifdef ANEURISM100170

#define INLET 1
#define OUTLET  2
#define OUTLET2 3
#define OUTLET3 4
#define OUTLET4 5
#define FLUIDINTERFACE 6
#define SOLIDINTERFACE 6

#define OUTERWALL 10
#define RING 0

#else                      // ifdef ANEURISM100170

#define OUTLET 3
#define INLET 2
#define FLUIDINTERFACE 1
#define SOLIDINTERFACE 1
#define INOUTEDGE 20

#define OUTERWALL 10
#define RING 0

#endif                     // ifdef ANEURISM100170

#else                  // ifdef ANEURISM
#ifdef AORTA
#define OUTLET 3
#define INLET 2
#define FLUIDINTERFACE 1
#define OUTERWALL 9
#define SOLIDINTERFACE 1
#define RING2 22
#define RING 2
#define RING3 3
#define RING4 4
#define RING5 5
#define RING6 6
#define RING7 7
#define RING8 8
#else

#ifdef BLUFF_CYLINDER
#define OUTLET 4
#else
#define OUTLET 3
#endif

#define INLET 2
#define FLUIDINTERFACE 1
#define SOLIDINTERFACE 1
#define INOUTEDGE 20

#define OUTERWALL 10
#define RING  2
#define RING2 3
#endif                 // ifdef AORTA
#endif                 // ifdef ANEURISM

namespace LifeV
{

typedef FSIOperator::fluid_raw_type fluid;
typedef FSIOperator::solid_raw_type solid;

FSIOperator::fluid_bchandler_type BCh_harmonicExtension(FSIOperator &_oper)
{

//     Debug(10000) << "SP harmonic extension\n";
//     fixedPoint *FPOper = dynamic_cast<fixedPoint *>(&_oper);

//     FPOper->setStructureDispToHarmonicExtension(_oper.lambdaFluid());

//    if (! _oper.isFluid() )
//        return FSIOperator::fluid_bchandler_type();

//    FPOper->bcvStructureDispToHarmonicExtension()->showMe(true,std::cout);

    // Boundary condition for the mesh
    Debug( 10000 ) << "Boundary condition for the harmonic extension\n";

    BCFunctionBase bcf(fZero);

    FSISolver::fluid_bchandler_type BCh_he(new FSIOperator::fluid_bchandler_raw_type );

    BCh_he->addBC("Top",   OUTLET,    Essential, Full, bcf,   3);
#ifdef ANEURISM100170
    BCh_he->addBC("Top",   OUTLET2,    Essential, Full, bcf,   3);
    BCh_he->addBC("Top",   OUTLET3,    Essential, Full, bcf,   3);
    BCh_he->addBC("Top",   OUTLET4,    Essential, Full, bcf,   3);
#else
#ifdef BLUFF_CYLINDER
#else
#ifdef AORTA
#else
    BCh_he->addBC("Edges", INOUTEDGE, Essential, Full, bcf,   3);
#endif
#endif
#endif
    BCh_he->addBC("Base",  INLET,     Essential, Full, bcf,   3);


    if (_oper.method() == "steklovPoincare")
    {
//         Debug(10000) << "SP harmonic extension\n";
//         steklovPoincare *SPOper = dynamic_cast<steklovPoincare *>(&_oper);
//         SPOper->setFluidInterfaceDisp((LifeV::Vector&) _oper.lambdaFluidRepeated());
//         BCh_he->addBC("Interface", 1, Essential, Full,
//                       *SPOper->bcvFluidInterfaceDisp(), 3);
    }
    else if (_oper.method() == "monolithic")
    {
        Debug(10000) << "EJ harmonic extension\n";
        Monolithic *EJOper = dynamic_cast<Monolithic *>(&_oper);
        EJOper->setStructureDispToHarmonicExtension(_oper.lambdaFluidRepeated());
        BCh_he->addBC("Interface", SOLIDINTERFACE, Essential, Full,
        *EJOper->bcvStructureDispToHarmonicExtension(), 3);
    }
   else if (_oper.method() == "fullMonolithic")
    {
        //        Monolithic *EJOper = dynamic_cast<Monolithic *>(&_oper);
        BCh_he->addBC("Interface", SOLIDINTERFACE, Essential, Full,
        bcf, 3);
    }
     else if (_oper.method() == "fixedPoint")
    {
        Debug(10000) << "FP harmonic extension\n";
        fixedPoint *FPOper = dynamic_cast<fixedPoint *>(&_oper);

        FPOper->setStructureDispToHarmonicExtension(_oper.lambdaFluidRepeated());
        BCh_he->addBC("Interface", FLUIDINTERFACE, Essential, Full,
                      *FPOper->bcvStructureDispToHarmonicExtension(), 3);
    }



    return BCh_he;
}


FSIOperator::fluid_bchandler_type BCh_monolithicFluid(FSIOperator &_oper)
{
    // Boundary conditions for the fluid velocity
    Debug( 10000 ) << "Boundary condition for the fluid\n";

    if (! _oper.isFluid() )
        return FSIOperator::fluid_bchandler_type();

    FSIOperator::fluid_bchandler_type BCh_fluid( new FSIOperator::fluid_bchandler_raw_type );

    BCFunctionBase bcf      (fZero);
#ifdef ANEURISM
    BCFunctionBase in_flow  (uInterpolated);
#else
#ifdef AORTA
    BCFunctionBase in_flow  (/*aortaPhisPress*/u2/*uInterpolated*/);
#else
    BCFunctionBase in_flow  (aortaPhisPress);
#endif
#endif
    BCFunctionBase out_flow (fZero);


#ifdef ANEURISM
    BCh_fluid->addBC("InFlow" , INLET,  Essential,   Full, in_flow, 3);
#else
#ifdef AORTA
    BCh_fluid->addBC("InFlow" , INLET,  Natural/*Essential*/, Normal, in_flow, 3);
#else
    BCh_fluid->addBC("InFlow" , INLET,  Natural,   Normal, in_flow/*, 2*/);
#endif
#endif
    BCh_fluid->addBC("OutFlow", OUTLET,  Natural,   Full, out_flow, 3);
#ifdef ANEURISM100170
    BCh_fluid->addBC("OutFlow", OUTLET2,  Natural,   Full, out_flow, 3);
    BCh_fluid->addBC("OutFlow", OUTLET3,  Natural,   Full, out_flow, 3);
    BCh_fluid->addBC("OutFlow", OUTLET4,  Natural,   Full, out_flow, 3);
#endif


    //    _oper.setHarmonicExtensionVelToFluid(_oper.veloFluidMesh());

#ifndef ANEURISM100170
    //    BCh_fluid->addBC("Edges",  INOUTEDGE, Essential, Full,
    //                     *_oper.bcvHarmonicExtensionVelToFluid(),  3);
#endif
    //    BCh_fluid->addBC("Interface",   FLUIDINTERFACE,  Essential, Full,
    //                     *_oper.bcvHarmonicExtensionVelToFluid(),  3);

    return BCh_fluid;
}

FSIOperator::fluid_bchandler_type BCh_monolithicFluidPrec(FSIOperator &_oper)
{
    // Boundary conditions for the fluid velocity
    Debug( 10000 ) << "Boundary condition for the fluid\n";

    if (! _oper.isFluid() )
        return FSIOperator::fluid_bchandler_type();

    FSIOperator::fluid_bchandler_type BCh_fluid( new FSIOperator::fluid_bchandler_raw_type );

    BCFunctionBase bcf      (fZero);
#ifdef ANEURISM
    BCFunctionBase in_flow  (uInterpolated);
#else
#ifdef AORTA
    BCFunctionBase in_flow  (/*u2*/uInterpolated);
#else
    BCFunctionBase in_flow  (u2);
#endif
#endif
    BCFunctionBase out_flow (fZero);


#ifdef ANEURISM
    BCh_fluid->addBC("InFlow" , INLET,  Essential,   Full, bcf, 2);
#else
    //    BCh_fluid->addBC("InFlow" , INLET,  Natural,   Full, in_flow, 3);
    //#endif
    //    BCh_fluid->addBC("OutFlow", OUTLET,  Natural,   Full, out_flow, 3);
    //#ifdef ANEURISM100170
    //    BCh_fluid->addBC("OutFlow", OUTLET2,  Natural,   Full, out_flow, 3);
    //    BCh_fluid->addBC("OutFlow", OUTLET3,  Natural,   Full, out_flow, 3);
    //    BCh_fluid->addBC("OutFlow", OUTLET4,  Natural,   Full, out_flow, 3);
#endif


    //    _oper.setHarmonicExtensionVelToFluid(_oper.veloFluidMesh());

    //#ifndef ANEURISM100170
    //    BCh_fluid->addBC("Edges",  INOUTEDGE, Essential, Full,
    //                     bcf,  3);
    //#endif
    //    BCh_fluid->addBC("Interface",   FLUIDINTERFACE,  Essential, Full,
    //                     *_oper.bcvHarmonicExtensionVelToFluid(),  3);
    //    if (_oper.method() == "monolithic")
    //    {
    //    BCh_fluid->addBC("Interface", FLUIDINTERFACE, Essential  , Full,
    //                     bcf, 3);
            //    }
    return BCh_fluid;
}


FSIOperator::solid_bchandler_type BCh_monolithicSolidPrec(FSIOperator &_oper)
{
    if (! _oper.isFluid() )
        return FSIOperator::solid_bchandler_type();

    // Boundary conditions for the fluid velocity
    Debug( 10000 ) << "Boundary condition for the linearized fluid\n";
    FSIOperator::solid_bchandler_type BCh_solid( new FSIOperator::fluid_bchandler_raw_type );

    BCFunctionBase bcf(fZero);
#ifdef ANEURISM
    //    BCFunctionBase in_flow  (uInterpolated);
#else
    //BCFunctionBase in_flow  (u2);
#endif
    BCFunctionBase in_flow  (fZero);
    BCFunctionBase out_flow (fZero);

    BCh_solid->addBC("Top",   RING, Essential, Full, bcf,  3);
#ifdef ANEURISM
#else
    BCh_solid->addBC("Base",  RING2, Essential, Full, bcf,  3);
#endif
    BCh_solid->addBC("OuterWall", OUTERWALL, Natural, Full, bcf,  3);
#ifndef ANEURISM100170
#ifdef BLUFF_CYLINDER
#else
#ifdef AORTA
#else
    BCh_solid->addBC("Edges", INOUTEDGE, Essential, Full, bcf,  3);
#endif
#endif
#endif
    //in BCh_fluid. Now it is set to 0 because the mesh displacement is zero in this part of the boundary

//    BCh_fluidLin->addBC("interface",  1,  Essential,   Full, bcf,     3);


//    if (_oper.method() == "monolithic")
//    {
//            BCh_solid->addBC("Interface", FLUIDINTERFACE, Essential  , Full,
//                                bcf, 3);
            //    }

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

    //BCh_solid->addBC("Top",   RING, Essential, Full, bcf,  3);//to uncomment
#ifdef ANEURISM
#else
    BCh_solid->addBC("Base",  RING2, Essential, Full, bcf,  3);
#endif
    BCh_solid->addBC("OuterWall", OUTERWALL, Natural, Full, bcf,  3);
#ifndef ANEURISM100170
#ifdef BLUFF_CYLINDER
#else
#ifdef AORTA
    BCh_solid->addBC("Base",  RING3, Essential, Full, bcf,  3);
    BCh_solid->addBC("Base",  RING4, Essential, Full, bcf,  3);
    BCh_solid->addBC("Base",  RING5, Essential, Full, bcf,  3);
    BCh_solid->addBC("Base",  RING6, Essential, Full, bcf,  3);
    BCh_solid->addBC("Base",  RING7, Essential, Full, bcf,  3);
    BCh_solid->addBC("Base",  RING8, Essential, Full, bcf,  3);
#else
    //BCh_solid->addBC("Edges", INOUTEDGE, Essential, Full, bcf,  3);//to uncomment
#endif
#endif
#endif


    return BCh_solid;
}


// FSIOperator::fluid_bchandler_type BCh_monolithicFluidLin(FSIOperator &_oper)
// {
//     if (! _oper.isFluid() )
//         return FSIOperator::fluid_bchandler_type();

//     // Boundary conditions for the fluid velocity
//     Debug( 10000 ) << "Boundary condition for the linearized fluid\n";
//     FSIOperator::fluid_bchandler_type BCh_fluidLin( new FSIOperator::fluid_bchandler_raw_type );

//     BCFunctionBase bcf(fZero);
//     BCFunctionBase in_flow(u2);

//     BCh_fluidLin->addBC("InFlow",  2,  Natural,   Full, bcf,     3);
//     BCh_fluidLin->addBC("outFlow", 3,  Natural,   Full, bcf,     3);
// #ifdef BLUFF_CYLINDER
// #else
//     BCh_fluidLin->addBC("Edges",  20,  Essential,   Full, bcf,     3);//this condition must be equal to the one
//     //in BCh_fluid. Now it is set to 0 because the mesh displacement is zero in this part of the boundary
// #endif

//     return BCh_fluidLin;
// }


// FSIOperator::solid_bchandler_type BCh_monolithicSolidLin(FSIOperator &_oper)
// {

//     if (! _oper.isSolid() )
//         return FSIOperator::solid_bchandler_type();

//     // Boundary conditions for the solid displacement
//     Debug( 10000 ) << "Boundary condition for the solid\n";
//     FSIOperator::solid_bchandler_type BCh_solid( new FSIOperator::solid_bchandler_raw_type );

//     BCFunctionBase bcf(fZero);

//     // offset still not set here (we set it in applyBoundaryConditions())

//     BCh_solid->addBC("Top",       3, Essential, Full, bcf,  3);
//     BCh_solid->addBC("Base",      2, Essential, Full, bcf,  3);
// #ifdef BLUFF_CYLINDER
// #else
//     BCh_solid->addBC("Edges",    20, Essential, Full, bcf,  3);
// #endif

// //     Debug(10000) << "SP harmonic extension\n";

//     return BCh_solid;
// }



}

#endif
