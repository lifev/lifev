/* -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-

  <short description here>

  Gilles Fourestey gilles.fourestey@epfl.ch

*/
/** \file boundaryConditions.hpp
*/

#ifndef BC_HPP
#define BC_HPP

#include <life/lifecore/life.hpp>
#include "ud_functions.hpp"
#include <life/lifefem/bcHandler.hpp>
#include <life/lifefem/bcFunction.hpp>


namespace LifeV
{

typedef FSI::fluid_raw_type fluid;
typedef FSI::solid_raw_type solid;

FSI::fluidBchandlerPtr_type BCh_harmonicExtension(FSI &_oper)
{

//     Debug(10000) << "SP harmonic extension\n";
//     fixedPoint *FPOper = dynamic_cast<fixedPoint *>(&_oper);

//     FPOper->setStructureDispToHarmonicExtension(_oper.lambdaFluid());

    if (! _oper.isFluid() )
        return FSI::fluidBchandlerPtr_type();

//    FPOper->bcvStructureDispToHarmonicExtension()->showMe(true,std::cout);

    // Boundary condition for the mesh
    Debug( 10000 ) << "Boundary condition for the harmonic extension\n";

    BCFunctionBase bcf(fZero);

    FSISolver::fluidBchandlerPtr_type BCh_he(new FSI::fluidBchandler_Type );

    BCh_he->addBC("Top",       3, Essential, Full, bcf,   3);
    BCh_he->addBC("Base",      2, Essential, Full, bcf,   3);
    BCh_he->addBC("Edges",    20, Essential, Full, bcf,   3);


    if (_oper.method() == "steklovPoincare")
    {
//         Debug(10000) << "SP harmonic extension\n";
//         steklovPoincare *SPOper = dynamic_cast<steklovPoincare *>(&_oper);
//         SPOper->setFluidInterfaceDisp((LifeV::Vector&) _oper.lambdaFluidRepeated());
//         BCh_he->addBC("Interface", 1, Essential, Full,
//                       *SPOper->bcvFluidInterfaceDisp(), 3);
    }
    else if (_oper.method() == "exactJacobian")
    {
        Debug(10000) << "EJ harmonic extension\n";
        exactJacobian *EJOper = dynamic_cast<exactJacobian *>(&_oper);
        EJOper->setStructureDispToHarmonicExtension(_oper.lambdaFluidRepeated());
        BCh_he->addBC("Interface", 1, Essential, Full,
                      *EJOper->bcvStructureDispToHarmonicExtension(), 3);
    }
    else if (_oper.method() == "fixedPoint")
    {
        Debug(10000) << "FP harmonic extension\n";
        fixedPoint *FPOper = dynamic_cast<fixedPoint *>(&_oper);

        FPOper->setStructureDispToHarmonicExtension(_oper.lambdaFluidRepeated());
        BCh_he->addBC("Interface", 1, Essential, Full,
                      *FPOper->bcvStructureDispToHarmonicExtension(), 3);
    }



    return BCh_he;
}


FSI::fluidBchandlerPtr_type BCh_fluid(FSI &_oper)
{
    // Boundary conditions for the fluid velocity
    Debug( 10000 ) << "Boundary condition for the fluid\n";

    if (! _oper.isFluid() )
        return FSI::fluidBchandlerPtr_type();

    FSI::fluidBchandlerPtr_type BCh_fluid( new FSI::fluidBchandler_Type );

    BCFunctionBase bcf      (fZero);
    BCFunctionBase in_flow  (u2);
    BCFunctionBase out_flow (fZero);


    BCh_fluid->addBC("InFlow" , 2,  Natural,   Full, in_flow, 3);
    BCh_fluid->addBC("OutFlow", 3,  Natural,   Full, out_flow, 3);

    _oper.setHarmonicExtensionVelToFluid(_oper.veloFluidMesh());

    BCh_fluid->addBC("Edges",  20, Essential, Full,
                     *_oper.bcvHarmonicExtensionVelToFluid(),  3);
    BCh_fluid->addBC("Interface",   1,  Essential, Full,
                     *_oper.bcvHarmonicExtensionVelToFluid(),  3);

    return BCh_fluid;
}


FSI::fluidBchandlerPtr_type BCh_fluidInv(FSI &_oper)
{

    if (! _oper.isFluid() )
        return FSI::fluidBchandlerPtr_type();

    // Boundary conditions for the fluid velocity
    Debug( 10000 ) << "Boundary condition for the inverse fluid\n";
    FSI::fluidBchandlerPtr_type BCh_fluidInv( new FSI::fluidBchandler_Type );

    BCFunctionBase bcf(fZero);
    BCFunctionBase in_flow(u2);

    BCh_fluidInv->addBC("InFlow", 2,  Natural,   Full, in_flow, 3);
    BCh_fluidInv->addBC("Edges",  20, Essential, Full, bcf,     3);

    return BCh_fluidInv;
}




FSI::fluidBchandlerPtr_type BCh_fluidLin(FSI &_oper)
{
    if (! _oper.isFluid() )
        return FSI::fluidBchandlerPtr_type();

    // Boundary conditions for the fluid velocity
    Debug( 10000 ) << "Boundary condition for the linearized fluid\n";
    FSI::fluidBchandlerPtr_type BCh_fluidLin( new FSI::fluidBchandler_Type );

    BCFunctionBase bcf(fZero);
    BCFunctionBase in_flow(u2);

    BCh_fluidLin->addBC("InFlow",  2,  Natural,   Full, bcf,     3);
    BCh_fluidLin->addBC("outFlow", 3,  Natural,   Full, bcf,     3);
    BCh_fluidLin->addBC("Edges",  20,  Essential,   Full, bcf,     3);//this condition must be equal to the one
    //in BCh_fluid. Now it is set to 0 because the mesh displacement is zero in this part of the boundary

//    BCh_fluidLin->addBC("interface",  1,  Essential,   Full, bcf,     3);


    if (_oper.method() == "steklovPoincare")
    {
//             steklovPoincare *SPOper = dynamic_cast<steklovPoincare *>(&_oper);
//             SPOper->setDerHarmonicExtensionVelToFluid(_oper.dw());
//             BCh_fluidLin->addBC("Wall"     , 1, Essential  , Full,
//                                 *SPOper->bcvDerHarmonicExtensionVelToFluid(), 3);
    }
    if (_oper.method() == "exactJacobian")
    {
        exactJacobian* EJOper = dynamic_cast<exactJacobian *>(&_oper);
        EJOper->setDerHarmonicExtensionVelToFluid(_oper.derVeloFluidMesh());
        BCh_fluidLin->addBC("Interface", 1, Essential  , Full,
                            *_oper.bcvDerHarmonicExtensionVelToFluid(), 3);
    }

    return BCh_fluidLin;
}



FSI::solidBchandlerPtr_type BCh_solid(FSI &_oper)
{

    if (! _oper.isSolid() )
        return FSI::solidBchandlerPtr_type();

    // Boundary conditions for the solid displacement
    Debug( 10000 ) << "Boundary condition for the solid\n";
    FSI::solidBchandlerPtr_type BCh_solid( new FSI::solidBchandler_Type );

    BCFunctionBase bcf(fZero);

    BCh_solid->addBC("Top",       3, Essential, Full, bcf,  3);
    BCh_solid->addBC("Base",      2, Essential, Full, bcf,  3);
    BCh_solid->addBC("Edges",    20, Essential, Full, bcf,  3);

//     Debug(10000) << "SP harmonic extension\n";

    if (_oper.method() == "steklovPoincare")
    {
//         steklovPoincare *SPOper = dynamic_cast<steklovPoincare *>(&_oper);
//         SPOper->setSolidInterfaceDisp((LifeV::Vector&) _oper.displacement());

//         BCh_solid->addBC("Interface", 1, Essential, Full,
//                          *SPOper->bcvSolidInterfaceDisp(), 3);
    }
    else if (_oper.method() == "exactJacobian")
    {
        exactJacobian  *EJOper = dynamic_cast<exactJacobian *>(&_oper);
        EJOper->setFluidLoadToStructure(_oper.sigmaSolidRepeated());

        BCh_solid->addBC("Interface", 1, Natural,   Full,
                         *EJOper->bcvFluidLoadToStructure(), 3);
    }
    else if (_oper.method() == "fixedPoint")
    {
        fixedPoint *FPOper = dynamic_cast<fixedPoint *>(&_oper);

        FPOper->setFluidLoadToStructure(_oper.sigmaSolidRepeated());

        BCh_solid->addBC("Interface", 1, Natural, Full,
                         *FPOper->bcvFluidLoadToStructure(), 3);
    }

    return BCh_solid;
}


FSI::solidBchandlerPtr_type BCh_solidLin(FSI &_oper)
{
    if (! _oper.isSolid() )
        return FSI::solidBchandlerPtr_type();

    // Boundary conditions for the solid displacement
    Debug( 10000 ) << "Boundary condition for the linear solid\n";
    FSI::solidBchandlerPtr_type BCh_solidLin( new FSI::solidBchandler_Type );

    BCFunctionBase bcf(fZero);

    BCh_solidLin->addBC("Top",       3, Essential, Full, bcf,  3);
    BCh_solidLin->addBC("Base",      2, Essential, Full, bcf,  3);
    BCh_solidLin->addBC("Edges",    20, Essential, Full, bcf,  3);

    if (_oper.method() == "steklovPoincare")
    {
//             steklovPoincare *SPOper = dynamic_cast<steklovPoincare *>(&_oper);
//             SPOper->setSolidLinInterfaceDisp((LifeV::Vector&) _oper.displacement());
//             BCh_solidLin->addBC("Interface", 1, Essential, Full,
//                                 *SPOper->bcvSolidLinInterfaceDisp(), 3);
    }
    else if (_oper.method() == "exactJacobian")
    {
        exactJacobian  *EJOper = dynamic_cast<exactJacobian *>(&_oper);
        EJOper->setDerFluidLoadToStructure(_oper.sigmaSolidRepeated());
        BCh_solidLin->addBC("Interface", 1, Natural,   Full,
                            *EJOper->bcvDerFluidLoadToStructure(), 3);
    }

    return BCh_solidLin;
}

FSI::solidBchandlerPtr_type BCh_solidInvLin(FSI &_oper)
{

    if (! _oper.isSolid() )
        return FSI::solidBchandlerPtr_type();

    // Boundary conditions for the solid displacement
    Debug( 10000 ) << "Boundary condition for the inverse linear solid\n";
    FSI::solidBchandlerPtr_type BCh_solidLinInv( new FSI::solidBchandler_Type );

    BCFunctionBase bcf(fZero);

    BCh_solidLinInv->addBC("Top",       3, Essential, Full, bcf,  3);
    BCh_solidLinInv->addBC("Base",      2, Essential, Full, bcf,  3);
    BCh_solidLinInv->addBC("Edges",    20, Essential, Full, bcf,  3);

//     if (_oper.method() == "steklovPoincare")
//     {
//         steklovPoincare *SPOper = dynamic_cast<steklovPoincare *>(&_oper);
//         SPOper->setSolidInvLinInterfaceStress((LifeV::Vector&) _oper.residual());

//         BCh_solidLinInv->addBC("Interface", 100, Natural, Full,
//                                *SPOper->bcvSolidInvLinInterfaceStress(), 3);
//     }
//     else
//     {
//         exactJacobian  *EJOper = dynamic_cast<exactJacobian *>(&_oper);
//         EJOper->setDerFluidLoadToStructure(_oper.fluid().residual());
//         BCh_solidLinInv->addBC("Interface", 1, Natural,   Full,
//                             *EJOper->bcvDerFluidLoadToStructure(), 3);
//     }

    return BCh_solidLinInv;
}



FSI::fluidBchandlerPtr_type BCh_reducedFluid(FSI &_oper)
{

    if (! _oper.isFluid() )
        return FSI::solidBchandlerPtr_type();

    Debug( 10000 ) << "Boundary condition for the reduced fluid\n";
    FSISolver::fluidBchandlerPtr_type BCh_reducedFluid(new FSI::fluidBchandler_Type );
    BCFunctionBase bcf(fZero);

    BCh_reducedFluid->addBC("Wall_Edges", 20, Essential, Scalar, bcf);
    BCh_reducedFluid->addBC("InFlow",      2, Essential, Scalar, bcf);
    BCh_reducedFluid->addBC("OutFlow",     3, Essential, Scalar, bcf);

//     steklovPoincare *SPOper = dynamic_cast<steklovPoincare *>(&_oper);
//     SPOper->setReducedFluidInterfaceAcc(SPOper->getReducedLinFluid()->dacc(), 2);

//     BCh_reducedFluid->addBC("Wall",        1, Natural,   Scalar,
//                             *SPOper->bcvReducedFluidInterfaceAcc());

    return BCh_reducedFluid;
}


FSI::fluidBchandlerPtr_type BCh_reducedFluidInv(FSI &_oper)
{
    if (! _oper.isFluid() )
        return FSI::solidBchandlerPtr_type();

    FSISolver::fluidBchandlerPtr_type BCh_reducedFluidInv( new FSI::fluidBchandler_Type );
    BCFunctionBase bcf(fZero);
    BCh_reducedFluidInv->addBC("Wall_Edges", 20, Essential, Scalar, bcf);
    BCh_reducedFluidInv->addBC("InFlow",      2, Essential, Scalar, bcf);
    BCh_reducedFluidInv->addBC("OutFlow",     3, Essential, Scalar, bcf);

    // the following not correct: this should be strong residual. However, strong residual
    // computation is false at the moment
//     _oper.setDerReducedFluidLoadToStructure(_oper.fluid().residual());

//     BCh_reducedFluidInv->addBC("Wall",        1, Essential, Scalar,//dr_wall);
//                         *_oper.bcvDerReducedFluidLoadToStructure());

    return BCh_reducedFluidInv;

}




}

#endif
