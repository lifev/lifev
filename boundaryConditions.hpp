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
//#include "life/lifesolver/exactJacobianBase.hpp"
#include "life/lifesolver/fixedPointBase.hpp"
//#include "life/lifesolver/Monolithic.hpp"


namespace LifeV
{

typedef FSIOperator::fluid_raw_type fluid;
typedef FSIOperator::solid_raw_type solid;

FSIOperator::fluid_bchandler_type BCh_harmonicExtension(FSIOperator &_oper)
{

//     Debug(10000) << "SP harmonic extension\n";
//     fixedPoint *FPOper = dynamic_cast<fixedPoint *>(&_oper);

//     FPOper->setStructureDispToHarmonicExtension(_oper.lambdaFluid());

    if (! _oper.isFluid() )
        return FSIOperator::fluid_bchandler_type();

//    FPOper->bcvStructureDispToHarmonicExtension()->showMe(true,std::cout);

    // Boundary condition for the mesh
    Debug( 10000 ) << "Boundary condition for the harmonic extension\n";

    BCFunctionBase bcf(fZero);

    FSISolver::fluid_bchandler_type BCh_he(new FSIOperator::fluid_bchandler_raw_type );

    BCh_he->addBC("Top",       3, Essential, Full, bcf,   3);
    BCh_he->addBC("Base",      2, Essential, Full, bcf,   3);
    BCh_he->addBC("Edges",    20, Essential, Full, bcf,   3);


    if (_oper.method() == "steklovPoincare")
    {
//         Debug(10000) << "SP harmonic extension\n";
//         steklovPoincare *SPOper = dynamic_cast<steklovPoincare *>(&_oper);
//         SPOper->setFluidInterfaceDisp((LifeV::Vector&) _oper.lambdaFluid());
//         BCh_he->addBC("Interface", 1, Essential, Full,
//                       *SPOper->bcvFluidInterfaceDisp(), 3);
    }
    else if (_oper.method() == "monolithic")
    {
                Debug(10000) << "EJ harmonic extension\n";
        Monolithic *EJOper = dynamic_cast<Monolithic *>(&_oper);
        //        exactJacobian *EJOper = dynamic_cast<exactJacobian *>(&_oper);

        EJOper->setStructureDispToHarmonicExtension(_oper.lambdaFluidRepeated());

        std::cout << "lambdaFluid: " << _oper.lambdaFluidRepeated().NormInf() << std::endl ;

        BCh_he->addBC("Interface", 1, Essential, Full,
        *EJOper->bcvStructureDispToHarmonicExtension(), 3);
    }
    else if (_oper.method() == "fixedPoint")
    {
        /*        Debug(10000) << "FP harmonic extension\n";
        fixedPoint *FPOper = dynamic_cast<fixedPoint *>(&_oper);

        FPOper->setStructureDispToHarmonicExtension(_oper.lambdaFluid());
        BCh_he->addBC("Interface", 1, Essential, Full,
        *FPOper->bcvStructureDispToHarmonicExtension(), 3);*/
    }



    return BCh_he;
}



FSIOperator::fluid_bchandler_type BCh_monolithicFluid(FSIOperator &_oper)
{
    // Boundary conditions for the fluid velocity
    Debug( 10000 ) << "Boundary condition for the fluid\n";

    if (! _oper.isFluid() )
        return FSIOperator::fluid_bchandler_type();

    FSISolver::fluid_bchandler_type BCh_fluid( new FSIOperator::fluid_bchandler_raw_type );

    BCFunctionBase bcf      (fZero);
    BCFunctionBase in_flow  (u2);
    BCFunctionBase out_flow (fZero);


    BCh_fluid->addBC("InFlow" , 2,  Natural,   Full, in_flow, 3);
    BCh_fluid->addBC("OutFlow", 3,  Natural,   Full, out_flow, 3);

    //    _oper.setHarmonicExtensionVelToFluid(_oper.veloFluidMesh());

    //    BCh_fluid->addBC("Edges",  20, Essential, Full,
    //                     *_oper.bcvHarmonicExtensionVelToFluid(),  3);
    //    BCh_fluid->addBC("Interface",   1,  Essential, Full,
    //                     *_oper.bcvHarmonicExtensionVelToFluid(),  3);

    return BCh_fluid;
}


FSIOperator::solid_bchandler_type BCh_monolithicSolid(FSIOperator &_oper)
{

    if (! _oper.isSolid() )
        return FSIOperator::solid_bchandler_type();

    // Boundary conditions for the solid displacement
    Debug( 10000 ) << "Boundary condition for the solid\n";
    FSIOperator::solid_bchandler_type BCh_solid( new FSIOperator::solid_bchandler_raw_type );

    BCFunctionBase bcf(fZero);

    // offset still not set here (we set it in applyBoundaryConditions())

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
    else if (_oper.method() == "monolithic")
    {
        //        Monolithic  *EJOper = dynamic_cast<Monolithic *>(&_oper);
        //        EJOper->setFluidLoadToStructure(_oper.sigmaSolid());

        //        BCh_solid->addBC("Interface", 1, Natural,   Full,
        //                   *EJOper->bcvFluidLoadToStructure(), 3);*/
    }
    else if (_oper.method() == "fixedPoint")
    {
        /*        fixedPoint *FPOper = dynamic_cast<fixedPoint *>(&_oper);

        FPOper->setFluidLoadToStructure(_oper.sigmaSolid());

        //        BCh_solid->addBC("Interface", 1, Natural, Full,
        //                         *FPOper->bcvFluidLoadToStructure(), 3);*/
    }

    return BCh_solid;
}



}

#endif
