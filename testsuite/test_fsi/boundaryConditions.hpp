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

namespace LifeV
{

// typedef VenantKirchhofSolver< RegionMesh3D_ALE<LinearTetra> >    solid_raw_type;
// typedef NavierStokesAleSolverPC< RegionMesh3D_ALE<LinearTetra> > fluid_raw_type;

// typedef solid_raw_type::bchandler_type        solid_bchandler_type;
// typedef solid_raw_type::bchandler_raw_type    solid_bchandler_raw_type;

// typedef fluid_raw_type::bchandler_type        fluid_bchandler_type;
// typedef fluid_raw_type::bchandler_raw_type    fluid_bchandler_raw_type;

FSIOperator::fluid_bchandler_type BCh_harmonicExtension()
{
    // Boundary condition for the mesh
    Debug( 10000 ) << "Boundary condition for the harmonic extension\n";

    BCFunctionBase bcf(fZero);
    FSISolver::fluid_bchandler_type BCh_hamonicExtension( new BCHandler );
    BCh_hamonicExtension->addBC("Top",       3, Essential, Full, bcf,   3);
    BCh_hamonicExtension->addBC("Base",      2, Essential, Full, bcf,   3);
    BCh_hamonicExtension->addBC("Edges",    20, Essential, Full, bcf,   3);

    return BCh_hamonicExtension;
}


FSIOperator::fluid_bchandler_type BCh_fluid()
{
    // Boundary conditions for the fluid velocity
    Debug( 10000 ) << "Boundary condition for the fluid\n";
    FSIOperator::fluid_bchandler_type BCh_fluid( new FSIOperator::fluid_bchandler_raw_type );

    BCFunctionBase bcf(fZero);
    BCFunctionBase in_flow(u2);

    BCh_fluid->addBC("InFlow", 2,  Natural,   Full, in_flow, 3);
    BCh_fluid->addBC("Edges",  20, Essential, Full, bcf,     3);

//     fsi.operFSI()->setBcvHarmonicExtensionVelToFluid(fsi.oper->fluid().wInterpolated());

//     M_BCh_fluid->addBC("Interface",   1,  Essential, Full,
//                        fsi.operFSI()->bcvHarmonicExtensionVelToFluid(),  3);
    return BCh_fluid;
}


FSIOperator::fluid_bchandler_type BCh_fluidInv()
{
    // Boundary conditions for the fluid velocity
    Debug( 10000 ) << "Boundary condition for the fluid\n";
    FSIOperator::fluid_bchandler_type BCh_fluidInv( new FSIOperator::fluid_bchandler_raw_type );

    BCFunctionBase bcf(fZero);
    BCFunctionBase in_flow(u2);

    BCh_fluidInv->addBC("InFlow", 2,  Natural,   Full, in_flow, 3);
    BCh_fluidInv->addBC("Edges",  20, Essential, Full, bcf,     3);

//     fsi.operFSI()->setBcvHarmonicExtensionVelToFluid(fsi.oper->fluid().wInterpolated());

//     M_BCh_fluid->addBC("Interface",   1,  Essential, Full,
//                        fsi.operFSI()->bcvHarmonicExtensionVelToFluid(),  3);
    return BCh_fluidInv;
}




FSIOperator::fluid_bchandler_type BCh_fluidLin()
{
    // Boundary conditions for the fluid velocity
    Debug( 10000 ) << "Boundary condition for the fluid\n";
    FSIOperator::fluid_bchandler_type BCh_fluidLin( new FSIOperator::fluid_bchandler_raw_type );

    BCFunctionBase bcf(fZero);
    BCFunctionBase in_flow(u2);

    BCh_fluidLin->addBC("InFlow", 2,  Natural,   Full, in_flow, 3);
    BCh_fluidLin->addBC("Edges",  20, Essential, Full, bcf,     3);

//     fsi.operFSI()->setBcvHarmonicExtensionVelToFluid(fsi.oper->fluid().wInterpolated());

//     M_BCh_fluid->addBC("Interface",   1,  Essential, Full,
//                        fsi.operFSI()->bcvHarmonicExtensionVelToFluid(),  3);
    return BCh_fluidLin;
}



FSIOperator::solid_bchandler_type BCh_solid()
{
    // Boundary conditions for the solid displacement
    Debug( 10000 ) << "Boundary condition for the solid\n";
    FSIOperator::solid_bchandler_type BCh_solid( new FSIOperator::solid_bchandler_raw_type );
    BCFunctionBase bcf(fZero);

    BCh_solid->addBC("Top",       3, Essential, Full, bcf,  3);
    BCh_solid->addBC("Base",      2, Essential, Full, bcf,  3);

    return BCh_solid;
}


FSIOperator::fluid_bchandler_type BCh_reducedFluid()
{
    FSISolver::fluid_bchandler_type BCh_reducedFluid(new FSIOperator::fluid_bchandler_raw_type );
    BCFunctionBase bcf(fZero);
//    BCh_reducedFluid->addBC("Wall",        1, Natural,   Scalar, da_wall);
    BCh_reducedFluid->addBC("Wall_Edges", 20, Essential, Scalar, bcf);
    BCh_reducedFluid->addBC("InFlow",      2, Essential, Scalar, bcf);
    BCh_reducedFluid->addBC("OutFlow",     3, Essential, Scalar, bcf);

    return BCh_reducedFluid;
}


FSIOperator::fluid_bchandler_type BCh_reducedFluidInv()
{
    FSISolver::fluid_bchandler_type BCh_reducedFluidInv( new FSIOperator::fluid_bchandler_raw_type );
    BCFunctionBase bcf(fZero);
//    BCh_reducedFluidInv->addBC("Wall",        1, Essential, Scalar, dr_wall);
    BCh_reducedFluidInv->addBC("Wall_Edges", 20, Essential, Scalar, bcf);
    BCh_reducedFluidInv->addBC("InFlow",      2, Essential, Scalar, bcf);
    BCh_reducedFluidInv->addBC("OutFlow",     3, Essential, Scalar, bcf);

    return BCh_reducedFluidInv;

}




}

#endif
