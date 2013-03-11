/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s):  A. Fumagalli  <alessio.fumagalli@mail.polimi.it>
       Date: 2010-07-29

  Copyright (C) 2010 EPFL, Politecnico di Milano

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA
*/
/**
   \file impes.cpp
   \author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
   \date 2010-07-29
 */

// ===================================================
//! Includes
// ===================================================

#include "impes.hpp"
#include "user_fun.hpp"

// ===================================================
//! Namespaces & define
// ===================================================

using namespace LifeV;

enum BCNAME
{

    FLUX0            = 0,
    INLETPRESSURE1   = 1,
    INLETPRESSURE2   = 2,
    OUTLETPRESSURE   = 3,
    FLUX1            = 4

    //    BACK   = 1,
    //    FRONT  = 2,
    //    LEFT   = 3,
    //    RIGHT  = 4,
    //    BOTTOM = 5,
    //    TOP    = 6
};

namespace dataProblem
{

// Standard functions
Real UOne ( const Real& /* t */,
            const Real& /* x */,
            const Real& /* y */,
            const Real& /* z */,
            const ID&   /* icomp */)
{
    return 1.;
}

Real UZero ( const Real& /* t */,
             const Real& /* x */,
             const Real& /* y */,
             const Real& /* z */,
             const ID&   /* icomp */)
{
    return 0.;
}

}
// ===================================================
//! Private Members
// ===================================================

struct impes::Private
{
    Private() {}

    // Policy for scalar functions
    typedef boost::function < Real ( const Real&, const Real&,
                                     const Real&, const Real&, const ID& ) >
    fct_type;

    // Policy for vector function
    typedef boost::function < Vector ( const Real&, const Real&,
                                       const Real&, const Real&,
                                       const std::vector<Real>& ) >
    Vfct_type;


    // Policy for matrix function
    typedef boost::function < Matrix ( const Real&, const Real&,
                                       const Real&, const Real&,
                                       const std::vector<Real>& ) >
    Mfct_type;

    std::string data_file_name;

    // General section
    std::string discretization_section;

    // Section for the Darcy solver
    std::string discretization_section_darcy;

    // Section for the hyperbolic solver
    std::string discretization_section_hyperbolic;

    // Section for the non-linear and transient Darcy solver
    std::string discretization_section_darcy_nonlin_trans;

    boost::shared_ptr<Epetra_Comm>   comm;

    // Function Types

    fct_type getUOne()
    {
        fct_type f;
        f = boost::bind ( &dataProblem::UOne, _1, _2, _3, _4, _5 );
        return f;
    }

    fct_type getUZero()
    {
        fct_type f;
        f = boost::bind ( &dataProblem::UZero, _1, _2, _3, _4, _5 );
        return f;
    }

    fct_type getPressureSource()
    {
        fct_type f;
        f = boost::bind ( &dataProblem::pressureSource, _1, _2, _3, _4, _5 );
        return f;
    }

    Mfct_type getPressurePermeability()
    {
        Mfct_type m;
        m = boost::bind ( &dataProblem::pressurePermeability, _1, _2, _3, _4, _5 );
        return m;
    }

    fct_type getSaturationInitialCondition()
    {
        fct_type f;
        f = boost::bind ( &dataProblem::saturationInitialCondition, _1, _2, _3, _4, _5 );
        return f;
    }

    Vfct_type getSaturationPhysicalFlux()
    {
        Vfct_type v;
        v = boost::bind ( &dataProblem::saturationPhysicalFlux, _1, _2, _3, _4, _5 );
        return v;
    }

    Vfct_type getSaturationFirstDerivativePhysicalFlux()
    {
        Vfct_type v;
        v = boost::bind ( &dataProblem::saturationFirstDerivativePhysicalFlux, _1, _2, _3, _4, _5 );
        return v;
    }

    fct_type getSaturationSource()
    {
        fct_type f;
        f = boost::bind ( &dataProblem::saturationSource, _1, _2, _3, _4, _5 );
        return f;
    }

    Mfct_type getSaturationPermeability()
    {
        Mfct_type mnl;
        mnl = boost::bind ( &dataProblem::saturationPermeability, _1, _2, _3, _4, _5 );
        return mnl;
    }

    fct_type getSaturationMass()
    {
        fct_type f;
        f = boost::bind ( &dataProblem::saturationMass, _1, _2, _3, _4, _5 );
        return f;
    }

};

// ===================================================
//! Constructors
// ===================================================

impes::impes ( int argc,
               char** argv )
    : Members ( new Private )
{
    GetPot command_line (argc, argv);
    const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile ( data_file_name );

    Members->data_file_name = data_file_name;
    Members->discretization_section = "impes";
    Members->discretization_section_darcy = Members->discretization_section + "/darcy";
    Members->discretization_section_hyperbolic = Members->discretization_section + "/hyperbolic";
    Members->discretization_section_darcy_nonlin_trans =  Members->discretization_section + "/darcy_transient_non_linear";

#ifdef EPETRA_MPI
    std::cout << "Epetra Initialization" << std::endl;
    Members->comm.reset ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
    int ntasks;
    MPI_Comm_size (MPI_COMM_WORLD, &ntasks);
#else
    Members->comm.reset ( new Epetra_SerialComm() );
#endif

}

// ===================================================
//! Methods
// ===================================================

Real
impes::run()
{
    typedef RegionMesh<LinearTetra>                                RegionMesh;
    typedef SolverAztecOO                                            solver_type;
    typedef DarcySolver< RegionMesh, solver_type >                   ds;
    typedef DarcySolverTransientNonLinear< RegionMesh, solver_type > dstnl;
    typedef HyperbolicSolver< RegionMesh, solver_type >              hyper;
    typedef ds::vector_Type                                          vector_type;
    typedef boost::shared_ptr<vector_type>                           vector_ptrtype;
    typedef FESpace< RegionMesh, MapEpetra >                         feSpace_Type;
    typedef boost::shared_ptr<feSpace_Type>                          feSpacePtr_Type;

    LifeChrono chronoTotal;
    LifeChrono chronoReadAndPartitionMesh;
    LifeChrono chronoBoundaryCondition;
    LifeChrono chronoFiniteElementSpace;
    LifeChrono chronoProblem;
    LifeChrono chronoProcess;
    LifeChrono chronoError;

    // Start chronoTotal for measure the total time for the computation.
    chronoTotal.start();

    // Reading from data file.
    GetPot dataFile ( Members->data_file_name.c_str() );

    // Create the leader process, i.e. the process with MyPID equal to zero.
    bool isLeader = ( Members->comm->MyPID() == 0 );

    //
    // The IMPES Solver.
    //

    if ( isLeader )
    {
        std::cout << "The IMPES solver" << std::endl << std::flush;
    }

    // Start chronoReadAndPartitionMesh for measure the total time for the creation of the local meshes.
    chronoReadAndPartitionMesh.start();

    // Create the data file for the pressure equation.
    DarcyData<RegionMesh> dataPressure;

    // Create the data file for the hyperbolic solver in the saturation equation.
    HyperbolicData<RegionMesh> dataSaturationHyperbolic;

    // Create the data file for the non-linear and transient Darcy solver in the saturation equation.
    DarcyData<RegionMesh> dataSaturationDarcyNonLinear;

    // Set up the data for the pressure equation.
    dataPressure.setup ( dataFile, Members->discretization_section_darcy );

    // Set up the data for the hyperbolic solver.
    dataSaturationHyperbolic.setup ( dataFile, Members->discretization_section_hyperbolic );

    // Set up the data for the non-linear and transient Darcy solver.
    dataSaturationDarcyNonLinear.setup ( dataFile, Members->discretization_section_darcy_nonlin_trans );

    // Create the mesh file handler.
    MeshData meshData;

    // Set up the mesh file.
    meshData.setup ( dataFile,  Members->discretization_section + "/space_discretization");

    // Create the mesh.
    boost::shared_ptr<RegionMesh> fullMeshPtr ( new RegionMesh ( * ( Members->comm ) ) );

    // Set up the mesh.
    readMesh ( *fullMeshPtr, meshData );

    // Partition the mesh using ParMetis.
    boost::shared_ptr<RegionMesh> meshPtr;
    {
        MeshPartitioner< RegionMesh >  meshPart ( fullMeshPtr, Members->comm );
        meshPtr = meshPart.meshPartition();
    }

    // Stop chronoReadAndPartitionMesh.
    chronoReadAndPartitionMesh.stop();

    // The leader process print chronoReadAndPartitionMesh.
    if ( isLeader )
        std::cout << "Time for read and partition the mesh " <<
                  chronoReadAndPartitionMesh.diff() << std::endl << std::flush;

    // Create the boundary conditions.

    // Start chronoBoundaryCondition for measure the total time for create the boundary conditions.
    chronoBoundaryCondition.start();

    BCFunctionBase pressureDirichletBDfun1,
                   pressureDirichletBDfun2,
                   pressureDirichletBDfun3,
                   pressureNeumannBDfun;

    BCFunctionBase saturationDirichletBDfun1,
                   saturationDirichletBDfun2,
                   saturationDirichletBDfun3,
                   saturationNeumannBDfun;

    //BCFunctionRobin pressureRobinBDfun, saturationRobinBDfun;

    // Set pressure bounday data.
    pressureDirichletBDfun1.setFunction ( dataProblem::pressureDirichlet1 );
    pressureDirichletBDfun2.setFunction ( dataProblem::pressureDirichlet2 );
    pressureDirichletBDfun3.setFunction ( dataProblem::pressureDirichlet3 );

    pressureNeumannBDfun.setFunction    ( dataProblem::pressureNeumann );

    // dp/dn = first_parameter + second_parameter * p.
    //pressureRobinBDfun.setFunctions_Robin ( dataProblem::pressureRobin,
    //                                         Members->getUOne() );

    // Set pressure bounday data.
    saturationDirichletBDfun1.setFunction ( dataProblem::saturationDirichlet1 );
    saturationDirichletBDfun2.setFunction ( dataProblem::saturationDirichlet2 );
    saturationDirichletBDfun3.setFunction ( dataProblem::saturationDirichlet3 );

    saturationNeumannBDfun.setFunction    ( dataProblem::saturationNeumann );
    // dp/dn = first_parameter + second_parameter * p.
    //pressureRobinBDfun.setFunctions_Robin ( dataProblem::saturationRobin,
    //                                      Members->getUOne() );

    // Boundary condition handler for the pressure equation.
    BCHandler bcPressure;

    // Boundary condition handler for the saturation equation.
    BCHandler bcSaturation;

    bcPressure.addBC ( "Top",            FLUX0,          Natural,   Full,    pressureNeumannBDfun, 0 );
    bcPressure.addBC ( "Top2",           FLUX1,          Natural,   Full,    pressureNeumannBDfun, 0 );
    bcPressure.addBC ( "InletPressure",  INLETPRESSURE1, Essential, Scalar,  pressureDirichletBDfun1  );
    bcPressure.addBC ( "InletPressure1", INLETPRESSURE2, Essential, Scalar,  pressureDirichletBDfun3  );
    bcPressure.addBC ( "OutletPressure", OUTLETPRESSURE, Essential, Scalar,  pressureDirichletBDfun2  );

    bcSaturation.addBC ( "Top",            FLUX0,          Natural,   Full,    saturationNeumannBDfun, 0 );
    bcSaturation.addBC ( "Top2",           FLUX1,          Natural,   Full,    saturationNeumannBDfun, 0 );
    bcSaturation.addBC ( "InletPressure",  INLETPRESSURE1, Essential, Scalar,  saturationDirichletBDfun1  );
    bcSaturation.addBC ( "InletPressure1", INLETPRESSURE2, Essential, Scalar,  saturationDirichletBDfun3  );
    bcSaturation.addBC ( "OutletPressure", OUTLETPRESSURE, Essential, Scalar,  saturationDirichletBDfun2  );

    // Setting the boundary condition for the pressure equation.
    //bcPressure.addBC(   "Top",    TOP,    Essential,  Scalar,  pressureDirichletBDfun1 );
    //bcPressure.addBC("Bottom", BOTTOM,    Essential,  Scalar,  pressureDirichletBDfun1 );
    //bcPressure.addBC(  "Left",   LEFT,    Essential,  Scalar,  pressureDirichletBDfun1 );
    //bcPressure.addBC( "Right",  RIGHT,    Essential,  Scalar,  pressureDirichletBDfun1 );
    //bcPressure.addBC( "Front",  FRONT,    Essential,  Scalar,  pressureDirichletBDfun1 );
    //bcPressure.addBC(  "Back",   BACK,    Essential,  Scalar,  pressureDirichletBDfun1 );

    // Setting the boundary condition for the saturation equation.
    //bcSaturation.addBC(   "Top",    TOP,    Essential,  Scalar,  saturationDirichletBDfun1 );
    //bcSaturation.addBC("Bottom", BOTTOM,    Essential,  Scalar,  saturationDirichletBDfun1 );
    //bcSaturation.addBC(  "Left",   LEFT,    Essential,  Scalar,  saturationDirichletBDfun1 );
    //bcSaturation.addBC( "Right",  RIGHT,    Essential,  Scalar,  saturationDirichletBDfun1 );
    //bcSaturation.addBC( "Front",  FRONT,    Essential,  Scalar,  saturationDirichletBDfun1 );
    //bcSaturation.addBC(  "Back",   BACK,    Essential,  Scalar,  saturationDirichletBDfun1 );

    // Stop chronoBoundaryCondition.
    chronoBoundaryCondition.stop();

    // The leader process print chronoBoundaryCondition.
    if ( isLeader )
    {
        std::cout << "Time for create the boundary conditions handler " <<
                  chronoBoundaryCondition.diff() << std::endl << std::flush;

    }

    // Create the solution spaces.

    // Start chronoFiniteElementSpace for measure the total time for create the finite element spaces.
    chronoFiniteElementSpace.start();

    // We impose directly the compatibily condition on the FE spaces.

    // Parameters for the pressure equation.
    const ReferenceFE*    pressure_refFE_primal ( static_cast<ReferenceFE*> (NULL) );
    const QuadratureRule* pressure_qR_primal    ( static_cast<QuadratureRule*> (NULL) );
    const QuadratureRule* pressure_bdQr_primal  ( static_cast<QuadratureRule*> (NULL) );

    pressure_refFE_primal = &feTetraP0;
    pressure_qR_primal    = &quadRuleTetra15pt;
    pressure_bdQr_primal  = &quadRuleTria4pt;

    // Dual solution parameters.
    const ReferenceFE*    pressure_refFE_dual ( static_cast<ReferenceFE*> (NULL) );
    const QuadratureRule* pressure_qR_dual    ( static_cast<QuadratureRule*> (NULL) );
    const QuadratureRule* pressure_bdQr_dual  ( static_cast<QuadratureRule*> (NULL) );

    pressure_refFE_dual = &feTetraRT0;
    pressure_qR_dual    = &quadRuleTetra15pt;
    pressure_bdQr_dual  = &quadRuleTria4pt;

    // Interpolate of dual solution parameters.
    const ReferenceFE*    pressure_refFE_dualInterpolate ( static_cast<ReferenceFE*> (NULL) );
    const QuadratureRule* pressure_qR_dualInterpolate    ( static_cast<QuadratureRule*> (NULL) );
    const QuadratureRule* pressure_bdQr_dualInterpolate  ( static_cast<QuadratureRule*> (NULL) );

    pressure_refFE_dualInterpolate = &feTetraP0;
    pressure_qR_dualInterpolate    = &quadRuleTetra15pt;
    pressure_bdQr_dualInterpolate  = &quadRuleTria4pt;

    // Hybrid solution parameters.
    // hybrid.
    const ReferenceFE*    pressure_refFE_hybrid ( static_cast<ReferenceFE*> (NULL) );
    const QuadratureRule* pressure_qR_hybrid    ( static_cast<QuadratureRule*> (NULL) );
    const QuadratureRule* pressure_bdQr_hybrid  ( static_cast<QuadratureRule*> (NULL) );

    pressure_refFE_hybrid = &feTetraRT0Hyb;
    pressure_qR_hybrid    = &quadRuleTetra15pt;
    pressure_bdQr_hybrid  = &quadRuleTria4pt;

    // dual dot outward unit normal.
    const ReferenceFE*    pressure_refFE_VdotN ( static_cast<ReferenceFE*> (NULL) );
    const QuadratureRule* pressure_qR_VdotN    ( static_cast<QuadratureRule*> (NULL) );
    const QuadratureRule* pressure_bdQr_VdotN  ( static_cast<QuadratureRule*> (NULL) );

    pressure_refFE_VdotN = &feTetraRT0VdotNHyb;
    pressure_qR_VdotN    = &quadRuleTetra15pt;
    pressure_bdQr_VdotN  = &quadRuleTria4pt;


    // Finite element space of the primal variable.
    feSpacePtr_Type pressure_p_FESpacePtr ( new feSpace_Type ( meshPtr, *pressure_refFE_primal, *pressure_qR_primal,
                                                               *pressure_bdQr_primal, 1, Members->comm ) );

    // Finite element space of the dual variable.
    feSpacePtr_Type pressure_u_FESpacePtr ( new feSpace_Type ( meshPtr, *pressure_refFE_dual, *pressure_qR_dual,
                                                               *pressure_bdQr_dual, 1, Members->comm ) );

    // Finite element space of the interpolation of dual variable.
    feSpacePtr_Type pressure_uInterpolate_FESpacePtr ( new feSpace_Type ( meshPtr, *pressure_refFE_dualInterpolate,
                                                                          *pressure_qR_dualInterpolate,
                                                                          *pressure_bdQr_dualInterpolate, 3, Members->comm ) );

    // Vector for the interpolated dual solution.
    vector_ptrtype pressure_dualInterpolated ( new vector_type ( pressure_uInterpolate_FESpacePtr->map(), Repeated ) );

    // Finite element space of the hybrid variable.
    FESpace< RegionMesh, MapEpetra > pressure_hybrid_FESpace ( meshPtr, *pressure_refFE_hybrid, *pressure_qR_hybrid,
                                                               *pressure_bdQr_hybrid, 1, Members->comm );

    // Finite element space of the  outward unit normal variable.
    FESpace< RegionMesh, MapEpetra > pressure_VdotN_FESpace ( meshPtr, *pressure_refFE_VdotN, *pressure_qR_VdotN,
                                                              *pressure_bdQr_VdotN, 1, Members->comm );

    // Parameters for the saturation equation.

    // Hyperbolic parameters.
    const ReferenceFE*    saturation_hyperbolic_refFE ( static_cast<ReferenceFE*> (NULL) );
    const QuadratureRule* saturation_hyperbolic_qR    ( static_cast<QuadratureRule*> (NULL) );
    const QuadratureRule* saturation_hyperbolic_bdQr  ( static_cast<QuadratureRule*> (NULL) );

    saturation_hyperbolic_refFE = pressure_refFE_primal;
    saturation_hyperbolic_qR    = pressure_qR_primal;
    saturation_hyperbolic_bdQr  = &quadRuleTria1pt;

    // Finite element space.
    feSpacePtr_Type saturation_hyperbolic_FESpacePtr ( new feSpace_Type ( meshPtr, *saturation_hyperbolic_refFE,
                                                                          *saturation_hyperbolic_qR,
                                                                          *saturation_hyperbolic_bdQr, 1, Members->comm ) );

    // Non-linear and transient Darcy parameters.

    // Primal solution parameters.
    const ReferenceFE*    saturation_darcy_refFE_primal ( static_cast<ReferenceFE*> (NULL) );
    const QuadratureRule* saturation_darcy_qR_primal    ( static_cast<QuadratureRule*> (NULL) );
    const QuadratureRule* saturation_darcy_bdQr_primal  ( static_cast<QuadratureRule*> (NULL) );

    saturation_darcy_refFE_primal = pressure_refFE_primal;
    saturation_darcy_qR_primal    = pressure_qR_primal;
    saturation_darcy_bdQr_primal  = pressure_bdQr_dual;

    // Dual solution parameters.
    const ReferenceFE*    saturation_darcy_refFE_dual ( static_cast<ReferenceFE*> (NULL) );
    const QuadratureRule* saturation_darcy_qR_dual    ( static_cast<QuadratureRule*> (NULL) );
    const QuadratureRule* saturation_darcy_bdQr_dual  ( static_cast<QuadratureRule*> (NULL) );

    saturation_darcy_refFE_dual = pressure_refFE_dual;
    saturation_darcy_qR_dual    = pressure_qR_dual;
    saturation_darcy_bdQr_dual  = pressure_bdQr_dual;

    // Hybrid solution parameters.
    // hybrid.
    const ReferenceFE*    saturation_darcy_refFE_hybrid ( static_cast<ReferenceFE*> (NULL) );
    const QuadratureRule* saturation_darcy_qR_hybrid    ( static_cast<QuadratureRule*> (NULL) );
    const QuadratureRule* saturation_darcy_bdQr_hybrid  ( static_cast<QuadratureRule*> (NULL) );

    saturation_darcy_refFE_hybrid = pressure_refFE_hybrid;
    saturation_darcy_qR_hybrid    = pressure_qR_hybrid;
    saturation_darcy_bdQr_hybrid  = pressure_bdQr_hybrid;

    // dual dot outward unit normal.
    const ReferenceFE*    saturation_darcy_refFE_VdotN ( static_cast<ReferenceFE*> (NULL) );
    const QuadratureRule* saturation_darcy_qR_VdotN    ( static_cast<QuadratureRule*> (NULL) );
    const QuadratureRule* saturation_darcy_bdQr_VdotN  ( static_cast<QuadratureRule*> (NULL) );

    saturation_darcy_refFE_VdotN = pressure_refFE_VdotN;
    saturation_darcy_qR_VdotN    = pressure_qR_VdotN;
    saturation_darcy_bdQr_VdotN  = pressure_bdQr_VdotN;


    // Finite element space of the primal variable.
    FESpace< RegionMesh, MapEpetra > saturation_darcy_p_FESpace ( meshPtr, *saturation_darcy_refFE_primal,
                                                                  *saturation_darcy_qR_primal,
                                                                  *saturation_darcy_bdQr_primal, 1, Members->comm );

    // Finite element space of the dual variable.
    FESpace< RegionMesh, MapEpetra > saturation_darcy_u_FESpace ( meshPtr, *saturation_darcy_refFE_dual,
                                                                  *saturation_darcy_qR_dual,
                                                                  *saturation_darcy_bdQr_dual, 1, Members->comm );

    // Finite element space of the hybrid variable.
    FESpace< RegionMesh, MapEpetra > saturation_darcy_hybrid_FESpace ( meshPtr, *saturation_darcy_refFE_hybrid,
                                                                       *saturation_darcy_qR_hybrid,
                                                                       *saturation_darcy_bdQr_hybrid, 1, Members->comm );

    // Finite element space of the  outward unit normal variable.
    FESpace< RegionMesh, MapEpetra > saturation_darcy_VdotN_FESpace ( meshPtr, *saturation_darcy_refFE_VdotN,
                                                                      *saturation_darcy_qR_VdotN,
                                                                      *saturation_darcy_bdQr_VdotN, 1, Members->comm );


    // Stop chronoFiniteElementSpace.
    chronoFiniteElementSpace.stop();

    // The leader process print chronoFiniteElementSpace.
    if ( isLeader )
        std::cout << "Time for create the finite element spaces " <<
                  chronoFiniteElementSpace.diff() << std::endl << std::flush;

    // Start chronoProblem for measure the total time for create the problem.
    chronoProblem.start();

    // Instantiation of the pressure equation solver, i.e. the Darcy solver.
    ds pressureSolver ( dataPressure, *pressure_p_FESpacePtr, *pressure_u_FESpacePtr,
                        pressure_hybrid_FESpace, pressure_VdotN_FESpace, Members->comm );

    // Instantiation of the saturation equation solver.

    // Instantiation of the hyperbolic solver.
    hyper saturationHyperbolicSolver ( dataSaturationHyperbolic, *saturation_hyperbolic_FESpacePtr, Members->comm );

    // Instantiation of the non-linear and transient Darcy solver.
    dstnl saturationDarcySolver ( dataSaturationDarcyNonLinear, saturation_darcy_p_FESpace,
                                  saturation_darcy_u_FESpace, saturation_darcy_hybrid_FESpace,
                                  saturation_darcy_VdotN_FESpace, Members->comm );

    // Stop chronoProblem.
    chronoProblem.stop();

    // The leader process print chronoProblem.
    pressureSolver.getDisplayer().leaderPrint ( "Time for create the problem ",
                                                chronoProblem.diff(), "\n" );

    // Process the problem.

    // Start chronoProcess for measure the total time for the simulation.
    chronoProcess.start();

    // Set up for the pressure equation.

    // Set up phase for the linear solver for the pressure equation.
    pressureSolver.setup();

    // Create the inverse permeability.
    inversePermeability < RegionMesh > invPermPress ( Members->getPressurePermeability(),
                                                      saturation_darcy_p_FESpace );

    // Set the saturation dependence for the permeability tensor.
    invPermPress.setField ( saturationDarcySolver.primalSolution() );

    // Set the inverse of the permeability in the pressure equation.
    pressureSolver.setInversePermeability ( invPermPress );

    // Set the source term in the pressure equation.
    pressureSolver.setSource ( Members->getPressureSource() );

    // Set the boudary conditions.
    pressureSolver.setBC ( bcPressure );

    // Set up for the saturation equation.

    // Set up for the hyperbolic solver in the saturation equation.
    saturationHyperbolicSolver.setup();

    // Create the numerical flux
    GodunovNumericalFlux < RegionMesh > numericalFlux ( Members->getSaturationPhysicalFlux(),
                                                        Members->getSaturationFirstDerivativePhysicalFlux(),
                                                        *pressure_uInterpolate_FESpacePtr,
                                                        dataFile,
                                                        "impes/hyperbolic/numerical_flux/" );

    // Set the dependence field
    numericalFlux.setExternalField ( pressure_dualInterpolated );

    // Set the numerical flux usign the physical flux
    saturationHyperbolicSolver.setNumericalFlux ( numericalFlux );

    // Set the porosity
    saturationHyperbolicSolver.setMassTerm ( Members->getSaturationMass() );

    // Set the boudary conditions.
    saturationHyperbolicSolver.setBoundaryCondition ( bcSaturation );

    // Set up for the non-linear and transient Darcy solver in the saturation equation.

    // Set up phase for the linear solver.
    saturationDarcySolver.setup();

    // Create the inverse permeability.
    inversePermeability < RegionMesh > invPermSat ( Members->getSaturationPermeability(),
                                                    saturation_darcy_p_FESpace );

    // Set the inverse of the permeability.
    saturationDarcySolver.setInversePermeability ( invPermSat );

    // Set the initial solution.
    saturationDarcySolver.setInitialPrimal ( Members->getSaturationInitialCondition() );

    // Set the source term.
    saturationDarcySolver.setSource ( Members->getSaturationSource() );

    // Set the boudary conditions.
    saturationDarcySolver.setBC ( bcSaturation );

    // Set the exporter for the solution.
    boost::shared_ptr< Exporter< RegionMesh > > exporter;

    // Shared pointer used in the exporter for the pressure in the pressure equation.
    vector_ptrtype pressureExporter;

    // Shared pointer  in the exporter for the saturation in the saturation equation.
    vector_ptrtype saturationExporter;

    // Type of the exporter.
    std::string const exporterType =  dataFile ( "exporter/type", "ensight");

    // Choose the exporter.
#ifdef HAVE_HDF5
    if ( exporterType.compare ("hdf5") == 0 )
    {
        exporter.reset ( new ExporterHDF5< RegionMesh > ( dataFile, dataFile ( "exporter/file_name", "PressureSaturation"  ) ) );

        // Set directory where to save the solution.
        exporter->setPostDir ( dataFile ( "exporter/folder", "./" ) );

        exporter->setMeshProcId ( meshPtr, Members->comm->MyPID() );
    }
    else
#endif
    {
        if ( exporterType.compare ("none") == 0 )
        {
            exporter.reset ( new ExporterEmpty< RegionMesh > ( dataFile, dataFile ( "exporter/file_name", "PressureSaturation"  ) ) );

            // Set directory where to save the solution.
            exporter->setPostDir ( dataFile ( "exporter/folder", "./" ) );

            exporter->setMeshProcId ( meshPtr, Members->comm->MyPID() );
        }
        else
        {
            exporter.reset ( new ExporterEnsight< RegionMesh > ( dataFile, dataFile ( "exporter/file_name", "PressureSaturation"  ) ) );

            // Set directory where to save the solution.
            exporter->setPostDir ( dataFile ( "exporter/folder", "./" ) );

            exporter->setMeshProcId ( meshPtr, Members->comm->MyPID() );
        }
    }

    // Set the exporter pressure pointer.
    pressureExporter.reset ( new vector_type ( *pressureSolver.primalSolution(), exporter->mapType() ) );

    // Add the pressure variable to the exporter.
    exporter->addVariable ( ExporterData<RegionMesh>::ScalarField, "Pressure",
                            pressure_p_FESpacePtr, pressureExporter,
                            static_cast<UInt> ( 0 ),
                            ExporterData< RegionMesh >::UnsteadyRegime,
                            ExporterData< RegionMesh >::Cell );

    // Set the exporter saturation pointer.
    saturationExporter.reset ( new vector_type ( *saturationHyperbolicSolver.solution(), exporter->mapType() ) );

    // Add the solution to the exporter.
    exporter->addVariable ( ExporterData<RegionMesh>::ScalarField, "Saturation",
                            saturation_hyperbolic_FESpacePtr, saturationExporter,
                            static_cast<UInt> ( 0 ),
                            ExporterData< RegionMesh >::UnsteadyRegime,
                            ExporterData< RegionMesh >::Cell );

    // Display the total number of unknowns.
    pressureSolver.getDisplayer().leaderPrint ( "Number of unknowns : ",
                                                2.*pressure_hybrid_FESpace.map().map (Unique)->NumGlobalElements() +
                                                saturation_hyperbolic_FESpacePtr->map().map (Unique)->NumGlobalElements(), "\n" );

    // Solve the problem.

    // Save the initial primal.

    // Copy the saturation initial solution to the exporter.
    *saturationExporter = *saturationDarcySolver.primalSolution();

    // Save the initial solution into the exporter.
    exporter->postProcess ( dataSaturationHyperbolic.dataTime()->initialTime() );

    // Update the time for the simulation
    dataSaturationDarcyNonLinear.dataTime()->updateTime();

    // Outher loop for the simulation, it starts from \Delta t and end in N \Delta t = T.
    for ( ; dataSaturationDarcyNonLinear.dataTime()->canAdvance(); dataSaturationDarcyNonLinear.dataTime()->updateTime() )
    {

        // The leader process prints the temporal data.
        if ( saturationDarcySolver.getDisplayer().isLeader() )
        {
            dataSaturationDarcyNonLinear.dataTime()->showMe();
        }

        // Solve the pressure equation.

        // Build the linear system and the right hand side.
        pressureSolver.buildSystem();

        // Solve the linear system.
        pressureSolver.solve();

        // Post process of the global pressure and total velocity.
        pressureSolver.computePrimalAndDual();

        // Interpolate the Darcy velocity
        *pressure_dualInterpolated = pressure_uInterpolate_FESpacePtr->feToFEInterpolate ( *pressure_u_FESpacePtr,
                                     *pressureSolver.dualSolution() );

        // Solve the saturation equation.

        // Solve the hyperbolic part of the saturation equation.

        // Set the initial condition for the inner loop
        saturationHyperbolicSolver.setSolution ( saturationDarcySolver.primalSolution()  );

        // Set the time parameters for the hyperbolic part of the saturation equation.

        // Set the initial time.
        dataSaturationHyperbolic.dataTime()->setInitialTime ( dataSaturationDarcyNonLinear.dataTime()->time() );

        // Set the current time as initial time.
        dataSaturationHyperbolic.dataTime()->setTime ( dataSaturationDarcyNonLinear.dataTime()->time() );

        // Define the inner time step.
        Real innerTimeStep ( 0. );

        // Set the end time.
        dataSaturationHyperbolic.dataTime()->setEndTime ( dataSaturationDarcyNonLinear.dataTime()->nextTime() );

        // Define if the current time is the last time step.
        bool isLastTimeStep ( false );

        // Inner loop for the simulation, it starts form N \Delta t and end in ( N + 1 ) \Delta t.
        for ( ; dataSaturationHyperbolic.dataTime()->canAdvance() && !isLastTimeStep; dataSaturationHyperbolic.dataTime()->updateTime() )
        {

            // The leader process prints the temporal data for the inner loop.
            saturationHyperbolicSolver.getDisplayer().leaderPrint ( "Inner loop for sub-temporal iteration for the hyperbolic equation.\n" );

            // Compute the new time step according to the CFL condition.
            //innerTimeStep = saturationHyperbolicSolver.CFL();

            // Check if the time step is consistent, i.e. if innerTimeStep + currentTime < endTime.
            if ( dataSaturationHyperbolic.dataTime()->isLastTimeStep() )
            {
                // Compute the last time step.
                innerTimeStep = dataSaturationHyperbolic.dataTime()->leftTime();

                // This is the last time step in the simulation
                isLastTimeStep = true;
            }

            // Set the new time step in the dataHyperbolic.
            //dataSaturationHyperbolic.dataTime()->setTimeStep( innerTimeStep );

            // The leader process prints the temporal data for the inner loop.
            if ( saturationHyperbolicSolver.getDisplayer().isLeader() )
            {
                dataSaturationHyperbolic.dataTime()->showMe();
            }

            // solve one step of the hyperbolic problem.
            saturationHyperbolicSolver.solveOneTimeStep();

        }

        // Solve the parabolic part of the saturation equation.

        // Set the new previous time step solution.
        saturationDarcySolver.setPrimalSolution ( saturationHyperbolicSolver.solution() );
        saturationDarcySolver.updatePrimalOldSolution();

        // Start the fixed point simulation.
        saturationDarcySolver.fixedPointScheme();

        // Save the solution.

        // Copy the pressure solution to the exporter.
        *pressureExporter = *pressureSolver.primalSolution();

        // Copy the saturation solution to the exporter.
        *saturationExporter = *saturationDarcySolver.primalSolution();

        // Save the solution into the exporter.
        exporter->postProcess ( dataSaturationDarcyNonLinear.dataTime()->time() );
    }


    // Stop chronoProcess.
    chronoProcess.stop();

    // The leader process print chronoProcess.
    pressureSolver.getDisplayer().leaderPrint ( "Time for process ",
                                                chronoProcess.diff(), "\n" );

    // Stop chronoTotal.
    chronoTotal.stop();

    // The leader process print chronoTotal.
    pressureSolver.getDisplayer().leaderPrint ( "Total time for the computation ",
                                                chronoTotal.diff(), "\n" );

    return 0.;

}
