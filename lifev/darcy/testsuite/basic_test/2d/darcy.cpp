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
/**
   @file darcy.cpp
   @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
   @date 2012-06-13
 */

// ===================================================
//! Includes
// ===================================================

#include "darcy.hpp"
#include "user_fun.hpp"

// ===================================================
//! Namespaces & define
// ===================================================

using namespace LifeV;

// ===================================================
//!              Standard functions
// ===================================================

Real UOne ( const Real& /* t */,
            const Real& /* x */,
            const Real& /* y */,
            const Real& /* z */,
            const ID&   /* icomp */)
{
    return 1.;
}

// ===================================================
//!                  Private Members
// ===================================================

struct darcy_nonlinear::Private
{
    Private() {}

    // Policy for scalar functions
    typedef boost::function < Real ( const Real&, const Real&,
                                     const Real&, const Real&, const ID& ) >
    fct_type;

    std::string data_file_name;
    std::string xml_file_name;
    std::string discretization_section;

    boost::shared_ptr<Epetra_Comm>   comm;

    // Function Types

    fct_type getUOne ( )
    {
        fct_type f;
        f = boost::bind ( &UOne, _1, _2, _3, _4, _5 );
        return f;
    }

    fct_type getAnalyticalSolution ( )
    {
        fct_type f;
        f = boost::bind ( &dataProblem::analyticalSolution, _1, _2, _3, _4, _5 );
        return f;
    }

    fct_type getAnalyticalFlux ( )
    {
        fct_type f;
        f = boost::bind ( &dataProblem::analyticalFlux, _1, _2, _3, _4, _5 );
        return f;
    }

};

// ===================================================
//!                  Constructors
// ===================================================

darcy_nonlinear::darcy_nonlinear ( int argc, char** argv )
    : Members ( new Private )
{
    GetPot command_line (argc, argv);
    Members->data_file_name = command_line.follow ("data", 2, "-f", "--file");
    Members->xml_file_name = command_line.follow ("parameterList.xml", "--xml");

    GetPot dataFile ( Members->data_file_name );

    Members->discretization_section = "darcy";

#ifdef EPETRA_MPI
    Members->comm.reset ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
    Members->comm.reset ( new Epetra_SerialComm() );
#endif

}

// ===================================================
//!                      Methods
// ===================================================

Real
darcy_nonlinear::run()
{
    using namespace dataProblem;

    // Life chonos
    LifeChrono chronoTotal;
    LifeChrono chronoReadAndPartitionMesh;
    LifeChrono chronoBoundaryCondition;
    LifeChrono chronoFiniteElementSpace;
    LifeChrono chronoProblem;
    LifeChrono chronoProcess;
    LifeChrono chronoError;

    // Displayer to print on screen
    boost::shared_ptr < Displayer > displayer ( new Displayer ( Members->comm ) );

    // Start chronoTotal for measure the total time for the computation
    chronoTotal.start();

    // Reading from data file
    GetPot dataFile ( Members->data_file_name.c_str() );

    //
    // The Darcy Solver
    //

    displayer->leaderPrint ( "The Darcy solver\n" );

    // Start chronoReadAndPartitionMesh for measure the total time for the creation of the local meshes
    chronoReadAndPartitionMesh.start();

    // Create the data file
    darcyDataPtr_Type darcyData ( new darcyData_Type );

    // Set up the data
    darcyData->setup ( dataFile );

    // Create a Teuchos parameter list for the linear algebra.
    darcyData_Type::paramListPtr_Type linearAlgebraList = Teuchos::rcp ( new darcyData_Type::paramList_Type );
    linearAlgebraList = Teuchos::getParametersFromXmlFile ( Members->xml_file_name );

    // Set the parameter list into the data darcy.
    darcyData->setLinearAlgebraList ( linearAlgebraList );

    // Create the mesh file handler
    MeshData meshData;

    // Set up the mesh file
    meshData.setup ( dataFile,  Members->discretization_section + "/space_discretization");

    // Create the the mesh
    regionMeshPtr_Type fullMeshPtr ( new regionMesh_Type ( Members->comm ) );

    // Select if the mesh is structured or not
    if ( meshData.meshType() != "structured" )
    {
        // Set up the mesh
        readMesh ( *fullMeshPtr, meshData );
    }
    else
    {
        // Section of the structured mesh
        const std::string structuredSection = Members->discretization_section + "/space_discretization/";

        // Set up the structured mesh
        regularMesh2D ( *fullMeshPtr, 0,
                        dataFile ( ( structuredSection + "nx" ).data(), 4 ),
                        dataFile ( ( structuredSection + "ny" ).data(), 4 ),
                        dataFile ( ( structuredSection + "verbose" ).data(), false ),
                        dataFile ( ( structuredSection + "lx" ).data(), 1. ),
                        dataFile ( ( structuredSection + "ly" ).data(), 1. ) );
    }

    // Create the local mesh ( local scope used to delete the meshPart object )
    regionMeshPtr_Type meshPtr;
    {
        // Create the partitioner
        MeshPartitioner< regionMesh_Type >  meshPart;

        // Partition the mesh using ParMetis
        meshPart.doPartition ( fullMeshPtr, Members->comm );

        // Get the local mesh
        meshPtr = meshPart.meshPartition();
    }

    // Stop chronoReadAndPartitionMesh
    chronoReadAndPartitionMesh.stop();

    // The leader process print chronoReadAndPartitionMesh
    displayer->leaderPrint ( "Time for read and partition the mesh ",
                             chronoReadAndPartitionMesh.diff(), "\n" );

    // Create the boundary conditions

    // Start chronoBoundaryCondition for measure the total time for create the boundary conditions
    chronoBoundaryCondition.start();

    bcHandlerPtr_Type bcDarcy ( new bcHandler_Type );

    setBoundaryConditions ( bcDarcy );

    // Stop chronoBoundaryCondition
    chronoBoundaryCondition.stop();

    // The leader process print chronoBoundaryCondition
    displayer->leaderPrint ( "Time for create the boundary conditions handler ",
                             chronoBoundaryCondition.diff(), "\n" );

    // Create the solution spaces

    // Start chronoFiniteElementSpace for measure the total time for create the finite element spaces
    chronoFiniteElementSpace.start();

    // Primal solution parameters
    const ReferenceFE* refFE_primal ( static_cast<ReferenceFE*> (NULL) );
    const QuadratureRule* qR_primal ( static_cast<QuadratureRule*> (NULL) );
    const QuadratureRule* bdQr_primal ( static_cast<QuadratureRule*> (NULL) );

    refFE_primal = &feTriaP0;
    qR_primal = &quadRuleTria4pt;
    bdQr_primal = &quadRuleSeg1pt;

    // Dual solution parameters
    const ReferenceFE* refFE_dual ( static_cast<ReferenceFE*> (NULL) );
    const QuadratureRule* qR_dual ( static_cast<QuadratureRule*> (NULL) );
    const QuadratureRule* bdQr_dual ( static_cast<QuadratureRule*> (NULL) );

    refFE_dual = &feTriaRT0;
    qR_dual = &quadRuleTria4pt;
    bdQr_dual = &quadRuleSeg1pt;

    // Interpolate of dual solution parameters
    const ReferenceFE* refFE_dualInterpolate ( static_cast<ReferenceFE*> (NULL) );
    const QuadratureRule* qR_dualInterpolate ( static_cast<QuadratureRule*> (NULL) );
    const QuadratureRule* bdQr_dualInterpolate ( static_cast<QuadratureRule*> (NULL) );

    refFE_dualInterpolate = &feTriaP0;
    qR_dualInterpolate = &quadRuleTria4pt;
    bdQr_dualInterpolate = &quadRuleSeg1pt;

    // Hybrid solution parameters
    // hybrid
    const ReferenceFE* refFE_hybrid ( static_cast<ReferenceFE*> (NULL) );
    const QuadratureRule* qR_hybrid ( static_cast<QuadratureRule*> (NULL) );
    const QuadratureRule* bdQr_hybrid ( static_cast<QuadratureRule*> (NULL) );

    refFE_hybrid = &feTriaRT0Hyb;
    qR_hybrid = &quadRuleTria4pt;
    bdQr_hybrid = &quadRuleSeg1pt;

    // Finite element space of the primal variable
    fESpacePtr_Type p_FESpacePtr ( new fESpace_Type ( meshPtr, *refFE_primal, *qR_primal,
                                                      *bdQr_primal, 1, Members->comm ) );

    // Finite element space of the dual variable
    fESpacePtr_Type u_FESpacePtr ( new fESpace_Type ( meshPtr, *refFE_dual, *qR_dual,
                                                      *bdQr_dual, 1, Members->comm ) );

    // Finite element space of the hybrid variable
    fESpacePtr_Type hybrid_FESpacePtr ( new fESpace_Type ( meshPtr, *refFE_hybrid, *qR_hybrid,
                                                           *bdQr_hybrid, 1, Members->comm ) );

    // Finite element space of the interpolation of dual variable
    fESpacePtr_Type uInterpolate_FESpacePtr ( new fESpace_Type ( meshPtr, *refFE_dualInterpolate, *qR_dualInterpolate,
                                                                 *bdQr_dualInterpolate, 2, Members->comm ) );

    // Stop chronoFiniteElementSpace
    chronoFiniteElementSpace.stop();

    // The leader process print chronoFiniteElementSpace
    displayer->leaderPrint ( "Time for create the finite element spaces ",
                             chronoFiniteElementSpace.diff(), "\n" );

    // Start chronoProblem for measure the total time for create the problem
    chronoProblem.start();

    // Instantiation of the DarcySolver class
    darcySolver_Type darcySolver;

    // Stop chronoProblem
    chronoProblem.stop();

    // The leader process print chronoProblem
    displayer->leaderPrint ( "Time for create the problem ", chronoProblem.diff(), "\n" );

    // Process the problem

    // Start chronoProcess for measure the total time for the simulation
    chronoProcess.start ();

    // Set the data for the solver.
    darcySolver.setData ( darcyData );

    // Set the displayer.
    darcySolver.setDisplayer ( displayer );

    // Setup phase for the linear solver
    darcySolver.setup ();

    // Create the fields for the solver

    // Scalar field for primal variable
    scalarFieldPtr_Type primalField ( new scalarField_Type ( p_FESpacePtr ) );

    // Scalar field for dual variable
    scalarFieldPtr_Type dualField ( new scalarField_Type ( u_FESpacePtr ) );

    // Scalar field for hybrid variable
    scalarFieldPtr_Type hybridField ( new scalarField_Type ( hybrid_FESpacePtr ) );

    // Vector for the interpolated dual solution
    vectorFieldPtr_Type dualInterpolated ( new vectorField_Type ( uInterpolate_FESpacePtr, Repeated ) );

    // Set the field for the solver

    // Set the h
    darcySolver.setHybridField ( hybridField );

    // Set the primal field
    darcySolver.setPrimalField ( primalField );

    // Set the dual field
    darcySolver.setDualField ( dualField );

    // Set the source term
    scalarFctPtr_Type scalarSourceFct ( new scalarSource );
    darcySolver.setScalarSource ( scalarSourceFct );

    // Set the vector source term
    vectorFctPtr_Type vectorSourceFct ( new vectorSource );
    darcySolver.setVectorSource ( vectorSourceFct );

    // Set the reaction term
    scalarFctPtr_Type reactionTermFct ( new reactionTerm );
    darcySolver.setReactionTerm ( reactionTermFct );

    // Set the inverse of the permeability
    matrixFctPtr_Type inversePermeabilityFct ( new inversePermeability );
    darcySolver.setInversePermeability ( inversePermeabilityFct );

    // Set the initial primal variable
    scalarFctPtr_Type initialPrimalFct ( new initialCondition );
    darcySolver.setInitialPrimal ( initialPrimalFct );

    // Set the mass function
    scalarFctPtr_Type massFct ( new massFunction );
    darcySolver.setMass ( massFct );

    // Set the boudary conditions
    darcySolver.setBoundaryConditions ( bcDarcy );

    // Set the exporter for the solution
    boost::shared_ptr< Exporter< regionMesh_Type > > exporter;

    // Shared pointer used in the exporter for the primal solution
    vectorPtr_Type primalExporter;

    // Shared pointer used in the exporter for the dual solution
    vectorPtr_Type dualExporter;

    // Type of the exporter
    std::string const exporterType = dataFile ( "exporter/type", "none" );

    // The name of the file
    const std::string exporterFileName = dataFile ( "exporter/file_name", "PressureVelocity" );

    // Choose the exporter
#ifdef HAVE_HDF5
    if ( exporterType.compare ("hdf5") == 0 )
    {
        exporter.reset ( new ExporterHDF5 < regionMesh_Type > ( dataFile, exporterFileName ) );
    }
    else
#endif
    {
        exporter.reset ( new ExporterEmpty < regionMesh_Type > ( dataFile, exporterFileName ) );
    }

    // Set directory where to save the solution
    exporter->setPostDir ( dataFile ( "exporter/folder", "./" ) );

    exporter->setMeshProcId ( meshPtr, Members->comm->MyPID() );

    // Set the exporter primal pointer
    primalExporter.reset ( new vector_Type ( primalField->getVector(), exporter->mapType() ) );

    // Add the primal variable to the exporter
    exporter->addVariable ( ExporterData < regionMesh_Type >::ScalarField,
                            dataFile ( "exporter/name_primal", "Pressure" ),
                            p_FESpacePtr,
                            primalExporter,
                            static_cast<UInt> ( 0 ),
                            ExporterData < regionMesh_Type >::UnsteadyRegime,
                            ExporterData < regionMesh_Type >::Cell );

    // Set the exporter dual pointer
    dualExporter.reset ( new vector_Type ( dualInterpolated->getVector(), exporter->mapType() ) );

    // Add the variable to the exporter
    exporter->addVariable ( ExporterData < regionMesh_Type >::VectorField,
                            dataFile ( "exporter/name_dual", "Velocity" ),
                            uInterpolate_FESpacePtr,
                            dualExporter,
                            static_cast<UInt> ( 0 ),
                            ExporterData < regionMesh_Type >::UnsteadyRegime,
                            ExporterData < regionMesh_Type >::Cell );

    // Display the total number of unknowns
    displayer->leaderPrint ( "Number of unknowns : ",
                             hybrid_FESpacePtr->map().map ( Unique )->NumGlobalElements(), "\n" );

    // Export the partitioning
    exporter->exportPID ( meshPtr, Members->comm );

    // Copy the initial primal to the exporter
    *primalExporter = primalField->getVector ();

    // Save the initial primal solution into the exporter
    exporter->postProcess ( darcyData->dataTimePtr()->initialTime() );

    // Define if the current time is the last time step.
    bool isLastTimeStep ( false );

    // Define the end time
    const Real endTime = darcyData->dataTimePtr()->endTime();

    // Define the tolerance for the elapsed time
    const Real tolTime = 1e-10;

    // Take the current time.
    Real currentTime = darcyData->dataTimePtr()->time();

    // A loop for the simulation
    while ( darcyData->dataTimePtr()->time() < endTime && !isLastTimeStep )
    {
        // Take the left time.
        const Real leftTime = endTime - currentTime;

        // Check if the time step is consistent
        if ( darcyData->dataTimePtr()->timeStep() > leftTime + tolTime )
        {
            // Compute the last time step.
            darcyData->dataTimePtr()->setTimeStep ( leftTime );

            // This is the last time step in the simulation
            isLastTimeStep = true;
        }

        // Advance the current time of \Delta t.
        darcyData->dataTimePtr()->updateTime();

        // Take the current time.
        currentTime = darcyData->dataTimePtr()->time();

        // The leader process prints the temporal data.
        if ( displayer->isLeader() )
        {
            darcyData->dataTimePtr()->showMe();
        }

        // Solve the problem.
        darcySolver.solve ();

        // Save the solutions.

        // Copy the primal solution to the exporter.
        *primalExporter = primalField->getVector ();

        // Interpolate the dual vector field spammed as Raviart-Thomas into a P0 vector field.
        dualInterpolated->getVector() = uInterpolate_FESpacePtr->feToFEInterpolate (
                                            *u_FESpacePtr, dualField->getVector() );

        // Copy the dual interpolated solution to the exporter.
        *dualExporter = dualInterpolated->getVector();

        // Save the solution into the exporter.
        exporter->postProcess ( currentTime );

    }

    // Stop chronoProcess
    chronoProcess.stop();

    // The leader process print chronoProcess
    displayer->leaderPrint ( "Time for process ", chronoProcess.diff(), "\n" );

    // Compute the errors
    // For non time dependences problem the time where the errors are computed is useless,
    // but thanks to common interface we handle both cases.

    // Start chronoError for measure the total time for computing the errors.
    chronoError.start();

    // Compute the error L2 norms
    Real primalL2Norm (0), exactPrimalL2Norm (0), primalL2Error (0), primalL2RelativeError (0);
    Real dualL2Norm (0), exactDualL2Norm (0), dualL2Error (0), dualL2RelativeError (0);

    // Norms and errors for the pressure
    displayer->leaderPrint ( "\nPRESSURE ERROR\n" );

    // Compute the L2 norm for the primal solution
    primalL2Norm = p_FESpacePtr->l2Norm ( primalField->getVector () );

    // Display the L2 norm for the primal solution
    displayer->leaderPrint ( " L2 norm of primal unknown:            ", primalL2Norm, "\n" );

    // Compute the L2 norm for the analytical primal
    exactPrimalL2Norm = p_FESpacePtr->l2NormFunction ( Members->getAnalyticalSolution(),
                                                       darcyData->dataTimePtr()->endTime() );

    // Display the L2 norm for the analytical primal
    displayer->leaderPrint ( " L2 norm of primal exact:              ", exactPrimalL2Norm, "\n" );

    // Compute the L2 error for the primal solution
    primalL2Error = p_FESpacePtr->l2ErrorWeighted ( Members->getAnalyticalSolution(),
                                                    primalField->getVector(),
                                                    Members->getUOne(),
                                                    darcyData->dataTimePtr()->endTime() );

    // Display the L2 error for the primal solution
    displayer->leaderPrint ( " L2 error of primal unknown:           ", primalL2Error, "\n" );

    // Compute the L2 realative error for the primal solution
    primalL2RelativeError = primalL2Error / exactPrimalL2Norm;

    // Display the L2 relative error for the primal solution
    displayer->leaderPrint ( " L2 relative error of primal unknown:  ", primalL2RelativeError, "\n" );

    // Norms and errors for the interpolated dual
    displayer->leaderPrint ( "\nINTERPOLATED DARCY VELOCITY ERROR\n" );

    // Compute the L2 norm for the interpolated dual solution
    dualL2Norm = uInterpolate_FESpacePtr->l2Norm ( dualInterpolated->getVector() );

    // Display the L2 norm for the interpolated dual solution
    displayer->leaderPrint ( " L2 norm of dual unknown:              ", dualL2Norm, "\n" );

    // Compute the L2 norm for the analytical dual
    exactDualL2Norm = uInterpolate_FESpacePtr->l2NormFunction ( Members->getAnalyticalFlux(),
                                                                darcyData->dataTimePtr()->endTime() );

    // Display the L2 norm for the analytical dual
    displayer->leaderPrint ( " L2 norm of dual exact:                ", exactDualL2Norm, "\n" );

    // Compute the L2 error for the dual solution
    dualL2Error = uInterpolate_FESpacePtr->l2Error ( Members->getAnalyticalFlux(),
                                                     dualInterpolated->getVector(),
                                                     darcyData->dataTimePtr()->endTime(),
                                                     NULL );

    // Display the L2 error for the dual solution
    displayer->leaderPrint ( " L2 error of dual unknown:             ", dualL2Error, "\n" );

    // Compute the L2 relative error for the dual solution
    dualL2RelativeError = dualL2Error / exactDualL2Norm;

    // Display the L2 relative error for the dual solution
    displayer->leaderPrint ( " L2 relative error of Dual unknown:    ", dualL2RelativeError, "\n" );

    // Stop chronoError
    chronoError.stop();

    // The leader process print chronoError
    displayer->leaderPrint ( "Time for compute errors ", chronoError.diff(), "\n" );

    // Stop chronoTotal
    chronoTotal.stop();

    // The leader process print chronoTotal
    displayer->leaderPrint ( "Total time for the computation ", chronoTotal.diff(), "\n" );

    // Return the error, needed for the succes/failure of the test
    return primalL2Error;

} // run
