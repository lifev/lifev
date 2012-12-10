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
   \file hyperbolic.cpp
   \author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
   \date 2010-07-29
 */

/*!
  Simple 3D hyperbolic test with Dirichlet, Neumann and Robin Boundary condition.
*/

// ===================================================
//! Includes
// ===================================================

#include "hyperbolic.hpp"
// ===================================================
//! Namespaces & define
// ===================================================

using namespace LifeV;

enum BCNAME
{
    BACK   = 1,
    FRONT  = 2,
    LEFT   = 3,
    RIGHT  = 4,
    BOTTOM = 5,
    TOP    = 6
};

// ===================================================
//! User functions
// ===================================================
namespace dataProblem
{
const Real Pi = 3.14159265358979323846264338328;

// Analytical solution
Real analyticalSolution( const Real& t,
                         const Real& x,
                         const Real& y,
                         const Real& /* z */,
                         const ID&   /* ic */)
{
    //Real u1 = 5;
    //if ( ( x <= u1*(t+0.25)/2 && x >= u1*(t - 0.15)/2 ) && ( y >= 0.25 && y <= 0.75 )  && ( z >= 0.25 && z <= 0.75 ) )
    ///   return 1;
    //else
    //   return 0;
    return 20.*std::exp( -( std::pow( x*std::cos(2.*Pi*t) + y*std::sin(2.*Pi*t), 2.) + std::pow( - x*std::sin(2.*Pi*t) + y*std::cos(2.*Pi*t), 2. ) ) );
}

// Physical flux function
Vector physicalFlux( const Real& /* t */,
                     const Real& x,
                     const Real& y,
                     const Real& /* z */,
                     const std::vector<Real>& u )
{
    Vector physicalFluxVector( static_cast<UInt>(3) );

    // First row
    Real Entry0 = y*u[0]*u[1];//u[0]*u[0]*u[1];

    // Second row
    Real Entry1 = x*u[0]*u[2];//u[2];

    // Third row
    Real Entry2 = u[3];

    physicalFluxVector( static_cast<UInt>(0) ) = Entry0;
    physicalFluxVector( static_cast<UInt>(1) ) = Entry1;
    physicalFluxVector( static_cast<UInt>(2) ) = Entry2;

    return physicalFluxVector;
}

// First derivative in u of the physical flux function
Vector firstDerivativePhysicalFlux( const Real& /* t */,
                                    const Real& x,
                                    const Real& y,
                                    const Real& /* z */,
                                    const std::vector<Real>& u )
{
    Vector firstDerivativePhysicalFluxVector( static_cast<UInt>(3) );

    // First row
    Real Entry0 = y*u[1];//2.*u[0]*u[1];

    // Second row
    Real Entry1 = x*u[2];//0.;

    // Third row
    Real Entry2 = 0.;

    firstDerivativePhysicalFluxVector( static_cast<UInt>(0) ) = Entry0;
    firstDerivativePhysicalFluxVector( static_cast<UInt>(1) ) = Entry1;
    firstDerivativePhysicalFluxVector( static_cast<UInt>(2) ) = Entry2;

    return firstDerivativePhysicalFluxVector;
}

// Initial time condition
Real dual( const Real& /* t */,
           const Real& /* x */,
           const Real& /* y */,
           const Real& /* z */,
           const ID&   ic )
{
    switch ( ic )
    {
    case 0:
        return -2.*Pi;//5
        break;
    case 1:
        return 2.*Pi; //0
        break;
    case 2:
        return 0.;
        break;
    }
    return 0.;
}

// Initial time condition
Real initialCondition( const Real& /* t */,
                       const Real& x,
                       const Real& y,
                       const Real& z,
                       const ID&   ic )
{
    return analyticalSolution( 0., x, y, z, ic );
}

// Mass function
Real mass( const Real& /* t */,
           const Real& /* x */,
           const Real& /* y */,
           const Real& /* z */,
           const ID&   /* ic */ )
{
    return 1.;
}

// Boundary condition of Dirichlet
Real dirichlet( const Real& t,
                const Real& x,
                const Real& y,
                const Real& z,
                const ID&   ic )
{
    return analyticalSolution( t, x, y, z, ic );
}


// Source term
Real source_in( const Real& /* t */,
                const Real& /* x */,
                const Real& /* y */,
                const Real& /* z */,
                const ID&   /* icomp */)
{
    return 0.;
}

// Standard functions
Real UOne( const Real& /* t */,
           const Real& /* x */,
           const Real& /* y */,
           const Real& /* z */,
           const ID&   /* icomp */)
{
    return 1.;
}

Real UZero( const Real& /* t */,
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

struct hyperbolic::Private
{
    Private() {}

    // Policy for scalar functions
    typedef boost::function<Real ( const Real&, const Real&,
                                   const Real&, const Real&, const ID& )>
    fct_type;

    // Policy for the flux function
    typedef boost::function<Vector ( const Real&, const Real&,
                                     const Real&, const Real&,
                                     const std::vector<Real>& )>
    vectorFct_type;

    std::string    data_file_name;
    std::string    discretization_section;

    boost::shared_ptr<Epetra_Comm>   comm;

    // Function Types

    fct_type getUOne()
    {
        fct_type f;
        f = boost::bind( &dataProblem::UOne, _1, _2, _3, _4, _5 );
        return f;
    }

    fct_type getUZero()
    {
        fct_type f;
        f = boost::bind( &dataProblem::UZero, _1, _2, _3, _4, _5 );
        return f;
    }

    fct_type getAnalyticalSolution()
    {
        fct_type f;
        f = boost::bind( &dataProblem::analyticalSolution, _1, _2, _3, _4, _5 );
        return f;
    }

    vectorFct_type getPhysicalFlux()
    {
        vectorFct_type f;
        f = boost::bind( &dataProblem::physicalFlux, _1, _2, _3, _4, _5 );
        return f;
    }

    vectorFct_type getFirstDerivativePhysicalFlux()
    {
        vectorFct_type f;
        f = boost::bind( &dataProblem::firstDerivativePhysicalFlux, _1, _2, _3, _4, _5 );
        return f;
    }

    fct_type getSource ( )
    {
        fct_type f;
        f = boost::bind( &dataProblem::source_in, _1, _2, _3, _4, _5 );
        return f;
    }

    fct_type getInitialCondition ( )
    {
        fct_type f;
        f = boost::bind( &dataProblem::initialCondition, _1, _2, _3, _4, _5 );
        return f;
    }

    fct_type getMass ( )
    {
        fct_type f;
        f = boost::bind( & dataProblem::initialCondition, _1, _2, _3, _4, _5 );
        return f;
    }

    fct_type getDual ( )
    {
        fct_type f;
        f = boost::bind( & dataProblem::dual, _1, _2, _3, _4, _5 );
        return f;
    }

};

// ===================================================
//! Constructors
// ===================================================

hyperbolic::hyperbolic( int argc,
                        char** argv )
        : Members( new Private )
{
    GetPot command_line(argc, argv);
    const string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );

    Members->data_file_name = data_file_name;
    Members->discretization_section = "hyperbolic";

#ifdef EPETRA_MPI
    std::cout << "Epetra Initialization" << std::endl;
    Members->comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );
#else
    Members->comm.reset( new Epetra_SerialComm() );
#endif

}

// ===================================================
//! Methods
// ===================================================

Real
hyperbolic::run()
{
    typedef RegionMesh<LinearTetra>                     RegionMesh;
    typedef SolverAztecOO                               solver_type;
    typedef HyperbolicSolver< RegionMesh, solver_type > hyper;
    typedef hyper::vector_Type                          vector_type;
    typedef boost::shared_ptr<vector_type>              vector_ptrtype;
    typedef FESpace< RegionMesh, MapEpetra >            feSpace_Type;
    typedef boost::shared_ptr< feSpace_Type >           feSpacePtr_Type;

    LifeChrono chronoTotal;
    LifeChrono chronoReadAndPartitionMesh;
    LifeChrono chronoBoundaryCondition;
    LifeChrono chronoFiniteElementSpace;
    LifeChrono chronoProblem;
    LifeChrono chronoProcess;
    LifeChrono chronoTimeStep;
    LifeChrono chronoError;

    // Start chronoTotal for measure the total time for the computation
    chronoTotal.start();

    // Reading from data file
    GetPot dataFile( Members->data_file_name.c_str() );

    // Create the leader process, i.e. the process with MyPID equal to zero
    bool isLeader = ( Members->comm->MyPID() == 0 );

    //
    // The Hyperbolic Solver
    //

    if ( isLeader )
        std::cout << "The hyperbolic solver" << std::endl << std::flush;

    // Start chronoReadAndPartitionMesh for measure the total time for the creation of the local meshes
    chronoReadAndPartitionMesh.start();

    // Create the data file
    HyperbolicData<RegionMesh> dataHyperbolic;

    // Set up the data
    dataHyperbolic.setup( dataFile );

    // Create the mesh file handler
    MeshData meshData;

    // Set up the mesh file
    meshData.setup( dataFile,  Members->discretization_section + "/space_discretization");

    // Create the mesh
    boost::shared_ptr<RegionMesh> fullMeshPtr( new RegionMesh( Members->comm ) );

    // Set up the mesh
    readMesh( *fullMeshPtr, meshData );

    // Partition the mesh using ParMetis
    boost::shared_ptr<RegionMesh> meshPtr;
    MeshPartitioner<RegionMesh>::GhostEntityDataMap_Type ghostDataMap;
    {
        MeshPartitioner< RegionMesh >  meshPart( fullMeshPtr, Members->comm );
        meshPtr = meshPart.meshPartition();
        ghostDataMap = meshPart.ghostDataMap();
    }

    // Stop chronoReadAndPartitionMesh
    chronoReadAndPartitionMesh.stop();

    // The leader process print chronoReadAndPartitionMesh
    if ( isLeader )
        std::cout << "Time for read and partition the mesh " <<
                  chronoReadAndPartitionMesh.diff() << std::endl << std::flush;

    // Create the boundary conditions

    // Start chronoBoundaryCondition for measure the total time for create the boundary conditions
    chronoBoundaryCondition.start();

    BCFunctionBase dirichletBDfun;

    dirichletBDfun.setFunction( dataProblem::dirichlet );

    BCHandler bcHyperbolic;

    bcHyperbolic.addBC(   "Top",    TOP,    Essential,  Scalar,  dirichletBDfun  );
    bcHyperbolic.addBC("Bottom", BOTTOM,    Essential,  Scalar,  dirichletBDfun  );
    bcHyperbolic.addBC(  "Left",   LEFT,    Essential,  Scalar,  dirichletBDfun  );
    bcHyperbolic.addBC( "Right",  RIGHT,    Essential,  Scalar,  dirichletBDfun  );
    bcHyperbolic.addBC( "Front",  FRONT,    Essential,  Scalar,  dirichletBDfun  );
    bcHyperbolic.addBC(  "Back",   BACK,    Essential,  Scalar,  dirichletBDfun  );

    // Stop chronoBoundaryCondition
    chronoBoundaryCondition.stop();

    // The leader process print chronoBoundaryCondition
    if ( isLeader )
    {
        std::cout << "Time for create the boundary conditions handler " <<
                  chronoBoundaryCondition.diff() << std::endl << std::flush;

    }

    // Create the solution spaces

    // Start chronoFiniteElementSpace for measure the total time for create the finite element spaces
    chronoFiniteElementSpace.start();

    // Primal solution parameters
    const ReferenceFE*    refFE ( static_cast<ReferenceFE*>(NULL) );
    const QuadratureRule* qR    ( static_cast<QuadratureRule*>(NULL) );
    const QuadratureRule* bdQr  ( static_cast<QuadratureRule*>(NULL) );

    refFE = &feTetraP0;
    qR    = &quadRuleTetra15pt;
    bdQr  = &quadRuleTria1pt;

    // Interpolate of dual solution parameters.
    const ReferenceFE*    pressure_refFE_dualInterpolate ( static_cast<ReferenceFE*>(NULL) );
    const QuadratureRule* pressure_qR_dualInterpolate    ( static_cast<QuadratureRule*>(NULL) );
    const QuadratureRule* pressure_bdQr_dualInterpolate  ( static_cast<QuadratureRule*>(NULL) );

    pressure_refFE_dualInterpolate = &feTetraP0;
    pressure_qR_dualInterpolate    = &quadRuleTetra15pt;
    pressure_bdQr_dualInterpolate  = &quadRuleTria4pt;

    // Finite element space of the interpolation of dual variable.
    FESpace< RegionMesh, MapEpetra > pressure_uInterpolate_FESpace( meshPtr, *pressure_refFE_dualInterpolate, *pressure_qR_dualInterpolate,
                                                                    *pressure_bdQr_dualInterpolate, 3, Members->comm );

    // Vector for the interpolated dual solution.
    vector_ptrtype pressure_dualInterpolated( new vector_type ( pressure_uInterpolate_FESpace.map(), Repeated ) );

    pressure_uInterpolate_FESpace.interpolate( static_cast<FESpace< RegionMesh, MapEpetra >::function_Type>( dataProblem::dual ), *pressure_dualInterpolated, 0 );

    // Finite element space
    feSpacePtr_Type feSpacePtr( new feSpace_Type( meshPtr,
                                                  *refFE,
                                                  *qR,
                                                  *bdQr,
                                                  1,
                                                  Members->comm ) );

    // Stop chronoFiniteElementSpace
    chronoFiniteElementSpace.stop();

    // The leader process print chronoFiniteElementSpace
    if ( isLeader )
        std::cout << "Time for create the finite element spaces " <<
                  chronoFiniteElementSpace.diff() << std::endl << std::flush;

    // Start chronoProblem for measure the total time for create the problem
    chronoProblem.start();

    // Instantiation of the HyperbolicSolver class
    hyper hyperbolicSolver ( dataHyperbolic,
                             *feSpacePtr,
                             Members->comm );

    // Stop chronoProblem
    chronoProblem.stop();

    // The leader process print chronoProblem
    hyperbolicSolver.getDisplayer().leaderPrint( "Time for create the problem ",
                                                 chronoProblem.diff(), "\n" );

    // Process the problem

    // Start chronoProcess for measure the total time for the simulation
    chronoProcess.start();

    // Setup phase
    hyperbolicSolver.setup();

    // Set the source term
    hyperbolicSolver.setSourceTerm( Members->getSource() );

    // Set the initial solution
    hyperbolicSolver.setInitialSolution( Members->getInitialCondition() );

    // Set the mass function
    hyperbolicSolver.setMassTerm( Members->getMass() );

    // Create the numerical flux.
    GodunovNumericalFlux < RegionMesh > numericalFlux ( Members->getPhysicalFlux(),
                                                        Members->getFirstDerivativePhysicalFlux(),
                                                        pressure_uInterpolate_FESpace,
                                                        dataFile,
                                                        "hyperbolic/numerical_flux/" );

    // Set the dependence field
    numericalFlux.setExternalField ( pressure_dualInterpolated );

    // Set the numerical flux usign the physical flux
    hyperbolicSolver.setNumericalFlux( numericalFlux );

    // Set the boudary conditions
    hyperbolicSolver.setBoundaryCondition( bcHyperbolic );

    // Set the exporter for the solution
    boost::shared_ptr< Exporter< RegionMesh > > exporter;

    // Shared pointer used in the exporter for the solution
    vector_ptrtype exporterSolution;

    // Type of the exporter
    std::string const exporterType =  dataFile( "exporter/type", "ensight");

    // Choose the exporter
#ifdef HAVE_HDF5
    if ( exporterType.compare("hdf5") == 0 )
    {
        exporter.reset( new ExporterHDF5< RegionMesh > ( dataFile, dataFile( "exporter/file_name", "Concentration" ) ) );

        // Set directory where to save the solution
        exporter->setPostDir( dataFile( "exporter/folder", "./" ) );

        exporter->setMeshProcId( meshPtr, Members->comm->MyPID() );
    }
    else
#endif
    {
        if ( exporterType.compare("none") == 0 )
        {
            exporter.reset( new ExporterEmpty< RegionMesh > ( dataFile, dataFile( "exporter/file_name", "Concentration" ) ) );

            // Set directory where to save the solution
            exporter->setPostDir( dataFile( "exporter/folder", "./" ) );

            exporter->setMeshProcId( meshPtr, Members->comm->MyPID() );
        }
        else
        {
            exporter.reset( new ExporterEnsight< RegionMesh > ( dataFile, dataFile( "exporter/file_name", "Concentration" ) ) );

            // Set directory where to save the solution
            exporter->setPostDir( dataFile( "exporter/folder", "./" ) );

            exporter->setMeshProcId( meshPtr, Members->comm->MyPID() );
        }
    }

    // Export the partitioning
    exporter->exportPID( meshPtr, Members->comm );

    // export the flags set on the mesh
    exporter->exportFlags( meshPtr, Members->comm );

    // Set the exporter solution
    exporterSolution.reset( new vector_type ( *hyperbolicSolver.solution(),
                                              exporter->mapType() ) );

    // Add the solution to the exporter
    exporter->addVariable( ExporterData<RegionMesh>::ScalarField,
                           "Concentration", feSpacePtr,
                           exporterSolution,
                           static_cast<UInt>( 0 ),
                           ExporterData<RegionMesh>::UnsteadyRegime,
                           ExporterData<RegionMesh>::Cell );

    // Display the total number of unknowns
    hyperbolicSolver.getDisplayer().leaderPrint( "Number of unknowns : ",
                                                 feSpacePtr->map().map(Unique)->NumGlobalElements(), "\n" );

    // Solve the problem

    // Save the initial primal

    // Copy the initial solution to the exporter
    *exporterSolution = *hyperbolicSolver.solution();

    // Save the initial solution into the exporter
    exporter->postProcess( dataHyperbolic.dataTime()->initialTime() );

    // Changing time step for the simulation
    Real timeStep(0.);

    // Flag for the last time step that does not coincide with the last advance
    bool isLastTimeStep( false );

    // A loop for the simulation, it starts from \Delta t and end in N \Delta t = T
    while ( dataHyperbolic.dataTime()->canAdvance() && !isLastTimeStep )
    {

        // Start chronoTimeStep for measure the time for the current time step
        chronoTimeStep.start();

        // update ghost values from neighboring processes
        hyperbolicSolver.updateGhostValues( ghostDataMap );

        // Check if the time step is consistent, i.e. if innerTimeStep + currentTime < endTime.
        if ( dataHyperbolic.dataTime()->isLastTimeStep() )
        {
            // Compute the last time step.
            timeStep = dataHyperbolic.dataTime()->leftTime();

            // This is the last time step in the simulation
            isLastTimeStep = true;
        }
        else
        {
            // Compute the new time step according to the CFL condition.
            timeStep = hyperbolicSolver.CFL();
        }

        // Set the last time step for the simulation.
        dataHyperbolic.dataTime()->setTimeStep( timeStep );

        // Update time
        dataHyperbolic.dataTime()->updateTime();

        // The leader process prints the temporal data.
        if ( hyperbolicSolver.getDisplayer().isLeader() )
        {
            dataHyperbolic.dataTime()->showMe();
        }

        // solve one step of the hyperbolic problem.
        hyperbolicSolver.solveOneTimeStep();

        // Save the solution

        // Copy the solution to the exporter
        *exporterSolution = *hyperbolicSolver.solution();

        // update the total time
        dataHyperbolic.dataTime()->updateTime();

        // Save the solution into the exporter
        exporter->postProcess( dataHyperbolic.dataTime()->time() );

        // Stop chronoTimeStep
        chronoTimeStep.stop();

        // The leader process print chronoTimeStep
        hyperbolicSolver.getDisplayer().leaderPrint( "Time for current time step ",
                                                     chronoTimeStep.diff(), "\n" );

    }

    // Stop chronoProcess
    chronoProcess.stop();

    // The leader process print chronoProcess
    hyperbolicSolver.getDisplayer().leaderPrint( "Time for process ",
                                                 chronoProcess.diff(), "\n" );

    // Compute the errors

    // Start chronoError for measure the total time for computing the errors.
    chronoError.start();

    // Compute the error L2 norms
    Real L2Norm(0), exactL2Norm(0), L2Error(0), L2RelativeError(0);

    // Norms and errors for the pressure
    hyperbolicSolver.getDisplayer().leaderPrint( "\nERROR\n" );

    // Compute the L2 norm for the solution
    L2Norm = feSpacePtr->l2Norm( *hyperbolicSolver.solution() );

    // Display the L2 norm for the solution
    hyperbolicSolver.getDisplayer().leaderPrint( " L2 norm of solution:            ",
                                                 L2Norm, "\n" );

    // Compute the L2 norm for the analytical solution
    exactL2Norm = feSpacePtr->l2NormFunction( Members->getAnalyticalSolution(),
                                              dataHyperbolic.dataTime()->endTime() );

    // Display the L2 norm for the analytical solution
    hyperbolicSolver.getDisplayer().leaderPrint( " L2 norm of exact solution:      ",
                                                 exactL2Norm, "\n" );

    // Compute the L2 error for the solution
    L2Error = feSpacePtr->l2ErrorWeighted( Members->getAnalyticalSolution(),
                                           *hyperbolicSolver.solution(),
                                           Members->getUOne(),
                                           dataHyperbolic.dataTime()->endTime() );

    // Display the L2 error for the solution
    hyperbolicSolver.getDisplayer().leaderPrint( " L2 error:                       ",
                                                 L2Error, "\n" );

    // Compute the L2 realative error for the solution
    L2RelativeError = L2Error / L2Norm;

    // Display the L2 relative error for the solution
    hyperbolicSolver.getDisplayer().leaderPrint( " L2 relative error:              ",
                                                 L2RelativeError, "\n" );

    // Stop chronoError
    chronoError.stop();

    // The leader process print chronoError
    hyperbolicSolver.getDisplayer().leaderPrint( "Time for compute errors ",
                                                 chronoError.diff(), "\n" );

    // Stop chronoTotal
    chronoTotal.stop();

    // The leader process print chronoTotal
    hyperbolicSolver.getDisplayer().leaderPrint( "Total time for the computation ",
                                                 chronoTotal.diff(), "\n" );

    // Return the error, needed for the succes/failure of the test
    return L2Error;

}
