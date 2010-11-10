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
#include <math.h>
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

// Analytical solution
Real analyticalSolution( const Real& t,
                         const Real& x,
                         const Real& y,
                         const Real& z,
                         const ID& /*ic*/)
{
    return 20.*exp( -( pow( x*cos(2.*Pi*t) + y*sin(2.*Pi*t), 2.) + pow( - x*sin(2.*Pi*t) + y*cos(2.*Pi*t), 2. ) ) );
}

// Physical flux function
Vector physicalFlux( const Real& /*t*/,
                     const Real& x,
                     const Real& y,
                     const Real& z,
                     const Real& u )
{
    Vector physicalFluxVector( static_cast<UInt>(3) );

    // First row
    Real Entry0 = -2.*Pi*y*u;

    // Second row
    Real Entry1 = 2.*Pi*x*u;

    // Third row
    Real Entry2 = 0.;

    physicalFluxVector( static_cast<UInt>(0) ) = Entry0;
    physicalFluxVector( static_cast<UInt>(1) ) = Entry1;
    physicalFluxVector( static_cast<UInt>(2) ) = Entry2;

    return physicalFluxVector;
}

// First derivative in u of the physical flux function
Vector firstDerivativePhysicalFlux( const Real& /*t*/,
                     const Real& x,
                     const Real& y,
                     const Real& z,
                     const Real& u )
{
    Vector firstDerivativePhysicalFluxVector( static_cast<UInt>(3) );

    // First row
    Real Entry0 = -2.*Pi*y;

    // Second row
    Real Entry1 = 2.*Pi*x;

    // Third row
    Real Entry2 = 0.;

    firstDerivativePhysicalFluxVector( static_cast<UInt>(0) ) = Entry0;
    firstDerivativePhysicalFluxVector( static_cast<UInt>(1) ) = Entry1;
    firstDerivativePhysicalFluxVector( static_cast<UInt>(2) ) = Entry2;

    return firstDerivativePhysicalFluxVector;
}

// Initial time condition
Real initialCondition( const Real& /*t*/,
                       const Real& x,
                       const Real& y,
                       const Real& z,
                       const ID&   ic )
{
      return analyticalSolution( 0., x, y, z, ic );
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
Real source_in( const Real& /*t*/,
                const Real& x,
                const Real& y,
                const Real& /*z*/,
                const ID&  /*icomp*/)
{
    return -2.*x*x - 4.*y*y - 8.*x*y;
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
    typedef boost::function<Vector ( const Real&, const Real&, const Real&,
                                     const Real&, const Real& )>
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

};

// ===================================================
//! Constructors
// ===================================================

hyperbolic::hyperbolic( int argc,
                        char** argv,
                        LifeV::AboutData const& /*ad*/,
                        LifeV::po::options_description const& /*od*/ ): Members( new Private )
{
    GetPot command_line(argc, argv);
    const string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );

    Members->data_file_name = data_file_name;
    Members->discretization_section = "hyperbolic";

	#ifdef EPETRA_MPI
        std::cout << "Epetra Initialization" << std::endl;
		Members->comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );
        int ntasks;
        MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
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
    typedef RegionMesh3D<LinearTetra>                   RegionMesh;
    typedef SolverTrilinos                              solver_type;
    typedef HyperbolicSolver< RegionMesh, solver_type > hyper;
    typedef hyper::vector_type                          vector_type;
    typedef boost::shared_ptr<vector_type>              vector_ptrtype;

    Chrono chronoTotal;
    Chrono chronoReadAndPartitionMesh;
    Chrono chronoBoundaryCondition;
    Chrono chronoFiniteElementSpace;
    Chrono chronoProblem;
    Chrono chronoProcess;
    Chrono chronoError;

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
    DataHyperbolic<RegionMesh> dataHyperbolic;

    // Set up the data
    dataHyperbolic.setup( dataFile );

    // Create the mesh file handler
    DataMesh dataMesh;

    // Set up the mesh file
    dataMesh.setup( dataFile,  Members->discretization_section + "/space_discretization");

    // Create the mesh
    boost::shared_ptr<RegionMesh> fullMeshPtr( new RegionMesh );

    // Set up the mesh
    readMesh( *fullMeshPtr, dataMesh );

    // Partition the mesh using ParMetis
    partitionMesh< RegionMesh >  meshPart( fullMeshPtr, Members->comm );

    // Stop chronoReadAndPartitionMesh
    chronoReadAndPartitionMesh.stop();

    // The leader process print chronoReadAndPartitionMesh
    if( isLeader )
        std::cout << "Time for read and partition the mesh " <<
                     chronoReadAndPartitionMesh.diff() << std::endl << std::flush;

    // Create the boundary conditions

    // Start chronoBoundaryCondition for measure the total time for create the boundary conditions
    chronoBoundaryCondition.start();

    BCFunctionBase dirichletBDfun;

    dirichletBDfun.setFunction( dataProblem::dirichlet );

	BCHandler bcHyperbolic( 6 );

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
    const RefFE*    refFE ( static_cast<RefFE*>(NULL) );
    const QuadRule* qR    ( static_cast<QuadRule*>(NULL) );
    const QuadRule* bdQr  ( static_cast<QuadRule*>(NULL) );

    refFE = &feTetraP0;
    qR    = &quadRuleTetra15pt;
    bdQr  = &quadRuleTria1pt;

    // Finite element space
    FESpace< RegionMesh, EpetraMap > fESpace( meshPart,
                                              *refFE,
                                              *qR,
                                              *bdQr,
                                              1,
                                              Members->comm );

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
                             fESpace,
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

    // Set the numerical flux usign the physical flux
    hyperbolicSolver.setNumericalFlux( Members->getPhysicalFlux(),
                                       Members->getFirstDerivativePhysicalFlux() );

    // Set the boudary conditions
    hyperbolicSolver.setBC( bcHyperbolic );

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
        exporter.reset( new Hdf5exporter< RegionMesh > ( dataFile, "Concentration" ) );

        // Set directory where to save the solution
        exporter->setDirectory( dataFile( "exporter/folder", "./" ) );

        exporter->setMeshProcId( meshPart.mesh(), Members->comm->MyPID() );
    }
    else
#endif
    {
        if ( exporterType.compare("none") == 0 )
        {
            exporter.reset( new NoExport< RegionMesh > ( dataFile, "Concentration" ) );

            // Set directory where to save the solution
            exporter->setDirectory( dataFile( "exporter/folder", "./" ) );

            exporter->setMeshProcId( meshPart.mesh(), Members->comm->MyPID() );
        }
        else
        {
            exporter.reset( new Ensight< RegionMesh > ( dataFile, "Concentration" ) );

            // Set directory where to save the solution
            exporter->setDirectory( dataFile( "exporter/folder", "./" ) );

            exporter->setMeshProcId( meshPart.mesh(), Members->comm->MyPID() );
        }
    }

    // Set the exporter solution
    exporterSolution.reset( new vector_type ( *hyperbolicSolver.solution(),
                                              exporter->mapType() ) );

    // Add the solution to the exporter
    exporter->addVariable( ExporterData::Scalar,
                           "Concentration",
                           exporterSolution,
                           static_cast<UInt>( 0 ),
                           static_cast<UInt>( fESpace.dof().numTotalDof() ),
                           static_cast<UInt>( 0 ),
                           ExporterData::Cell );

    // Display the total number of unknowns
    hyperbolicSolver.getDisplayer().leaderPrint( "Number of unknowns : ",
                                                 fESpace.map().getMap(Unique)->NumGlobalElements(), "\n" );

    // Solve the problem

    // Save the initial primal

    // Copy the initial solution to the exporter
    *exporterSolution = *hyperbolicSolver.solution();

    // Save the initial solution into the exporter
    exporter->postProcess( dataHyperbolic.dataTime()->getInitialTime() );

    Real timeStep(0.);

    // A loop for the simulation, it starts from \Delta t and end in N \Delta t = T
    while( !dataHyperbolic.dataTime()->isLastTimeStep() )
    {

        // Compute the new time step according to the CFL condition.
        timeStep = hyperbolicSolver.CFL();

        // Set the new time step in the dataHyperbolic.
        dataHyperbolic.dataTime()->setTimeStep( timeStep );

        // Advance the current time of \Delta t.
        dataHyperbolic.dataTime()->updateTime();

        // The leader process prints the temporal data.
        if ( hyperbolicSolver.getDisplayer().isLeader() )
        {
            dataHyperbolic.dataTime()->showMe();
        }

        // solve one step of the hyperbolic problem.
        hyperbolicSolver.solveOneStep();

        // Save the solution

        // Copy the solution to the exporter
        *exporterSolution = *hyperbolicSolver.solution();

        // Save the solution into the exporter
        exporter->postProcess( dataHyperbolic.dataTime()->getTime() );

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
    L2Norm = fESpace.L2Norm( *hyperbolicSolver.solution() );

    // Display the L2 norm for the solution
    hyperbolicSolver.getDisplayer().leaderPrint( " L2 norm of solution:            ",
                                                 L2Norm, "\n" );

    // Compute the L2 norm for the analytical solution
    exactL2Norm = fESpace.L2NormFunction( Members->getAnalyticalSolution(),
                                          dataHyperbolic.dataTime()->getEndTime() );

    // Display the L2 norm for the analytical solution
    hyperbolicSolver.getDisplayer().leaderPrint( " L2 norm of exact solution:              ",
                                                 exactL2Norm, "\n" );

    // Compute the L2 error for the solution
    L2Error = fESpace.L2ErrorWeighted( Members->getAnalyticalSolution(),
                                       *hyperbolicSolver.solution(),
                                       Members->getUOne(),
                                       dataHyperbolic.dataTime()->getEndTime() );

    // Display the L2 error for the solution
    hyperbolicSolver.getDisplayer().leaderPrint( " L2 error:           ",
                                                 L2Error, "\n" );

    // Compute the L2 realative error for the solution
    L2RelativeError = L2Error / L2Norm;

    // Display the L2 relative error for the solution
    hyperbolicSolver.getDisplayer().leaderPrint( " L2 relative error:  ",
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
