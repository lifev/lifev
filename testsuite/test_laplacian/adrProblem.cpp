/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s):  L. Iapichino  <laura.iapichino@epfl.ch>
              C. Malossi    <cristiano.malossi@epfl.ch>
              A. Manzoni    <andrea.manzoni@epfl.ch>
              T. Passerini  <tiziano@mathcs.emory.edu>

  Copyright (C) 2009-2010 EPFL, Emory University

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
   \file adrProblem.cpp
   \author L. Iapichino <laura.iapichino@epfl.ch>
           C. Malossi <cristiano.malossi@epfl.ch>
           A. Manzoni <andrea.manzoni@epfl.ch>
           T. Passerini <tiziano@mathcs.emory.edu>
 */

/* ========================================================

Simple Laplacian test with Dirichlet Boundary condition

Solve the problem

               - \Delta u = f

               u = 0 on the boundary


 3D: with the source term f = 12 \pi^2 sin(2 \pi x) sin(2 \pi y) sin (2 \pi z) on a cube
 2D: with the source term f = 8 \pi^2 sin(2 \pi x) sin(2 \pi y) on a square

 the rhs is computed as rhs = Mass_Matrix * f_iterpolated


 More generally this test can solve the problem:

               - \nu \Delta u + \beta \nabla u + \sigma u = f

               u = g on the boundary

 being \nu and \sigma constants defined in the data file and \beta interpolated.

 */



// ===================================================
//! Includes
// ===================================================
#include <Epetra_Comm.h>

#include <life/lifemesh/dataMesh.hpp>
#include <life/lifemesh/partitionMesh.hpp>

#include <life/lifefilters/ensight.hpp>

#include <life/lifesolver/AdvectionDiffusionReactionSolver.hpp>

#include "analyticalSol.hpp"
#include "adrProblem.hpp"

#define POSTPROCESS 1


// ===================================================
//! Namespaces & Define
// ===================================================
using namespace LifeV;

#ifdef TWODIM

typedef RegionMesh2D<LinearTriangle> mesh_type;

const std::string discretization_section="adr/space_discretization2D";

#elif defined THREEDIM

typedef RegionMesh3D<LinearTetra> mesh_type;

const std::string discretization_section="adr/space_discretization3D";

#endif


// ===================================================
//! Private Members
// ===================================================
struct ADRProblem::Private
{
    Private() {}

    std::string    data_file_name;
    boost::shared_ptr<Epetra_Comm>   commPtr;
    Real squareFunction (const Real& t) { return t*t; }

};





// ===================================================
//! Constructors
// ===================================================
ADRProblem::ADRProblem( int argc,
                        char** argv,
                        LifeV::AboutData const& /*ad*/,
                        LifeV::po::options_description const& /*od*/ ): Members( new Private )
{
    GetPot command_line(argc, argv);
    const string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );

    Members->data_file_name = data_file_name;

#ifdef EPETRA_MPI
    std::cout << "  problem - Epetra Initialization" << std::endl;
    Members->commPtr.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );
#else
    Members->commPtr.reset( new Epetra_SerialComm() );
#endif
}





// ===================================================
//! Methods
// ===================================================
void
ADRProblem::run()
{
    // ===================================================
    // Type definitions
    // ===================================================
    typedef SolverTrilinos                                     linearSolver_type;
    //typedef SolverAmesos                                     linearSolver_type;
    typedef ADRSolver< mesh_type, linearSolver_type >          solver_type;
    typedef solver_type::vector_type                           vector_type;
    typedef boost::shared_ptr<vector_type>                     vector_ptr_type;
    typedef solver_type::timeIntegrator_type                   timeIntegrator_type;
    typedef solver_type::timeIntegrator_ptr_type               timeIntegrator_ptr_type;
    typedef boost::shared_ptr<mesh_type>                       mesh_ptr_type;
    typedef solver_type::data_type                             data_type;
    typedef solver_type::data_ptr_type                         data_ptr_type;
    typedef data_type::dataTime_type                           dataTime_type;
    typedef data_type::dataTime_ptr_type                       dataTime_ptr_type;
    typedef solver_type::fespace_type                          fespace_type;
    typedef solver_type::fespace_ptr_type                      fespace_ptr_type;
    typedef AnalyticalSol                                      solution_type;

    // ===================================================
    // Data structures declarations
    // ===================================================

    // data file descriptor
    GetPot dataFile( Members->data_file_name.c_str() );

    // Controller for verbose level
    bool verbose = ( Members->commPtr->MyPID() == 0 );

    // A chronometer
    Chrono chrono;

    // User-defined parameters for the mesh
    DataMesh dataMesh;

    // The mesh object
    mesh_ptr_type meshPtr( new mesh_type() );

    // User-defined parameters for the time manager
    dataTime_ptr_type dataTimePtr( new dataTime_type() );

    // User-defined parameters for the solver
    data_ptr_type dataADRPtr( new data_type() );

    // The solver
    solver_type adrSolver;

    // The time integrator for unsteady problems
    timeIntegrator_ptr_type timeIntegratorBDFPtr( new timeIntegrator_type() );

    // The list of choices for the analytic solution type
    DataStringList solutionTypeList("Solution type");
    solutionTypeList.add( "steadyPolynomial", ADR_STEADY_POLYNOMIAL, "polynomial solution, steady problem" );
    solutionTypeList.add( "unsteadyPolynomial", ADR_UNSTEADY_POLYNOMIAL, "polynomial solution, unsteady problem" );

    // An instance of AnalyticalSol class
    solution_type uExact;
    // A shortcut to store the time at which to evaluate time dependent functions
    Real timeEval;
    // Error L2 and H1 Norms
    Real H1_Error, H1_RelError, L2_Error, L2_RelError;


    // ===================================================
    // Load user-defined specs
    // ===================================================
    dataMesh.setup( dataFile, discretization_section );
    dataTimePtr->setup( dataFile, "adr/time_discretization" );
    dataADRPtr->setup( dataFile, "adr" );
    adrSolver.setupLinearSolver( dataFile, "adr" );
    solution_type::setup( dataFile, "adr" );

    ADRProblemSolution solutionType = ADRProblemSolution(
            solutionTypeList.value( dataFile( "adr/problem/solutionType", "steadyPolynomial") ) );
    solution_type::problemSolution = solutionType;
    bool steadyProblem( solutionType == ADR_STEADY_POLYNOMIAL );

    if(verbose)
    {
        dataMesh.showMe();
        solution_type::showMe();
        dataTimePtr->showMe();
        dataADRPtr->showMe();
        // timeIntegratorBDFPtr->showMe();
        if( steadyProblem )
            std::cout << "  problem- STEADY problem" << std::endl;
        else
            std::cout << "  problem- UNSTEADY problem" << std::endl;
    }

    // ===================================================
    // Mesh construction and partitioning
    // ===================================================
    readMesh(*meshPtr, dataMesh);

    // partition the mesh across processes
    partitionMesh< mesh_type> meshPart(meshPtr, Members->commPtr);


    // ===================================================
    // Definition of the FE spaces for the problem under study
    // ===================================================
    // Solution Space
    if (verbose) std::cout << "  problem - Building the solution FE space ... " << std::flush;

    fespace_ptr_type uFESpacePtr( new fespace_type( meshPart, dataADRPtr->solFEType(),
                                                    dataADRPtr->solutionFieldDimension(),
                                                    Members->commPtr ) );

    if (verbose) std::cout << "ok." << std::endl;


    // ===================================================
    // Manager for the postprocess
    // ===================================================
#if POSTPROCESS
    boost::shared_ptr< Exporter<mesh_type> > exporter;
    exporter.reset( new Ensight<mesh_type> ( dataFile, meshPart.mesh(),
                                             "adrProblem", Members->commPtr->MyPID()) );
#endif


    // ===================================================
    // Boundary conditions definition
    // ===================================================
    std::string dirichletList = dataFile( "adr/problem/dirichletList", "" );
    std::list<UInt> dirichletMarkers; parseList( dirichletList, dirichletMarkers );
    std::string neumannList = dataFile( "adr/problem/neumannList", "" );
    std::list<UInt> neumannMarkers; parseList( neumannList, neumannMarkers );

    BCHandler::BCHints hint(BCHandler::HINT_BC_NONE);
    if( !neumannMarkers.size() )
    {
        if( dirichletMarkers.size() ){
            if(verbose)
                std::cout << "  problem - Warning: only Dirichlet boundary conditions have been imposed!" << std::endl;
            hint = BCHandler::HINT_BC_ONLY_ESSENTIAL;
        }
        else
            if(verbose)
                std::cout << "  problem - Warning: NO boundary conditions have been imposed!" << std::endl;
    }

    BCHandler bcH( 0, hint );
    BCFunctionBase uDirichlet( solution_type::u_ex );
    BCFunctionBase uNeumann( solution_type::fNeumann );

    for (std::list<UInt>::const_iterator it = dirichletMarkers.begin();
            it != dirichletMarkers.end(); ++it)
    {
        // std::cout << "\nDirichlet BC on section " << *it << std::flush;
        bcH.addBC( "Dirichlet", *it, Essential, Full, uDirichlet, dataADRPtr->solutionFieldDimension() );
    }
    for (std::list<UInt>::const_iterator it = neumannMarkers.begin();
            it != neumannMarkers.end(); ++it)
    {
        // std::cout << "\nNeumann BC on section " << *it << std::flush;
        bcH.addBC( "Neumann", *it, Natural, Full, uNeumann, dataADRPtr->solutionFieldDimension() );
    }
    std::cout << std::endl;


    // ===================================================
    // Set up the problem
    // ===================================================
    // Classes managing the time advance and integration
    dataADRPtr->setDataTimePtr( dataTimePtr );
    timeIntegratorBDFPtr->setup( dataTimePtr->getBDF_order() );

    // post processing setup
#if POSTPROCESS
    vector_ptr_type computed_solution( new vector_type( uFESpacePtr->map(), Repeated ) );
    vector_ptr_type exact_solution( new vector_type( uFESpacePtr->map(), Repeated ) );

    exporter->addVariable( ExporterData::Scalar, "computedSolution", computed_solution,
                           UInt(0), UInt(uFESpacePtr->dof().numTotalDof()));
    exporter->addVariable( ExporterData::Scalar, "exactSolution", exact_solution,
                           UInt(0), UInt(uFESpacePtr->dof().numTotalDof()));
#endif

    // Initialize the solver internal data structures
    adrSolver.setDataPtr( dataADRPtr );
    adrSolver.setUFESpacePtr( uFESpacePtr );
    adrSolver.setup();

    // Initialize the solution
    adrSolver.initialize( solution_type::u_ex );
    vector_type u0( adrSolver.solution() );
    timeIntegratorBDFPtr->initialize_unk( solution_type::u_ex, u0,
                                          *uFESpacePtr, dataTimePtr->getInitialTime(),
                                          dataTimePtr->getTimeStep() );
    // timeIntegratorBDFPtr->initialize_unk( adrSolver.solution(), dataTimePtr->getTimeStep() );

    timeEval = dataTimePtr->getInitialTime();
    L2_Error = uFESpacePtr->L2Error(uExact, adrSolver.solution(), timeEval, &L2_RelError);
    H1_Error = uFESpacePtr->H1Error(uExact, adrSolver.solution(), timeEval, &H1_RelError);

    if (verbose)
        std::cout << "At time " << timeEval
        << ",\nError Norm L2: " << L2_Error <<
        "\nRelative Error Norm L2: " << L2_RelError<<
        "\nError Norm H1: " << H1_Error <<
        "\nRelative Error Norm H1: " << H1_RelError<<std::endl;


    // Set advection and source terms
    adrSolver.setAdvectionField(solution_type::fAdvection);
    adrSolver.setSourceTerm(solution_type::fSource);

    // ===================================================
    // Set up the problem
    // ===================================================
    chrono.start();

    // assembly of time-constant matrices
    adrSolver.computeConstantMatrices();

    // enter the time loop if the problem is time dependent
    for( dataTimePtr->updateTime(); dataTimePtr->canAdvance(); dataTimePtr->updateTime() )
    {
        if (verbose) std::cout << "\n  problem- We are now at time ... "
                << dataTimePtr->getTime() << std::endl;

        // updating the time dependent coefficients
        dataADRPtr->setDiffusionCoefficient( Members->squareFunction( dataTimePtr->getTime() ) );
        solution_type::diffusionCoeff = dataADRPtr->diffusionCoefficient();

        // updating the system
        if( steadyProblem )
            adrSolver.updateMatrix();
        else
            adrSolver.updateMatrix( timeIntegratorBDFPtr );

        // computing the rhs
        if( steadyProblem )
            adrSolver.updateRHS();
        else
            adrSolver.updateRHS( timeIntegratorBDFPtr );

        // solve the linear system
        adrSolver.iterate(bcH);

        timeIntegratorBDFPtr->shift_right( adrSolver.solution() );

#if POSTPROCESS
        *computed_solution = adrSolver.solution();
        uFESpacePtr->interpolate( solution_type::u_ex, *exact_solution,
                                  dataTimePtr->getTime() );

        exporter->postProcess( dataTimePtr->getTime() );
#endif

        if( steadyProblem ) break;
    }
    chrono.stop();

    if (verbose) std::cout << "\n \n  problem- Total time = " << chrono.diff() << std::endl << std::endl;

    timeEval = dataTimePtr->getTime() - dataTimePtr->getTimeStep();

    L2_Error = uFESpacePtr->L2Error(uExact, *computed_solution, timeEval, &L2_RelError);
    H1_Error = uFESpacePtr->H1Error(uExact, *computed_solution, timeEval, &H1_RelError);

    if (verbose)
        std::cout << "At time " << timeEval
        << ",\nError Norm L2: " << L2_Error <<
        "\nRelative Error Norm L2: " << L2_RelError<<
        "\nError Norm H1: " << H1_Error <<
        "\nRelative Error Norm H1: " << H1_RelError<<std::endl;

}
