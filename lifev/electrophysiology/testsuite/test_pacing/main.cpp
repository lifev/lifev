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

/*!
    @file
    @brief Test for using the pacing protocol

    In this test we use show an example of the use of the pacing protocols.

    @date 03 - 2014
    @author Simone Rossi <simone.rossi@epfl.ch>

    @contributor
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
 */

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

//We use this to transform the mesh from [0,1]x[0,1] to [0,5]x[0,5]
#include <lifev/core/mesh/MeshUtility.hpp>

//We will use the ETAMonodomainSolver
#include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>

//We will use the Aliev Panfilov model for simplicity
#include <lifev/electrophysiology/solver/IonicModels/IonicAlievPanfilov.hpp>

//Include LifeV
#include <lifev/core/LifeV.hpp>

//Include utilities for the pacing protocols
#include <lifev/electrophysiology/stimulus/StimulusPacingProtocol.hpp>

//this is useful to save the out in a separate folder
#include <sys/stat.h>

using namespace LifeV;

//Norm of the solution in order to check if the test failed
#define solutionNorm 38.2973164904655

//Forward declaration of the function defining the fiber direction.
//See the end of the file for the implementation
Real fiberDistribution ( const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& z, const ID&   id);

// Test pacing
Int main ( Int argc, char** argv )
{

    //! Initializing Epetra communicator
    MPI_Init (&argc, &argv);
    std::shared_ptr<Epetra_Comm>  Comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "% using MPI" << std::endl;
    }

    //*********************************************//
    // creating output folder
    //*********************************************//
    GetPot commandLine ( argc, argv );
    std::string problemFolder = commandLine.follow ( "Output", 2, "-o", "--output" );
    // Create the problem folder
    if ( problemFolder.compare ("./") )
    {
        problemFolder += "/";

        if ( Comm->MyPID() == 0 )
        {
            mkdir ( problemFolder.c_str(), 0777 );
        }
    }

    //********************************************//
    // Some typedefs                              //
    //********************************************//

    typedef RegionMesh<LinearTetra>                         mesh_Type;

    typedef std::function < Real (const Real& /*t*/,
                                    const Real &   x,
                                    const Real &   y,
                                    const Real& /*z*/,
                                    const ID&   /*i*/ ) >   function_Type;

    typedef ElectroETAMonodomainSolver< mesh_Type, IonicAlievPanfilov > monodomainSolver_Type;
    typedef std::shared_ptr< monodomainSolver_Type >                 monodomainSolverPtr_Type;

    typedef VectorEpetra                                               vector_Type;
    typedef std::shared_ptr<vector_Type>                             vectorPtr_Type;

    //********************************************//
    // Import parameters from an xml list. Use    //
    // Teuchos to create a list from a given file //
    // in the execution directory.                //
    //********************************************//

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Importing parameters list...";
    }
    Teuchos::ParameterList monodomainList = * ( Teuchos::getParametersFromXmlFile ( "MonodomainSolverParamList.xml" ) );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }

    //********************************************//
    // Creates a new model object representing the//
    // model from Aliev and Panfilov 1996.  The   //
    // model input are the parameters. Pass  the  //
    // parameter list in the constructor          //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Building Constructor for Aliev Panfilov model ... ";
    }
    std::shared_ptr<IonicAlievPanfilov>  model ( new IonicAlievPanfilov() );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }


    //********************************************//
    // In the parameter list we need to specify   //
    // the mesh name and the mesh path.           //
    //********************************************//

    std::string meshName = monodomainList.get ("mesh_name", "lid16.mesh");
    std::string meshPath = monodomainList.get ("mesh_path", "./");

    //********************************************//
    // We need the GetPot datafile for to setup   //
    // the preconditioner.                        //
    //********************************************//

    GetPot dataFile (argc, argv);

    //********************************************//
    // We create the monodomain solver            //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Building Monodomain Solvers... ";
    }

    monodomainSolverPtr_Type monodomain ( new monodomainSolver_Type ( meshName, meshPath, dataFile, model ) );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << " Splitting solver done... ";
    }

    //********************************************//
    // We transform the mesh to get a larger      //
    // domain.                                    //
    //********************************************//
    std::vector<Real> scale (3, 1.0);
    scale[2] = 5.0;
    scale[1] = 5.0;
    std::vector<Real> rotate (3, 0.0);
    std::vector<Real> translate (3, 0.0);
    MeshUtility::MeshTransformer<mesh_Type> transformer (* (monodomain -> localMeshPtr() ) );
    transformer.transformMesh (scale, rotate, translate);


    //********************************************//
    // Initialize the solution                    //
    //********************************************//
    monodomain -> setInitialConditions();

    //********************************************//
    // Setting up the monodomain solver parameters//
    //********************************************//
    monodomain -> setParameters ( monodomainList );

    //********************************************//
    // fiber direction                            //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nSetting fibers:  " ;
    }
    // This is used only to initialize the pointer in the solver
    // could be replaced with a method "initialize fibers"
    // where we just create the pointer to a vectorial vector.
    // The method generates a default constant fiber direction (0,0,1).
    // At least, in case of misuse, the fibers will not be zero
    // and the signal may propagate.
    monodomain -> setupFibers();
    function_Type fibreFunction = &fiberDistribution;
    ElectrophysiologyUtility::setFibersFromFunction (monodomain -> fiberPtr(), monodomain -> localMeshPtr(), fibreFunction);
    ElectrophysiologyUtility::normalize (* (monodomain -> fiberPtr() ) );

    //********************************************//
    // Set some noise in the fiber.               //
    // Set to false by default.                   //
    //********************************************//
    bool randomNoise = monodomainList.get ("noise", false);
    if (randomNoise)
    {
        std::vector<bool> component (3, false);
        component[0] = true;
        Real magnitude = monodomainList.get ("noise_magnitude", 1e-3);
        ElectrophysiologyUtility::addNoiseToFibers (* (monodomain -> fiberPtr() ), magnitude, component );
    }
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Done! \n" ;
    }

    //********************************************//
    // Saving Fiber direction to file             //
    //********************************************//
    monodomain -> exportFiberDirection (problemFolder);

    //********************************************//
    // Create the global matrix: mass + stiffness //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nSetup operators:  " ;
    }
    monodomain -> setupLumpedMassMatrix();
    monodomain -> setupStiffnessMatrix();
    monodomain -> setupGlobalMatrix();
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Done! \n" ;
    }

    //********************************************//
    // Creating exporters to save the solution    //
    //********************************************//
    ExporterHDF5< RegionMesh <LinearTetra> > exporter;
    monodomain -> setupExporter ( exporter, monodomainList.get ("OutputFile", "Solution"), problemFolder );
    monodomain -> exportSolution ( exporter, 0);

    //********************************************//
    // Solving the system                         //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nstart solving:  " ;
    }

    Real dt = monodomain -> timeStep();
    //Uncomment for proper use
    Real TF = monodomainList.get ("endTime", 48.0);
    Int iter = monodomainList.get ("saveStep", 1.0) / dt;


    //********************************************//
    // Define the pacing protocol                 //
    //********************************************//
    StimulusPacingProtocol pacing;

    //********************************************//
    // Import parameters for the pacing protocol  //
    //********************************************//

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Importing pacing protocol parameters ...";
    }
    std::string pacingListFileName = monodomainList.get ("pacing_parameter_list_filename", "PacingParameters.xml");
    Teuchos::ParameterList pacingList = * ( Teuchos::getParametersFromXmlFile ( pacingListFileName ) );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }

    pacing.setParameters (pacingList);
    pacing.setTimeStep ( monodomain -> timeStep() );

    //********************************************//
    // Loop over time solving with L-ICI          //
    //********************************************//
    int loop = 0;
    for (Real t = monodomain -> initialTime(); t < TF;)
    {
        //********************************************//
        // Set the applied current from the pacing    //
        // protocol                                   //
        //********************************************//
        monodomain -> setAppliedCurrentFromElectroStimulus ( pacing, t);
        loop++;
        t += dt;

        //********************************************//
        // Solve the monodomain equations.            //
        //********************************************//

        monodomain -> solveOneStepGatingVariablesFE();
        monodomain -> solveOneICIStep();

        //********************************************//
        // Exporting the solution.                    //
        //********************************************//
        if (loop % iter == 0 )
        {
            if ( Comm->MyPID() == 0 )
            {
                std::cout << "\ntime = " << t;
            }
            exporter.postProcess (t);
        }
    }

    //********************************************//
    // Close the exporter                         //
    //********************************************//

    exporter.closeFile();


    //********************************************//
    // Check if the test failed                   //
    //********************************************//
    Real newSolutionNorm = monodomain -> potentialPtr() -> norm2();

    monodomain.reset();
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Finalize();


    Real err = std::abs (newSolutionNorm - solutionNorm) / std::abs (solutionNorm);
    if ( err > 1e-8 )
    {
        std::cout << "\nTest Failed: " <<  err << "\n" << "\nSolution Norm: " <<  newSolutionNorm << "\n";
        return EXIT_FAILURE; // Norm of solution did not match
    }
    else
    {
        return EXIT_SUCCESS;
    }
}

#undef solutionNorm

//Definition of the fiber direction
Real fiberDistribution ( const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& z, const ID&   id)
{
    Real y0 = 2.5;
    Real z0 = 2.5;
    Real r = std::sqrt ( (y - y0) * (y - y0) + (z - z0) * (z - z0) );

    switch ( id )
    {
        case 0:
            return 0.0;
        case 1:
            if (r > 1e-16)
            {
                return (z - z0) / r;
            }
            else
            {
                return 0.0;
            }
        case 2:
            if (r > 1e-16)
            {
                return (y0 - y) / r;
            }
            else
            {
                return 1.0;
            }
        default:
            return 0.0;
    }
}
