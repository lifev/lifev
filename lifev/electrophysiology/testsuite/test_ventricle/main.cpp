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
    @brief Test on an idealized ventricle

    In this test I solve the monodomain on the idealized ventricle
    setting with the minimal model and SVI.
    To initiate the excitation I set the initial condition of
    the potential to one on the endocardium.
    The fiber field is loaded from file.

    Note that the given mesh is super coarse.
    Using SVI the propagation will be super fast.
    Please refine the mesh and run this test again.
    You can create a fiber field using the test_fibers

    @date 04 - 2014
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

//We will use the ETAMonodomainSolver
#include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>

//We will use the Aliev Panfilov model for simplicity
#include <lifev/electrophysiology/solver/IonicModels/IonicMinimalModel.hpp>

//Include LifeV
#include <lifev/core/LifeV.hpp>

//Include utilities for importing the fiber field
#include <lifev/electrophysiology/util/HeartUtility.hpp>

//this is useful to save the out in a separate folder
#include <sys/stat.h>

using namespace LifeV;

//Norm of the solution in order to check if the test failed
#define solutionNorm 48.66815714672445381

// Test pacing
Int main ( Int argc, char** argv )
{

    //! Initializing Epetra communicator
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm>  Comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    if ( Comm->MyPID() == 0 )
    {
        cout << "% using MPI" << endl;
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

    typedef RegionMesh<LinearTetra>                         			mesh_Type;

    typedef IonicMinimalModel											ionicModel_Type;

    typedef ElectroETAMonodomainSolver< mesh_Type, ionicModel_Type > 	monodomainSolver_Type;

    typedef boost::shared_ptr< monodomainSolver_Type >                  monodomainSolverPtr_Type;

    typedef VectorEpetra                                               	vector_Type;

    typedef boost::shared_ptr<vector_Type>                             	vectorPtr_Type;

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
        std::cout << " Done!" << endl;
    }

    //********************************************//
    // Creates a new model object representing the//
    // model from Aliev and Panfilov 1996.  The   //
    // model input are the parameters. Pass  the  //
    // parameter list in the constructor          //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Building Constructor for Minimal Model with parameters ... ";
    }
    boost::shared_ptr<ionicModel_Type>  model ( new ionicModel_Type() );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }


    //********************************************//
    // In the parameter list we need to specify   //
    // the mesh name and the mesh path.           //
    //********************************************//

    std::string meshName = monodomainList.get ("mesh_name", "idealized.mesh");
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
    // Initialize the solution                    //
    //********************************************//
    monodomain -> setInitialConditions();
	ElectrophysiologyUtility::setValueOnBoundary( *(monodomain -> potentialPtr() ), monodomain -> fullMeshPtr(), 1.0, 5 );


    //********************************************//
    // Setting up the monodomain solver parameters//
    //********************************************//
    monodomain -> setParameters ( monodomainList );

    //********************************************//
    // fiber direction                            //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        cout << "\nSetting fibers:  " ;
    }
    // This is used only to initialize the pointer in the solver
    // could be replaced with a method "initialize fibers"
    // where we just create the pointer to a vectorial vector.
    // The method generates a default constant fiber direction (0,0,1).
    // At least, in case of misuse, the fibers will not be zero
    // and the signal may propagate.
    monodomain -> setupFibers();
    std::string fiberFile(monodomainList.get ("solid_fiber_file", ""));
    std::string fiberField(monodomainList.get ("solid_fiber_field", ""));
    ElectrophysiologyUtility::importVectorField(monodomain -> fiberPtr(),  fiberFile,  fiberField, monodomain -> localMeshPtr() );

    if ( Comm->MyPID() == 0 )
    {
        cout << "Done! \n" ;
    }

    //********************************************//
    // Saving Fiber direction to file             //
    //********************************************//
    monodomain -> exportFiberDirection(problemFolder);

    //********************************************//
    // Create the global matrix: mass + stiffness //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        cout << "\nSetup operators:  " ;
    }
    monodomain -> setupLumpedMassMatrix();
    monodomain -> setupStiffnessMatrix();
    monodomain -> setupGlobalMatrix();
    if ( Comm->MyPID() == 0 )
    {
        cout << "Done! \n" ;
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
    // Loop over time solving with L-ICI          //
    //********************************************//
	int loop = 0;
	for (Real t = monodomain -> initialTime(); t < TF;)
    {
	    //********************************************//
	    // Set the applied current from the pacing    //
		// protocol                                   //
	    //********************************************//
		loop++;
        t += dt;

	    //********************************************//
	    // Solve the monodomain equations.            //
	    //********************************************//

        monodomain -> solveOneStepGatingVariablesFE();
        monodomain -> solveOneSVIStep();

        //********************************************//
		// Exporting the solution.                    //
		//********************************************//
        if(loop % iter == 0 )
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


    Real err = std::abs (newSolutionNorm - solutionNorm) / std::abs(solutionNorm);
    if ( err > 1e-8 )
    {
    	std::cout << std::setprecision(20) << "\nTest Failed: " <<  err <<"\n" << "\nSolution Norm: " <<  newSolutionNorm << "\n";
        return EXIT_FAILURE; // Norm of solution did not match
    }
    else
    {
        return EXIT_SUCCESS;
    }
}

#undef solutionNorm


