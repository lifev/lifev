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
    @brief Test to restart a simulation from an hdf5 exported solution

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

//this is useful to save the out in a separate folder
#include <sys/stat.h>

using namespace LifeV;

//Forward declaration of the function defining the fiber direction.
//See the end of the file for the implementation
Real cut (const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& /*z*/, const ID& /*i*/);
Real initialCondition ( const Real& /*t*/, const Real& x, const Real& /*y*/, const Real& /*z*/, const ID&   /*id*/);

//Forward declaration of the function defining the fiber direction.
//See the end of the file for the implementation
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

    typedef RegionMesh<LinearTetra>                         mesh_Type;

    typedef boost::function < Real (const Real& /*t*/,
                                    const Real &   x,
                                    const Real &   y,
                                    const Real& /*z*/,
                                    const ID&   /*i*/ ) >   function_Type;

    typedef ElectroETAMonodomainSolver< mesh_Type, IonicAlievPanfilov > monodomainSolver_Type;
    typedef boost::shared_ptr< monodomainSolver_Type >                 monodomainSolverPtr_Type;

    typedef VectorEpetra                                               vector_Type;
    typedef boost::shared_ptr<vector_Type>                             vectorPtr_Type;

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
    boost::shared_ptr<IonicAlievPanfilov>  model ( new IonicAlievPanfilov() );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
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
    // Computing the norm of the precomputed at   //
    // final time                                 //
    //********************************************//
    std::string prefix = monodomainList.get ("importPrefix", "ciccia");
    std::string dir = monodomainList.get ("importDir", "./");
    monodomain -> importSolution (prefix, dir, 100.0);
    Real solutionNorm = monodomain -> potentialPtr() -> norm2();

    //********************************************//
    // Setting up the initial condition form      //
    // importing the old solution                 //
    //********************************************//

    if ( Comm->MyPID() == 0 )
    {
        cout << "\nInitializing potential:  " ;
    }

    Real initialTime = monodomainList.get ("importTime", 0.0);
    monodomain -> importSolution (prefix, dir, initialTime);

    if ( Comm->MyPID() == 0 )
    {
        cout << "Done! \n" ;
    }

    //********************************************//
    // Setting up the monodomain solver           //
    //********************************************//

    monodomain -> setParameters ( monodomainList );

    //********************************************//
    // fiber direction                            //
    //********************************************//

    if ( Comm->MyPID() == 0 )
    {
        cout << "\nImporting fibers:  " ;
    }

    VectorSmall<3> fibers;
    fibers[0] = monodomainList.get ("fibers_X", 0.0);
    fibers[1] = monodomainList.get ("fibers_Y", 0.0);
    fibers[2] = monodomainList.get ("fibers_Z", 1.0);

    monodomain -> setupFibers (fibers);
    if ( Comm->MyPID() == 0 )
    {
        cout << "Done! \n" ;
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
    monodomain -> exportSolution ( exporter, initialTime);

    //********************************************//
    // Solving the system                         //
    //********************************************//

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nstart solving:  " ;
    }

    Real dt = monodomain -> timeStep();
    Real cutTime = monodomainList.get ("cutTime", 150.0);
    //Uncomment for proper use
    Real TF =  monodomain -> endTime();
    //Real TF = 100.0;
    Int iter = monodomainList.get ("saveStep", 1.0) / dt;


    //********************************************//
    // Defining the cut for generating the spiral //
    //********************************************//
    vectorPtr_Type spiral ( new vector_Type ( ( monodomain -> globalSolution().at (0) ) -> map() ) );
    function_Type f = &cut;
    monodomain -> feSpacePtr() -> interpolate (
        static_cast<FESpace<RegionMesh<LinearTetra>, MapEpetra>::function_Type> (f),
        *spiral, 0.0);


    //********************************************//
    // Loop over time solving with L-ICI          //
    //********************************************//
    int loop = 0;
    for (Real t = initialTime; t < (TF - dt * 1e-4);) // the -dt*1e-4 is needed or you do an additional iteration
    {
        loop++;
        t += dt;

        if (t >= cutTime && t <= cutTime + dt)
        {
            * (monodomain -> potentialPtr() ) *= *spiral;
        }
        monodomain -> solveOneStepGatingVariablesFE();
        monodomain -> solveOneICIStep();
        if (loop % iter == 0 )
        {
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
    std::cout << std::setprecision (20) << "\nError: " <<  err << "\nSolution Norm: " <<  newSolutionNorm << "\n";
    std::cout << std::setprecision (20) << "\nImported solution Norm: " <<  solutionNorm << "\n";
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



// Function to cut the spiral
Real cut (const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& /*z*/, const ID& /*i*/)
{
    if ( y >= 2.5)
    {
        return 1.0;
    }
    else
    {
        return 0.0;
    }
}

//Initial condition used for the precomputed solution
Real initialCondition ( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& z, const ID&   /*id*/)
{
    if ( z <= 0.4 )
    {
        return 1.0;
    }
    else
    {
        return 0;
    }
}


