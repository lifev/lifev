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
    @brief Electrophysiology benchmark proposed in Niederer et al. 2011

    This test is the basic test to start using the electrophysiology module.
    We consider a parallelepiped excited in a small region in one corner.
    The excitation wave propagates and excites the whole domain.
    We check the time it takes for the wave to activate the whole domain.
    We provide two coarse meshes with this test.
    Be aware that with such mesh sizes we are far from convergence.

    The parameter list is an xml file. To work properly it should be called
    MonodomainSolverParamList.xml
    as the provided one, or change the name of the file (see below in the code.)

    Run it using e.g.:
    mpirun -n 2 Electrophysiology_benchmark -o OutputFolder


    WARNING: The provided mesh is coarse!!! You should not expect nice results
    using such mesh. A finer mesh can be obtained just using gmsh and using
    the refine by splitting function. The required mesh would be too large
    and we decided to provide only a coarse one.

    NOTICE: The value of the membrane capacitance has a big influence on the
    propagation. For example the TenTusscher model has Cm = 2 and final activation
    time (in this benchmark) about twice the one computed with Luo-Rudy I or the Minimal
    Model. Imposing Cm = 1, then, these ionic models give similar results in terms of
    activation times.
    Be aware of what you are doing! Don't use the default parameters as black boxes!!!

    @date 01-2013
    @author Simone Rossi <simone.rossi@epfl.ch>

    @contributor
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
 */

// ---------------------------------------------------------------
// In order to solve the monodomain model we need
// to include the ETAMonodmomainSolver
// ---------------------------------------------------------------

#include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>

// ---------------------------------------------------------------
//  We created a separate file where we collect all the functions
// needed to run the benchmark. In particular this utility file
// contains the following functions:
// - chooseIonicModel,
//          to be chosen among:
//           - AlievPanfilov
//           - LuoRudyI
//           - TenTusscher06
//           - HodgkinHuxley
//           - NoblePurkinje
//           - MinimalModel
//           - Fox (tested with timestep 0.0025 ms, RushLarsen method not implemented)
//
// - pacingProtocolMM,
//           to pace with the minimal model and AlievPanfilov
// - pacingProtocolHH,
//           to pace with the Hodgkin-Huxley model
// - pacingProtocol,
//           to pace with the other models
// - setStimulus,
//           to actually set the above pacings depending on the chosen ionic model
// ---------------------------------------------------------------

#include <lifev/electrophysiology/testsuite/test_benchmark/benchmarkUtility.hpp>

// ---------------------------------------------------------------
//  This include is necessary in order to save the output of the
// simulation on a specified folder. To specify the output folder
// you should call "-o Folder", e.g.:
//
// mpirun -n 2 Electrophysiology_benchmark -o OutputFolder
// ---------------------------------------------------------------

#include <sys/stat.h>

// ---------------------------------------------------------------
// As usual, we work in the LifeV namespace.
// ---------------------------------------------------------------

using namespace LifeV;

// ---------------------------------------------------------------
//  The test is checked only with the AlievPanfilov and it compares
// the final activation time with a precomputed one.
// This final time was computed using gcc on Ubuntu/Linux
// On Mac Os X 10.9 Mavercik with clang 503.0.38 the test has final time 40.
// ---------------------------------------------------------------

#define finalActivationTime  53.1

// ---------------------------------------------------------------
//  We start the programm by the definition of the communicator.
// We assume you will be using mpi even for a single process.
// ---------------------------------------------------------------

Int main ( Int argc, char** argv )
{

    //! Initializing Epetra communicator
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm>  Comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "% using MPI" << std::endl;
    }

    // ---------------------------------------------------------------
    //  We create the output folder where we save the solution.
    // The folder can be specified appending "-o FolderName" when
    // in the command  line when executing the test. If nothing is
    // specified we create by default a folder called "Ouput".
    // ---------------------------------------------------------------

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

    // ---------------------------------------------------------------
    //  For convenience we create some typedefs for the mesh,
    // the vectors, the matrices, the functions (used to impose the
    // external stimulus), the ionic model and the monodomain solver.
    // ---------------------------------------------------------------

    typedef RegionMesh<LinearTetra>                                     mesh_Type;

    typedef VectorEpetra                                                vector_Type;
    typedef boost::shared_ptr<vector_Type>                              vectorPtr_Type;

    typedef MatrixEpetra<Real>                                          matrix_Type;
    typedef boost::shared_ptr<matrix_Type>                              matrixPtr_Type;

    typedef boost::function < Real (const Real& /*t*/,
                                    const Real &   x,
                                    const Real &   y,
                                    const Real& /*z*/,
                                    const ID&   /*i*/ ) >               function_Type;

    typedef ElectroIonicModel                                           ionicModel_Type;
    typedef boost::shared_ptr<ionicModel_Type>                          ionicModelPtr_Type;

    typedef ElectroETAMonodomainSolver< mesh_Type, ionicModel_Type >    monodomainSolver_Type;
    typedef boost::shared_ptr< monodomainSolver_Type >                  monodomainSolverPtr_Type;

    // ---------------------------------------------------------------
    //  We set up a chronometer, to check how long it takes to setup
    // the monodomain problem with different methods. The chronometer
    // will be check the time only on the processor 0.
    // ---------------------------------------------------------------

    LifeChrono chronoinitialsettings;

    if ( Comm->MyPID() == 0 )
    {
        chronoinitialsettings.start();
    }

    // ---------------------------------------------------------------
    //  We read the xml parameter list file using Teuchos.
    // ---------------------------------------------------------------

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Importing parameters list...";
    }
    Teuchos::ParameterList monodomainList = * ( Teuchos::getParametersFromXmlFile ( "MonodomainSolverParamList.xml" ) );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }

    // ---------------------------------------------------------------
    //  From the parameter list we read the ionic model. We use the
    // function defined in the benchmarkUtility.hpp file to choose
    // between the models.
    //           - AlievPanfilov
    //           - LuoRudyI
    //           - TenTusscher06
    //           - HodgkinHuxley
    //           - NoblePurkinje
    //           - MinimalModel
    //           - Fox (tested with timestep 0.0025 ms)
    // These are not the only model available. All the ionic models can
    // be found in the folder solver/IonicModels/ of this module.
    // The function chooseIonicModel returns the value for which we
    // consider tissue activation. For example in adimensional models,
    // such as the Aliev-Panfilov or the MinimalModel, this value
    // is set to 0.95. The same values is also used for the other
    // cardiac models. For the Hodgkin-Huxley model instead we set
    // this threshold value equal to ten.
    // ---------------------------------------------------------------

    std::string ionic_model ( monodomainList.get ("ionic_model", "minimalModel") );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nIonic_Model:" << ionic_model;
    }
    ionicModelPtr_Type  model;
    Real activationThreshold = BenchmarkUtility::chooseIonicModel (model, ionic_model, *Comm );


    // ---------------------------------------------------------------
    //  We read from the parameter list the name of the mesh file
    // and the path to the mesh. Please use absolute paths.
    // Make sure you specify the mesh name in the xml file, otherwise
    // the test will fail.
    // ---------------------------------------------------------------

    std::string meshName = monodomainList.get ("mesh_name", "");
    std::string meshPath = monodomainList.get ("mesh_path", "./");

    // ---------------------------------------------------------------
    //  We create a GetPot datafile. Although the module do not require
    // the use of datafile, the object is required in order to
    // setup the preconditioners.
    // ---------------------------------------------------------------

    GetPot dataFile  (argc, argv);

    // ---------------------------------------------------------------
    // The standard way to create the monodomain solver is to call the
    // constructur with the following arguments:
    // - std::string meshName (the name of the mesh file)
    // - std::string meshPath (the path to the mesh)
    // - GetPot dataFile      (as explained above, used only to setup the preconditioners)
    // - ionicModelPtr_Type   (a shared_ptr to the ionic model)
    // The monodomain solver class is templetized over the IonicModel class
    // as it does not make sense to run a monodomain simulation without
    // an ionic model.
    // For the bistable equation consider using the ADRassembler available in LifeV.
    // ---------------------------------------------------------------

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nBuilding Monodomain Solver ... \n";
    }

    monodomainSolverPtr_Type solver ( new monodomainSolver_Type ( meshName, meshPath, dataFile , model ) );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\t...  solver created!";
    }

    // ---------------------------------------------------------------
    // We initialize the potential and the applied current to zero.
    // Each ionic model has its own resting conditions.
    // The method setInitialConditions, initialize all the variables
    // in the ionic model with this default values.
    // If different initial conditions are needed it is possible
    // to use setters and getters.
    // ---------------------------------------------------------------

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nInitializing potential and gating variables ...  " ;
    }
    solver -> setInitialConditions();
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Done!" ;
    }

    // ---------------------------------------------------------------
    // In order to set the parameters of the modomain solver
    // we call the method setParameters which take as argument
    // the Teuchos parameter list we have imported at the beginning
    // This method sets:
    // - the time step                  (default: 0.01 ms)
    // - the initial and ending time    (default: 0 and 100 ms)
    // - the surface to volume ratio    (default: 2400.0)
    // - the conductivity coefficients  (default: 0.001 isotropic)
    // - the order of the elements      (default: P1. Not tested with higher order elements)
    // - the lumping of the mass matrix (default: false)
    // Take care that the default values do not make sense in particular
    // if the surface to volume ratio is of the order of 1e3 then the conductivity should
    // be around 1e-1 - 1e0.
    // ---------------------------------------------------------------

    solver -> setParameters ( monodomainList );

    // ---------------------------------------------------------------
    // Cardiac tissue is typically model as transversely isotropic.
    // We need to define the preferred direction (or fiber direction).
    // In this simple case, we specify a uniform fiber field described
    // by the vector (0, 0, 1). The setupFiber method setup
    // an EpetraVector  with the given direction. In general
    // the fiber field must be computed using the rule-based algorithm
    // available in the module.
    // ---------------------------------------------------------------

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nSetting fibers ...  " ;
    }

    VectorSmall<3> fibers;
    fibers[0] = 0.0;
    fibers[1] = 0.0;
    fibers[2] = 1.0;
    solver -> setupFibers ( fibers );

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Done!" ;
    }

    // ---------------------------------------------------------------
    // Let's choose the method to solve the monodomain
    // model:
    // - 1st order operator splitting (2nd order available but still experimental)
    // - L-ICI method (or Full lumping, all mass matrices are lumped - behavior equivalent to operator splitting)
    // - ICI method (or Half-Lumping, only the mass matrix relative to the time dependent term is lumped)
    // - SVI method (the most computationally expensive method)
    // ---------------------------------------------------------------

    std::string solutionMethod = monodomainList.get ("solutionMethod", "splitting");

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nSolving the monodomain using " << solutionMethod;
    }

    // ---------------------------------------------------------------
    // The monodomain is typically solved using an IMEX method, where
    // the linear part is treated implicitely and the nonlinear
    // explicitely. Therefore the operator is constant and has
    // two contributions: a mass matrix and a stiffness matrix.
    // Here we builde the two matrices and we add them together in
    // a global matrix.
    //
    // If we choose to lump the mass matrix (recommended), the lumping
    // is achieved by using nodal integration.
    // If we want to solve the system using the so called ICI method
    // ( sometimes called half-lumping method ) we need both the full
    // mass matrix and the lumped one.
    // So if we want to lump the mass matrix we first create an auxiliary
    // full mass matrix that we will use to solve with the ICI method.
    // ---------------------------------------------------------------

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nSetup operators:  \n" ;
    }

    // Read if you are going to use a lumped mass matrix;
    // This is set by the xml file when you create the monodomain
    bool lumpedMass = solver -> lumpedMassMatrix();

    //safe check! You should not do this error.
    // If you are using L-ICI, you should set setLumperMassMatrix to false
    // If you do not set anything in the xml is set to false by default
    if (solutionMethod == "L-ICI" && !lumpedMass)
    {
        std::cout << "===============================================\n";
        std::cout << "You are using L-ICI without lumping! You can't!\n";
        std::cout << "This time I'm fixing it for you ... be careful.\n";
        std::cout << "===============================================\n";
        solver -> setLumpedMassMatrix (true);
    }

    //We create a pointer to store a full mass matrix
    matrixPtr_Type fullMass;

    //if we are using ICI then we need to compute the fullMass matrix even
    // if we are using lumping
    if ( lumpedMass && solutionMethod == "ICI")
    {
        solver -> setLumpedMassMatrix (false);
        solver -> setupMassMatrix();
        fullMass.reset (new matrix_Type ( * (solver -> massMatrixPtr() ) ) );
        solver -> setFullMassMatrixPtr (fullMass);
        solver -> setLumpedMassMatrix (lumpedMass);
    }

    //Build the solver mass matrix
    solver -> setupMassMatrix();

    //safety check!!! In this tesst ICI is solved using the solveOneICIStepWithFullMass()
    //Therefore we need to set the fullMass pointer. If we are using lumping that we
    // already set it before, but if you are not lumping (not recommended choice)
    // then the full mass matrix is equal to the solver mass matrix.
    if ( !lumpedMass && solutionMethod == "ICI")
    {
        solver -> setFullMassMatrixPtr (solver -> massMatrixPtr() );
    }

    //Building the stiffness matrix and the global matrix (stifness and mass)
    solver -> setupStiffnessMatrix ();
    solver -> setupGlobalMatrix();

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Done!" ;
    }

    // ---------------------------------------------------------------
    // The exporter is used to save the solution of the simulation.
    // It's an external object and not part of the solver. Therefore
    // we need to create it ourself. On the other hand the monodomain
    // solver can set it up in order to save all the solution solved
    // by the solver, which may vary depending on the ionic model.
    // We pass the exporter, the name of the file where we want to
    // save the solution, and the folder where we want it to be saved.
    // We export the initial conditions at time 0.
    // We always use HDF5 ouput.
    // ---------------------------------------------------------------

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nSetting up the exporter ... " ;
    }
    ExporterHDF5< RegionMesh <LinearTetra> > exporter;
    solver -> setupExporter ( exporter, monodomainList.get ("OutputFile", "Solution") , problemFolder);
    if ( Comm->MyPID() == 0 )
    {
        std::cout << " exporting initial solution ... " ;
    }
    solver -> exportSolution ( exporter, 0);

    // ---------------------------------------------------------------
    // We want to save the activation times in the domains.
    // Therefore, we create a vector which is initialized with the value -1.
    // At every timestep, we will check if the nodes in the mesh have
    // been activated, that is we check if the value of the potential
    // is bigger than a given threshold (which was defined at the beninning
    // when choosing the ionic model).
    // Moreover, we want to export the activation time. We therefore create
    // another HDF5 exporter to save the activation times on a separate file.
    // ---------------------------------------------------------------

    vectorPtr_Type activationTimeVector ( new vector_Type ( solver -> potentialPtr() -> map() ) );
    *activationTimeVector = -1.0;

    ExporterHDF5< RegionMesh <LinearTetra> > activationTimeExporter;
    activationTimeExporter.setMeshProcId (solver -> localMeshPtr(), solver -> commPtr() ->MyPID() );
    activationTimeExporter.addVariable (ExporterData<mesh_Type>::ScalarField, "Activation Time",
                                        solver -> feSpacePtr(), activationTimeVector, UInt (0) );
    activationTimeExporter.setPrefix ("ActivationTime");
    activationTimeExporter.setPostDir (problemFolder);

    // ---------------------------------------------------------------
    // We are ready to solve the monodomain model. We will not save the
    // solution at every timestep. We put in the xml file the timestep
    // between each save. By default, we save evry 1ms.
    // We will keep track of the time used for solving the system and therefore
    // we initialize some variables to record these values.
    // ---------------------------------------------------------------

    Real dt (solver -> timeStep() );
    Int iter = monodomainList.get ("saveStep", 1.0) / dt;
    Int subiter = monodomainList.get ("subiter", 10);
    Int k (0);

    Real timeReac = 0.0;
    Real timeDiff = 0.0;
    Real timeReacDiff = 0.0;
    LifeChrono chrono;

    // ---------------------------------------------------------------
    // The external stimulus is given as a boost function.
    // We use the function in the benchamarkUtility.hpp to set up
    // the boost function depending on the chosen ionic model.
    // Then, we call the setAppliedCurrentFromFunction method of the
    // monodomain solver wich sets the external stimulus in the solver.
    // The value 0.0 is the time.
    // ---------------------------------------------------------------

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nSetting up the external stimulus ...  " ;
    }
    function_Type stimulus;
    BenchmarkUtility::setStimulus (stimulus, ionic_model);
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Done!" ;
    }


    // ---------------------------------------------------------------
    // We are done with the setup of the solver. We show how long it
    // took.
    // ---------------------------------------------------------------

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nMonodomain solver setup done in " << chronoinitialsettings.diff() << " s.";
    }

    // ---------------------------------------------------------------
    //  We perform the for loop over time and we start solving the
    // monodomain model.
    // ---------------------------------------------------------------

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nstart solving:  " ;
    }

    for ( Real t = solver -> initialTime(); t < solver -> endTime(); )
    {

        // ---------------------------------------------------------------
        // The first step is to update the applied current with respect
        // to time. We call again the setAppliedCurrentFromFunction
        // method passing the time t.
        // ---------------------------------------------------------------

        solver -> setAppliedCurrentFromFunction ( stimulus, t );
        // ---------------------------------------------------------------
        // Next we consider case by case the different solution methods:
        // ---------------------------------------------------------------

        // ---------------------------------------------------------------
        // Solving with first order operator splitting:
        // ---------------------------------------------------------------

        if ( solutionMethod == "splitting" )
        {

            // ---------------------------------------------------------------
            // We first solve the system of ODEs in the ionic model.
            // For some models it is possible to solve them with the
            // Rush-Larsen method, otherwise we use simple forward Euler.
            // The Rush-Larsen method is directly implemented in the ionic models,
            // while the forward Euler scheme is implemented in the MonodomainSolver.
            // We also save the time needed to compute the reaction step.
            // ---------------------------------------------------------------

            chrono.reset();
            chrono.start();
            if (ionic_model != "MinimalModel" && ionic_model != "AlievPanfilov" && ionic_model != "Fox")
            {
                solver->solveOneReactionStepRL();
            }
            else
            {
                for (int j = 0; j < subiter; j++)
                {
                    solver->solveOneReactionStepFE (subiter);
                }
            }
            chrono.stop();

            timeReac += chrono.globalDiff ( *Comm );

            // ---------------------------------------------------------------
            // We solve the diffusion step. First we need to update the rhs,
            // and then we solve the linear system. Eventually we save the time
            // needed to compute the diffusion step.
            // ---------------------------------------------------------------

            (*solver->rhsPtrUnique() ) *= 0.0;
            solver->updateRhs();

            chrono.reset();
            chrono.start();
            solver->solveOneDiffusionStepBE();
            chrono.stop();
            timeDiff += chrono.globalDiff ( *Comm );
        }

        // ---------------------------------------------------------------
        // Solving with Lumped Ionic Current Interpolation (L-ICI)
        // ---------------------------------------------------------------

        else if ( solutionMethod == "L-ICI" )
        {

            // ---------------------------------------------------------------
            // We first solve the system of ODEs in the ionic model.
            // For some models it is possible to solve them with the
            // Rush-Larsen method, otherwise we use simple forward Euler.
            // The Rush-Larsen method is directly implemented in the ionic models,
            // while the forward Euler scheme is implemented in the MonodomainSolver.
            // We also save the time needed to compute the solution.
            // ---------------------------------------------------------------

            chrono.reset();
            chrono.start();
            if (ionic_model != "MinimalModel" && ionic_model != "AlievPanfilov" && ionic_model != "Fox")
            {
                solver -> solveOneStepGatingVariablesRL();
            }
            else
            {
                solver -> solveOneStepGatingVariablesFE();
            }


            // ---------------------------------------------------------------
            // The solution with L-ICI consist in calling directly the solve ICI
            // method using a lumped mass matrix.
            // ---------------------------------------------------------------

            solver -> solveOneICIStep();
            chrono.stop();
            timeReacDiff += chrono.globalDiff ( *Comm );
        }
        // ---------------------------------------------------------------
        // Solving with Lumped Ionic Current Interpolation (L-ICI)
        // ---------------------------------------------------------------

        else if ( solutionMethod == "ICI" )
        {
            // ---------------------------------------------------------------
            // We first solve the system of ODEs in the ionic model.
            // For some models it is possible to solve them with the
            // Rush-Larsen method, otherwise we use simple forward Euler.
            // The Rush-Larsen method is directly implemented in the ionic models,
            // while the forward Euler scheme is implemented in the MonodomainSolver.
            // We also save the time needed to compute the solution.
            // ---------------------------------------------------------------

            chrono.reset();
            chrono.start();
            if (ionic_model != "MinimalModel" && ionic_model != "AlievPanfilov" && ionic_model != "Fox")
            {
                solver -> solveOneStepGatingVariablesRL();
            }
            else
            {
                solver -> solveOneStepGatingVariablesFE();
            }

            // ---------------------------------------------------------------
            // The solution with L-ICI consist in calling directly the solve ICI
            // method using a lumped mass matrix for the time dependent terms
            // while using the full mass matrix for the reaction part.
            // Therefore if we have lumped the mass matrix, we pass the full mass
            // matrix as argument  in the solveOneICIStep method.
            // ---------------------------------------------------------------

            solver -> solveOneICIStepWithFullMass();
            chrono.stop();
            timeReacDiff += chrono.globalDiff ( *Comm );
        }

        // ---------------------------------------------------------------
        // Solving with State Variable Interpolation (SVI)
        // ---------------------------------------------------------------

        else if ( solutionMethod == "SVI" )
        {

            // ---------------------------------------------------------------
            // We first solve the system of ODEs in the ionic model.
            // For some models it is possible to solve them with the
            // Rush-Larsen method, otherwise we use simple forward Euler.
            // The Rush-Larsen method is directly implemented in the ionic models,
            // while the forward Euler scheme is implemented in the MonodomainSolver.
            // We also save the time needed to compute the solution.
            // ---------------------------------------------------------------

            chrono.reset();
            chrono.start();
            if (ionic_model != "MinimalModel" && ionic_model != "AlievPanfilov" && ionic_model != "Fox")
            {
                solver -> solveOneStepGatingVariablesRL();
            }
            else
            {
                solver -> solveOneStepGatingVariablesFE();
            }


            // ---------------------------------------------------------------
            // We call the SolveOneSVIStep method of the solver.
            // ---------------------------------------------------------------

            solver -> solveOneSVIStep();
            chrono.stop();
            timeReacDiff += chrono.globalDiff ( *Comm );
        }

        // ---------------------------------------------------------------
        // We update the iteration number k, and the time.
        // ---------------------------------------------------------------

        k++;
        t = t + dt;

        // ---------------------------------------------------------------
        // We  save the activation time in the vector  (*activationTimeVector)
        // ---------------------------------------------------------------

        solver -> registerActivationTime (*activationTimeVector, t, activationThreshold);

        // ---------------------------------------------------------------
        // If it's time to save the solution we export using the exportSolution method
        // ---------------------------------------------------------------

        if ( k % iter == 0 )
        {
            if ( Comm->MyPID() == 0 )
            {
                std::cout << "\nTime : " << t;
            }
            solver -> exportSolution (exporter, t);
        }

    }

    // ---------------------------------------------------------------
    // We close the solution exporter. Then we export the activation
    // times and we close the relative exporter.
    // ---------------------------------------------------------------

    exporter.closeFile();
    activationTimeExporter.postProcess (0);
    activationTimeExporter.closeFile();

    // ---------------------------------------------------------------
    // We also export the fiber direction.
    // ---------------------------------------------------------------

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nExporting fibers ...  ";
    }
    solver -> exportFiberDirection (problemFolder);

    // ---------------------------------------------------------------
    // We destroy the solver
    // ---------------------------------------------------------------

    solver.reset();

    // ---------------------------------------------------------------
    // We show show long it took to solve the problem and we thank
    // for using the monodomain solver
    // ---------------------------------------------------------------

    if ( Comm->MyPID() == 0 )
    {
        chronoinitialsettings.stop();
        std::cout << "\n\n\nTotal lapsed time : " << chronoinitialsettings.diff() << std::endl;
        if ( solutionMethod == "splitting" )
        {
            std::cout << "Diffusion time : " << timeDiff << std::endl;
            std::cout << "Reaction time : " << timeReac << std::endl;
        }
        else
        {
            std::cout << "Solution time : " << timeReacDiff << std::endl;
        }

        std::cout << "\n\nThank you for using ETA_MonodomainSolver.\nI hope to meet you again soon!\n All the best for your simulation :P\n  " ;
    }

    // ---------------------------------------------------------------
    // Before ending we test if the test has succeeded.
    // We compute the last activation
    // time and we compare it with the precomputed value.
    // ---------------------------------------------------------------
    Real fullActivationTime = activationTimeVector -> maxValue();
    activationTimeVector.reset();
    Real returnValue;

    Real err = std::abs (fullActivationTime - finalActivationTime) / std::abs (finalActivationTime);

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nError: " <<  err << "\n" << "\nActivation time: " <<  fullActivationTime << "\n";
    }

    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Finalize();
    if ( err > 1e-12 )
    {
        if ( Comm->MyPID() == 0 )
        {
            std::cout << "\nTest failed!\n";
        }
        returnValue = EXIT_FAILURE; // Norm of solution did not match
    }
    else
    {
        returnValue = EXIT_SUCCESS;
    }



    return ( returnValue );
}

