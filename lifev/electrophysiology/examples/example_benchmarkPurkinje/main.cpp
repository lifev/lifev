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
    @brief Electrophysiology benchmark with Noble-Purkinje excitation model

    @date 01−2013
    @author Simone Rossi <simone.rossi@epfl.ch>
    @author Simone Palamara <palamara.simone@gmail.com>

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



#include <fstream>
#include <string>

#include <lifev/core/array/VectorSmall.hpp>

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/electrophysiology/solver/IonicModels/IonicMinimalModel.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicLuoRudyI.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicTenTusscher06.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicHodgkinHuxley.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicNoblePurkinje.hpp>
#include <lifev/electrophysiology/util/ElectrophysiologyUtility.hpp>
#include <lifev/core/LifeV.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <sys/stat.h>

using namespace LifeV;

Real PacingProtocolMM ( const Real& t, const Real& x, const Real& y, const Real& z, const ID&   /*id*/)
{

    Real pacingSite_X = 0.0;
    Real pacingSite_Y = 0.0;
    Real pacingSite_Z = 0.0;
    Real stimulusRadius = 0.15;
    Real stimulusValue = 10;

    Real returnValue;

    if ( std::abs ( x - pacingSite_X ) <= stimulusRadius
            &&
            std::abs ( z - pacingSite_Z ) <= stimulusRadius
            &&
            std::abs ( y - pacingSite_Y ) <= stimulusRadius
            &&
            t <= 2)
    {
        returnValue = stimulusValue;
    }
    else
    {
        returnValue = 0.;
    }

    return returnValue;
}

Real PacingProtocolHH ( const Real& t, const Real& x, const Real& y, const Real& z, const ID&   /*id*/)
{

    Real pacingSite_X = 0.0;
    Real pacingSite_Y = 0.0;
    Real pacingSite_Z = 0.0;
    Real stimulusRadius = 0.15;
    Real stimulusValue = 500.;

    Real returnValue;

    if ( std::abs ( x - pacingSite_X ) <= stimulusRadius
            &&
            std::abs ( z - pacingSite_Z ) <= stimulusRadius
            &&
            std::abs ( y - pacingSite_Y ) <= stimulusRadius
            &&
            t <= 2)
    {
        returnValue = stimulusValue;
    }
    else
    {
        returnValue = 0.;
    }

    return returnValue;
}

Real PacingProtocol ( const Real& t, const Real& x, const Real& y, const Real& z, const ID&   /*id*/)
{

    Real pacingSite_X = 0.0;
    Real pacingSite_Y = 0.0;
    Real pacingSite_Z = 0.0;
    Real stimulusRadius = 0.15;
    Real stimulusValue = 50;

    Real returnValue;

    if ( std::abs ( x - pacingSite_X ) <= stimulusRadius
            &&
            std::abs ( z - pacingSite_Z ) <= stimulusRadius
            &&
            std::abs ( y - pacingSite_Y ) <= stimulusRadius
            &&
            t <= 2)
    {
        returnValue = stimulusValue;
    }
    else
    {
        returnValue = 0.;
    }

    return returnValue;
}


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
    // Starts the chronometer.                    //
    //********************************************//

    typedef RegionMesh<LinearTetra>                         mesh_Type;
    typedef std::function < Real (const Real& /*t*/,
                                    const Real &   x,
                                    const Real &   y,
                                    const Real& /*z*/,
                                    const ID&   /*i*/ ) >   function_Type;

    typedef ElectroIonicModel                                        ionicModel_Type;
    typedef std::shared_ptr<ionicModel_Type>                       ionicModelPtr_Type;
    typedef ElectroETAMonodomainSolver< mesh_Type, ionicModel_Type > monodomainSolver_Type;
    typedef std::shared_ptr< monodomainSolver_Type >               monodomainSolverPtr_Type;

    typedef VectorEpetra                                             vector_Type;
    typedef std::shared_ptr<vector_Type>                           vectorPtr_Type;


    LifeChrono chronoinitialsettings;

    if ( Comm->MyPID() == 0 )
    {
        chronoinitialsettings.start();
    }

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

    std::string ionic_model ( monodomainList.get ("ionic_model", "minimalModel") );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nIonic_Model:" << ionic_model;
    }



    //********************************************//
    // Creates a new model object representing the//
    // model from Aliev and Panfilov 1996.  The   //
    // model input are the parameters. Pass  the  //
    // parameter list in the constructor          //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nBuilding Constructor for " << ionic_model << " Model with parameters ... ";
    }
    ionicModelPtr_Type  model;
    std::shared_ptr<IonicNoblePurkinje> excitationModel (new IonicNoblePurkinje() );

    if ( ionic_model == "LuoRudyI" )
    {
        model.reset ( new IonicLuoRudyI() );
    }
    if ( ionic_model == "TenTusscher06")
    {
        model.reset (new IonicTenTusscher06() );
    }
    if ( ionic_model == "HodgkinHuxley")
    {
        model.reset (new IonicHodgkinHuxley() );
    }
    if ( ionic_model == "NoblePurkinje")
    {
        model.reset (new IonicNoblePurkinje() );
    }
    if ( ionic_model == "MinimalModel")
    {
        model.reset ( new IonicMinimalModel() );
    }

    if ( Comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }

    model -> showMe();


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
    GetPot command_line (argc, argv);
    const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);

    //********************************************//
    // We create three solvers to solve with:     //
    // 1) Operator Splitting method               //
    // 2) Ionic Current Interpolation             //
    // 3) State Variable Interpolation            //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Building Monodomain Solvers... ";
    }

    monodomainSolverPtr_Type solver ( new monodomainSolver_Type ( meshName, meshPath, dataFile, model ) );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << " solver done... ";
    }




    //********************************************//
    // Setting up the initial condition form      //
    // a given function.                          //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nInitializing potential and gating variables:  " ;
    }

    // Initial pacing
    solver -> initializePotential();
    solver -> initializeAppliedCurrent();
    solver -> setInitialConditions();

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nDone.  " << std::flush ;
    }

    //********************************************//
    // MESH STUFF      //
    //********************************************//

    std::vector<Real> junction (3, 0.0);
    //Real Radius = 0.1;
    UInt numberOfSources (4);
    UInt numberLocalSources (numberOfSources / Comm->NumProc() );
    if ( Comm->MyPID() == 0 )
    {
        numberLocalSources += numberOfSources % Comm->NumProc();
    }
    UInt flag (300);
    std::vector<ID> containerPointsGivenFlag, containerPointsWithAppliedCurrent;
    std::vector<Real> activationTime;
    double deltaT (2.0);
    ElectrophysiologyUtility::allIdsPointsWithGivenFlag<mesh_Type> (containerPointsGivenFlag, flag, solver -> appliedCurrentPtr(),  solver -> fullMeshPtr() );
    ElectrophysiologyUtility::randomNPointsInSetAndNeighborhood<mesh_Type> (containerPointsGivenFlag, containerPointsWithAppliedCurrent, activationTime, deltaT, numberLocalSources,  * (solver -> fullMeshPtr() ), Comm );



    //********************************************//
    // Setting up the time data                   //
    //********************************************//
    solver -> setParameters ( monodomainList );

    //********************************************//
    // Create a fiber direction                   //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nSetting fibers:  " ;
    }

    VectorSmall<3> fibers;
    fibers[0] = 0.0;
    fibers[1] = 0.0;
    fibers[2] = 1.0;
    solver -> setupFibers ( fibers );

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Done! \n" ;
    }

    //********************************************//
    // Create the global matrix: mass + stiffness //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nSetup operators:  " ;
    }

    bool lumpedMass = monodomainList.get ("LumpedMass", true);
    if ( lumpedMass)
    {
        solver -> setupLumpedMassMatrix();
    }
    else
    {
        solver -> setupMassMatrix();
    }


    solver -> setupStiffnessMatrix ( solver -> diffusionTensor() );
    solver -> setupGlobalMatrix();
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Done! \n" ;
    }

    //********************************************//
    // Creating exporters to save the solution    //
    //********************************************//
    ExporterHDF5< RegionMesh <LinearTetra> > exporter;

    solver -> setupExporter ( exporter, monodomainList.get ("OutputFile", "Solution") );
    exporter.setPostDir (problemFolder);
    solver -> exportSolution ( exporter, 0);

    //********************************************//
    // Activation time                            //
    //********************************************//
    vectorPtr_Type activationTimeVector ( new vector_Type ( solver -> potentialPtr() -> map() ) );
    *activationTimeVector = -1.0;

    ExporterHDF5< RegionMesh <LinearTetra> > activationTimeExporter;
    activationTimeExporter.setMeshProcId (solver -> localMeshPtr(), solver -> commPtr() ->MyPID() );
    activationTimeExporter.addVariable (ExporterData<mesh_Type>::ScalarField, "Activation Time",
                                        solver -> feSpacePtr(), activationTimeVector, UInt (0) );
    activationTimeExporter.setPrefix ("ActivationTime");
    activationTimeExporter.setPostDir (problemFolder);

    //********************************************//
    // Solving the system                         //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nstart solving:  " ;
    }

    Real dt = monodomainList.get ("timeStep", 0.1);
    Real TF = monodomainList.get ("endTime", 150.0);
    Int iter = monodomainList.get ("saveStep", 1.0) / dt;
    Int k (0);

    Real timeReac = 0.0;
    Real timeDiff = 0.0;
    Real timeReacDiff = 0.0;
    LifeChrono chrono;

    std::string solutionMethod = monodomainList.get ("solutionMethod", "splitting");

    std::vector<Real> statesExcitationModel (excitationModel->Size(), 0);
    std::vector<Real> rhsExcitationModel (excitationModel->Size(), 0);
    std::vector<Real> valueAppliedCurrent (TF / dt, 0);
    std::vector<UInt> shiftVector (TF / dt, -1);
    excitationModel->initialize (statesExcitationModel);
    statesExcitationModel[0] = -70.0;
    statesExcitationModel[1] = excitationModel->mInf (-70.0);
    statesExcitationModel[2] = excitationModel->nInf (-70.0);
    statesExcitationModel[3] = excitationModel->hInf (-70.0);
    int localIndexTime (0);
    for ( Real t = 0.0; t < TF; )
    {
        //        solver -> setAppliedCurrentFromFunction ( stimulus, t );


        excitationModel->computeGatingVariablesWithRushLarsen (statesExcitationModel, dt);
        statesExcitationModel.at (0) = statesExcitationModel.at (0)  + dt * (excitationModel->computeLocalPotentialRhs ( statesExcitationModel) );
        valueAppliedCurrent[localIndexTime] = excitationModel->Itotal();
        localIndexTime++;
        //ElectrophysiologyUtility::appliedCurrentPointsWithinRadius<mesh_Type>(junction,Radius,solver -> appliedCurrentPtr(),excitationModel->Itotal()/85.7,solver -> fullMeshPtr() );


        //ElectrophysiologyUtility::applyCurrentGivenSetOfPoints<ID>(containerPointsWithAppliedCurrent,solver -> appliedCurrentPtr(), excitationModel->Itotal() );//
        ElectrophysiologyUtility::applyCurrentGivenSetOfPoints<ID> (containerPointsWithAppliedCurrent, activationTime, solver -> appliedCurrentPtr(), valueAppliedCurrent, shiftVector, t);
        if ( solutionMethod == "splitting" )
        {
            chrono.reset();
            chrono.start();
            if (ionic_model != "MinimalModel" && ionic_model != "HodgkinHuxley")
            {
                solver->solveOneReactionStepRL();
            }
            else
            {
                solver->solveOneReactionStepFE();
            }
            chrono.stop();

            timeReac += chrono.globalDiff ( *Comm );

            (*solver->rhsPtrUnique() ) *= 0.0;
            solver->updateRhs();

            chrono.reset();
            chrono.start();
            solver->solveOneDiffusionStepBE();
            chrono.stop();
            timeDiff += chrono.globalDiff ( *Comm );
        }
        else if ( solutionMethod == "ICI" )
        {
            chrono.reset();
            chrono.start();
            if (ionic_model != "MinimalModel" && ionic_model != "HodgkinHuxley")
            {
                solver -> solveOneStepGatingVariablesRL();
            }
            else
            {
                solver -> solveOneStepGatingVariablesFE();
            }
            solver -> solveOneICIStep();
            chrono.stop();
            timeReacDiff += chrono.globalDiff ( *Comm );
        }
        else if ( solutionMethod == "SVI" )
        {
            chrono.reset();
            chrono.start();
            if (ionic_model != "MinimalModel" && ionic_model != "HodgkinHuxley")
            {
                solver -> solveOneStepGatingVariablesRL();
            }
            else
            {
                solver -> solveOneStepGatingVariablesFE();
            }
            solver -> solveOneSVIStep();
            chrono.stop();
            timeReacDiff += chrono.globalDiff ( *Comm );
        }

        //register activation time
        k++;
        t = t + dt;
        solver -> registerActivationTime (*activationTimeVector, t, 0.95);

        if ( k % iter == 0 )
        {
            solver -> exportSolution (exporter, t);
        }
        if ( Comm->MyPID() == 0 )
        {
            std::cout << "\n\n\nActual time : " << t << std::endl << std::endl << std::endl;
        }
    }

    Real normSolution = ( ( solver -> globalSolution().at (0) )->norm2() );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "2-norm of potential solution: " << normSolution << std::endl;
    }

    exporter.closeFile();
    activationTimeExporter.postProcess (0);
    activationTimeExporter.closeFile();

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Exporting fibers: " << std::endl;
    }

    //********************************************//
    // Saving Fiber direction to file             //
    //********************************************//
    solver -> exportFiberDirection();


    if ( Comm->MyPID() == 0 )
    {
        chronoinitialsettings.stop();
        std::cout << "\n\n\nTotal lapsed time : " << chronoinitialsettings.diff() << std::endl;
        if ( solutionMethod == "splitting" )
        {
            std::cout << "Diffusion time : " << timeDiff << std::endl;
            std::cout << "Reaction time : " << timeReac << std::endl;
        }
        else if ( solutionMethod == "ICI" )
        {
            std::cout << "Solution time : " << timeReacDiff << std::endl;
        }
        else if ( solutionMethod == "SVI" )
        {
            std::cout << "Solution time : " << timeReacDiff << std::endl;
        }

        std::cout << "\n\nThank you for using ETA_MonodomainSolver.\nI hope to meet you again soon!\n All the best for your simulation :P\n  " ;
    }
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Finalize();

    return ( EXIT_SUCCESS );
}
