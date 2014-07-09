//HEADER
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
    @brief Test for ElectroETAbidomainSolver and IonicMinimalModel

    This example make use of a bidomain solver which should not
    be considered as a solver! If you would like to implement
    a bidomain solver, please make write it from scratch.

    @date 08 - 013
    @author Toni Lassila <toni.lassila@epfl.ch>

    @contributor
    @mantainer Toni Lassila <toni.lassila@epfl.ch>
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
#include <lifev/electrophysiology/examples/example_bidomain/ElectroETABidomainSolver.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/electrophysiology/solver/IonicModels/IonicMinimalModel.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/LifeV.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"


using namespace LifeV;

Real smoothing (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{

    Real pacingSite_X = 0.0;
    Real pacingSite_Y = 0.0;
    Real pacingSite_Z = 0.0;
    Real stimulusRadius = 0.15;

    if ( std::abs ( x - pacingSite_X ) <= stimulusRadius && std::abs ( z - pacingSite_Z ) <= stimulusRadius && std::abs ( y - pacingSite_Y ) <= stimulusRadius)
    {
        return 1.0;
    }
    else
    {
        return 0.0;
    }
}

Real PacingProtocol ( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID&   /*id*/)
{

    Real pacingSite_X = 0.0;
    Real pacingSite_Y = 0.0;
    Real pacingSite_Z = 0.0;
    Real stimulusRadius = 0.15;
    Real stimulusValue = 10.0;

    Real returnValue1;

    if ( std::abs ( x - pacingSite_X ) <= stimulusRadius && std::abs ( z - pacingSite_Z ) <= stimulusRadius && std::abs ( y - pacingSite_Y ) <= stimulusRadius)
    {
        returnValue1 = stimulusValue;
    }
    else
    {
        returnValue1 = 0.;
    }

    return returnValue1;
}

static Real bcDirichletZero (const Real&, const Real&, const Real&, const Real&, const LifeV::ID&)
{
    return 0.000;
}

static Real bcDirichletOne (const Real&, const Real&, const Real&, const Real&, const LifeV::ID&)
{
    return 1;
}


Int main ( Int argc, char** argv )
{

    //! Initializing Epetra communicator
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm>  Comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "% using MPI" << std::endl;
    }

    //********************************************//
    // Starts the chronometer.                    //
    //********************************************//
    //  LifeChrono chronoinitialsettings;
    //  chronoinitialsettings.start();

    typedef RegionMesh<LinearTetra>                         mesh_Type;
    typedef boost::function < Real (const Real& /*t*/,
                                    const Real &   x,
                                    const Real &   y,
                                    const Real& /*z*/,
                                    const ID&   /*i*/ ) >   function_Type;

    typedef ElectroETABidomainSolver< mesh_Type, IonicMinimalModel > bidomainSolver_Type;
    typedef boost::shared_ptr< bidomainSolver_Type >                 bidomainSolverPtr_Type;
    typedef VectorEpetra                                               vector_Type;
    typedef boost::shared_ptr<vector_Type>                             vectorPtr_Type;

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
    Teuchos::ParameterList bidomainList = * ( Teuchos::getParametersFromXmlFile ( "BidomainSolverParamList.xml" ) );
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
        std::cout << "Building Constructor for Minimal Model with parameters ... ";
    }
    boost::shared_ptr<IonicMinimalModel>  model ( new IonicMinimalModel() );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }


    //********************************************//
    // In the parameter list we need to specify   //
    // the mesh name and the mesh path.           //
    //********************************************//
    std::string meshName = bidomainList.get ("mesh_name", "lid16.mesh");
    std::string meshPath = bidomainList.get ("mesh_path", "./");

    //********************************************//
    // We need the GetPot datafile for to setup   //
    // the preconditioner.                        //
    //********************************************//
    GetPot dataFile (argc, argv);

    //********************************************//
    // We create three solvers to solve with:     //
    // 1) Operator Splitting method               //
    // 2) Ionic Current Interpolation             //
    // 3) State Variable Interpolation            //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Building bidomain Solvers... ";
    }

    bidomainSolverPtr_Type splitting ( new bidomainSolver_Type ( meshName, meshPath, dataFile, model ) );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << " Splitting solver done... ";
    }


    //********************************************//
    // Setting up the initial condition form      //
    // a given function.                          //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nInitializing potential:  " ;
    }

    // Initial pacing
    function_Type pacing = &PacingProtocol;
    splitting -> initializePotential();

    //HeartUtility::setValueOnBoundary( *(splitting -> potentialTransPtr() ), splitting -> fullMeshPtr(), 1.0, 6 );

    //function_Type f = &smoothing;
    //vectorPtr_Type smoother( new vector_Type( splitting -> potentialTransPtr() -> map() ) );
    //splitting -> feSpacePtr() -> interpolate ( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type > ( f ), *smoother , 0);
    //(*smoother) *= *(splitting -> potentialTransPtr() );
    //splitting -> setPotentialTransPtr(smoother);

    //setting up initial conditions
    * ( splitting -> globalSolution().at (1) ) = 1.0;
    * ( splitting -> globalSolution().at (2) ) = 1.0;
    * ( splitting -> globalSolution().at (3) ) = 0.021553043080281;

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Done! \n" ;
    }

    //********************************************//
    // Setting up the time data                   //
    //********************************************//
    splitting -> setParameters ( bidomainList );

    //********************************************//
    // Settung up the Boundary condition          //
    //********************************************//

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "-- Reading bc ..";
    }

    boost::shared_ptr<LifeV::BCHandler> bcs (new LifeV::BCHandler() );

    LifeV::BCFunctionBase zero (bcDirichletZero);

    std::vector<LifeV::ID> compx (1, 0), compy (1, 1), compz (1, 2);
    bcs->addBC ("boundaryDirichletZero", 600, LifeV::Essential, LifeV::Full, zero, 3);

    //bcs->addBC("boundaryNeumannBase", 99, LifeV::Natural, LifeV::Full, one,1);
    //partition mesh
    bcs->bcUpdate ( *splitting->localMeshPtr(), splitting->feSpacePtr()->feBd(), splitting->feSpacePtr()->dof() );

    if ( Comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }


    //********************************************//
    // Create a fiber direction                   //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nImporting fibers:  " ;
    }

    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > Space3D
    ( new FESpace< mesh_Type, MapEpetra > ( splitting -> localMeshPtr(), "P1", 3, splitting -> commPtr() ) );

    boost::shared_ptr<VectorEpetra> fiber ( new VectorEpetra ( Space3D -> map() ) );
    std::string nm = bidomainList.get ("fiber_file", "FiberDirection") ;
    ElectrophysiologyUtility::setupFibers ( *fiber, 0.0, 0.0, 1.0 );
    splitting -> setFiberPtr (fiber);

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

    //splitting -> setupLumpedMassMatrix();
    splitting -> setupMassMatrix();
    splitting -> setupStiffnessMatrix();
    splitting -> setupGlobalMatrix();
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Done! \n" ;
    }

    //********************************************//
    // Creating exporters to save the solution    //
    //********************************************//
    ExporterHDF5< RegionMesh <LinearTetra> > exporterSplitting;

    splitting -> setupExporter ( exporterSplitting, bidomainList.get ("OutputFile", "Splitting") );

    splitting -> exportSolution ( exporterSplitting, 0);

    //********************************************//
    // Solving the system                         //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "\nstart solving:  " ;
    }

    bidomainSolver_Type::vectorPtr_Type dtVec ( new VectorEpetra ( splitting->feSpacePtr() -> map(), LifeV::Unique ) );
    ExporterHDF5<mesh_Type> Exp;
    Exp.setMeshProcId ( splitting -> localMeshPtr(), splitting -> commPtr() -> MyPID() );
    Exp.setPrefix (bidomainList.get ("OutputTimeSteps", "TimeSteps") );
    Exp.addVariable ( ExporterData<mesh_Type>::ScalarField,  "dt", splitting->feSpacePtr(), dtVec, UInt (0) );

    Real dt = bidomainList.get ("timeStep", 0.1);
    Real TF = bidomainList.get ("endTime", 150.0);
    Int iter = bidomainList.get ("saveStep", 1.0) / dt;
    Int meth = bidomainList.get ("meth", 1 );
    Real dt_min = dt / 50.0;
    Int k (0), j (0);
    Int nodes;

    Real timeReac = 0.0;
    Real timeDiff = 0.0;
    LifeChrono chrono;

    Real stimulusStart = 0.0;
    Real stimulusStop  = 2.0;

    for ( Real t = 0.0; t < TF; )
    {

        if (  (t >= stimulusStart ) &&   (t <=  stimulusStop + dt) )
        {
            splitting -> setAppliedCurrentFromFunctionIntra ( pacing );
            //   splitting->appliedCurrentIntraPtr()->showMe();
        }
        else
        {
            splitting -> initializeAppliedCurrentIntra();
        }

        chrono.reset();
        if (meth == 1)
        {
            chrono.start();
            //            splitting->solveOneReactionStepROS3P(dtVec, dt_min);
            chrono.stop();
        }
        else
        {
            chrono.start();
            splitting->solveOneReactionStepFE( );
            chrono.stop();
        }

        timeReac += chrono.globalDiff ( *Comm );

        (*splitting->rhsPtrUnique() ) *= 0.0;
        splitting->updateRhs();

        chrono.reset();
        //bcManage ( *splitting->stiffnessMatrixPtr() , *splitting->rhsPtrUnique(), *splitting->localMeshPtr(), splitting->feSpacePtr()->dof(), *bcs,  splitting->feSpacePtr()->feBd(), 1.0, t );
        //splitting->rhsPtrUnique()->spy("rhs.dat");
        //splitting->stiffnessMatrixPtr()->spy("Stiffness");
        chrono.start();
        splitting->solveOneDiffusionStepBE();
        chrono.stop();
        timeDiff += chrono.globalDiff ( *Comm );

        if ( k % iter == 0 )
        {
            splitting -> exportSolution (exporterSplitting, t);
            Exp.postProcess (t);
        }

        nodes = dtVec->epetraVector().MyLength();
        j = dtVec->blockMap().GID (0);
        dt_min = (*dtVec) [j];
        for (int i = 1; i < nodes; i++)
        {
            j = dtVec->blockMap().GID (i);
            if (dt_min > (*dtVec) [j])
            {
                dt_min = (*dtVec) [j];
            }
        }

        k++;

        t = t + dt;

        if ( Comm->MyPID() == 0 )
        {
            std::cout << "\n\n\nActual time : " << t << std::endl << std::endl << std::endl;
        }

    }

    Real normSolution = ( ( splitting -> globalSolution().at (0) )->norm2() );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "2-norm of potential solution: " << normSolution << std::endl;
    }

    exporterSplitting.closeFile();
    Exp.closeFile();


    //********************************************//
    // Saving Fiber direction to file             //
    //********************************************//
    splitting -> exportFiberDirection();


    if ( Comm->MyPID() == 0 )
    {
        chronoinitialsettings.stop();
        std::cout << "\n\n\nTotal lapsed time : " << chronoinitialsettings.diff() << std::endl;
        std::cout << "Diffusion time : " << timeDiff << std::endl;
        std::cout << "Reaction time : " << timeReac << std::endl;
        std::cout << "\nThank you for using ETA_bidomainSolver.\nI hope to meet you again soon!\n All the best for your simulation :P\n  " ;
    }
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Finalize();

    Real returnValue;

    if (std::abs (normSolution - 14.182) > 1e-4 )
    {
        returnValue = EXIT_FAILURE; // Norm of solution did not match
    }
    else
    {
        returnValue = EXIT_SUCCESS;
    }
    return ( returnValue );
}
