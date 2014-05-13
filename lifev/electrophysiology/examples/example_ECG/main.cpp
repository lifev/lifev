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
    @brief 0D test with the Negroni Lascano model of 1996.

    @date 01−2013
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

#include <lifev/electrophysiology/solver/IonicModels/IonicAlievPanfilov.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicMinimalModel.hpp>
#include <lifev/core/LifeV.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"
// ---------------------------------------------------------------
// In order to use the ETA framework, a special version of the
// FESpace structure must be used. It is called ETFESpace and
// has basically the same role as the FESpace.
// ---------------------------------------------------------------

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

//--------------------------------------------------------
// For the pseudo- ECG
//--------------------------------------------------------
#include <lifev/electrophysiology/examples/example_ECG/Norm.hpp>
#include <lifev/core/solver/ADRAssembler.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>

// ---------------------------------------------------------------
// The most important file to include is the Integrate.hpp file
// which contains all the definitions required to perform the
// different integrations.
// ---------------------------------------------------------------

//#include <lifev/eta/expression/Integrate.hpp>
//
//#include <lifev/eta/expression/ExpressionDot.hpp>


using std::cout;
using std::endl;
using namespace LifeV;

// Choice of the fibers direction : ||.||=1
Real Fibers (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& z, const ID& i)
{
    Teuchos::ParameterList monodomainList = * ( Teuchos::getParametersFromXmlFile ( "MonodomainSolverParamList.xml" ) );
    Real zmin = 0.0;
    Real zmax = monodomainList.get ("domain_Z", 1. ); // 1.0;
    Real L = zmax - zmin;
    Real thetamin = M_PI / 3.;
    Real thetamax = -M_PI / 3.;

    Real ztheta = (L - z) / L;
    Real theta = (thetamax - thetamin) * ztheta + thetamin;

    switch (i)
    {
        case 0:
            return  std::cos (theta); // x_fib; //y/N; //std::sqrt (2.0) / 2.0; //
            break;
        case 1:
            return  std::sin (theta); // y_fib; //-x/N; //std::sqrt (2.0) / 2.0; //
            break;
        case 2:
            return 0.0;
            break;
        default:
            return 0.0;
            break;
    }
}

Real Stimulus2 (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    Teuchos::ParameterList monodomainList = * ( Teuchos::getParametersFromXmlFile ( "MonodomainSolverParamList.xml" ) );
    Real pacingSite_X = monodomainList.get ("pacingSite_X", 0.);
    Real pacingSite_Y = monodomainList.get ("pacingSite_Y", 0.);
    Real pacingSite_Z = monodomainList.get ("pacingSite_Z", 0. );
    Real stimulusRadius = 0.1; // monodomainList.get ("stimulusRadius", 0.1);

    if (  ( ( x - pacingSite_X ) * ( x - pacingSite_X ) +  ( y - pacingSite_Y ) * ( y - pacingSite_Y ) +  ( z - pacingSite_Z ) * ( z - pacingSite_Z )  )
            <= ( stimulusRadius * stimulusRadius ) )
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
    Teuchos::ParameterList monodomainList = * ( Teuchos::getParametersFromXmlFile ( "MonodomainSolverParamList.xml" ) );
    Real pacingSite_X = monodomainList.get ("pacingSite_X", 0.);
    Real pacingSite_Y = monodomainList.get ("pacingSite_Y", 0.);
    Real pacingSite_Z = monodomainList.get ("pacingSite_Z", 1.);
    Real stimulusRadius = monodomainList.get ("stimulusRadius", 0.1);
    Real stimulusValue = monodomainList.get ("stimulusValue", 1.);

    Real returnValue1;

    // --- Pacing protocol parameters -----------------------------------

    std::vector<double> returnPeriods;
    std::vector<double> returnStimulusTime;

    if (  ( ( x - pacingSite_X ) * ( x - pacingSite_X ) +  ( y - pacingSite_Y ) * ( y - pacingSite_Y ) +  ( z - pacingSite_Z ) * ( z - pacingSite_Z )  )
            <= ( stimulusRadius * stimulusRadius ) )
    {
        returnValue1 = stimulusValue;
    }
    else
    {
        returnValue1 = 0.;
    }

    Real returnValue = returnValue1;
    return returnValue;
}


Int main ( Int argc, char** argv )
{
    //! Initializing Epetra communicator
    MPI_Init (&argc, &argv);
    {
        boost::shared_ptr<Epetra_Comm>  Comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
        if ( Comm->MyPID() == 0 )
        {
            cout << "% using MPI" << endl;
        }

        //********************************************//
        // Starts the chronometer.                    //
        //********************************************//
        //  LifeChrono chronoinitialsettings;
        //  chronoinitialsettings.start();

        typedef RegionMesh<LinearTetra>            mesh_Type;
        typedef boost::shared_ptr<mesh_Type>       meshPtr_Type;
        typedef boost::shared_ptr<VectorEpetra>    vectorPtr_Type;
        typedef FESpace< mesh_Type, MapEpetra >    feSpace_Type;
        typedef boost::shared_ptr<feSpace_Type>    feSpacePtr_Type;
        typedef boost::function < Real (const Real& /*t*/,
                                        const Real &   x,
                                        const Real &   y,
                                        const Real& /*z*/,
                                        const ID&   /*i*/ ) >   function_Type;
        typedef VectorEpetra                        vector_Type;
        typedef MatrixEpetra<Real>                  matrix_Type;
        typedef LifeV::Preconditioner               basePrec_Type;
        typedef boost::shared_ptr<basePrec_Type>    basePrecPtr_Type;
        typedef LifeV::PreconditionerML             prec_Type;
        typedef boost::shared_ptr<prec_Type>        precPtr_Type;
        typedef ElectroETAMonodomainSolver< mesh_Type, IonicMinimalModel >        monodomainSolver_Type;
        typedef boost::shared_ptr< monodomainSolver_Type >                        monodomainSolverPtr_Type;

        typedef boost::shared_ptr< LifeV::Exporter<LifeV::RegionMesh<LifeV::LinearTetra> > >    filterPtr_Type;
        typedef LifeV::ExporterHDF5< RegionMesh<LinearTetra> >                                  hdf5Filter_Type;
        typedef boost::shared_ptr<hdf5Filter_Type>                                              hdf5FilterPtr_Type;
        typedef ExporterData<mesh_Type>                           exporterData_Type;
        typedef Exporter< mesh_Type >                             IOFile_Type;
        typedef boost::shared_ptr< IOFile_Type >                  IOFilePtr_Type;
        typedef ExporterHDF5< mesh_Type >                         hdf5IOFile_Type;

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
            std::cout << "Building Constructor for Aliev Panfilov Model with parameters ... ";
        }
        boost::shared_ptr<IonicMinimalModel>  model ( new IonicMinimalModel() );
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
        GetPot command_line (argc, argv);
        const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
        GetPot dataFile (data_file_name);

        // +-----------------------------------------------+
        // |               Loading the mesh                |
        // +-----------------------------------------------+
        if ( Comm->MyPID() == 0 )
        {
            std::cout << std::endl << "[Loading the mesh]" << std::endl;
        }

        meshPtr_Type fullMeshPtr ( new mesh_Type ( Comm ) );

        std::vector<Real> meshDim (3, 0);
        meshDim[0] =  monodomainList.get ("meshDim_X", 10 );
        meshDim[1] =  monodomainList.get ("meshDim_Y", 10 );
        meshDim[2] =  monodomainList.get ("meshDim_Z", 10 );
        std::vector<Real> domain (3, 0);
        domain[0] =  monodomainList.get ("domain_X", 1. );
        domain[1] =  monodomainList.get ("domain_Y", 1. );
        domain[2] =  monodomainList.get ("domain_Z", 1. );

        // Building the mesh from the source
        regularMesh3D ( *fullMeshPtr,
                        1,
                        meshDim[0], meshDim[1], meshDim[2],
                        false,
                        domain[0], domain[1], domain[2],
                        0.0, 0.0, 0.0 );

        if ( Comm->MyPID() == 0 )
        {
            std::cout << "Mesh size  : " << MeshUtility::MeshStatistics::computeSize ( *fullMeshPtr ).maxH << std::endl;
        }
        if ( Comm->MyPID() == 0 )
        {
            std::cout << "Partitioning the mesh ... " << std::endl;
        }
        meshPtr_Type meshPtr;
        {
            MeshPartitioner< mesh_Type > meshPart ( fullMeshPtr, Comm );
            meshPtr = meshPart.meshPartition();
        }
        fullMeshPtr.reset(); //Freeing the global mesh to save memory

        //********************************************//
        // We create the solvers to solve with:     //
        // Operator Splitting method               //
        //********************************************//
        if ( Comm->MyPID() == 0 )
        {
            std::cout << "Building Monodomain Solvers... ";
        }

        monodomainSolverPtr_Type splitting ( new monodomainSolver_Type ( dataFile, model, meshPtr ) );
        const feSpacePtr_Type FESpacePtr =  splitting->feSpacePtr(); //FE Space

        if ( Comm->MyPID() == 0 )
        {
            std::cout << " Splitting solver done... ";
        }

        function_Type pacing = &PacingProtocol;
        function_Type f = &Stimulus2;

        //    //***********************************************//
        //    //              RESTART Protocol                 //
        //    //***********************************************//
        //
        //    string filenameStart =  monodomainList.get ("StartFile", "MMStartSplitting" );
        //    std::string sol ( filenameStart ); // name of the file from which we want to restart
        //    string iterationString =  monodomainList.get ("StartIteration", "00050" );
        //
        //
        //    // Setting up the importer
        //    filterPtr_Type importer ( new hdf5Filter_Type( ) );
        //    importer->setMeshProcId ( splitting -> localMeshPtr(), Comm -> MyPID() );
        //    importer->setPrefix ( sol );
        //
        //    // Import the value of the potential and gating variable
        //    monodomainSolver_Type::vectorPtr_Type newSol;
        //    newSol.reset(new VectorEpetra ( splitting -> feSpacePtr() -> map(), LifeV::Unique ));
        //
        //    exporterData_Type ImportDataPotential (exporterData_Type::ScalarField, ("Variable0." + iterationString) , splitting -> feSpacePtr(),
        //                                              newSol, UInt (0), exporterData_Type::UnsteadyRegime);
        //
        //    monodomainSolver_Type::vectorPtr_Type newGating1;
        //    newGating1.reset(new VectorEpetra ( splitting -> feSpacePtr() -> map(), LifeV::Unique ));
        //
        //    exporterData_Type ImportDataGating1 (exporterData_Type::ScalarField, ("Variable1." + iterationString), splitting -> feSpacePtr(),
        //                                                  newGating1, UInt (0), exporterData_Type::UnsteadyRegime);
        //
        //    monodomainSolver_Type::vectorPtr_Type newGating2;
        //    newGating2.reset(new VectorEpetra ( splitting -> feSpacePtr() -> map(), LifeV::Unique ));
        //
        //    exporterData_Type ImportDataGating2 (exporterData_Type::ScalarField, ( "Variable2." + iterationString), splitting -> feSpacePtr(),
        //                                                  newGating2, UInt (0), exporterData_Type::UnsteadyRegime);
        //
        //    monodomainSolver_Type::vectorPtr_Type newGating3;
        //    newGating3.reset(new VectorEpetra ( splitting -> feSpacePtr() -> map(), LifeV::Unique ));
        //
        //    exporterData_Type ImportDataGating3 (exporterData_Type::ScalarField, ( "Variable3." + iterationString), splitting -> feSpacePtr(),
        //                                                  newGating3, UInt (0), exporterData_Type::UnsteadyRegime);
        //
        //    importer->readVariable ( ImportDataPotential );
        //    importer->readVariable( ImportDataGating1 );
        //    importer->readVariable( ImportDataGating2 );
        //    importer->readVariable( ImportDataGating3 );
        //
        //    //********************************************//
        //    // Setting up the initial condition form      //
        //    // a given function.                          //
        //    //********************************************//
        //    if ( Comm->MyPID() == 0 )
        //    {
        //        cout << "\nInitializing potential:  " ;
        //    }
        //
        //    splitting-> setPotentialPtr(newSol);
        //    * ( splitting -> globalSolution().at (1) ) = *newGating1;
        //    * ( splitting -> globalSolution().at (2) ) = *newGating2;
        //    * ( splitting -> globalSolution().at (3) ) = *newGating3;

        //*****************************************************************************************************
        // END OF RESTART
        //*****************************************************************************************************

        //********************************************//
        // Setting up the initial condition form      //
        // a given function.                          //
        //********************************************//
        if ( Comm->MyPID() == 0 )
        {
            cout << "\nInitializing potential:  " ;
        }
        //Compute the potential at t0

        splitting -> setPotentialFromFunction ( f );
        //setting up initial conditions
        * ( splitting -> globalSolution().at (1) ) = 1.0;
        * ( splitting -> globalSolution().at (2) ) = 1.0;
        * ( splitting -> globalSolution().at (3) ) = 0.021553043080281;

        // APD calculation variables
        Int sz = 0;
        sz = (* (splitting->globalSolution().at (0) ) ).size();
        Real threshold = monodomainList.get ("threshold", 0.1);
        Real trep = 0.;
        vector_Type tact = (* (splitting->globalSolution().at (0) ) );
        tact *= 0;
        vector_Type apd = tact;
        vector_Type delta_apd = tact;
        vector_Type previouspotential = (* (splitting->globalSolution().at (0) ) );



        if ( Comm->MyPID() == 0 )
        {
            cout << "Done! \n" ;
        }

        //********************************************//
        // Setting up the pacing protocol             //
        //********************************************//
        bool PacingProtocol = monodomainList.get ("pacingProtocol", false);
        std::vector<double> returnStimulusTime;
        std::vector<double> returnPeriods;

        Real stimulusStart = monodomainList.get ("stimulusStart", 0.);
        Real stimulusStop = monodomainList.get ("stimulusStop", 0.05);
        Int stimulusNumber;
        Int NumberPacingPeriods;

        int i (0);

        if (PacingProtocol)
        {
            Real pacingPeriod = monodomainList.get ("pacingPeriod", 500.);
            Real pacingPeriodMin = monodomainList.get ("pacingPeriodMin", 400.);
            Real pacingDelta = monodomainList.get ("pacingDelta", 0.);
            stimulusNumber = monodomainList.get ("stimulusNumber", 1);
            NumberPacingPeriods = (pacingPeriod - pacingPeriodMin) / pacingDelta;
            //--- Pacing method
            if ( pacingDelta > 0 )    // IF pacing
            {
                for (int k = 0; k <= NumberPacingPeriods - 1; k++ )
                {
                    for (i = stimulusNumber * k; i <= stimulusNumber * (k + 1) - 1; i++)
                    {
                        returnPeriods.push_back ( pacingPeriod - k * pacingDelta );
                        if ( i == 0 )
                        {
                            returnStimulusTime.push_back ( returnPeriods[0] );
                        }
                        else
                        {
                            returnStimulusTime.push_back ( returnStimulusTime[i - 1] + returnPeriods[i] );
                        }
                    }
                }
            }
        }
        else
        {
            returnPeriods.push_back ( monodomainList.get ("stimulus0", 500.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus1", 450.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus2", 425.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus3", 400.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus4", 375.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus5", 350.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus6", 325.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus7", 300.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus8", 275.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus9", 275.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus10", 275.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus11", 275.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus12", 450.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus13", 400.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus14", 375.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus15", 350.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus16", 325.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus17", 300.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus18", 275.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus19", 275.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus20", 275.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus21", 275.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus22", 450.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus23", 400.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus24", 375.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus25", 350.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus26", 325.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus27", 300.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus28", 275.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus29", 275.) );
            returnPeriods.push_back ( monodomainList.get ("stimulus30", 275.) );

            NumberPacingPeriods = monodomainList.get ("NbStimulusPeriod", 12);

            stimulusNumber = 1;
            for (int k = 0; k <= NumberPacingPeriods - 1; k++ )
            {
                if (k == 0)
                {
                    returnStimulusTime.push_back ( returnPeriods[0] );
                }
                else
                {
                    returnStimulusTime.push_back ( returnStimulusTime[k - 1] + returnPeriods[k] );
                }
            }
        }
        //*******************************************//
        // Setting up the pseudo-ECG                 //
        //*******************************************//
        Real ecg_position_X = monodomainList.get ("ecg_position_X", 1.);
        Real ecg_position_Y = monodomainList.get ("ecg_position_Y", 1.);
        Real ecg_position_Z = monodomainList.get ("ecg_position_Z", 0.5);
        vector_Type ecgDistance = (* (splitting->globalSolution().at (0) ) ) ;
        ecgDistance *= 0.;
        FESpacePtr->interpolate ( static_cast<function_Type> ( Norm::f ), ecgDistance, 0.0 );
        Norm::setPosition ( ecg_position_X , ecg_position_Y , ecg_position_Z ); // Set electrode position

        Real pseudoEcgReal (0.);
        Real Global_pseudoEcgReal (0.);
        vector_Type solutionLaplacian = ( (* (splitting->globalSolution().at (0) ) ) );
        vector_Type pseudoEcgVec = ( (* (splitting->globalSolution().at (0) ) ) );

        // Discrete Laplacian matrix
        // setting up the assembler
        ADRAssembler<mesh_Type, matrix_Type, vector_Type> adrAssembler;
        adrAssembler.setup ( FESpacePtr, FESpacePtr );
        // define the matrices
        boost::shared_ptr<matrix_Type> systemMatrixL ( new matrix_Type ( FESpacePtr->map() ) );
        boost::shared_ptr<matrix_Type> systemMatrixM ( new matrix_Type ( FESpacePtr->map() ) );
        boost::shared_ptr<vector_Type> rhs_Laplacian ( new vector_Type ( FESpacePtr->map() ) );
        boost::shared_ptr<vector_Type> pseudoEcgVec_ptr ( new vector_Type ( FESpacePtr->map() ) );

        // fill the matrix
        adrAssembler.addDiffusion ( systemMatrixL, -1.0 );
        adrAssembler.addMass ( systemMatrixM, 1.0 );
        // closed
        systemMatrixL->globalAssemble();
        systemMatrixM->globalAssemble();

        // uncomment to check the matrices with MATLAB
        //    matrix_Type LaplacianMatrix ( *systemMatrixL );
        //    matrix_Type MassMatrix ( *systemMatrixM );
        //    LaplacianMatrix.spy("matriceL_check");
        //    MassMatrix.spy("matriceM_check");
        //********************************************//
        // Setting up the time data                   //
        //********************************************//
        splitting -> setParameters ( monodomainList );

        //********************************************//
        // Create a fiber direction                   //
        //********************************************//
        VectorSmall<3> fibers;
        fibers[0] =  monodomainList.get ("fiber_X", std::sqrt (2.0) / 2.0 );
        fibers[1] =  monodomainList.get ("fiber_Y", std::sqrt (2.0) / 2.0 );
        fibers[2] =  monodomainList.get ("fiber_Z", 0.0 );

        splitting ->setupFibers (fibers);

        //********************************************//
        // Create the global matrix: mass + stiffness //
        //********************************************//
        splitting -> setupLumpedMassMatrix();
        splitting -> setupStiffnessMatrix();
        splitting -> setupGlobalMatrix();

        //********************************************//
        // Creating exporters to save the solution    //
        //********************************************//
        ExporterHDF5< RegionMesh <LinearTetra> > exporterSplitting;
        string filenameSplitting =  monodomainList.get ("OutputFile", "MinMod" );
        filenameSplitting += "Splitting";
        splitting -> setupPotentialExporter ( exporterSplitting, filenameSplitting );

        ExporterHDF5< RegionMesh <LinearTetra> > exporterSplittingRestart;
        string filenameSplittingLast =  monodomainList.get ("OutputFile", "MinMod" );
        filenameSplittingLast += "RESTART";
        splitting -> setupExporter ( exporterSplittingRestart, filenameSplittingLast );

        vectorPtr_Type APDptr ( new vector_Type (apd, Repeated ) );
        exporterSplitting.addVariable ( ExporterData<mesh_Type>::ScalarField,  "apd", FESpacePtr,
                                        APDptr, UInt (0) );

        vectorPtr_Type DELTA_APDptr ( new vector_Type (delta_apd, Repeated ) );
        exporterSplitting.addVariable ( ExporterData<mesh_Type>::ScalarField,  "delta_apd", FESpacePtr,
                                        DELTA_APDptr, UInt (0) );

        *APDptr = apd;
        *DELTA_APDptr = delta_apd;

        std::ofstream output  ("ecg_output.txt");

        //**************************************************//
        // Solver initialization for the discrete Laplacian //
        //**************************************************//
        if (  Comm->MyPID() == 0 )
        {
            std::cout << std::endl << "[Solvers initialization]" << std::endl;
        }
        prec_Type* precRawPtr;
        basePrecPtr_Type precPtr;
        precRawPtr = new prec_Type;
        precRawPtr->setDataFromGetPot ( dataFile, "prec" );
        precPtr.reset ( precRawPtr );
        if (  Comm->MyPID() == 0 )
        {
            std::cout << "Setting up LinearSolver (Belos)... " << std::flush;
        }
        Teuchos::RCP< Teuchos::ParameterList > belosList2 = Teuchos::rcp ( new Teuchos::ParameterList );
        belosList2 = Teuchos::getParametersFromXmlFile ( "SolverParamList2.xml" );
        LinearSolver linearSolver2;
        linearSolver2.setCommunicator ( Comm );
        linearSolver2.setParameters ( *belosList2 );
        linearSolver2.setPreconditioner ( precPtr );
        if (  Comm->MyPID() == 0 )
        {
            std::cout << "done" << std::endl;
        }
        linearSolver2.showMe();

        //********************************************//
        // Solving the system                         //
        //********************************************//
        if ( Comm->MyPID() == 0 )
        {
            cout << "\nstart solving:  " ;
        }

        Real Savedt = monodomainList.get ("saveStep", 0.1);
        Real timeStep = monodomainList.get ("timeStep", 0.01);
        Real endTime = monodomainList.get ("endTime", 10.);
        Real initialTime = monodomainList.get ("initialTime", 0.);

        vectorPtr_Type previousPotential0Ptr ( new vector_Type ( FESpacePtr->map() ) );
        * (previousPotential0Ptr) = * (splitting->globalSolution().at (0) );
        int control = 0;
        int iter ( (Savedt / timeStep) + 1e-9);
        int iterRestart ( (500 / timeStep) + 1e-9);
        int nbTimeStep (1);
        int k (0);
        if (endTime > timeStep)
        {
            for (Real t = initialTime; t < endTime;)
            {

                // APD calculation
                previouspotential = (* (splitting->globalSolution().at (0) ) );
                //--------------------------------------
                // ECG initialization
                pseudoEcgReal = 0.;
                pseudoEcgVec_ptr.reset ( new vector_Type ( FESpacePtr->map(), Unique ) );
                rhs_Laplacian.reset ( new vector_Type ( FESpacePtr->map(), Unique ) );
                //-----------------------------------------------------------------

                control = 0;

                t += timeStep;
                for (i = 0; i <= NumberPacingPeriods * stimulusNumber - 1; i++)
                {
                    if ( control < 1 )
                    {
                        if ( (t >= (returnStimulusTime[i] + stimulusStart) && t <= (returnStimulusTime[i] + stimulusStop + timeStep) ) )
                        {
                            splitting -> setAppliedCurrentFromFunction ( pacing );
                            std::cout << "stim_time " << t << std::endl;
                            control = control + 1;
                        }
                        else
                        {
                            splitting -> initializeAppliedCurrent();
                        }
                    }
                }
                k++;

                if (nbTimeStep == 1)
                {
                    splitting->solveOneReactionStepFE();
                    (* (splitting->rhsPtrUnique() ) ) *= 0;
                    splitting->updateRhs();
                    splitting->solveOneDiffusionStepBE();
                    splitting->exportSolution (exporterSplitting, t);
                }
                else
                {
                    * (previousPotential0Ptr) = * (splitting->globalSolution().at (0) );
                    splitting->solveOneReactionStepFE (2);
                    (* (splitting->rhsPtrUnique() ) ) *= 0;
                    splitting->updateRhs();
                    splitting->solveOneDiffusionStepBDF2 (previousPotential0Ptr);
                    splitting->solveOneReactionStepFE (2);
                    if (k % iter == 0)
                    {
                        splitting->exportSolution (exporterSplitting, t);
                    }
                    if (k % iterRestart == 0)
                    {
                        splitting->exportSolution (exporterSplittingRestart, t);
                    }
                }
                nbTimeStep++;

                // ECG : discrete laplacian of the solution
                (*rhs_Laplacian) = (*systemMatrixL) * (* (splitting->globalSolution().at (0) ) );

                linearSolver2.setOperator ( systemMatrixM );
                linearSolver2.setRightHandSide ( rhs_Laplacian );
                linearSolver2.solve ( pseudoEcgVec_ptr );

                pseudoEcgVec = (*pseudoEcgVec_ptr) / ecgDistance;

                //      // APD calculation
                for (int i = 0; i <= sz - 1; i++)
                {
                    if ( (* (splitting->globalSolution().at (0) ) ).isGlobalIDPresent (i) )
                    {

                        if ( ( previouspotential[i] < threshold ) && ( (* (splitting->globalSolution().at (0) ) ) [i] >= threshold ) )
                        {
                            tact[i] = t - ( (-threshold + previouspotential[i]) / ( (* (splitting->globalSolution().at (0) ) ) [i] - previouspotential[i]) ) * timeStep;
                        }
                        else if ( ( previouspotential[i] >= threshold ) && ( (* (splitting->globalSolution().at (0) ) ) [i] < threshold ) )
                        {
                            trep = t - ( (-threshold + previouspotential[i]) / ( (* (splitting->globalSolution().at (0) ) ) [i] - previouspotential[i]) ) * timeStep;
                            delta_apd[i] = (trep - tact[i]) - apd[i];
                            apd[i] = trep - tact[i];
                        }

                        //          // Pseudo-ECG summation
                        pseudoEcgReal += pseudoEcgVec[i];
                    }
                }
                *APDptr = apd;
                *DELTA_APDptr = delta_apd;

                MPI_Allreduce (&pseudoEcgReal, &Global_pseudoEcgReal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // rapporte a une variable connue de tous les procs

                if (  Comm->MyPID() == 0 )
                {
                    output << Global_pseudoEcgReal << "\n";
                }
            }
        }
        splitting->exportSolution (exporterSplittingRestart, endTime);
        exporterSplitting.closeFile();

        //********************************************//
        // Saving Fiber direction to file             //
        //********************************************//
        splitting -> exportFiberDirection();

        if ( Comm->MyPID() == 0 )
        {
            cout << "\nThank you for using ETA_MonodomainSolver.\nI hope to meet you again soon!\n All the best for your simulation :P\n  " ;
        }
        MPI_Barrier (MPI_COMM_WORLD);
    }
    MPI_Finalize();
    return ( EXIT_SUCCESS );
}
