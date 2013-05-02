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
  @brief Cardiac Electrophysiology Test
  @author Lucia Mirabella <lucia.mirabella@mail.polimi.it> and Mauro Perego <mauro.perego@polimi.it>
  @date 11-2007
  @contributors Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>, Simone Rossi <simone.rossi@epfl.ch>
  @last update 11-2010
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

#include "heart.hpp"

using namespace LifeV;

typedef RegionMesh<LinearTetra> mesh_Type;

//! Identifiers for heart boundaries
const Int EPICARDIUM    = 40;
const Int ENDOCARDIUM   = 60;
const Int TRUNC_SEC     = 50;

Real zero_scalar ( const Real& /* t */,
                   const Real& /* x */,
                   const Real& /* y */,
                   const Real& /* z */,
                   const ID& /* i */ )
{
    return 0.;
}

// ===================================================
//! Constructors
// ===================================================

Heart::Heart ( Int argc,
               char** argv )
{
    GetPot command_line (argc, argv);
    const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);

    //! Pointer to access functors
    M_heart_fct.reset (new HeartFunctors ( dataFile) );
    ion_model = dataFile ("electric/physics/ion_model", 1);
    M_heart_fct->M_comm.reset (new Epetra_MpiComm ( MPI_COMM_WORLD ) );

    if (!M_heart_fct->M_comm->MyPID() )
    {
        std::cout << "My PID = " << M_heart_fct->M_comm->MyPID() << std::endl;
    }
}


// ===================================================
//! Methods
// ===================================================

void
Heart::run()
{
    typedef FESpace< mesh_Type, MapEpetra > feSpace_Type;
    typedef boost::shared_ptr<feSpace_Type> feSpacePtr_Type;

    LifeChrono chronoinitialsettings;
    LifeChrono chronototaliterations;
    chronoinitialsettings.start();
    Real normu;
    Real meanu;
    Real minu;

    //! Construction of data classes

#ifdef MONODOMAIN
    HeartMonodomainData _data (M_heart_fct);
#else
    HeartBidomainData _data (M_heart_fct);
#endif
    HeartIonicData _dataIonic (M_heart_fct->M_dataFile);

    MeshData meshData;
    meshData.setup (M_heart_fct->M_dataFile, "electric/space_discretization");
    boost::shared_ptr<mesh_Type > fullMeshPtr ( new mesh_Type ( M_heart_fct->M_comm ) );
    readMesh (*fullMeshPtr, meshData);
    bool verbose = (M_heart_fct->M_comm->MyPID() == 0);

    //! Boundary conditions handler and function
    BCFunctionBase uZero ( zero_scalar );
    BCHandler bcH;
    bcH.addBC ( "Endo",      ENDOCARDIUM,    Natural,    Full,   uZero,  1 );
    bcH.addBC ( "Epi",       EPICARDIUM,     Natural,    Full,   uZero,  1 );
    bcH.addBC ( "Trunc",     TRUNC_SEC,      Natural,    Full,   uZero,  1 );

    const ReferenceFE*    refFE_w;
    const QuadratureRule* qR_w;
    const QuadratureRule* bdQr_w;

    const ReferenceFE*    refFE_u;
    const QuadratureRule* qR_u;
    const QuadratureRule* bdQr_u;


    //! Construction of the partitioned mesh
    boost::shared_ptr<mesh_Type> localMeshPtr;
    {
        MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, M_heart_fct->M_comm);
        localMeshPtr = meshPart.meshPartition();
    }
    std::string uOrder =  M_heart_fct->M_dataFile ( "electric/space_discretization/u_order", "P1");

    //! Initialization of the FE type and quadrature rules for both the variables
    if ( uOrder.compare ("P1") == 0 )
    {
        if (verbose)
        {
            std::cout << "P1 potential " << std::flush;
        }
        refFE_u = &feTetraP1;
        qR_u    = &quadRuleTetra15pt;
        bdQr_u  = &quadRuleTria3pt;
    }
    else
    {
        cout << "\n " << uOrder << " finite element not implemented yet \n";
        exit (1);
    }

    std::string wOrder =  M_heart_fct->M_dataFile ( "electric/space_discretization/w_order", "P1");
    if ( wOrder.compare ("P1") == 0 )
    {
        if (verbose)
        {
            std::cout << "P1 recovery variable " << std::flush;
        }
        refFE_w = &feTetraP1;
        qR_w    = &quadRuleTetra4pt;
        bdQr_w  = &quadRuleTria3pt;
    }
    else
    {
        cout << "\n " << wOrder << " finite element not implemented yet \n";
        exit (1);
    }

    //! Construction of the FE spaces
    if (verbose)
    {
        std::cout << "Building the potential FE space ... " << std::flush;
    }

    feSpacePtr_Type uFESpacePtr ( new feSpace_Type ( localMeshPtr,
                                                     *refFE_u,
                                                     *qR_u,
                                                     *bdQr_u,
                                                     1,
                                                     M_heart_fct->M_comm) );

#ifdef BIDOMAIN
    feSpacePtr_Type _FESpacePtr ( new feSpace_Type (localMeshPtr,
                                                    *refFE_u,
                                                    *qR_u,
                                                    *bdQr_u,
                                                    2,
                                                    M_heart_fct->M_comm) );
#endif
    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }
    if (verbose)
    {
        std::cout << "Building the recovery variable FE space ... " << std::flush;
    }
    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    UInt totalUDof  = uFESpacePtr->map().map (Unique)->NumGlobalElements();
    if (verbose)
    {
        std::cout << "Total Potential DOF = " << totalUDof << std::endl;
    }
    if (verbose)
    {
        std::cout << "Calling the electric model constructor ... ";
    }

#ifdef MONODOMAIN
    HeartMonodomainSolver< mesh_Type > electricModel (_data, *uFESpacePtr, bcH, M_heart_fct->M_comm);
#else
    HeartBidomainSolver< mesh_Type > electricModel (_data, *_FESpacePtr, *uFESpacePtr, bcH, M_heart_fct->M_comm);
#endif

    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }
    MapEpetra fullMap (electricModel.getMap() );
    vector_Type rhs ( fullMap);
    electricModel.setup ( M_heart_fct->M_dataFile );
    std::cout << "setup ok" << std::endl;

    if (verbose)
    {
        std::cout << "Calling the ionic model constructor ... ";
    }
    boost::shared_ptr< HeartIonicSolver< mesh_Type > > ionicModel;
    if (ion_model == 1)
    {
        if (verbose)
        {
            std::cout << "Ion Model = Rogers-McCulloch" << std::endl << std::flush;
        }
        ionicModel.reset (new RogersMcCulloch< mesh_Type > (_dataIonic,
                                                            *localMeshPtr,
                                                            *uFESpacePtr,
                                                            *M_heart_fct->M_comm) );
    }
    else if (ion_model == 2)
    {
        if (verbose)
        {
            std::cout << "Ion Model = Luo-Rudy" << std::endl << std::flush;
        }
        ionicModel.reset (new LuoRudy< mesh_Type > (_dataIonic,
                                                    *localMeshPtr,
                                                    *uFESpacePtr,
                                                    *M_heart_fct->M_comm) );
    }
    else if (ion_model == 3)
    {
        if (verbose)
        {
            std::cout << "Ion Model = Mitchell-Schaeffer" << std::endl << std::flush;
        }
        ionicModel.reset (new MitchellSchaeffer< mesh_Type > (_dataIonic,
                                                              *localMeshPtr,
                                                              *uFESpacePtr,
                                                              *M_heart_fct->M_comm) );
    }

#ifdef MONODOMAIN
    electricModel.initialize ( M_heart_fct->initialScalar() );
#else
    electricModel.initialize ( M_heart_fct->initialScalar(),
                               M_heart_fct->zeroScalar() );
#endif

    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    ionicModel->initialize( );

    //! Building time-independent part of the system
    electricModel.buildSystem( );
    std::cout << "buildsystem ok" << std::endl;
    //! Initialization
    Real dt     = _data.timeStep();
    Real t0     = 0;
    Real tFinal = _data.endTime();
    MPI_Barrier (MPI_COMM_WORLD);

    if (verbose)
    {
        std::cout << "Setting the initial solution ... " << std::endl << std::endl;
    }
    _data.setTime (t0);
    electricModel.resetPreconditioner();
    if (verbose)
    {
        std::cout << " ok " << std::endl;
    }

    //! Setting generic Exporter postprocessing
    boost::shared_ptr< Exporter<mesh_Type > > exporter;
    std::string const exporterType =  M_heart_fct->M_dataFile ( "exporter/type", "ensight");
#ifdef HAVE_HDF5
    if (exporterType.compare ("hdf5") == 0)
    {
        exporter.reset ( new ExporterHDF5<mesh_Type > ( M_heart_fct->M_dataFile,
                                                        "heart" ) );
        exporter->setPostDir ( "./" ); // This is a test to see if M_post_dir is working
        exporter->setMeshProcId ( localMeshPtr, M_heart_fct->M_comm->MyPID() );
    }
    else
#endif
    {
        if (exporterType.compare ("none") == 0)
        {
            exporter.reset ( new ExporterEmpty<mesh_Type > ( M_heart_fct->M_dataFile,
                                                             localMeshPtr,
                                                             "heart",
                                                             M_heart_fct->M_comm->MyPID() ) );
        }
        else
        {
            exporter.reset ( new ExporterEnsight<mesh_Type > ( M_heart_fct->M_dataFile,
                                                               localMeshPtr,
                                                               "heart",
                                                               M_heart_fct->M_comm->MyPID() ) );
        }
    }


    vectorPtr_Type Uptr ( new vector_Type (electricModel.solutionTransmembranePotential(), Repeated ) );

    exporter->addVariable ( ExporterData<mesh_Type>::ScalarField,  "potential", uFESpacePtr,
                            Uptr, UInt (0) );

#ifdef BIDOMAIN
    vectorPtr_Type Ueptr ( new vector_Type (electricModel.solutionExtraPotential(), Repeated ) );
    exporter->addVariable ( ExporterData<mesh_Type>::ScalarField,  "potential_e", _FESpacePtr,
                            Ueptr, UInt (0) );
#endif

    vectorPtr_Type Fptr ( new vector_Type (electricModel.fiberVector(), Repeated ) );

    if (_data.hasFibers() )
        exporter->addVariable ( ExporterData<mesh_Type>::VectorField,
                                "fibers",
                                uFESpacePtr,
                                Fptr,
                                UInt (0) );
    exporter->postProcess ( 0 );

    MPI_Barrier (MPI_COMM_WORLD);
    chronoinitialsettings.stop();

    //! Temporal loop
    LifeChrono chrono;
    Int iter = 1;
    chronototaliterations.start();
    for ( Real time = t0 + dt ; time <= tFinal + dt / 2.; time += dt, iter++)
    {
        _data.setTime (time);
        if (verbose)
        {
            std::cout << std::endl;
            std::cout << "We are now at time " << _data.time() << " s. " << std::endl;
            std::cout << std::endl;
        }
        chrono.start();
        MPI_Barrier (MPI_COMM_WORLD);
        ionicModel->solveIonicModel ( electricModel.solutionTransmembranePotential(), _data.timeStep() );
        rhs *= 0;
        computeRhs ( rhs, electricModel, ionicModel, _data );
        electricModel.updatePDESystem ( rhs );
        electricModel.PDEiterate ( bcH );
        normu = electricModel.solutionTransmembranePotential().norm2();
        electricModel.solutionTransmembranePotential().epetraVector().MeanValue (&meanu);
        electricModel.solutionTransmembranePotential().epetraVector().MaxValue (&minu);
        if (verbose)
        {
            std::cout << "norm u " << normu << std::endl;
            std::cout << "mean u " << meanu << std::endl;
            std::cout << "max u " << minu << std::endl << std::flush;
        }

        *Uptr = electricModel.solutionTransmembranePotential();
#ifdef BIDOMAIN
        *Ueptr = electricModel.solutionExtraPotential();
#endif

        exporter->postProcess ( time );
        MPI_Barrier (MPI_COMM_WORLD);
        chrono.stop();
        if (verbose)
        {
            std::cout << "Total iteration time " << chrono.diff() << " s." << std::endl;
        }
        chronototaliterations.stop();
    }

    if (verbose)
    {
        std::cout << "Total iterations time " << chronototaliterations.diff() << " s." << std::endl;
    }
    if (verbose)
    {
        std::cout << "Total initial settings time " << chronoinitialsettings.diff() << " s." << std::endl;
    }
    if (verbose)
    {
        std::cout << "Total execution time " << chronoinitialsettings.diff() + chronototaliterations.diff() << " s." << std::endl;
    }
}

#ifdef MONODOMAIN
void Heart::computeRhs ( vector_Type& rhs,
                         HeartMonodomainSolver< mesh_Type >& electricModel,
                         boost::shared_ptr< HeartIonicSolver< mesh_Type > > ionicModel,
                         HeartMonodomainData& data )
{
    bool verbose = (M_heart_fct->M_comm->MyPID() == 0);
    if (verbose)
    {
        std::cout << "  f-  Computing Rhs ...        " << "\n" << std::flush;
    }
    LifeChrono chrono;
    chrono.start();

    //! u, w with repeated map
    vector_Type uVecRep (electricModel.solutionTransmembranePotential(), Repeated);
    ionicModel->updateRepeated();
    VectorElemental elvec_Iapp ( electricModel.potentialFESpace().fe().nbFEDof(), 2 ),
                    elvec_u ( electricModel.potentialFESpace().fe().nbFEDof(), 1 ),
                    elvec_Iion ( electricModel.potentialFESpace().fe().nbFEDof(), 1 );

    for (UInt iVol = 0; iVol < electricModel.potentialFESpace().mesh()->numVolumes(); ++iVol)
    {
        electricModel.potentialFESpace().fe().updateJacQuadPt ( electricModel.potentialFESpace().mesh()->volumeList ( iVol ) );
        elvec_Iapp.zero();
        elvec_u.zero();
        elvec_Iion.zero();

        UInt eleIDu = electricModel.potentialFESpace().fe().currentLocalId();
        UInt nbNode = ( UInt ) electricModel.potentialFESpace().fe().nbFEDof();

        //! Filling local elvec_u with potential values in the nodes
        for ( UInt iNode = 0 ; iNode < nbNode ; iNode++ )
        {
            Int  ig = electricModel.potentialFESpace().dof().localToGlobalMap ( eleIDu, iNode );
            elvec_u.vec() [ iNode ] = uVecRep[ig];
        }

        ionicModel->updateElementSolution (eleIDu);
        ionicModel->computeIonicCurrent (data.membraneCapacitance(), elvec_Iion, elvec_u, electricModel.potentialFESpace() );

        //! Computing the current source of the righthand side, repeated
        source (M_heart_fct->stimulus(),
                elvec_Iapp,
                electricModel.potentialFESpace().fe(),
                data.time(),
                0);
        source (M_heart_fct->stimulus(),
                elvec_Iapp,
                electricModel.potentialFESpace().fe(),
                data.time(),
                1);

        //! Assembling the righthand side
        for ( UInt i = 0 ; i < electricModel.potentialFESpace().fe().nbFEDof() ; i++ )
        {
            Int  ig = electricModel.potentialFESpace().dof().localToGlobalMap ( eleIDu, i );
            rhs.sumIntoGlobalValues (ig, (data.conductivityRatio() * elvec_Iapp.vec() [i] +
                                          elvec_Iapp.vec() [i + nbNode]) /
                                     (1 + data.conductivityRatio() ) + data.volumeSurfaceRatio() * elvec_Iion.vec() [i] );
        }
    }
    rhs.globalAssemble();
    Real coeff = data.volumeSurfaceRatio() * data.membraneCapacitance() / data.timeStep();
    vector_Type tmpvec (electricModel.solutionTransmembranePotential() );
    tmpvec *= coeff;
    rhs += electricModel.massMatrix() * tmpvec;
    MPI_Barrier (MPI_COMM_WORLD);
    chrono.stop();
    if (verbose)
    {
        std::cout << "done in " << chrono.diff() << " s." << std::endl;
    }
}
#else
void Heart::computeRhs ( vector_Type& rhs,
                         HeartBidomainSolver< mesh_Type >& electricModel,
                         boost::shared_ptr< HeartIonicSolver< mesh_Type > > ionicModel,
                         HeartBidomainData& data )
{
    bool verbose = (M_heart_fct->M_comm->MyPID() == 0);
    if (verbose)
    {
        std::cout << "  f-  Computing Rhs ...        " << "\n" << std::flush;
    }
    LifeChrono chrono;
    chrono.start();

    //! u, w with repeated map
    vector_Type uVecRep (electricModel.solutionTransmembranePotential(), Repeated);
    ionicModel->updateRepeated();

    VectorElemental elvec_Iapp ( electricModel.potentialFESpace().fe().nbFEDof(), 2 ),
                    elvec_u ( electricModel.potentialFESpace().fe().nbFEDof(), 1 ),
                    elvec_Iion ( electricModel.potentialFESpace().fe().nbFEDof(), 1 );
    for (UInt iVol = 0; iVol < electricModel.potentialFESpace().mesh()->numVolumes(); ++iVol)
    {
        electricModel.potentialFESpace().fe().updateJacQuadPt ( electricModel.potentialFESpace().mesh()->volumeList ( iVol ) );
        elvec_u.zero();
        elvec_Iion.zero();
        elvec_Iapp.zero();

        UInt eleIDu = electricModel.potentialFESpace().fe().currentLocalId();
        UInt nbNode = ( UInt ) electricModel.potentialFESpace().fe().nbFEDof();
        for ( UInt iNode = 0 ; iNode < nbNode ; iNode++ )
        {
            Int ig = electricModel.potentialFESpace().dof().localToGlobalMap ( eleIDu, iNode );
            elvec_u.vec() [ iNode ] = uVecRep[ig];
        }

        UInt eleID = electricModel.potentialFESpace().fe().currentLocalId();
        ionicModel->updateElementSolution (eleID);
        ionicModel->computeIonicCurrent (data.membraneCapacitance(), elvec_Iion, elvec_u, electricModel.potentialFESpace() );

        //! Computing Iapp
        source (M_heart_fct->stimulus(),
                elvec_Iapp,
                electricModel.potentialFESpace().fe(),
                data.time(), 0);
        source (M_heart_fct->stimulus(),
                elvec_Iapp,
                electricModel.potentialFESpace().fe(),
                data.time(),
                1);
        UInt totalUDof  = electricModel.potentialFESpace().map().map (Unique)->NumGlobalElements();

        for ( UInt iNode = 0 ; iNode < nbNode ; iNode++ )
        {
            Int ig = electricModel.potentialFESpace().dof().localToGlobalMap ( eleIDu, iNode );
            rhs.sumIntoGlobalValues (ig, elvec_Iapp.vec() [iNode] +
                                     data.volumeSurfaceRatio() * elvec_Iion.vec() [iNode] );
            rhs.sumIntoGlobalValues (ig + totalUDof,
                                     -elvec_Iapp.vec() [iNode + nbNode] -
                                     data.volumeSurfaceRatio() * elvec_Iion.vec() [iNode] );
        }
    }
    rhs.globalAssemble();

    rhs += electricModel.matrMass() * data.volumeSurfaceRatio() *
           data.membraneCapacitance() * electricModel.BDFIntraExtraPotential().time_der (data.timeStep() );

    MPI_Barrier (MPI_COMM_WORLD);

    chrono.stop();
    if (verbose)
    {
        std::cout << "done in " << chrono.diff() << " s." << std::endl;
    }
}
#endif
