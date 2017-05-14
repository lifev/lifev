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
    @brief

    @author Claudia Colciago <claudia.colciago@epfl.ch>
    @date 08-11-2012
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


#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/navier_stokes/fem/TimeAdvanceBDFNavierStokes.hpp>

#include <lifev/core/fem/BCManage.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEnsight.hpp>

#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/navier_stokes/solver/OseenData.hpp>
#include <lifev/navier_stokes/solver/StabilizationIP.hpp>


#include "ETRobinMembraneSolver.hpp"
#include "ud_functions.hpp"

#include <iostream>

//#define FLUX 1

using namespace LifeV;


const Real PI = 3.141592653589793;


struct ETRobinMembraneSolver::Private
{
    typedef std::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > fct_type;

    std::string    data_file_name;

    // Data for the fluid
    Real         Re;
    Real         nu;         /**< viscosity (in m^2/s) */
    Real         H;          /**< height and width of the domain (in m) */
    Real         D;          /**< diameter of the cylinder (in m) */
    Real         density;
    Real         R;          //radius


    //Data for the membrnae
    Real rhos;
    Real Hs;
    Real ni;
    Real E;

    //Discretization choices
    UInt transpirationOrder;
    std::string stabilization;
    std::string initial_sol;
    Real numLM;
    bool useFlowRate;

    std::shared_ptr<Epetra_Comm>   comm;

};

ETRobinMembraneSolver::ETRobinMembraneSolver ( int argc, char** argv )
    :
    M_d ( new Private )
{
    GetPot command_line (argc, argv);
    std::string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile ( data_file_name );

    M_d->data_file_name  = data_file_name;

    M_d->Re              = dataFile ( "fluid/physics/Re", 1. );
    M_d->nu              = dataFile ( "fluid/physics/viscosity", 1. );
    M_d->H               = dataFile ( "fluid/physics/H", 20. );
    M_d->D               = dataFile ( "fluid/physics/D", 1. );
    M_d->R               = M_d->D / 2;
    M_d->density         = dataFile ( "fluid/physics/density", 1. );

    M_d->rhos            = dataFile ( "membrane/physics/density_sol", 1. );
    M_d->E               = dataFile ( "membrane/physics/young", 1. );
    M_d->ni              = dataFile ( "membrane/physics/poisson", 1. );
    M_d->Hs              = dataFile ( "membrane/physics/wall_thickness", 1. );


    M_d->initial_sol = (std::string) dataFile ( "fluid/problem/initial_sol", "none");
    M_d->stabilization = (std::string) dataFile ( "fluid/problem/stabilization", "none");
    M_d->numLM       = dataFile ( "fluid/problem/numLM"      , 0    );
    M_d->useFlowRate = 0;
    if ( M_d->numLM > 0 )
    {
        M_d->useFlowRate = 1 ;
    }
    M_d->transpirationOrder       = dataFile ( "fluid/problem/transpiration_order"      , 0     );

#ifdef EPETRA_MPI
    std::cout << "mpi initialization ... " << std::endl;

    int ntasks = 0;
    M_d->comm.reset ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
    if (!M_d->comm->MyPID() )
    {
        std::cout << "My PID            = " << M_d->comm->MyPID() << " out of " << ntasks << " running." << std::endl;
        std::cout << "fluid density     = " << M_d->density       << std::endl
                  << "nu                = " << M_d->nu            << std::endl
                  << "poisson           = " << M_d->ni            << std::endl
                  << "initial solution  = " << M_d->initial_sol   << std::endl
                  << "Young modulus     = " << M_d->E             << std::endl
                  << "Wall thickness    = " << M_d->Hs            << std::endl;

    }

#else
    M_d->comm.reset ( new Epetra_SerialComm() );
#endif

}

void
ETRobinMembraneSolver::run()

{
    bool verbose = (M_d->comm->MyPID() == 0);

    //------------Creating file to store data------------------------------

    GetPot dataFile ( M_d->data_file_name );

    // NS
    Real NSSupg (dataFile ("fluid/space_discretization/supg_coeff", 1e-3) ); // 1e-3
    Real NSPspg (dataFile ("fluid/space_discretization/pspg_coeff", 1e-3) ); // 1e-3
    Real NSdivdiv (dataFile ("fluid/space_discretization/dividiv_coeff", 5e-2) );
    Real OUTLET (dataFile ("fluid/problem/boundary_flags/outlet", 3) );
    Real INLET (dataFile ("fluid/problem/boundary_flags/inlet", 2) );
    Real WALL (dataFile ("fluid/problem/boundary_flags/wall", 100) );
    Real RING (dataFile ("fluid/problem/boundary_flags/ring_in", 20) );
    Real RING2 (dataFile ("fluid/problem/boundary_flags/ring_out", 30) );

    std::shared_ptr<OseenData> oseenData (new OseenData() );
    oseenData->setup ( dataFile );

    MeshData meshData;
    meshData.setup (dataFile, "fluid/space_discretization");

    Real LameI = (/*M_d->Hs * */M_d->E * M_d->ni ) / ( ( 1 - M_d->ni * M_d->ni ) );
    Real LameII = /*M_d->Hs * */M_d->E / ( 2 * ( 1 + M_d->ni ) );

    if (verbose)
    {
        std::cout << " Outlet : " << OUTLET << std::endl;
        std::cout << " Wall : " << WALL << std::endl;
    }
    if (verbose)
    {
        std::cout << " LameI : " << LameI << std::endl;
        std::cout << " LameII : " << LameII << std::endl;
    }

    //-----------------------------------------------------------------------

    //----------Creating Mesh and FESpaces-----------------------------------

    std::shared_ptr< mesh_type > fullMeshPtr (new mesh_type);
    readMesh (*fullMeshPtr, meshData);

    MeshPartitioner< mesh_type >   meshPart (fullMeshPtr, M_d->comm);

    std::string uOrder =  dataFile ( "fluid/space_discretization/vel_order", "P1");
    std::string pOrder =  dataFile ( "fluid/space_discretization/press_order", "P1");

    if (verbose)
    {
        std::cout << "Building the FE spaces ... " << std::flush;
    }

    M_uFESpace.reset ( new FESpace< mesh_type, MapEpetra > ( meshPart, uOrder, 3,  M_d->comm ) );
    M_uCompFESpace.reset ( new FESpace< mesh_type, MapEpetra > ( meshPart, uOrder, 1,  M_d->comm ) );
    M_pFESpace.reset ( new FESpace< mesh_type, MapEpetra > ( meshPart, pOrder, 1,  M_d->comm ) );

    M_ETuFESpace.reset ( new ETFESpace< mesh_type, MapEpetra, 3, 3 > ( meshPart, & (M_uFESpace->refFE() ), M_d->comm ) );
    M_ETpFESpace.reset ( new ETFESpace< mesh_type, MapEpetra, 3, 1 > ( meshPart, & (M_pFESpace->refFE() ), M_d->comm ) );

    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    UInt totalVelDof   = M_uFESpace->map().map ( Unique )->NumGlobalElements();
    UInt totalPressDof = M_pFESpace->map().map ( Unique )->NumGlobalElements();

    if (verbose)
    {
        std::cout << "Total Velocity DOF = " << totalVelDof << std::endl;
    }
    if (verbose)
    {
        std::cout << "Total Pressure DOF = " << totalPressDof << std::endl;
    }
    if (verbose)
    {
        std::cout << "Total FLux DOF = " << M_d->numLM << std::endl;
    }

    MapEpetra fluidMap ( M_uFESpace->map() );
    MapEpetra fluxMap ( M_d->numLM, M_d->comm );
    MapEpetra fullMap ( M_uFESpace->map() + M_pFESpace->map() + fluxMap );

    DOFInterface3Dto3D interfaceDOF ( M_uFESpace->refFE() , M_uFESpace->dof() );
    interfaceDOF.update ( *M_uFESpace->mesh() , WALL , *M_uFESpace->mesh() , WALL , 0);
    createInterfaceMap ( interfaceDOF.localDofMap() , M_uFESpace->dof() );

    //------------------FLUID PROBLEM--------------------------------

    if (verbose)
    {
        std::cout << "Calling the solver constructors ... ";
    }

    SolverAztecOO NSSolver;
    NSSolver.setCommunicator (M_d->comm);
    NSSolver.setDataFromGetPot (dataFile, "solver");
    NSSolver.setupPreconditioner (dataFile, "prec");

    if (verbose)
    {
        std::cout << "done." << std::endl;
    }

    //---------------------------------------------------------

    //--------------------Creating BDF objects---------------------------

    if (verbose)
    {
        std::cout << "Calling the TimeAdvance constructors ... ";
    }

    Real dt     = oseenData->dataTime()->timeStep();
    Real t0     = oseenData->dataTime()->initialTime();
    Real tFinal = oseenData->dataTime()->endTime();

    TimeAdvanceBDFNavierStokes<vector_type> fluidTimeAdvance;
    fluidTimeAdvance.setup (oseenData->dataTimeAdvance()->orderBDF() );

    TimeAdvanceBDF<vector_type> dispTimeAdvance;
    dispTimeAdvance.setup (oseenData->dataTimeAdvance()->orderBDF() );

    //TimeAdvance for the fluid
    fluidTimeAdvance.bdfVelocity().setTimeStep (oseenData->dataTime()->timeStep() );
    if (verbose)
    {
        fluidTimeAdvance.showMe();
    }

    //TimeAdvance for the displacement
    dispTimeAdvance.setTimeStep (oseenData->dataTime()->timeStep() );
    if (verbose)
    {
        dispTimeAdvance.showMe();
    }

    Real alpha = fluidTimeAdvance.bdfVelocity().coefficientFirstDerivative (t0);

    if (verbose)
    {
        std::cout << " done. \n";
    }

    //---------------------------------------------------------

    //--------------Boundary Conditions-------------------------

    if (verbose)
    {
        std::cout << "Calling the BCHandler constructor ... ";
    }

    //creating BCHandler objects
    BCHandler bcHFluid;
    BCFunctionBase vel_in (flatNormalVelInlet /*linearVelInletCylinder*/);
    BCFunctionBase flux_in (linearInletCylinder);
    BCFunctionBase uZero ( fZero );
    BCFunctionBase press_out ( p0 );

    if (M_d->useFlowRate == true )
    {
        bcHFluid.addBC ("InFlow" ,   INLET,    Flux,                Normal, flux_in    );
        bcHFluid.addBC ("InFlow_2" , OUTLET,   Natural,             Full,   uZero,   3 );

        //In order to have a well defined problem you need to set conditions on boundaries of boundaries
        bcHFluid.addBC ("Ring" ,     RING,     EssentialEdges,   Full,   uZero,   3 );
        bcHFluid.addBC ("Ring3" ,    RING2,    EssentialEdges,   Full,   uZero,   3 );

        bcHFluid.setOffset ("InFlow", totalVelDof + totalPressDof);

    }
    else
    {

        bcHFluid.addBC ("InFlow" ,   INLET,    Essential,           Full,   vel_in, 3 );
        bcHFluid.addBC ("InFlow_2" , OUTLET,   Natural,             Normal,   press_out );

        //In order to have a well defined problem you need to set conditions on boundaries of boundaries
        bcHFluid.addBC ("Ring" ,     RING,     EssentialEdges,   Full,   uZero,   3 );
        bcHFluid.addBC ("Ring3" ,    RING2,    EssentialEdges,   Full,   uZero,   3 );

    }

    if (verbose)
    {
        std::cout << " BC done \n";
    }

    //----------------------------------------------------------------------------

    vector_block_type NSSolution (M_ETuFESpace->map() | M_ETpFESpace->map() | fluxMap, Unique);
    vector_type velocitySolution (M_ETuFESpace->map(), Repeated);
    vector_type dispSolution (M_ETuFESpace->map(), Repeated);

    NSSolution *= 0.0;

    //-------------Creating Exporter Object------------------------------------------------

    if (verbose)
    {
        std::cout << "Calling the Exporter constructor ... ";
    }

    std::string const exporterFileName    =  dataFile ( "exporter/filename", "cylinder");
    UInt  exportEach = dataFile ("exporter/each", 1);

    ExporterHDF5< mesh_type > exporter ( dataFile, meshPart.meshPartition(), exporterFileName, M_d->comm->MyPID() );

    vectorPtr_type velAndPressureExporter ( new vector_type (NSSolution, exporter.mapType() ) );
    vectorPtr_type dispExporter (new vector_type (dispSolution, Repeated ) );

    exporter.addVariable ( ExporterData< mesh_type >::VectorField, "f-velocity",
                           M_uFESpace, velAndPressureExporter, UInt (0) );

    exporter.addVariable ( ExporterData< mesh_type >::ScalarField, "f-pressure",
                           M_pFESpace, velAndPressureExporter, UInt (3 * M_uFESpace->dof().numTotalDof() ) );

    exporter.addVariable ( ExporterData< mesh_type >::VectorField, "s-displacement",
                           M_uFESpace, dispExporter, UInt (0) ) ;
    if (verbose)
    {
        std::cout << " done. \n";
    }

    //-----------------------------------------------------------------------------

    //--------------------------Initialization--------------------------------------



    if (M_d->initial_sol == "restart")
    {
        ASSERT (oseenData->dataTimeAdvance()->orderBDF() == 1, "Restart works only for BDF1");

        if (verbose)
        {
            std::cout << std::endl;
        }
        if (verbose)
        {
            std::cout << "Restoring the previous solution ... " << std::endl;
        }

        std::string filename = dataFile ("importer/filename", "cylinder");

        LifeV::ExporterHDF5<mesh_type> importer ( dataFile, filename );
        importer.setMeshProcId ( M_uFESpace->mesh(), M_d->comm->MyPID() );

        vectorPtr_type velAndPressureImporter ( new vector_type (NSSolution, importer.mapType() ) );
        vectorPtr_type dispImporter (new vector_type (dispSolution, importer.mapType() ) );

        importer.addVariable ( ExporterData<mesh_type>::VectorField,
                               "f-velocity",
                               M_uFESpace,
                               velAndPressureImporter,
                               UInt ( 0 ) );

        importer.addVariable ( ExporterData<mesh_type>::ScalarField,
                               "f-pressure",
                               M_pFESpace,
                               velAndPressureImporter,
                               3 * M_uFESpace->dof().numTotalDof() );

        importer.addVariable ( ExporterData<mesh_type>::VectorField,
                               "s-displacement",
                               M_uFESpace,
                               dispImporter,
                               UInt ( 0 ) );

        exporter.setTimeIndex ( importer.importFromTime (t0) );

        Real norm = velAndPressureImporter->norm2();
        *velAndPressureExporter = *velAndPressureImporter;
        *dispExporter = *dispImporter;
        if (verbose)
        {
            std::cout << "   f- restart solution norm = " << norm << std::endl;
        }


    }

    velocitySolution.subset ( *velAndPressureExporter );
    dispSolution = *dispExporter;

    fluidTimeAdvance.bdfVelocity().setInitialCondition ( * (velAndPressureExporter) );
    dispTimeAdvance.setInitialCondition ( dispSolution );

    fluidTimeAdvance.bdfVelocity().updateRHSContribution ( oseenData->dataTime()->timeStep() );
    dispTimeAdvance.updateRHSContribution (oseenData->dataTime()->timeStep() );

    if (verbose)
    {
        std::cout << " done \n";
    }


    //-------------------------------------------------------------------------------------

    std::shared_ptr<matrix_block_type> NSMatrixConstant (new matrix_block_type ( M_ETuFESpace->map() | M_ETpFESpace->map() | fluxMap ) );
    *NSMatrixConstant *= 0.0;
    std::shared_ptr<matrix_block_type> NSMatrixSteadyStokes (new matrix_block_type ( M_ETuFESpace->map() | M_ETpFESpace->map() | fluxMap ) );
    *NSMatrixSteadyStokes *= 0.0;

    //----------------------Temporal Loop----------------------------------

    if (verbose)
    {
        std::cout << std::endl;
    }
    if (verbose)
    {
        std::cout << " ### Simulation times ### " << std::endl;
    }
    if (verbose)
    {
        std::cout << " From " << t0 << " to " << tFinal << std::endl;
    }
    if (verbose)
    {
        std::cout << " Time step: " << dt << std::endl;
    }
    if (verbose)
    {
        std::cout << std::endl;
    }

    Real currentTime (t0);
    UInt niter (0);

    *NSMatrixConstant *= 0.0;

    QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria4pt) );
    vector_type robinExt ( M_pFESpace->map(), Unique);
    M_pFESpace->interpolate ( uZero , robinExt, 0);
    vector_type hSolid (M_pFESpace->map(), Unique);
    hSolid += M_d->Hs; // * robinExt;
    //robinExt *= 0.0;

    {

        using namespace ExpressionAssembly;

        integrate (
            elements (M_ETuFESpace->mesh() ), // Mesh

            M_uFESpace->qr(), // QR

            M_ETuFESpace,
            M_ETuFESpace,

            // Viscous Term
            M_d->nu * dot (grad (phi_i) , grad (phi_j) )

        )
                >> NSMatrixSteadyStokes->block (0, 0);



        integrate (
            elements (M_ETuFESpace->mesh() ), // Mesh

            M_uFESpace->qr(), // QR

            M_ETuFESpace,
            M_ETpFESpace,

            value (-1.0) *phi_j * div (phi_i)


        )
                >> NSMatrixSteadyStokes->block (0, 1);

        integrate (
            elements (M_ETuFESpace->mesh() ), // Mesh

            M_uFESpace->qr(), // QR

            M_ETpFESpace,
            M_ETuFESpace,

            value (1.0) *phi_i * div (phi_j)

        )
                >> NSMatrixSteadyStokes->block (1, 0);

        NSMatrixSteadyStokes->globalAssemble();

        integrate (
            elements (M_ETuFESpace->mesh() ), // Mesh

            M_uFESpace->qr(), // QR

            M_ETuFESpace,
            M_ETuFESpace,

            // Intertial Term
            M_d->density
            * ( value (alpha / dt) * dot (phi_i, phi_j) )

        )
                >> NSMatrixConstant->block (0, 0);


        integrate ( boundary (M_ETuFESpace->mesh(), WALL),
                    myBDQR,

                    M_ETuFESpace,
                    M_ETuFESpace,

                    //Boundary Mass
                    ( value (M_d->rhos /*M_d->Hs*/  * alpha / dt) * value ( M_ETpFESpace, hSolid ) +  value ( M_ETpFESpace, robinExt ) * ( 0.1 + dt ) ) * dot ( phi_j, phi_i )

                  )
                >> NSMatrixConstant->block (0, 0);

    }


    MatrixSmall<3, 3> Eye;
    Eye *= 0.0;
    Eye[0][0] = 1;
    Eye[1][1] = 1;
    Eye[2][2] = 1;

    {
        using namespace ::LifeV::ExpressionAssembly;


        integrate ( boundary (M_ETuFESpace->mesh(), WALL),
                    myBDQR,

                    M_ETuFESpace,
                    M_ETuFESpace,

                    //Boundary Stiffness
                    ( dt / alpha ) *
                    2  * LameII * value ( M_ETpFESpace, hSolid )   *
                    0.5 * dot ( ( grad (phi_j) + (-1) * grad (phi_j) * outerProduct ( Nface, Nface ) )
                                + transpose (grad (phi_j) + (-1) * grad (phi_j) * outerProduct ( Nface, Nface ) ),
                                ( grad (phi_i) + ( (-1) * grad (phi_i) * outerProduct ( Nface, Nface ) ) ) ) +

                    ( dt / alpha ) * value ( M_ETpFESpace, hSolid ) *
                    LameI  * dot ( value ( Eye ) , ( grad (phi_j) + (-1) *  grad (phi_j) * outerProduct ( Nface, Nface ) ) )
                    * dot ( value ( Eye ) ,  ( grad (phi_i) + (-1) * grad (phi_i) * outerProduct ( Nface, Nface ) )  )

                  )
                >> NSMatrixConstant->block (0, 0);

    }

    *NSMatrixConstant += *NSMatrixSteadyStokes;
    NSMatrixConstant->globalAssemble();

    details::StabilizationIP<mesh_type, DOF> M_ipStabilization;
    M_ipStabilization.setFeSpaceVelocity ( *M_uFESpace );
    M_ipStabilization.setViscosity ( M_d->nu );

    // Parameters from J. Michalik
    M_ipStabilization.setGammaBeta ( dataFile ( "ipstab/gammaBeta", 1.0) );
    M_ipStabilization.setGammaDiv  ( dataFile ( "ipstab/gammaDiv", 0.2) );
    M_ipStabilization.setGammaPress ( dataFile ( "ipstab/gammaPress", 0.5) );

    vector_type steadyResidual ( fullMap, Unique );

    //     if (M_d->initial_sol == "steady")
    //     {

    // #ifdef FLUX
    //         std::shared_ptr<matrix_block_type> NSMatrix (new matrix_block_type ( M_ETuFESpace->map() | M_ETpFESpace->map() | fluxMap ) );
    //         *NSMatrix *= 0.0;

    // #else

    //         std::shared_ptr<matrix_block_type> NSMatrix (new matrix_block_type ( M_ETuFESpace->map() | M_ETpFESpace->map() ) );
    //         *NSMatrix *= 0.0;
    // #endif

    //     vector_type NSRhsUnique ( fullMap, Unique );

    //     std::shared_ptr<matrix_type> stabMatrix (new matrix_type ( fullMap ));
    //     M_ipStabilization.apply ( *stabMatrix, velocitySolution, false );
    //     stabMatrix->globalAssemble();
    //     *NSMatrix += *NSMatrixSteadyStokes;
    //     *NSMatrix += *stabMatrix;

    //     NSMatrix->globalAssemble();

    //         bcHFluid.bcUpdate ( *meshPart.meshPartition(), M_uFESpace->feBd(), M_uFESpace->dof() );
    //        // bcManage (*NSMatrix, NSRhsUnique,
    //        //            *M_uFESpace->mesh(), M_uFESpace->dof(),
    //        //            bcHFluid, M_uFESpace->feBd(), 1.0, currentTime);

    //     std::shared_ptr<matrix_type> NSMatrixDiri( new matrix_type( fullMap ) );
    //     *NSMatrixDiri += *NSMatrix;
    //     NSMatrixDiri->globalAssemble();

    //     BCHandler bcHInit;
    //     bcHInit.addBC ("Wall" ,   WALL,    Essential,           Full,   uZero, 3 );
    //     bcHInit.bcUpdate ( *meshPart.meshPartition(), M_uFESpace->feBd(), M_uFESpace->dof() );
    //         // bcManage (*NSMatrixDiri, NSRhsUnique,
    //         //           *M_uFESpace->mesh(), M_uFESpace->dof(),
    //         //           bcHInit, M_uFESpace->feBd(), 1.0, currentTime);

    //         NSSolver.setMatrix (*NSMatrixDiri);
    //         NSSolver.solveSystem (NSRhsUnique, NSSolution, NSMatrixDiri);

    //     *velAndPressureExporter = NSSolution;

    //     steadyResidual = (-1) * (*NSMatrix * NSSolution);

    //     vector_type dispRhsInit( M_ETuFESpace->map(), Unique );
    //     dispRhsInit.subset( steadyResidual );

    //     vector_type dispSolutionInit( M_ETuFESpace->map(), Unique );
    //         std::shared_ptr<matrix_type> DispMatrixInit (new matrix_type ( M_ETuFESpace->map() ) );
    //         *DispMatrixInit *= 0.0;
    //      DispMatrixInit->insertOneDiagonal();
    //      DispMatrixInit->insertValueDiagonal( -1 , (*M_interfaceMap) );

    //      {
    //          using namespace ::LifeV::ExpressionAssembly;


    //          integrate ( boundary (M_ETuFESpace->mesh(), WALL),
    //              myBDQR,

    //              M_ETuFESpace,
    //              M_ETuFESpace,

    //              //Boundary Stiffness
    //              //Boundary Mass
    //              value( M_ETpFESpace, robinExt ) * dot ( phi_j, phi_i )

    //              + 2  * LameII  *
    //              0.5 * dot ( ( grad (phi_j) + (-1) * grad (phi_j) * outerProduct ( Nface, Nface ) )
    //                      + transpose (grad (phi_j) + (-1) * grad (phi_j) * outerProduct ( Nface, Nface ) ),
    //                      ( grad (phi_i) + ( (-1) * grad (phi_i) * outerProduct ( Nface, Nface ) ) ) ) +

    //              LameI  * dot ( value ( Eye ) , ( grad (phi_j) + (-1) *  grad (phi_j) * outerProduct ( Nface, Nface ) ) )
    //              * dot ( value ( Eye ) ,  ( grad (phi_i) + (-1) * grad (phi_i) * outerProduct ( Nface, Nface ) )  )

    //              )
    //                 >> *DispMatrixInit;

    //      }

    //     DispMatrixInit->globalAssemble();

    //     DispMatrixInit->spy("dispMatrix");

    //     BCHandler bcHDiri;
    //     bcHDiri.addBC ("Inlet" ,   INLET,    Essential,           Full,   uZero, 3 );
    //     bcHDiri.addBC ("Outlet" ,   OUTLET,    Essential,           Full,   uZero, 3 );
    //     bcHDiri.addBC ("Ring" ,   RING2,    EssentialEdges,           Full,   uZero, 3 );
    //     bcHDiri.addBC ("Ring" ,   RING,    EssentialEdges,           Full,   uZero, 3 );
    //     bcHDiri.bcUpdate ( *meshPart.meshPartition(), M_uFESpace->feBd(), M_uFESpace->dof() );
    //        // bcManage (*DispMatrixInit, dispRhsInit,
    //        //            *M_uFESpace->mesh(), M_uFESpace->dof(),
    //        //            bcHDiri, M_uFESpace->feBd(), 1.0, currentTime);

    //     NSSolver.setMatrix (*DispMatrixInit);
    //     NSSolver.solveSystem (dispRhsInit, dispSolutionInit , DispMatrixInit);

    //     *dispExporter = dispSolutionInit;

    //     dispSolution = dispSolutionInit;

    //     dispTimeAdvance.shiftRight ( dispSolution );

    //        dispTimeAdvance.updateRHSContribution (oseenData->dataTime()->timeStep() );

    //     velocitySolution.subset (NSSolution);

    //     fluidTimeAdvance.bdfVelocity().shiftRight ( NSSolution );

    //     fluidTimeAdvance.bdfVelocity().updateRHSContribution ( oseenData->dataTime()->timeStep() );
    //     }

    exporter.postProcess (t0);


    while ( currentTime < tFinal)
    {
        LifeChrono ChronoIteration;
        ChronoIteration.start();

        currentTime += dt;
        niter += 1;

        if (verbose)
        {
            std::cout << std::endl;
        }
        if (verbose)
        {
            std::cout << "----------------------------" << std::endl;
        }
        if (verbose)
        {
            std::cout << " Time : " << currentTime << std::endl;
        }
        if (verbose)
        {
            std::cout << " Iter : " << niter << std::endl;
        }
        if (verbose)
        {
            std::cout << "----------------------------" << std::endl;
        }
        if (verbose)
        {
            std::cout << std::endl;
        }

        vector_type velocityExtrapolated (velocitySolution, Repeated);
        fluidTimeAdvance.bdfVelocity().extrapolation ( velocityExtrapolated );
        vector_type velocityBdfRHS (velocitySolution, Repeated);
        velocityBdfRHS = fluidTimeAdvance.bdfVelocity().rhsContributionFirstDerivative();

        vector_type dispExtrapolated (dispSolution, Repeated);
        dispTimeAdvance.extrapolation ( dispExtrapolated );
        vector_type dispBdfRHS (fluidMap, Repeated);
        dispBdfRHS = dispTimeAdvance.rhsContributionFirstDerivative();
        vector_type dtDisp (fluidMap, Repeated);
        dtDisp = ( alpha / dt ) * dispSolution - dispTimeAdvance.rhsContributionFirstDerivative();

        std::shared_ptr<matrix_block_type> NSMatrix (new matrix_block_type ( M_ETuFESpace->map() | M_ETpFESpace->map() | fluxMap ) );
        *NSMatrix *= 0.0;

#define DIVDIV_TEST value(NSdivdiv) * h_K * div(phi_i)

#define SUPG_TEST value(NSSupg) * h_K * (grad(phi_i)*eval(normalize,value(M_ETuFESpace,velocitySolution)))

#define PSPG_TEST value(NSPspg) * h_K * grad(phi_i)

#define TRANSP_GRADGRAD grad(phi_j) * ( grad(M_ETuFESpace, dispExtrapolated) ) * ( Eye + (-1) * outerProduct(Nface, Nface) )

        if (verbose)
        {
            std::cout << "[Navier-Stokes] Assembling the matrix ... " << std::flush;
        }

        LifeChrono ChronoItem;
        ChronoItem.start();

        {
            std::shared_ptr<NormalizeFct> normalize (new NormalizeFct);
            using namespace ::LifeV::ExpressionAssembly;

            integrate (
                elements (M_ETuFESpace->mesh() ), // Mesh
                M_uFESpace->qr(), // QR

                M_ETuFESpace,
                M_ETuFESpace,

                // Advection Term
                M_d->density
                * ( dot (grad (phi_j) * value (M_ETuFESpace, velocityExtrapolated) , phi_i) )

            )
                    >> NSMatrix->block (0, 0);


            if (M_d->stabilization == "sup")
            {

                integrate (
                    elements (M_ETuFESpace->mesh() ), // Mesh
                    M_uFESpace->qr(), // QR

                    M_ETuFESpace,
                    M_ETuFESpace,

                    // SUPG
                    dot (
                        grad (phi_j) *value (M_ETuFESpace, velocityExtrapolated) * M_d->density
                        + value (alpha / dt) * M_d->density * phi_j
                        , SUPG_TEST)

                    // Div div
                    + DIVDIV_TEST
                    *div (phi_j)


                )
                        >> NSMatrix->block (0, 0);


                integrate (
                    elements (M_ETuFESpace->mesh() ), // Mesh
                    M_uFESpace->qr(), // QR

                    M_ETuFESpace,
                    M_ETpFESpace,

                    // SUPG
                    dot ( grad (phi_j) , SUPG_TEST)

                )
                        >> NSMatrix->block (0, 1);


                integrate (
                    elements (M_ETuFESpace->mesh() ), // Mesh
                    M_uFESpace->qr(), // QR

                    M_ETpFESpace,
                    M_ETuFESpace,

                    // PSPG
                    dot (
                        grad (phi_j) *value (M_ETuFESpace, velocityExtrapolated) * value ( M_d->density ) + value (1.0 / dt) * value ( M_d->density ) * phi_j
                        , PSPG_TEST)
                )
                        >> NSMatrix->block (1, 0);


                integrate (
                    elements (M_ETuFESpace->mesh() ), // Mesh

                    M_uFESpace->qr(), // QR

                    M_ETpFESpace,
                    M_ETpFESpace,

                    // PSPG
                    dot (
                        grad (phi_j)
                        , PSPG_TEST)

                )
                        >> NSMatrix->block (1, 1);
            }

        }

        if ( M_d->transpirationOrder == 1 )
        {

            {
                using namespace ::LifeV::ExpressionAssembly;


                integrate ( boundary (M_ETuFESpace->mesh(), WALL),
                            myBDQR,

                            M_ETuFESpace,
                            M_ETuFESpace,

                            //Boundary Stiffness
                            ( dt / alpha ) *
                            2  * LameII * value ( M_ETpFESpace, hSolid ) *
                            0.5 * dot ( TRANSP_GRADGRAD
                                        + transpose ( TRANSP_GRADGRAD ),
                                        ( grad (phi_i) + ( (-1) * grad (phi_i) * outerProduct ( Nface, Nface ) ) ) ) +

                            ( dt / alpha ) *
                            LameI * value ( M_ETpFESpace, hSolid ) * dot ( value ( Eye ) , TRANSP_GRADGRAD )
                            * dot ( value ( Eye ) ,  ( grad (phi_i) + (-1) * grad (phi_i) * outerProduct ( Nface, Nface ) )  )

                          )
                        >> NSMatrix->block (0, 0);


                integrate ( boundary (M_ETuFESpace->mesh(), WALL),
                            myBDQR,

                            M_ETuFESpace,
                            M_ETuFESpace,

                            //Boundary Mass
                            value (M_d->rhos* /*M_d->Hs*/ alpha / dt) * value ( M_ETpFESpace, hSolid) * dot ( grad (phi_j) * value ( M_ETuFESpace, dispExtrapolated) , phi_i )
                            + dot ( grad (phi_j) * value ( M_ETuFESpace, dtDisp ) , phi_i )
                            + value ( M_ETpFESpace, robinExt ) * ( 0.1 + dt / alpha )  * dot ( grad (phi_j) * value ( M_ETuFESpace, dispExtrapolated) , phi_i )

                          )
                        >> NSMatrix->block (0, 0);



            }

        }


        if (M_d->stabilization == "ip")
        {

            std::shared_ptr<matrix_type> stabMatrix (new matrix_type ( fullMap ) );
            M_ipStabilization.apply ( *stabMatrix, velocityExtrapolated, false );
            stabMatrix->globalAssemble();
            *NSMatrix += *stabMatrix;
        }



        ChronoItem.stop();
        if (verbose)
        {
            std::cout << ChronoItem.diff() << " s" << std::endl;
        }


        if (verbose)
        {
            std::cout << "[Navier-Stokes] Adding constant parts ... " << std::flush;
        }
        ChronoItem.start();

        *NSMatrix += *NSMatrixConstant;

        ChronoItem.stop();
        if (verbose)
        {
            std::cout << ChronoItem.diff() << " s" << std::endl;
        }


        if (verbose)
        {
            std::cout << "[Navier-Stokes] Assembling the rhs ... " << std::flush;
        }
        ChronoItem.start();

        vector_block_type NSRhs ( M_ETuFESpace->map() | M_ETpFESpace->map() | fluxMap, Repeated );
        NSRhs *= 0.0;

        {
            std::shared_ptr<NormalizeFct> normalize (new NormalizeFct);

            using namespace ExpressionAssembly;

            integrate (
                elements (M_ETuFESpace->mesh() ), // Mesh

                M_uFESpace->qr(), // QR

                M_ETuFESpace,

                //Inertial Term
                M_d->density * dot (value (M_ETuFESpace, velocityBdfRHS), phi_i )

            )
                    >> NSRhs.block (0);

            if (M_d->stabilization == "sup")
            {

                integrate (
                    elements (M_ETuFESpace->mesh() ), // Mesh

                    M_uFESpace->qr(), // QR

                    M_ETuFESpace,

                    // SUPG
                    dot (
                        M_d->density * value (M_ETuFESpace, velocityBdfRHS)
                        , SUPG_TEST)

                )
                        >> NSRhs.block (0);

                integrate (
                    elements (M_ETuFESpace->mesh() ), // Mesh

                    M_uFESpace->qr(), // QR

                    M_ETpFESpace,


                    // PSPG
                    dot (
                        M_d->density  * value (M_ETuFESpace, velocityBdfRHS)
                        , PSPG_TEST)


                )
                        >> NSRhs.block (1);
            }

        }

        ChronoItem.stop();
        if (verbose)
        {
            std::cout << ChronoItem.diff() << " s" << std::endl;
        }

        if (verbose)
        {
            std::cout << "[Navier-Stokes] Boundary Intergrals in the rhs ... " << std::flush;
        }
        ChronoItem.start();

        {
            using namespace ExpressionAssembly;

            integrate ( boundary (M_ETuFESpace->mesh(), WALL),
                        myBDQR,

                        M_ETuFESpace,

                        //BOUNDARY STIFFNESS
                        value (-1) * ( dt / alpha ) *
                        2  * LameII * value ( M_ETpFESpace, hSolid ) *
                        0.5 * dot ( ( grad (M_ETuFESpace, dispBdfRHS ) + (-1) * grad (M_ETuFESpace, dispBdfRHS ) * outerProduct ( Nface, Nface ) )
                                    + transpose (grad (M_ETuFESpace, dispBdfRHS  ) + (-1) * grad (M_ETuFESpace, dispBdfRHS  ) * outerProduct ( Nface, Nface ) ),
                                    ( grad (phi_i) + ( (-1) * grad (phi_i) * outerProduct ( Nface, Nface ) ) ) ) +

                        value (-1) * ( dt / alpha ) *
                        LameI * value ( M_ETpFESpace, hSolid ) * dot ( value ( Eye ) , ( grad (M_ETuFESpace, dispBdfRHS ) + (-1) *  grad (M_ETuFESpace, dispBdfRHS ) * outerProduct ( Nface, Nface ) ) )
                        * dot ( value ( Eye ) ,  ( grad (phi_i) + (-1) * grad (phi_i) * outerProduct ( Nface, Nface ) )  )
                      ) >> NSRhs.block (0);

            integrate (
                boundary (M_ETuFESpace->mesh(), WALL), // Mesh

                myBDQR, // QR

                M_ETuFESpace,

                //BOUNDARY MASS
                ( /*M_d->Hs*/value ( M_ETpFESpace, hSolid ) * M_d->rhos * alpha) * dot (value (M_ETuFESpace, velocityBdfRHS), phi_i )
                - dt * dot ( value ( M_ETpFESpace, robinExt) * value (M_ETuFESpace, dispBdfRHS ) , phi_i )

            )
                    >> NSRhs.block (0);

        }

        if ( M_d->transpirationOrder == 1 )
        {
            {
                using namespace ExpressionAssembly;

                integrate (
                    boundary (M_ETuFESpace->mesh(), WALL), // Mesh

                    myBDQR, // QR

                    M_ETuFESpace,

                    //BOUNDARY MASS
                    ( /*M_d->Hs*/value ( M_ETpFESpace, hSolid) * M_d->rhos * alpha) *  dot (grad (M_ETuFESpace, velocityBdfRHS) * value ( M_ETuFESpace, dispExtrapolated ), phi_i )


                )
                        >> NSRhs.block (0);
            }
        }

        ChronoItem.stop();
        if (verbose)
        {
            std::cout << ChronoItem.diff() << " s" << std::endl;
        }

        if (verbose)
        {
            std::cout << "[Navier-Stokes] Closing the matrix and the rhs ... " << std::flush;
        }

        ChronoItem.start();

        NSMatrix->globalAssemble();

        NSRhs.globalAssemble();
        vector_block_type NSRhsUnique ( NSRhs, Unique );
        //NSRhsUnique += (-1) * steadyResidual;
        //steadyResidual *= 0.0;

        ChronoItem.stop();
        if (verbose)
        {
            std::cout << ChronoItem.diff() << " s" << std::endl;
        }


        if (verbose)
        {
            std::cout << "[Navier-Stokes] Applying boundary conditions ... " << std::flush;
        }

        bcHFluid.bcUpdate ( *meshPart.meshPartition(), M_uFESpace->feBd(), M_uFESpace->dof() );

        std::shared_ptr<matrix_type> NSMatrixNoBlock (new matrix_type ( fullMap ) );
        *NSMatrixNoBlock += *NSMatrix;
        NSMatrixNoBlock->globalAssemble();
        bcManage (*NSMatrixNoBlock, NSRhsUnique,
                  *M_uFESpace->mesh(), M_uFESpace->dof(),
                  bcHFluid, M_uFESpace->feBd(), 1.0, currentTime);

        ChronoItem.stop();
        if (verbose)
        {
            std::cout << ChronoItem.diff() << " s" << std::endl;
        }

        if (verbose)
        {
            std::cout << "[Navier-Stokes] Solving the system " << std::endl;
        }

        NSSolver.setMatrix (*NSMatrixNoBlock);

        NSSolver.solveSystem (NSRhsUnique, NSSolution, NSMatrixNoBlock);

        fluidTimeAdvance.bdfVelocity().shiftRight ( NSSolution );

        velocitySolution.subset (NSSolution);

        vector_type dispInterface (M_interfaceMap, Repeated);

        dispInterface = dt / alpha * ( velocitySolution + dispBdfRHS );

        dispSolution *= 0.0;

        dispSolution.subset (dispInterface, *M_interfaceMap, 0, 0);

        dispTimeAdvance.shiftRight ( dispSolution );

        fluidTimeAdvance.bdfVelocity().updateRHSContribution ( oseenData->dataTime()->timeStep() );
        dispTimeAdvance.updateRHSContribution (oseenData->dataTime()->timeStep() );

        *velAndPressureExporter = NSSolution;
        *dispExporter = dispSolution;

        //dispExporter->spy ("dispExporter");

        if (niter % exportEach == 0)
        {
            exporter.postProcess (currentTime);
        }

        ChronoIteration.stop();
        if (verbose)
        {
            std::cout << std::endl << " Total iteration time : " << ChronoIteration.diff() << " s" << std::endl;
        }

    } // end time loop


    exporter.closeFile();

}


void ETRobinMembraneSolver::createInterfaceMap ( std::map<ID, ID> const& locDofMap, const DOF& dof )
{
    Displayer disp (M_d->comm);
    disp.leaderPrint ("Building the Interface Map ...             ");

    std::vector<int> dofInterfaceFluid;

    typedef std::map<ID, ID>::const_iterator iterator_Type;

    //std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
    dofInterfaceFluid.reserve ( locDofMap.size() );

    for (UInt dim = 0; dim < nDimensions; ++dim)
        for ( iterator_Type i = locDofMap.begin(); i != locDofMap.end(); ++i )
        {
            dofInterfaceFluid.push_back (i->second + dim * dof.numTotalDof() );    // in solid numerotation
        }

    int* pointerToDofs (0);
    if (dofInterfaceFluid.size() > 0)
    {
        pointerToDofs = &dofInterfaceFluid[0];
    }

    M_interfaceMap.reset ( new MapEpetra ( -1,
                                           static_cast<int> (dofInterfaceFluid.size() ),
                                           pointerToDofs,
                                           M_d->comm ) );
    disp.leaderPrint ("done\n");
    M_d->comm->Barrier();

}
