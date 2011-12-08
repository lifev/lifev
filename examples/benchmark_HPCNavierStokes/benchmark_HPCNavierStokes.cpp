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
    @brief Benchmark to test the performance of LifeV for HPC

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 2011-06-09
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

#include <life/lifecore/LifeV.hpp>
#include <life/lifemesh/RegionMesh3DStructured.hpp>
#include <life/lifemesh/MeshData.hpp>
#include <life/lifemesh/RegionMesh3D.hpp>
#include <life/lifemesh/MeshPartitioner.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/BCManage.hpp>
#include <life/lifearray/MatrixEpetra.hpp>
#include <life/lifesolver/OseenAssembler.hpp>
#include <life/lifealg/PreconditionerIfpack.hpp>
#include <life/lifealg/PreconditionerML.hpp>
#include <life/lifealg/SolverAztecOO.hpp>
#include <life/lifefilters/ExporterHDF5.hpp>
#include <life/lifefem/TimeAdvanceBDF.hpp>

using namespace LifeV;

namespace
{

enum DiffusionType{ViscousStress, StiffStrain};

typedef RegionMesh3D<LinearTetra>         mesh_type;
typedef MatrixEpetra<Real>                matrix_type;
typedef VectorEpetra                      vector_type;
typedef boost::shared_ptr<VectorEpetra>   vectorPtr_type;
typedef FESpace< mesh_type, MapEpetra >   fespace_type;
typedef boost::shared_ptr< fespace_type > fespacePtr_type;

// +-----------------------------------------------+
// | Data and functions for the boundary conditions|
// +-----------------------------------------------+
const Int INFLOW   = 1;
const Int OUTFLOW1 = 2;
const Int OUTFLOW2 = 3;
const Int OUTFLOW3 = 4;
const Int OUTFLOW4 = 5;
const Int WALL     = 6;

Real fZero( const Real& /* t */,
                  const Real& /* x */,
                  const Real& /* y */,
                  const Real& /* z */,
                  const ID& /* i */ )
{
    return 0.0;
}

/*
   fluxes as fourier interpolation of physiological values, taken from
   Baek H, Jayaraman MV, Richardson PD, Karniadakis GE.
   Flow instability and wall shear stress variation in intracranial aneurysms.
   J R Soc Interface. 2010 Jun 6;7(47):967-88. Epub 2009 Dec 18.

*/
Real aneurismFluxIn(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    // We change the flux for our geometry
    const Real pi         = 3.141592653589793;
    const Real area       = 0.0907122;
    const Real areaFactor = area/((0.6/2)*(0.6/2)*pi);

    // Unit conversion from ml/min to cm^3/s
    const Real unitFactor = 1./ 60.;

    // T is the period of the cardiac cycle
    const Real T          = 1.0;

    // a0 is the average VFR (the value is taken from Karniadakis p970)
    const Real a0         = 255;

    // Fourrier
    const Int M(7);
    const Real a[M] = {-0.152001,-0.111619, 0.043304,0.028871,0.002098,-0.027237,-0.000557};
    const Real b[M] = { 0.129013,-0.031435,-0.086106,0.028263,0.010177, 0.012160,-0.026303};

    Real flux(0);
    const Real xi(2*pi*t/T);

    flux = a0;
    Int k(1);
    for (; k<=M ; ++k)
        flux += a0*(a[k-1]*cos(k*xi) + b[k-1]*sin(k*xi));

    return - flux * areaFactor * unitFactor;
}

/*
 * This function imposes a flat profile for the inflow.
 * It is not the best choice as a boundary condition.
 * However this is satisfactory enough for a benchmark
 */
Real aneurismFluxInVectorial(const Real&  t, const Real& x, const Real& y, const Real& z, const ID& i)
{
    Real n1(-0.000019768882940);
    Real n2(-0.978289616345544);
    Real n3( 0.207242433298975);
    Real flux(aneurismFluxIn(t,x,y,z,i));
    Real area(0.0907122); // Computed with the triangle on the INLET boundary

    switch(i) {
    case 0:
        return n1*flux/area;
    case 1:
        return n2*flux/area;
    case 2:
        return n3*flux/area;
    default:
        return 0.0;
    }
}

}


Int
main( Int argc, char** argv )
{
    // +-----------------------------------------------+
    // |            Initialization of MPI              |
    // +-----------------------------------------------+
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm(new Epetra_MpiComm(MPI_COMM_WORLD));
    Int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#else
    boost::shared_ptr<Epetra_Comm> Comm(new Epetra_SerialComm);
#endif

    const bool verbose(Comm->MyPID()==0);
    if(verbose){
        std::cout
            << " +-----------------------------------------------+" << std::endl
            << " |          HPC Navier-Stokes benchmark          |" << std::endl
            << " +-----------------------------------------------+" << std::endl
            << std::endl
            << " +-----------------------------------------------+" << std::endl
            << " |           Author: Gwenol Grandperrin          |" << std::endl
            << " |             Date: 2011-06-28                  |" << std::endl
            << " +-----------------------------------------------+" << std::endl
            << std::endl;

        std::cout << "[Initilization of MPI]" << std::endl;
#ifdef HAVE_MPI
        std::cout << "Using MPI (" << nproc << " proc.)" << std::endl;
#else
        std::cout << "Using serial version" << std::endl;
#endif
    }

    // +-----------------------------------------------+
    // |               Loading the data                |
    // +-----------------------------------------------+
    if (verbose) std::cout << std::endl << "[Loading the data]" << std::endl;
    LifeChrono globalChrono;
    LifeChrono initChrono;
    LifeChrono iterChrono;

    globalChrono.start();
    initChrono.start();

    // ******* GetPot stuff ********
    GetPot command_line(argc,argv);
    const std::string dataFileName = command_line.follow("data", 2, "-f","--file");
    GetPot dataFile(dataFileName);
    // *****************************

    // Physical quantity
    const Real viscosity       = dataFile("problem/viscosity",0.035);
    const Real density         = dataFile("problem/density",1.0);

    // Time discretization
    const Real initialTime     = dataFile("problem/initial_time",0.0);
    const Real endTime         = dataFile("problem/end_time",1e-2);
    const Real timestep        = dataFile("problem/timestep",1e-3);

    if (verbose) std::cout << "Viscosity    : " << viscosity << std::endl;
    if (verbose) std::cout << "Density      : " << density << std::endl;
    if (verbose) std::cout << "Initial time : " << initialTime << std::endl;
    if (verbose) std::cout << "End time     : " << endTime << std::endl;
    if (verbose) std::cout << "Timestep     : " << timestep << std::endl;

    // Finite element
    std::string uOrder("P2");
    std::string pOrder("P1");

    // Numerical scheme
    const DiffusionType diffusionType = ViscousStress;
    UInt BDFOrder = 1;

    // Preconditioner
    std::string precName       = dataFile("prec/prectype","none");

    const bool reusePreconditioner = false;

    // +-----------------------------------------------+
    // |               Loading the mesh                |
    // +-----------------------------------------------+
    if (verbose) std::cout << std::endl << "[Loading the mesh]" << std::endl;

    boost::shared_ptr<RegionMesh3D<LinearTetra> > fullMeshPtr(new RegionMesh3D<LinearTetra>);

    // Building the mesh from the source
    MeshData meshData;
    meshData.setup(dataFile, "problem");
    if (verbose) std::cout << "Mesh source: file("
                           << meshData.meshDir() << meshData.meshFile() << ")" << std::endl;
    readMesh(*fullMeshPtr, meshData);

    if (verbose)
        std::cout << "Mesh size  : " <<
        MeshUtility::MeshStatistics::computeSize(*fullMeshPtr).maxH << std::endl;
    if (verbose) std::cout << "Partitioning the mesh ... " << std::endl;
    MeshPartitioner< RegionMesh3D<LinearTetra> >   meshPart(fullMeshPtr, Comm);
    fullMeshPtr.reset(); //Freeing the global mesh to save memory

    // +-----------------------------------------------+
    // |            Creating the FE spaces             |
    // +-----------------------------------------------+
    if (verbose) std::cout << std::endl << "[Creating the FE spaces]" << std::endl;

    if (verbose) std::cout << "FE for the velocity: " << uOrder << std::endl
                           << "FE for the pressure: " << pOrder << std::endl;

    if (verbose) std::cout << "Building the velocity FE space ... " << std::flush;
    fespacePtr_type uFESpace( new FESpace< mesh_type, MapEpetra >(meshPart,uOrder, nDimensions, Comm));
    if (verbose) std::cout << "ok." << std::endl;

    if (verbose) std::cout << "Building the pressure FE space ... " << std::flush;
    fespacePtr_type pFESpace( new FESpace< mesh_type, MapEpetra >(meshPart,pOrder, 1, Comm));
    if (verbose) std::cout << "ok." << std::endl;

    // Creation of the total map
    MapEpetra solutionMap(uFESpace->map()+pFESpace->map());

    // Pressure offset in the vector
    UInt pressureOffset = nDimensions * uFESpace->dof().numTotalDof();

    if (verbose) std::cout << "Total Velocity Dof: " << pressureOffset << std::endl;
    if (verbose) std::cout << "Total Pressure Dof: " << pFESpace->dof().numTotalDof() << std::endl;

    // +-----------------------------------------------+
    // |             Boundary conditions               |
    // +-----------------------------------------------+
    if (verbose) std::cout << std::endl << "[Boundary conditions]" << std::endl;
    boost::shared_ptr<BCHandler> bcHandler;
    bcHandler.reset(new BCHandler);
    BCFunctionBase uZero(fZero);
    BCFunctionBase uFlux(aneurismFluxInVectorial);

    if (verbose) std::cout << "Setting BC... " << std::flush;
    bcHandler->addBC( "Inflow"  , INFLOW  , Essential, Full, uFlux, 3 );
    bcHandler->addBC( "Outflow1", OUTFLOW1, Natural  , Full, uZero, 3 );
    bcHandler->addBC( "Outflow2", OUTFLOW2, Natural  , Full, uZero, 3 );
    bcHandler->addBC( "Outflow3", OUTFLOW3, Natural  , Full, uZero, 3 );
    bcHandler->addBC( "Outflow4", OUTFLOW4, Natural  , Full, uZero, 3 );
    bcHandler->addBC( "Wall"    , WALL    , Essential, Full, uZero, 3 );
    if (verbose) std::cout << "ok." << std::endl;

    // Update the BCHandler (internal data related to FE)
    bcHandler->bcUpdate( *meshPart.meshPartition(), uFESpace->feBd(), uFESpace->dof());

    // +-----------------------------------------------+
    // |              Matrices Assembly                |
    // +-----------------------------------------------+
    if (verbose) std::cout << std::endl << "[Matrices Assembly]" << std::endl;

    if (verbose) std::cout << "Setting up assembler... " << std::flush;
    OseenAssembler<mesh_type,matrix_type,vector_type> oseenAssembler;
    oseenAssembler.setup(uFESpace,pFESpace);
    if (verbose) std::cout << "done" << std::endl;

    if (verbose) std::cout << "Defining the matrices... " << std::flush;
    boost::shared_ptr<matrix_type> systemMatrix(new matrix_type( solutionMap ));
    *systemMatrix *=0.0;
    boost::shared_ptr<matrix_type> baseMatrix(new matrix_type( solutionMap ));
    *baseMatrix *=0.0;
    boost::shared_ptr<matrix_type> massMatrix(new matrix_type( solutionMap ));
    *massMatrix *=0.0;
    if (verbose) std::cout << "done" << std::endl;

    // Perform the assembly of the matrix
    switch(diffusionType)
    {
        case ViscousStress:
            if (verbose) std::cout << "Adding the viscous stress... " << std::flush;
            oseenAssembler.addViscousStress(*baseMatrix,viscosity/density);
            if (verbose) std::cout << "done" << std::endl;
            break;
        case StiffStrain:
            if (verbose) std::cout << "Adding the stiff strain... " << std::flush;
            oseenAssembler.addStiffStrain(*baseMatrix,viscosity/density);
            if (verbose) std::cout << "done" << std::endl;
            break;
        default:
            cerr << "[Error] Diffusion type unknown" << std::endl;
            exit(1);
            break;
    }

    if (verbose) std::cout << "Adding the gradient of the pressure... " << std::flush;
    oseenAssembler.addGradPressure(*baseMatrix);
    if (verbose) std::cout << "done" << std::endl;

    if (verbose) std::cout << "Adding the divergence free constraint... " << std::flush;
    oseenAssembler.addDivergence(*baseMatrix,-1.0);
    if (verbose) std::cout << "done" << std::endl;

    if (verbose) std::cout << "Adding the mass... " << std::flush;
    oseenAssembler.addMass(*massMatrix,1.0);
    if (verbose) std::cout << "done" << std::endl;

    if (verbose) std::cout << "Closing the matrices... " << std::flush;
    baseMatrix->globalAssemble();
    massMatrix->globalAssemble();
    if (verbose) std::cout << "done" << std::endl;

    // +-----------------------------------------------+
    // |            Solver initialization              |
    // +-----------------------------------------------+
    if (verbose) std::cout << std::endl << "[Solver initialization]" << std::endl;
    SolverAztecOO linearSolver;

    if (verbose) std::cout << "Setting up the solver... " << std::flush;
    linearSolver.setDataFromGetPot(dataFile,"solver");

    // Creating the preconditioner
    linearSolver.setupPreconditioner(dataFile,"prec");

    if (verbose) std::cout << "done" << std::endl;

    linearSolver.setCommunicator(Comm);

    // +-----------------------------------------------+
    // |       Initialization of the simulation        |
    // +-----------------------------------------------+
    if (verbose) std::cout<< std::endl << "[Initialization of the simulation]" << std::endl;
    if (verbose) std::cout << "Creation of vectors... " << std::flush;
    vectorPtr_type rhs;
    rhs.reset(new vector_type(solutionMap,Unique));

    vectorPtr_type beta;
    beta.reset(new vector_type(solutionMap,Repeated));

    vectorPtr_type velocity;
    velocity.reset(new vector_type(uFESpace->map(),Unique));

    vectorPtr_type pressure;
    pressure.reset(new vector_type(pFESpace->map(),Unique));

    vectorPtr_type solution;
    solution.reset(new vector_type(solutionMap,Unique));
    if (verbose) std::cout << "done" << std::endl;

    if (verbose) std::cout << "Computing the initial solution ... " << std::endl;

    // bdf object to store the previous solutions
    TimeAdvanceBDF<vector_type> bdf;
    bdf.setup(BDFOrder);
    Real currentTime = initialTime-timestep*BDFOrder;

    *velocity *= 0;
    *pressure *= 0;
    *solution *= 0;
    bdf.setInitialCondition( *solution );

    // Initial solution (obtained from a Stokes problem)
    linearSolver.setTolerance(1e-6);
    currentTime += timestep;
    for ( ; currentTime <=  initialTime + timestep/2.; currentTime += timestep)
    {

        if (verbose) std::cout << "Updating the system... " << std::flush;
        *rhs  *= 0.;
        rhs->globalAssemble();
        systemMatrix.reset(new matrix_type( solutionMap ));
        *systemMatrix += *baseMatrix;
        if (verbose) std::cout << "done" << std::endl;

        if (verbose) std::cout << "Applying BC... " << std::flush;
        bcManage(*systemMatrix,*rhs,*uFESpace->mesh(),uFESpace->dof(),*bcHandler,uFESpace->feBd(),1.0,currentTime);
        systemMatrix->globalAssemble();
        if (verbose) std::cout << "done" << std::endl;

        if (verbose) std::cout << "Solving the system... " << std::endl;
        *solution *= 0;
        linearSolver.setMatrix(*systemMatrix);
        linearSolver.solveSystem(*rhs,*solution,systemMatrix);

        // Updating bdf
        bdf.shiftRight( *solution );

    }

    linearSolver.resetPreconditioner();

    // +-----------------------------------------------+
    // |             Setting the exporter              |
    // +-----------------------------------------------+
    if (verbose) std::cout << "Defining the exporter... " << std::flush;
    ExporterHDF5<mesh_type> exporter ( dataFile, "OseenAssembler");
    exporter.setPostDir( "./" ); // This is a test to see if M_post_dir is working
    exporter.setMeshProcId( meshPart.meshPartition(), Comm->MyPID() );
    if (verbose) std::cout << "done" << std::endl;

    if (verbose) std::cout << "Updating the exporter... " << std::flush;
    exporter.addVariable( ExporterData<mesh_type>::VectorField, "velocity", uFESpace,
                          solution, UInt(0));
    exporter.addVariable( ExporterData<mesh_type>::ScalarField, "pressure", pFESpace,
                          solution, pressureOffset );
    if (verbose) std::cout << "done" << std::endl;

    if (verbose) std::cout << "Exporting solution at time t=" << initialTime << "... " << std::endl;
    exporter.postProcess(initialTime);

    initChrono.stop();
    if (verbose) std::cout << "Initialization time: " << initChrono.diff() << " s." << std::endl;

    // Reset the preconditioner
    linearSolver.resetPreconditioner();
    linearSolver.setReusePreconditioner(reusePreconditioner);
    linearSolver.setTolerance(1e-10);

    // +-----------------------------------------------+
    // |             Solving the problem               |
    // +-----------------------------------------------+
    if (verbose) std::cout<< std::endl << "[Solving the problem]" << std::endl;
    Int iter = 1;

    for ( ; currentTime <= endTime + timestep/2.; currentTime += timestep, iter++)
    {
        iterChrono.reset();
        iterChrono.start();

        if (verbose) std::cout << std::endl << "[t = "<< currentTime << " s.]" << std::endl;

        if (verbose) std::cout << "Updating the system... " << std::flush;
        bdf.updateRHSContribution( timestep );
        *rhs  = *massMatrix*bdf.rhsContributionFirstDerivative();

        systemMatrix.reset(new matrix_type( solutionMap ));
        Real alpha = bdf.coefficientFirstDerivative( 0 ) / timestep;
        *systemMatrix += *massMatrix*alpha;
        *systemMatrix += *baseMatrix;

        // SemiImplicit
        *beta = bdf.extrapolation(); // Extrapolation for the convective term
        oseenAssembler.addConvection(*systemMatrix,*beta);

        if (verbose) std::cout << "done" << std::endl;

        if (verbose) std::cout << "Applying BC... " << std::flush;
        bcManage(*systemMatrix,*rhs,*uFESpace->mesh(),uFESpace->dof(),*bcHandler,uFESpace->feBd(),1.0,currentTime);
        systemMatrix->globalAssemble();
        if (verbose) std::cout << "done" << std::endl;

        if (verbose) std::cout << "Solving the system... " << std::endl;
        *solution *= 0;
        linearSolver.setMatrix(*systemMatrix);
        linearSolver.solveSystem(*rhs,*solution,systemMatrix);

        // Updating the BDF scheme
        bdf.shiftRight( *solution );

        // Exporting the solution
        exporter.postProcess( currentTime );

        iterChrono.stop();
        if (verbose) std::cout << "Iteration time: " << iterChrono.diff() << " s." << std::endl;

        MPI_Barrier(MPI_COMM_WORLD);
    }

    // +-----------------------------------------------+
    // |            Ending the simulation              |
    // +-----------------------------------------------+
    exporter.closeFile();

    globalChrono.stop();
    if (verbose) std::cout << std::endl << "Total simulation time: " << globalChrono.diff() << " s." << std::endl;
    if (verbose) std::cout << std::endl << "[[END_SIMULATION]]" << std::endl;

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return( EXIT_SUCCESS );
}

