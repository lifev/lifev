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

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 25-03-2011
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
#include <life/lifefem/TimeAdvanceBDFNavierStokes.hpp>
#include <life/lifefem/TimeAdvanceBDF.hpp>

#include "EthierSteinmanUnsteady.hpp"

using namespace LifeV;

namespace
{
static bool regIF = (PRECFactory::instance().registerProduct( "Ifpack", &createIfpack ));
static bool regML = (PRECFactory::instance().registerProduct( "ML", &createML ));

enum DiffusionType{ViscousStress, StiffStrain};
enum MeshType{RegularMesh, File};
enum InitType{Interpolation, Projection};
enum ConvectionType{Explicit, SemiImplicit, KIO91};

typedef RegionMesh3D<LinearTetra> mesh_type;
typedef MatrixEpetra<Real> matrix_type;
typedef VectorEpetra vector_type;
typedef boost::shared_ptr<VectorEpetra> vectorPtr_type;
typedef FESpace< mesh_type, MapEpetra > fespace_type;
typedef boost::shared_ptr< fespace_type > fespacePtr_type;
}

void printErrors(const vector_type& solution, const Real& currentTime, fespacePtr_type uFESpace, fespacePtr_type pFESpace, bool verbose)
{
    vector_type velocity(uFESpace->map(),Repeated);
    vector_type pressure(pFESpace->map(),Repeated);
    if (verbose) std::cout << std::endl << "[Computed errors]" << std::endl;
    velocity.subset(solution);
    pressure.subset(solution, 3 * uFESpace->dof().numTotalDof());
    Real uRelativeError, pRelativeError, uL2Error, pL2Error;
    uL2Error = uFESpace->l2Error (EthierSteinmanUnsteady::uexact, velocity, currentTime, &uRelativeError );
    pL2Error = pFESpace->l20Error(EthierSteinmanUnsteady::pexact, pressure, currentTime, &pRelativeError );
    if (verbose) std::cout << "Velocity" << std::endl;
    if (verbose) std::cout << "  L2 error      : " << uL2Error << std::endl;
    if (verbose) std::cout << "  Relative error: " << uRelativeError << std::endl;
    if (verbose) std::cout << "Pressure" << std::endl;
    if (verbose) std::cout << "  L2 error      : " << pL2Error << std::endl;
    if (verbose) std::cout << "  Relative error: " << pRelativeError << std::endl;
}


int
main( int argc, char** argv )
{
    // +-----------------------------------------------+
    // |            Initialization of MPI              |
    // +-----------------------------------------------+
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm(new Epetra_MpiComm(MPI_COMM_WORLD));
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#else
    boost::shared_ptr<Epetra_Comm> Comm(new Epetra_SerialComm);
#endif

    const bool verbose(Comm->MyPID()==0);
    if(verbose){
        std::cout
            << " +-----------------------------------------------+" << std::endl
            << " |            OseenAssembler example             |" << std::endl
            << " +-----------------------------------------------+" << std::endl
            << std::endl
            << " +-----------------------------------------------+" << std::endl
            << " |           Author: Gwenol Grandperrin          |" << std::endl
            << " |             Date: 2010-03-25                  |" << std::endl
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

    // **** Stupid GetPot stuff ****
    GetPot command_line(argc,argv);
    const std::string dataFileName = command_line.follow("data", 2, "-f","--file");
    GetPot dataFile(dataFileName);
    // *****************************

    // Physical quantity
    const Real viscosity      = 0.01;
    const Real density        = 1.0;

    // Time discretization
    const Real initialTime    = 0.0;
    const Real endTime        = 1e-2;
    const Real timestep       = 1e-3;

    // Space discretization
    const UInt numDimensions  = 3;
    const MeshType meshSource = RegularMesh;
    const UInt numMeshElem    = 10;

    // Numerical scheme
    const DiffusionType diffusionType = ViscousStress;
          UInt BDFOrder = 3;
    const InitType initializationMethod = Interpolation;
    const ConvectionType convectionTerm = KIO91;

    // EthierSteinman data
    EthierSteinmanUnsteady::setA(1.0);
    EthierSteinmanUnsteady::setD(1.0);
    EthierSteinmanUnsteady::setViscosity(viscosity);
    EthierSteinmanUnsteady::setDensity(density);

    // +-----------------------------------------------+
    // |               Loading the mesh                |
    // +-----------------------------------------------+
    if (verbose) std::cout << std::endl << "[Loading the mesh]" << std::endl;

    boost::shared_ptr<RegionMesh3D<LinearTetra> > fullMeshPtr(new RegionMesh3D<LinearTetra>);

    // Building the mesh from the source
    if(meshSource == RegularMesh)
    {
        regularMesh3D( *fullMeshPtr,
                       1,
                       numMeshElem, numMeshElem, numMeshElem,
                       false,
                       2.0,   2.0,   2.0,
                       -1.0,  -1.0,  -1.0);

        if (verbose) std::cout << "Mesh source: regular mesh("
                               << numMeshElem << "x" << numMeshElem << "x" << numMeshElem << ")" << std::endl;
    }
    else if(meshSource == File)
    {
        MeshData meshData;
        meshData.setup(dataFile, "fluid/space_discretization");
        readMesh(*fullMeshPtr, meshData);

        if (verbose) std::cout << "Mesh source: file("
                               << meshData.meshDir() << meshData.meshFile() << ")" << std::endl;
    }
    else
    {
        if (verbose) std::cout << std::endl << "Error: Unknown source type for the mesh" << std::endl;
        exit(1);
    }

    if (verbose) std::cout << "Mesh size  : " << fullMeshPtr->maxH() << std::endl;
    if (verbose) std::cout << "Partitioning the mesh ... " << std::endl;
    MeshPartitioner< RegionMesh3D<LinearTetra> >   meshPart(fullMeshPtr, Comm);
    fullMeshPtr.reset(); //Freeing the global mesh to save memory

    // +-----------------------------------------------+
    // |            Creating the FE spaces             |
    // +-----------------------------------------------+
    if (verbose) std::cout << std::endl << "[Creating the FE spaces]" << std::endl;
    std::string uOrder("P2");
    std::string pOrder("P1");

    if (verbose) std::cout << "FE for the velocity: " << uOrder << std::endl
                           << "FE for the pressure: " << pOrder << std::endl;

    if (verbose) std::cout << "Building the velocity FE space ... " << std::flush;
    fespacePtr_type uFESpace( new FESpace< mesh_type, MapEpetra >(meshPart,uOrder, numDimensions, Comm));
    if (verbose) std::cout << "ok." << std::endl;

    if (verbose) std::cout << "Building the pressure FE space ... " << std::flush;
    fespacePtr_type pFESpace( new FESpace< mesh_type, MapEpetra >(meshPart,pOrder, 1, Comm));
    if (verbose) std::cout << "ok." << std::endl;

    // Creation of the total map
    MapEpetra solutionMap(uFESpace->map()+pFESpace->map());

    // Pressure offset in the vector
    UInt pressureOffset = numDimensions * uFESpace->dof().numTotalDof();

    if (verbose) std::cout << "Total Velocity Dof: " << pressureOffset << std::endl;
    if (verbose) std::cout << "Total Pressure Dof: " << pFESpace->dof().numTotalDof() << std::endl;

    // +-----------------------------------------------+
    // |             Boundary conditions               |
    // +-----------------------------------------------+
    if (verbose) std::cout << std::endl << "[Boundary conditions]" << std::endl;
    BCHandler bcHandler;
    BCFunctionBase uDirichlet( EthierSteinmanUnsteady::uexact );
    BCFunctionBase uNeumann( EthierSteinmanUnsteady::fNeumann );

    if (verbose) std::cout << "Setting Neumann BC... " << std::flush;
    bcHandler.addBC( "Flux", 1, Natural, Full, uNeumann, 3 );
    if (verbose) std::cout << "ok." << std::endl;

    if (verbose) std::cout << "Setting Dirichlet BC... " << std::flush;
    for (UInt iDirichlet(2);iDirichlet<=26;++iDirichlet)
    {
        bcHandler.addBC( "Wall", iDirichlet, Essential, Full, uDirichlet, 3 );
    }
    if (verbose) std::cout << "ok." << std::endl;

    // Update the BCHandler (internal data related to FE)
    bcHandler.bcUpdate( *meshPart.meshPartition(), uFESpace->feBd(), uFESpace->dof());

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
            oseenAssembler.addViscousStress(baseMatrix,viscosity/density);
            if (verbose) std::cout << "done" << std::endl;
            break;
        case StiffStrain:
            if (verbose) std::cout << "Adding the stiff strain... " << std::flush;
            oseenAssembler.addStiffStrain(baseMatrix,viscosity/density);
            if (verbose) std::cout << "done" << std::endl;
            break;
        default:
            cerr << "[Error] Diffusion type unknown" << std::endl;
            exit(1);
            break;
    }

    if (verbose) std::cout << "Adding the gradient of the pressure... " << std::flush;
    oseenAssembler.addGradPressure(baseMatrix);
    if (verbose) std::cout << "done" << std::endl;

    if (verbose) std::cout << "Adding the divergence free constraint... " << std::flush;
    oseenAssembler.addDivergence(baseMatrix,-1.0);
    if (verbose) std::cout << "done" << std::endl;

    if (verbose) std::cout << "Adding the mass... " << std::flush;
    oseenAssembler.addMass(massMatrix,1.0);
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
    if(convectionTerm == KIO91) BDFOrder = 3;

    // bdf object to store the previous solutions
    TimeAdvanceBDFNavierStokes<vector_type> bdf;
    bdf.setup(BDFOrder);
    TimeAdvanceBDF<vector_type> bdfConvection;
    bdfConvection.setup(BDFOrder);
    TimeAdvanceBDF<vector_type> bdfConvectionInit; // Just for KIO91
    bdfConvectionInit.setup(BDFOrder);
    Real currentTime = initialTime-timestep*BDFOrder;

    if(convectionTerm == KIO91)
    {
        uFESpace->interpolate(EthierSteinmanUnsteady::uexact,*velocity,currentTime);
        *solution *= 0;
        *solution = *velocity;
        *beta *= 0;
        oseenAssembler.addConvectionRhs(*beta,*solution);
        bdfConvection.setInitialCondition( *beta );
        bdfConvectionInit.setInitialCondition( *beta );

        if(initializationMethod == Projection)
        {
            for(UInt i(0);i<BDFOrder;++i)
            {
                uFESpace->interpolate(EthierSteinmanUnsteady::uexact,*velocity,currentTime-(3-i)*timestep);
                *solution = *velocity;
                *beta *= 0;
                oseenAssembler.addConvectionRhs(*beta,*solution);
                bdfConvectionInit.shiftRight( *beta );
            }
        }
    }

    uFESpace->interpolate(EthierSteinmanUnsteady::uexact,*velocity,currentTime);
    pFESpace->interpolate(EthierSteinmanUnsteady::pexact,*pressure,currentTime);
    *solution *= 0;
    *solution = *velocity;
    solution->add(*pressure,pressureOffset);
    bdf.bdfVelocity().setInitialCondition( *solution );

    // Initial solution (interpolation or projection)
    currentTime += timestep;
    for ( ; currentTime <=  initialTime + timestep/2.; currentTime += timestep)
    {
        *rhs  *= 0.;
        *beta *= 0.;
        *solution *= 0;

        uFESpace->interpolate(EthierSteinmanUnsteady::uexact,*velocity,currentTime);
        pFESpace->interpolate(EthierSteinmanUnsteady::pexact,*pressure,currentTime);
        *solution = *velocity;
        solution->add(*pressure,pressureOffset);

        if (initializationMethod == Projection)
        {
            uFESpace->interpolate(EthierSteinmanUnsteady::uderexact, *rhs, currentTime);
            //uFESpace.l2ScalarProduct(EthierSteinmanUnsteady::uderexact, rhs, currentTime);
            rhs->globalAssemble();
            *rhs *= -1.;
            *rhs = (*massMatrix)*(*rhs);

            if (verbose) std::cout << "Updating the system... " << std::flush;
            systemMatrix.reset(new matrix_type( solutionMap ));
            *systemMatrix += *baseMatrix;
            if(convectionTerm == SemiImplicit)
            {
                oseenAssembler.addConvection(systemMatrix,*solution);
            }
            else if(convectionTerm == Explicit)
            {
                oseenAssembler.addConvectionRhs(*rhs,*solution);
            }
            else if(convectionTerm == KIO91)
            {
                vector_type tmp(uFESpace->map(),Unique); // we do not want the pressure part
                *rhs -= bdfConvectionInit.extrapolation();
                *beta *= 0;
                oseenAssembler.addConvectionRhs(*beta,*solution);
                bdfConvectionInit.shiftRight(*beta);
            }

            if (verbose) std::cout << "done" << std::endl;

            if (verbose) std::cout << "Applying BC... " << std::flush;
            bcManage(*systemMatrix,*rhs,*uFESpace->mesh(),uFESpace->dof(),bcHandler,uFESpace->feBd(),1.0,currentTime);
            systemMatrix->globalAssemble();
            if (verbose) std::cout << "done" << std::endl;

            if (verbose) std::cout << "Solving the system... " << std::endl;
            *solution *= 0;
            linearSolver.setMatrix(*systemMatrix);
            linearSolver.solveSystem(*rhs,*solution,systemMatrix);
        }

        // Updating bdf
        bdf.bdfVelocity().shiftRight( *solution );

        if(convectionTerm == KIO91)
        {
            *beta *= 0;
            oseenAssembler.addConvectionRhs(*beta,*solution);
            bdfConvection.shiftRight(*beta);
        }

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
    exporter.addVariable( ExporterData::Vector, "velocity", solution,
                          UInt(0), uFESpace->dof().numTotalDof() );
    exporter.addVariable( ExporterData::Scalar, "pressure", solution,
                          pressureOffset,
                          UInt(pFESpace->dof().numTotalDof()) );
    if (verbose) std::cout << "done" << std::endl;

    if (verbose) std::cout << "Exporting solution at time t=" << initialTime << "... " << std::endl;
    exporter.postProcess(initialTime);

    initChrono.stop();
    if (verbose) std::cout << "Initialization time: " << initChrono.diff() << " s." << std::endl;

    // +-----------------------------------------------+
    // |             Computing the error               |
    // +-----------------------------------------------+
    printErrors(*solution,currentTime, uFESpace,pFESpace,verbose);

    // +-----------------------------------------------+
    // |             Solving the problem               |
    // +-----------------------------------------------+
    if (verbose) std::cout<< std::endl << "[Solving the problem]" << std::endl;
    int iter = 1;

    for ( ; currentTime <= endTime + timestep/2.; currentTime += timestep, iter++)
    {
        iterChrono.reset();
        iterChrono.start();

        if (verbose) std::cout << std::endl << "[t = "<< currentTime << " s.]" << std::endl;

        if (verbose) std::cout << "Updating the system... " << std::flush;
        bdf.bdfVelocity().updateRHSContribution( timestep );
        *rhs  = *massMatrix*bdf.bdfVelocity().rhsContributionFirstDerivative();

        systemMatrix.reset(new matrix_type( solutionMap ));
        double alpha = bdf.bdfVelocity().coefficientFirstDerivative( 0 ) / timestep;
        *systemMatrix += *massMatrix*alpha;
        *systemMatrix += *baseMatrix;

        if(convectionTerm == SemiImplicit)
        {
            *beta = bdf.bdfVelocity().extrapolation(); // Extrapolation for the convective term
            oseenAssembler.addConvection(systemMatrix,*beta);
        }
        else if(convectionTerm == Explicit)
        {
            oseenAssembler.addConvectionRhs(*rhs,*solution);
        }
        else if(convectionTerm == KIO91)
        {
            vector_type tmp(uFESpace->map(),Unique); // we do not want the pressure part
            *rhs -= bdfConvection.extrapolation();
        }
        if (verbose) std::cout << "done" << std::endl;

        if (verbose) std::cout << "Applying BC... " << std::flush;
        bcManage(*systemMatrix,*rhs,*uFESpace->mesh(),uFESpace->dof(),bcHandler,uFESpace->feBd(),1.0,currentTime);
        systemMatrix->globalAssemble();
        if (verbose) std::cout << "done" << std::endl;

        if (verbose) std::cout << "Solving the system... " << std::endl;
        *solution *= 0;
        linearSolver.setMatrix(*systemMatrix);
        linearSolver.solveSystem(*rhs,*solution,systemMatrix);

        // Updating the BDF scheme
        bdf.bdfVelocity().shiftRight( *solution );

        if(convectionTerm == KIO91)
        {
            *beta *= 0;
            oseenAssembler.addConvectionRhs(*beta,*solution);
            bdfConvection.shiftRight(*beta);
        }

        // Exporting the solution
        exporter.postProcess( currentTime );

        iterChrono.stop();
        if (verbose) std::cout << "Iteration time: " << iterChrono.diff() << " s." << std::endl;

        // +-----------------------------------------------+
        // |             Computing the error               |
        // +-----------------------------------------------+
        printErrors(*solution,currentTime, uFESpace,pFESpace,verbose);

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


