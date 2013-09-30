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


#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <lifev/core/algorithm/SolverAztecOO.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>
#include <lifev/navier_stokes/solver/OseenAssembler.hpp>
#include <lifev/navier_stokes/function/RossEthierSteinmanDec.hpp>

using namespace LifeV;

namespace
{

enum DiffusionType {ViscousStress, StiffStrain};
enum MeshType {RegularMesh, File};
enum InitType {Interpolation, Projection};
enum ConvectionType {Explicit, SemiImplicit, KIO91};

/*
 * Some references for the KIO91 scheme:
 *
 * Karniadakis, G.E., Israeli, M. and Orszag, S.A. (1991)
 * High-OrderSplitting Methods for the Incompressible Navier-Stokes Equations,
 * Journal of Computational Physics, 97, 414-443.
 *
 * Canuto, C., Hussaini, M.Y., Quarteroni, A., Zang, T.A. (2007),
 * Spectral Methods Evolution to Complex Geometries and Applications to Fluid Dynamics
 */

typedef RegionMesh<LinearTetra> mesh_Type;
typedef MatrixEpetra<Real> matrix_Type;
typedef VectorEpetra vector_Type;
typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;
typedef FESpace< mesh_Type, MapEpetra > fespace_Type;
typedef boost::shared_ptr< fespace_Type > fespacePtr_Type;
typedef LifeV::RossEthierSteinmanUnsteadyDec problem_Type;
}

Real fluxFunction (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 1;
}


void printErrors (const vector_Type& solution, const Real& currentTime, fespacePtr_Type uFESpace, fespacePtr_Type pFESpace, bool verbose)
{
    vector_Type velocity (uFESpace->map(), Repeated);
    vector_Type pressure (pFESpace->map(), Repeated);
    if (verbose)
    {
        std::cout << std::endl << "[Computed errors]" << std::endl;
    }
    velocity.subset (solution);
    pressure.subset (solution, 3 * uFESpace->dof().numTotalDof() );
    Real uRelativeError, pRelativeError, uL2Error, pL2Error;
    uL2Error = uFESpace->l2Error (problem_Type::uexact, velocity, currentTime, &uRelativeError );
    pL2Error = pFESpace->l20Error (problem_Type::pexact, pressure, currentTime, &pRelativeError );
    if (verbose)
    {
        std::cout << "Velocity" << std::endl;
    }
    if (verbose)
    {
        std::cout << "  L2 error      : " << uL2Error << std::endl;
    }
    if (verbose)
    {
        std::cout << "  Relative error: " << uRelativeError << std::endl;
    }
    if (verbose)
    {
        std::cout << "Pressure" << std::endl;
    }
    if (verbose)
    {
        std::cout << "  L2 error      : " << pL2Error << std::endl;
    }
    if (verbose)
    {
        std::cout << "  Relative error: " << pRelativeError << std::endl;
    }
}


int
main ( int argc, char** argv )
{
    // +-----------------------------------------------+
    // |            Initialization of MPI              |
    // +-----------------------------------------------+
#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
    int nproc;
    MPI_Comm_size (MPI_COMM_WORLD, &nproc);
#else
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif

    const bool verbose (Comm->MyPID() == 0);
    if (verbose)
    {
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
    if (verbose)
    {
        std::cout << std::endl << "[Loading the data]" << std::endl;
    }
    LifeChrono globalChrono;
    LifeChrono initChrono;
    LifeChrono iterChrono;

    globalChrono.start();
    initChrono.start();

    // **** Stupid GetPot stuff ****
    GetPot command_line (argc, argv);
    const std::string dataFileName = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (dataFileName);
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
    const UInt numMeshElem    = 3;

    // Numerical scheme
    const DiffusionType diffusionType = ViscousStress;
    UInt BDFOrder = 3;
    const InitType initializationMethod = Interpolation;
    const ConvectionType convectionTerm = KIO91;

    // Reading settings from file
    //dataFile;


    // EthierSteinman data
    problem_Type::setA (1.0);
    problem_Type::setD (1.0);
    problem_Type::setViscosity (viscosity);
    problem_Type::setDensity (density);

    if (diffusionType == StiffStrain)
    {
        problem_Type::setFlagStrain (1);
    }
    else
    {
        problem_Type::setFlagStrain (0);
    }

    // +-----------------------------------------------+
    // |               Loading the mesh                |
    // +-----------------------------------------------+
    if (verbose)
    {
        std::cout << std::endl << "[Loading the mesh]" << std::endl;
    }

    boost::shared_ptr<RegionMesh<LinearTetra> > fullMeshPtr ( new RegionMesh<LinearTetra> ( Comm ) );

    // Building the mesh from the source
    if (meshSource == RegularMesh)
    {
        regularMesh3D ( *fullMeshPtr,
                        1,
                        numMeshElem, numMeshElem, numMeshElem,
                        false,
                        2.0,   2.0,   2.0,
                        -1.0,  -1.0,  -1.0);

        if (verbose) std::cout << "Mesh source: regular mesh("
                                   << numMeshElem << "x" << numMeshElem << "x" << numMeshElem << ")" << std::endl;
    }
    else if (meshSource == File)
    {
        MeshData meshData;
        meshData.setup (dataFile, "fluid/space_discretization");
        readMesh (*fullMeshPtr, meshData);

        if (verbose) std::cout << "Mesh source: file("
                                   << meshData.meshDir() << meshData.meshFile() << ")" << std::endl;
    }
    else
    {
        if (verbose)
        {
            std::cout << std::endl << "Error: Unknown source type for the mesh" << std::endl;
        }
        exit (1);
    }

    if (verbose)
    {
        std::cout << "Mesh size  : " <<
                  MeshUtility::MeshStatistics::computeSize (*fullMeshPtr).maxH << std::endl;
    }
    if (verbose)
    {
        std::cout << "Partitioning the mesh ... " << std::endl;
    }
    boost::shared_ptr<RegionMesh<LinearTetra> > meshPtr;
    {
        MeshPartitioner< RegionMesh<LinearTetra> >   meshPart (fullMeshPtr, Comm);
        meshPtr = meshPart.meshPartition();
    }
    fullMeshPtr.reset(); //Freeing the global mesh to save memory

    // +-----------------------------------------------+
    // |            Creating the FE spaces             |
    // +-----------------------------------------------+
    if (verbose)
    {
        std::cout << std::endl << "[Creating the FE spaces]" << std::endl;
    }
    std::string uOrder ("P2");
    std::string pOrder ("P1");

    if (verbose) std::cout << "FE for the velocity: " << uOrder << std::endl
                               << "FE for the pressure: " << pOrder << std::endl;

    if (verbose)
    {
        std::cout << "Building the velocity FE space ... " << std::flush;
    }
    fespacePtr_Type uFESpace ( new fespace_Type ( meshPtr, uOrder, numDimensions, Comm ) );
    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    if (verbose)
    {
        std::cout << "Building the pressure FE space ... " << std::flush;
    }
    fespacePtr_Type pFESpace ( new fespace_Type ( meshPtr, pOrder, 1, Comm ) );
    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    // Creation of the total map
    MapEpetra solutionMap (uFESpace->map() + pFESpace->map() );

    // Pressure offset in the vector
    UInt pressureOffset = numDimensions * uFESpace->dof().numTotalDof();

    if (verbose)
    {
        std::cout << "Total Velocity Dof: " << pressureOffset << std::endl;
    }
    if (verbose)
    {
        std::cout << "Total Pressure Dof: " << pFESpace->dof().numTotalDof() << std::endl;
    }

    // +-----------------------------------------------+
    // |             Boundary conditions               |
    // +-----------------------------------------------+
    if (verbose)
    {
        std::cout << std::endl << "[Boundary conditions]" << std::endl;
    }
    BCHandler bcHandler;
    BCFunctionBase uDirichlet ( problem_Type::uexact );
    BCFunctionBase uNeumann ( problem_Type::fNeumann );

    if (verbose)
    {
        std::cout << "Setting Neumann BC... " << std::flush;
    }
    bcHandler.addBC ( "Flux", 1, Natural, Full, uNeumann, 3 );
    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    if (verbose)
    {
        std::cout << "Setting Dirichlet BC... " << std::flush;
    }
    for (UInt iDirichlet (2); iDirichlet <= 26; ++iDirichlet)
    {
        bcHandler.addBC ( "Wall", iDirichlet, Essential, Full, uDirichlet, 3 );
    }
    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    // Update the BCHandler (internal data related to FE)
    bcHandler.bcUpdate ( *meshPtr, uFESpace->feBd(), uFESpace->dof() );

    // +-----------------------------------------------+
    // |              Matrices Assembly                |
    // +-----------------------------------------------+
    if (verbose)
    {
        std::cout << std::endl << "[Matrices Assembly]" << std::endl;
    }

    if (verbose)
    {
        std::cout << "Setting up assembler... " << std::flush;
    }
    OseenAssembler<mesh_Type, matrix_Type, vector_Type> oseenAssembler;
    oseenAssembler.setup (uFESpace, pFESpace);
    if (verbose)
    {
        std::cout << "done" << std::endl;
    }

    if (verbose)
    {
        std::cout << "Defining the matrices... " << std::flush;
    }
    boost::shared_ptr<matrix_Type> systemMatrix (new matrix_Type ( solutionMap ) );
    *systemMatrix *= 0.0;
    boost::shared_ptr<matrix_Type> baseMatrix (new matrix_Type ( solutionMap ) );
    *baseMatrix *= 0.0;
    boost::shared_ptr<matrix_Type> massMatrix (new matrix_Type ( solutionMap ) );
    *massMatrix *= 0.0;
    if (verbose)
    {
        std::cout << "done" << std::endl;
    }

    // Perform the assembly of the matrix
    switch (diffusionType)
    {
        case ViscousStress:
            if (verbose)
            {
                std::cout << "Adding the viscous stress... " << std::flush;
            }
            oseenAssembler.addViscousStress (*baseMatrix, viscosity / density);
            if (verbose)
            {
                std::cout << "done" << std::endl;
            }
            break;
        case StiffStrain:
            if (verbose)
            {
                std::cout << "Adding the stiff strain... " << std::flush;
            }
            oseenAssembler.addStiffStrain (*baseMatrix, viscosity / density);
            if (verbose)
            {
                std::cout << "done" << std::endl;
            }
            break;
        default:
            cerr << "[Error] Diffusion type unknown" << std::endl;
            return ( EXIT_FAILURE );
            break;
    }

    if (verbose)
    {
        std::cout << "Adding the gradient of the pressure... " << std::flush;
    }
    oseenAssembler.addGradPressure (*baseMatrix);
    if (verbose)
    {
        std::cout << "done" << std::endl;
    }

    if (verbose)
    {
        std::cout << "Adding the divergence free constraint... " << std::flush;
    }
    oseenAssembler.addDivergence (*baseMatrix, -1.0);
    if (verbose)
    {
        std::cout << "done" << std::endl;
    }

    if (verbose)
    {
        std::cout << "Adding the mass... " << std::flush;
    }
    oseenAssembler.addMass (*massMatrix, 1.0);
    if (verbose)
    {
        std::cout << "done" << std::endl;
    }

    if (verbose)
    {
        std::cout << "Closing the matrices... " << std::flush;
    }
    baseMatrix->globalAssemble();
    massMatrix->globalAssemble();
    if (verbose)
    {
        std::cout << "done" << std::endl;
    }

    // +-----------------------------------------------+
    // |  Flux computation: vector initialization      |
    // +-----------------------------------------------+
    BCFunctionBase flow (fluxFunction);

    BCHandler fluxHandler;
    fluxHandler.addBC ("Flux" , 1,  Flux, Normal, flow);

    // Update the BCHandler (internal data related to FE)
    fluxHandler.bcUpdate ( *meshPtr, uFESpace->feBd(), uFESpace->dof() );

    vector_Type fluxVector (solutionMap);
    oseenAssembler.addFluxTerms (fluxVector, fluxHandler);


    // +-----------------------------------------------+
    // |            Solver initialization              |
    // +-----------------------------------------------+
    if (verbose)
    {
        std::cout << std::endl << "[Solver initialization]" << std::endl;
    }
    SolverAztecOO linearSolver;

    if (verbose)
    {
        std::cout << "Setting up the solver... " << std::flush;
    }
    linearSolver.setDataFromGetPot (dataFile, "solver");
    linearSolver.setupPreconditioner (dataFile, "prec");
    if (verbose)
    {
        std::cout << "done" << std::endl;
    }

    linearSolver.setCommunicator (Comm);

    // +-----------------------------------------------+
    // |       Initialization of the simulation        |
    // +-----------------------------------------------+
    if (verbose)
    {
        std::cout << std::endl << "[Initialization of the simulation]" << std::endl;
    }
    if (verbose)
    {
        std::cout << "Creation of vectors... " << std::flush;
    }
    vectorPtr_Type rhs;
    rhs.reset (new vector_Type (solutionMap, Unique) );

    vectorPtr_Type beta;
    beta.reset (new vector_Type (solutionMap, Repeated) );

    vector_Type convect (rhs->map() );

    vectorPtr_Type velocity;
    velocity.reset (new vector_Type (uFESpace->map(), Unique) );

    vectorPtr_Type pressure;
    pressure.reset (new vector_Type (pFESpace->map(), Unique) );

    vectorPtr_Type solution;
    solution.reset (new vector_Type (solutionMap, Unique) );
    if (verbose)
    {
        std::cout << "done" << std::endl;
    }

    if (verbose)
    {
        std::cout << "Computing the initial solution ... " << std::endl;
    }
    if (convectionTerm == KIO91)
    {
        BDFOrder = 3;
    }

    // bdf object to store the previous solutions
    TimeAdvanceBDF<vector_Type> bdf;
    bdf.setup (BDFOrder);
    TimeAdvanceBDF<vector_Type> bdfConvection;
    bdfConvection.setup (BDFOrder);
    TimeAdvanceBDF<vector_Type> bdfConvectionInit; // Just for KIO91
    bdfConvectionInit.setup (BDFOrder);
    Real currentTime = initialTime - timestep * BDFOrder;

    if (convectionTerm == KIO91)
    {
        uFESpace->interpolate ( static_cast<FESpace< mesh_Type, MapEpetra >::function_Type> ( problem_Type::uexact ), *velocity, currentTime );
        *solution *= 0;
        *solution = *velocity;
        *beta *= 0;
        oseenAssembler.addConvectionRhs (*beta, 1., *solution);
        bdfConvection.setInitialCondition ( *beta );
        bdfConvectionInit.setInitialCondition ( *beta );

        if (initializationMethod == Projection)
        {
            for (UInt i (0); i < BDFOrder; ++i)
            {
                uFESpace->interpolate ( static_cast<FESpace< mesh_Type, MapEpetra >::function_Type> ( problem_Type::uexact ), *velocity, currentTime - (3 - i) *timestep );
                *solution = *velocity;
                *beta *= 0;
                oseenAssembler.addConvectionRhs (*beta, 1., *solution);
                bdfConvectionInit.shiftRight ( *beta );
            }
        }
    }

    uFESpace->interpolate ( static_cast<FESpace< mesh_Type, MapEpetra >::function_Type> ( problem_Type::uexact ), *velocity, currentTime );
    pFESpace->interpolate ( static_cast<FESpace< mesh_Type, MapEpetra >::function_Type> ( problem_Type::pexact ), *pressure, currentTime );
    *solution *= 0;
    *solution = *velocity;
    solution->add (*pressure, pressureOffset);
    bdf.setInitialCondition ( *solution );

    // Compute initial flux trough face 1
    Real fluxThrou1 = fluxVector.dot (*solution);
    if (verbose) std::cout << "Flux through face 1 = "
                               << fluxThrou1 << std::endl;

    // Initial solution (interpolation or projection)
    currentTime += timestep;
    for ( ; currentTime <=  initialTime + timestep / 2.; currentTime += timestep)
    {
        *rhs  *= 0.;
        *beta *= 0.;
        *solution *= 0;

        uFESpace->interpolate ( static_cast<FESpace< mesh_Type, MapEpetra >::function_Type> ( problem_Type::uexact ), *velocity, currentTime );
        pFESpace->interpolate ( static_cast<FESpace< mesh_Type, MapEpetra >::function_Type> ( problem_Type::pexact ), *pressure, currentTime );
        *solution = *velocity;
        solution->add (*pressure, pressureOffset);

        if (initializationMethod == Projection)
        {
            oseenAssembler.addMassRhs (*rhs, problem_Type::uderexact , currentTime);
            *rhs *= -1.;


            if (verbose)
            {
                std::cout << "Updating the system... " << std::flush;
            }
            systemMatrix.reset (new matrix_Type ( solutionMap ) );
            *systemMatrix += *baseMatrix;
            if (convectionTerm == SemiImplicit)
            {
                oseenAssembler.addConvection (*systemMatrix, 1.0, *solution);
            }
            else if (convectionTerm == Explicit)
            {
                oseenAssembler.addConvectionRhs (*rhs, 1., *solution);
            }
            else if (convectionTerm == KIO91)
            {
                bdfConvectionInit.extrapolation (convect);
                *rhs -= convect;
                *beta *= 0;
                oseenAssembler.addConvectionRhs (*beta, 1., *solution);
                bdfConvectionInit.shiftRight (*beta);
            }

            if (verbose)
            {
                std::cout << "done" << std::endl;
            }

            if (verbose)
            {
                std::cout << "Applying BC... " << std::flush;
            }
            bcManage (*systemMatrix, *rhs, *uFESpace->mesh(), uFESpace->dof(), bcHandler, uFESpace->feBd(), 1.0, currentTime);
            systemMatrix->globalAssemble();
            if (verbose)
            {
                std::cout << "done" << std::endl;
            }

            if (verbose)
            {
                std::cout << "Solving the system... " << std::endl;
            }
            *solution *= 0;
            linearSolver.setMatrix (*systemMatrix);
            linearSolver.solveSystem (*rhs, *solution, systemMatrix);
        }

        // Updating bdf
        bdf.shiftRight ( *solution );

        if (convectionTerm == KIO91)
        {
            *beta *= 0;
            oseenAssembler.addConvectionRhs (*beta, 1., *solution);
            bdfConvection.shiftRight (*beta);
        }

    }

    linearSolver.resetPreconditioner();

    // +-----------------------------------------------+
    // |             Setting the exporter              |
    // +-----------------------------------------------+
    if (verbose)
    {
        std::cout << "Defining the exporter... " << std::flush;
    }
    ExporterHDF5<mesh_Type> exporter ( dataFile, "OseenAssembler");
    exporter.setPostDir ( "./" ); // This is a test to see if M_post_dir is working
    exporter.setMeshProcId ( meshPtr, Comm->MyPID() );
    if (verbose)
    {
        std::cout << "done" << std::endl;
    }

    if (verbose)
    {
        std::cout << "Updating the exporter... " << std::flush;
    }
    exporter.addVariable ( ExporterData<mesh_Type>::VectorField, "velocity", uFESpace,
                           solution, UInt (0) );
    exporter.addVariable ( ExporterData<mesh_Type>::ScalarField, "pressure", pFESpace,
                           solution, pressureOffset );
    if (verbose)
    {
        std::cout << "done" << std::endl;
    }

    if (verbose)
    {
        std::cout << "Exporting solution at time t=" << initialTime << "... " << std::endl;
    }
    exporter.postProcess (initialTime);

    initChrono.stop();
    if (verbose)
    {
        std::cout << "Initialization time: " << initChrono.diff() << " s." << std::endl;
    }

    // +-----------------------------------------------+
    // |             Computing the error               |
    // +-----------------------------------------------+
    printErrors (*solution, currentTime, uFESpace, pFESpace, verbose);

    // +-----------------------------------------------+
    // |             Solving the problem               |
    // +-----------------------------------------------+
    if (verbose)
    {
        std::cout << std::endl << "[Solving the problem]" << std::endl;
    }
    int iter = 1;

    for ( ; currentTime <= endTime + timestep / 2.; currentTime += timestep, iter++)
    {
        iterChrono.reset();
        iterChrono.start();

        if (verbose)
        {
            std::cout << std::endl << "[t = " << currentTime << " s.]" << std::endl;
        }

        if (verbose)
        {
            std::cout << "Updating the system... " << std::flush;
        }
        bdf.updateRHSContribution ( timestep );
        *rhs  = *massMatrix * bdf.rhsContributionFirstDerivative();

        systemMatrix.reset (new matrix_Type ( solutionMap ) );
        double alpha = bdf.coefficientFirstDerivative ( 0 ) / timestep;
        *systemMatrix += *massMatrix * alpha;
        *systemMatrix += *baseMatrix;

        if (convectionTerm == SemiImplicit)
        {
            //  *beta = bdf.extrapolation(); // Extrapolation for the convective term
            bdf.extrapolation ( *beta ); // Extrapolation for the convective term
            oseenAssembler.addConvection (*systemMatrix, 1.0, *beta);
        }
        else if (convectionTerm == Explicit)
        {
            oseenAssembler.addConvectionRhs (*rhs, 1., *solution);
        }
        else if (convectionTerm == KIO91)
        {
            //   *rhs -= bdfConvection.extrapolation();
            convect *= 0.;
            bdfConvection.extrapolation (convect);
            *rhs -= convect;
        }
        if (verbose)
        {
            std::cout << "done" << std::endl;
        }

        if (verbose)
        {
            std::cout << "Applying BC... " << std::flush;
        }
        bcManage (*systemMatrix, *rhs, *uFESpace->mesh(), uFESpace->dof(), bcHandler, uFESpace->feBd(), 1.0, currentTime);
        systemMatrix->globalAssemble();
        if (verbose)
        {
            std::cout << "done" << std::endl;
        }

        if (verbose)
        {
            std::cout << "Solving the system... " << std::endl;
        }
        *solution *= 0;
        linearSolver.setMatrix (*systemMatrix);
        linearSolver.solveSystem (*rhs, *solution, systemMatrix);

        // Updating the BDF scheme
        bdf.shiftRight ( *solution );

        if (convectionTerm == KIO91)
        {
            *beta *= 0;
            oseenAssembler.addConvectionRhs (*beta, 1., *solution);
            bdfConvection.shiftRight (*beta);
        }

        // Exporting the solution
        exporter.postProcess ( currentTime );

        iterChrono.stop();
        if (verbose)
        {
            std::cout << "Iteration time: " << iterChrono.diff() << " s." << std::endl;
        }

        // +-----------------------------------------------+
        // |             Computing the error               |
        // +-----------------------------------------------+
        printErrors (*solution, currentTime, uFESpace, pFESpace, verbose);

        MPI_Barrier (MPI_COMM_WORLD);
    }

    // +-----------------------------------------------+
    // |            Ending the simulation              |
    // +-----------------------------------------------+
    exporter.closeFile();

    globalChrono.stop();
    if (verbose)
    {
        std::cout << std::endl << "Total simulation time: " << globalChrono.diff() << " s." << std::endl;
    }
    if (verbose)
    {
        std::cout << std::endl << "[[END_SIMULATION]]" << std::endl;
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return ( EXIT_SUCCESS );
}


