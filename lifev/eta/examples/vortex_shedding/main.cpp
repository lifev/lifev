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

    @author Toni Lassila <toni.lassila@epfl.ch>
    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 23-11-2012
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

#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>

#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>

#include <lifev/core/fem/BCManage.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorBlockMonolithicEpetra.hpp>


#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <lifev/core/algorithm/SolverAztecOO.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>

#include <lifev/navier_stokes/solver/StabilizationIP.hpp>
#include <lifev/core/fem/Assembly.hpp>
#include <lifev/core/fem/AssemblyElemental.hpp>

using namespace LifeV;

namespace
{

enum DiffusionType {ViscousStress, StiffStrain};
enum MeshType {RegularMesh, File};
enum InitType {Interpolation, Projection};
enum ConvectionType {Explicit, SemiImplicit};
enum StabilizationType {None, VMS, IP};

typedef RegionMesh<LinearTetra> mesh_Type;
typedef MatrixEpetra<Real> matrix_Type;
typedef VectorEpetra vector_Type;
typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;

typedef MatrixEpetraStructured<Real> matrix_block_type;
typedef VectorBlockMonolithicEpetra vector_block_type;
typedef boost::shared_ptr<VectorBlockMonolithicEpetra> vector_blockPtr_type;

typedef FESpace< mesh_Type, MapEpetra > fespace_Type;
typedef boost::shared_ptr< fespace_Type > fespacePtr_Type;

typedef FESpace< mesh_Type, MapEpetra > fespace_Type;
typedef boost::shared_ptr< fespace_Type > fespacePtr_Type;

typedef ETFESpace< mesh_Type, MapEpetra, 3, 3 > etaUspace_Type;
typedef boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 3 > > etaUspacePtr_Type;

typedef ETFESpace< mesh_Type, MapEpetra, 3, 1 > etaPspace_Type;
typedef boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > > etaPspacePtr_Type;

}

Real fluxFunction (const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return (1 * t / 0.5 * (t < 0.5) + 1 * (t > 0.5) ); // Ramp function
}

Real zeroFunction (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0;
}

Real inflowFunction (const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    if (i == 0)
    {
        Real ux = 1; // flat velocity profile
        return (1 * ux * t / 0.5 * (t < 0.5) + 1 * ux * (t > 0.5) );
    }
    else
    {
        return 0;
    }

}

/* NormalizeFct */
class NormalizeFct
{
public:
    typedef VectorSmall<3> return_Type;

    return_Type operator() (const VectorSmall<3>& value)
    {
        Real norm (sqrt ( value[0]*value[0] + value[1]*value[1] + value[2]*value[2]) );

        if (norm > 0)
        {
            return value * (1.0 / norm);
        }
        return value;
    }

    NormalizeFct() {}
    NormalizeFct (const NormalizeFct&) {}
    ~NormalizeFct() {}
};


/* NormalizeFct */
class NormFct
{
public:
    typedef Real return_Type;

    return_Type operator() (const VectorSmall<3> value)
    {
        Real norm (sqrt ( value[0]*value[0] + value[1]*value[1] + value[2]*value[2]) );

        return norm;
    }

    NormFct() {}
    NormFct (const NormFct&) {}
    ~NormFct() {}
};

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
                << " |   Vortex shedding example w/ETA assembly      |" << std::endl
                << " +-----------------------------------------------+" << std::endl
                << std::endl
                << " +-----------------------------------------------+" << std::endl
                << " |           Author: Toni Lassila                |" << std::endl
                << " |             Date: 2012-11-13                  |" << std::endl
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

    // Physical quantity (corresponds to Re = 100)
    const Real viscosity      = 0.01 / 5;
    const Real density        = 1.0;

    // Time discretization
    const Real initialTime    = 0.0;
    const Real endTime        = 100.0;
    const Real timestep       = 1e-1;

    // Space discretization
    const UInt numDimensions  = 3;

    // Time discretization
    UInt BDFOrder = 2;

    // Stabilization method
    const StabilizationType stabilizationMethod = VMS;

    // Numerical scheme
    const InitType initializationMethod = Interpolation;
    const ConvectionType convectionTerm = SemiImplicit;

    // +-----------------------------------------------+
    // |               Loading the mesh                |
    // +-----------------------------------------------+
    if (verbose)
    {
        std::cout << std::endl << "[Loading the mesh]" << std::endl;
    }

    boost::shared_ptr<RegionMesh<LinearTetra> > fullMeshPtr (new RegionMesh<LinearTetra>);


    MeshData meshData;
    meshData.setup (dataFile, "fluid/space_discretization");
    readMesh (*fullMeshPtr, meshData);

    if (verbose) std::cout << "Mesh source: file("
                               << meshData.meshDir() << meshData.meshFile() << ")" << std::endl;

    if (verbose)
    {
        std::cout << "Mesh size  : " <<
                  MeshUtility::MeshStatistics::computeSize (*fullMeshPtr).maxH << std::endl;
    }
    if (verbose)
    {
        std::cout << "Partitioning the mesh ... " << std::endl;
    }
    MeshPartitioner< RegionMesh<LinearTetra> >   meshPart (fullMeshPtr, Comm);
    fullMeshPtr.reset(); //Freeing the global mesh to save memory

    // +-----------------------------------------------+
    // |            Creating the FE spaces             |
    // +-----------------------------------------------+
    if (verbose)
    {
        std::cout << std::endl << "[Creating the FE spaces]" << std::endl;
    }
    std::string uOrder ("P1");
    std::string pOrder ("P1");

    if (verbose) std::cout << "FE for the velocity: " << uOrder << std::endl
                               << "FE for the pressure: " << pOrder << std::endl;

    if (verbose)
    {
        std::cout << "Building the velocity FE space ... " << std::flush;
    }
    fespacePtr_Type uFESpace ( new FESpace< mesh_Type, MapEpetra > (meshPart, uOrder, numDimensions, Comm) );
    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    if (verbose)
    {
        std::cout << "Building the pressure FE space ... " << std::flush;
    }
    fespacePtr_Type pFESpace ( new FESpace< mesh_Type, MapEpetra > (meshPart, pOrder, 1, Comm) );
    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    // Define the ETA spaces for velocity and pressure
    etaUspacePtr_Type ETuFESpace ( new etaUspace_Type (meshPart, & (uFESpace->refFE() ), Comm) );
    etaPspacePtr_Type ETpFESpace ( new etaPspace_Type (meshPart, & (pFESpace->refFE() ), Comm) );

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
    BCFunctionBase uZero ( zeroFunction );
    BCFunctionBase uInflow ( inflowFunction );


    std::vector<LifeV::ID> zComp (1);
    zComp[0] = 2;

    if (verbose)
    {
        std::cout << "Setting Neumann BC... " << std::flush;
    }
    bcHandler.addBC ( "Outflow", 3, Natural, Full, uZero, 3 );
    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    if (verbose)
    {
        std::cout << "Setting Dirichlet BC... " << std::flush;
    }
    bcHandler.addBC ( "Inflow",   1, Essential, Full,      uInflow, 3 );
    bcHandler.addBC ( "Wall",     2, Essential, Full,      uInflow, 3 );
    bcHandler.addBC ( "Cube",     4, Essential, Full,      uZero,   3 );
    bcHandler.addBC ( "Symmetry", 5, Essential, Component, uZero,   zComp );

    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    // Update the BCHandler (internal data related to FE)
    bcHandler.bcUpdate ( *meshPart.meshPartition(), uFESpace->feBd(), uFESpace->dof() );

    // +-----------------------------------------------+
    // |              Matrices Assembly                |
    // +-----------------------------------------------+
    if (verbose)
    {
        std::cout << std::endl << "[Matrices Assembly]" << std::endl;
    }

    if (verbose)
    {
        std::cout << "Defining the matrices... " << std::flush;
    }

    // Initialize the full matrix, the nonconvective part, and the mass matrix (for time advancing)
    boost::shared_ptr<matrix_block_type> systemMatrix (new matrix_block_type ( ETuFESpace->map() | ETpFESpace->map() ) );
    boost::shared_ptr<matrix_block_type> baseMatrix (new matrix_block_type ( ETuFESpace->map() | ETpFESpace->map() ) );
    boost::shared_ptr<matrix_block_type> convMatrix (new matrix_block_type ( ETuFESpace->map() | ETpFESpace->map() ) );
    boost::shared_ptr<matrix_block_type> massMatrix (new matrix_block_type ( ETuFESpace->map() | ETpFESpace->map() ) );

    *systemMatrix *= 0.0;
    *baseMatrix   *= 0.0;
    *convMatrix   *= 0.0;
    *massMatrix   *= 0.0;

    if (verbose)
    {
        std::cout << "done" << std::endl;
    }

    // Perform the assembly of the base matrix with ETA

    {
        boost::shared_ptr<NormalizeFct> normalize (new NormalizeFct);

        using namespace ExpressionAssembly;

        // Use of stiff-strain formulation due to stabilization is mandatory
        integrate (
            elements (ETuFESpace->mesh() ), // Mesh
            uFESpace->qr(), // QR
            ETuFESpace,
            ETuFESpace,
            0.5 * viscosity * dot (grad (phi_i) + transpose (grad (phi_i) ), grad (phi_j) + transpose (grad (phi_j) ) )
        )
                >> baseMatrix->block (0, 0);

        integrate (
            elements (ETuFESpace->mesh() ), // Mesh
            uFESpace->qr(), // QR
            ETuFESpace,
            ETpFESpace,
            value (-1.0) * phi_j * div (phi_i)
        )
                >> baseMatrix->block (0, 1);


        integrate (
            elements (ETuFESpace->mesh() ), // Mesh
            uFESpace->qr(), // QR
            ETpFESpace,
            ETuFESpace,
            phi_i * div (phi_j)
        )
                >> baseMatrix->block (1, 0);

        integrate (
            elements (ETuFESpace->mesh() ), // Mesh
            uFESpace->qr(), // QR
            ETuFESpace,
            ETuFESpace,
            // NS
            dot (phi_i, phi_j)
        )
                >> massMatrix->block (0, 0);
    }

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
    convMatrix->globalAssemble();

    if (verbose)
    {
        std::cout << "done" << std::endl;
    }

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
    vectorPtr_Type prevRhs;

    rhs.reset (new vector_Type (solutionMap, Unique) );
    prevRhs.reset (new vector_Type (solutionMap, Unique) );

    vectorPtr_Type velocityExtrapolated;
    velocityExtrapolated.reset (new vector_Type (solutionMap, Repeated) );

    vector_Type convect (rhs->map() );

    vectorPtr_Type velocity;
    velocity.reset (new vector_Type (uFESpace->map(), Unique) );

    vectorPtr_Type pressure;
    vectorPtr_Type prevPressure;
    pressure.reset (new vector_Type (pFESpace->map(), Unique) );
    prevPressure.reset (new vector_Type (pFESpace->map(), Unique) );

    vectorPtr_Type solution;
    vectorPtr_Type prevSolution;
    vectorPtr_Type prevSolutionTimeDerivative;

    solution.reset (new vector_Type (solutionMap, Unique) );
    prevSolution.reset (new vector_Type (solutionMap, Unique) );
    prevSolutionTimeDerivative.reset (new vector_Type (solutionMap, Unique) );
    if (verbose)
    {
        std::cout << "done" << std::endl;
    }

    if (verbose)
    {
        std::cout << "Computing the initial solution ... " << std::endl;
    }

    // TimeAdvanceBDF object to store the previous solutions
    TimeAdvanceBDF<vector_Type> bdf;
    bdf.setup (BDFOrder);

    TimeAdvanceBDF<vector_Type> bdfConvection;
    bdfConvection.setup (BDFOrder);

    Real currentTime = initialTime - timestep * BDFOrder;

    // Start from zero velocity and ramp up
    *velocity *= 0;
    *pressure *= 0;

    *solution *= 0;
    *prevRhs  *= 0.;
    *prevSolutionTimeDerivative *= 0.;
    *prevSolution *= 0;

    bdf.setInitialCondition ( *solution );

    // Initial solution (interpolation or projection)
    currentTime += timestep;
    for ( ; currentTime <=  initialTime + timestep / 2.; currentTime += timestep)
    {
        *rhs  *= 0.;
        *velocityExtrapolated *= 0.;
        *solution *= 0;

        if (initializationMethod == Projection)
        {
            *rhs *= -1.;

            if (verbose)
            {
                std::cout << "Updating the system... " << std::flush;
            }
            systemMatrix.reset (new matrix_block_type ( solutionMap ) );
            *systemMatrix += *baseMatrix;

            // Assemble the convective term

            if (convectionTerm == SemiImplicit)
            {
                convMatrix.reset (new matrix_block_type ( solutionMap ) );

                // Perform the assembly of the convection matrix with ETA

                {
                    using namespace ExpressionAssembly;

                    integrate (
                        elements (ETuFESpace->mesh() ), // Mesh
                        uFESpace->qr(), // QR
                        ETuFESpace,
                        ETuFESpace,
                        dot (grad (phi_j) * value (ETuFESpace, *velocityExtrapolated), phi_i)
                    )
                            >> convMatrix->block (0, 0);
                }
                *systemMatrix += *convMatrix;
            }
            else if (convectionTerm == Explicit)
            {
                //oseenAssembler.addConvectionRhs(*rhs,1.,*solution);
                /*{
                    using namespace ExpressionAssembly;

                    integrate(
                                    elements(ETuFESpace->mesh()), // Mesh
                                    uFESpace->qr(), // QR
                                    ETuFESpace,
                                    ETuFESpace,
                                    dot( value(ETuFESpace,*solution), phi_i)
                    )
                    >> rhs->block(0);
                }*/
            }

            boost::shared_ptr<matrix_block_type> stabMatrix (new matrix_block_type ( ETuFESpace->map() | ETpFESpace->map() ) );
            *stabMatrix *= 0;
            stabMatrix->globalAssemble();
            *systemMatrix += *stabMatrix;

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

            // LinearSolver needs the monolithic matrix
            boost::shared_ptr<matrix_Type> systemMatrixNoBlock (new matrix_Type ( systemMatrix->matrixPtr() ) );
            linearSolver.setMatrix (*systemMatrix);
            linearSolver.solveSystem (*rhs, *solution, systemMatrixNoBlock);
        }

        // Updating bdf
        bdf.shiftRight ( *solution );

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
    exporter.setMeshProcId ( meshPart.meshPartition(), Comm->MyPID() );
    if (verbose)
    {
        std::cout << "done" << std::endl;
    }

    if (verbose)
    {
        std::cout << "Updating the exporter... " << std::flush;
    }
    exporter.addVariable ( ExporterData<mesh_Type>::VectorField, "velocity", uFESpace, solution, UInt (0) );
    exporter.addVariable ( ExporterData<mesh_Type>::ScalarField, "pressure", pFESpace, solution, pressureOffset );
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

        systemMatrix.reset (new matrix_block_type ( solutionMap ) );
        double alpha = bdf.coefficientFirstDerivative ( 0 ) / timestep;
        *systemMatrix += *massMatrix * alpha;
        *systemMatrix += *baseMatrix;

        if (convectionTerm == SemiImplicit)
        {
            bdf.extrapolation ( *velocityExtrapolated ); // Extrapolation for the convective term

            convMatrix.reset (new matrix_block_type ( solutionMap ) );

            // Perform the assembly of the convection matrix with ETA

            {
                using namespace ExpressionAssembly;

                integrate (
                    elements (ETuFESpace->mesh() ), // Mesh
                    uFESpace->qr(), // QR
                    ETuFESpace,
                    ETuFESpace,
                    dot (grad (phi_j) * value (ETuFESpace, *velocityExtrapolated), phi_i)
                )
                        >> convMatrix->block (0, 0);
            }
            convMatrix->globalAssemble();
            *systemMatrix += *convMatrix;
        }
        else if (convectionTerm == Explicit)
        {
            //oseenAssembler.addConvectionRhs(*rhs,1.,*solution);
            /*{
                using namespace ExpressionAssembly;

                integrate(
                                elements(ETuFESpace->mesh()), // Mesh
                                uFESpace->qr(), // QR
                                ETuFESpace,
                                ETuFESpace,
                                dot( value(ETuFESpace,*solution), phi_i)
                )
                >> rhs->block(0);
            }*/
        }

        // Assemble the stabilization terms

        // Stabilization parameters for SUPG/PSPG/LSIC/VMS/LES
        Real TauM (1e-3); // 1e-3
        Real TauC (5e-2); // 5e-2

        /* Stabilization Macros for the stabilization terms (later add the ones coming from VMS) */

#define SUPG_TEST    value(TauM) * h_K * (grad(phi_i) * value(ETuFESpace, *velocityExtrapolated))
#define VMS_TEST     value(TauM) * h_K * (transpose(grad(phi_i)) * value(ETuFESpace, *velocityExtrapolated))
#define PSPG_TEST    value(TauM) * h_K * grad(phi_i)
#define DIVDIV_TEST  value(TauC) * h_K * div(phi_i)

#define RES_MOMENTUM grad(phi_j) * value(ETuFESpace, *velocityExtrapolated) * density + value(alpha) * density * phi_j
#define RES_MASS     div(phi_j)
#define RES_PRESSURE grad(phi_j)

#define RESIDUAL_EXPLICIT value(TauM) * h_K * value(ETuFESpace, *residualVector)

        boost::shared_ptr<matrix_block_type> stabMatrix (new matrix_block_type ( ETuFESpace->map() | ETpFESpace->map() ) );
        boost::shared_ptr<matrix_block_type> residualMatrix (new matrix_block_type ( ETuFESpace->map() | ETpFESpace->map() ) );
        vector_block_type NSRhs ( ETuFESpace->map() | ETpFESpace->map(), Repeated );
        vectorPtr_Type residualVector;
        residualVector.reset (new vector_Type (solutionMap, Unique) );

        *stabMatrix     *= 0;
        *residualMatrix *= 0;
        *residualVector *= 0;
        NSRhs *= 0.0;

        if (stabilizationMethod == VMS)
        {

            boost::shared_ptr<NormalizeFct> normalize (new NormalizeFct);

            using namespace ExpressionAssembly;

            // Residual computation
            integrate (
                elements (ETuFESpace->mesh() ),
                uFESpace->qr(),
                ETuFESpace,
                ETuFESpace,
                dot (RES_MOMENTUM, phi_i)
            )
                    >> residualMatrix->block (0, 0);

            integrate (
                elements (ETuFESpace->mesh() ),
                uFESpace->qr(),
                ETuFESpace,
                ETpFESpace,
                dot (RES_PRESSURE, phi_i)
            )
                    >> residualMatrix->block (0, 1);

            // Assemble the system matrix used for the residual computation
            residualMatrix->globalAssemble();
            bcManage (*residualMatrix, *prevRhs, *uFESpace->mesh(), uFESpace->dof(), bcHandler, uFESpace->feBd(), 1.0, currentTime);

            // Compute the residual from the previous time step for the turbulence terms
            *residualVector = *prevRhs - (*residualMatrix * (*prevSolution) );

            // Stabilization, SUPG and DIV/DIV (1) and (2), VMS
            integrate (
                elements (ETuFESpace->mesh() ),
                uFESpace->qr(),
                ETuFESpace,
                ETuFESpace,
                RES_MASS * DIVDIV_TEST // OK
                + dot (RES_MOMENTUM, SUPG_TEST) // OK
                + dot (RES_MOMENTUM, VMS_TEST) // OK
            )
                    >> stabMatrix->block (0, 0);

            // Stabilization, SUPG (3), VMS
            integrate (
                elements (ETuFESpace->mesh() ),
                uFESpace->qr(),
                ETuFESpace,
                ETpFESpace,
                dot (RES_PRESSURE, SUPG_TEST) // OK
                + dot (RES_PRESSURE, VMS_TEST)
            )
                    >> stabMatrix->block (0, 1);


            // Stabilization, PSPG (4)
            integrate (
                elements (ETuFESpace->mesh() ),
                uFESpace->qr(),
                ETpFESpace,
                ETuFESpace,
                dot (RES_MOMENTUM, PSPG_TEST) // OK
            )
                    >> stabMatrix->block (1, 0);

            // Stabilization, PSPG (5)
            integrate (
                elements (ETuFESpace->mesh() ),
                uFESpace->qr(),
                ETpFESpace,
                ETpFESpace,
                dot (RES_PRESSURE, PSPG_TEST) // OK
            )
                    >> stabMatrix->block (1, 1);


            stabMatrix->globalAssemble();
            *systemMatrix += *stabMatrix;
        }
        else if (stabilizationMethod == IP)
        {
            details::StabilizationIP<mesh_Type, DOF> M_ipStabilization;
            M_ipStabilization.setFeSpaceVelocity ( *uFESpace );
            M_ipStabilization.setViscosity ( viscosity );

            // Parameters from J. Michalik
            M_ipStabilization.setGammaBeta ( 0.1 );
            M_ipStabilization.setGammaDiv  ( 0.1 );
            M_ipStabilization.setGammaPress ( 0.1 );

            M_ipStabilization.apply ( *stabMatrix, *velocityExtrapolated, false );
            stabMatrix->globalAssemble();
            *systemMatrix += *stabMatrix;
        }

        // Now we can assemble the consistency terms on the RHS

        if (verbose)
        {
            std::cout << "done" << std::endl;
        }

        if (stabilizationMethod == VMS)
        {
            boost::shared_ptr<NormalizeFct> normalize (new NormalizeFct);

            using namespace ExpressionAssembly;
            // RHS, consistency term for SUPG (6) and turbulence model (8)
            integrate (
                elements (ETuFESpace->mesh() ),
                uFESpace->qr(),
                ETuFESpace,
                dot (value (ETuFESpace, *rhs), SUPG_TEST)
                + dot (value (ETuFESpace, *rhs), VMS_TEST)
                + dot (grad (phi_i), outerProduct ( RESIDUAL_EXPLICIT, RESIDUAL_EXPLICIT ) ) // Explicit turbulence model
            )
                    >> NSRhs.block (0);

            // RHS, consistency term for PSPG (7)
            integrate (
                elements (ETuFESpace->mesh() ),
                pFESpace->qr(),
                ETpFESpace,
                dot (value (ETuFESpace, *rhs), PSPG_TEST)
            )
                    >> NSRhs.block (1);


            vector_block_type NSRhsUnique ( NSRhs, Unique );
            *rhs += NSRhsUnique;
        }
        if (verbose)
        {
            std::cout << "done" << std::endl;
        }

        // RHS has to assembled next

        if (verbose)
        {
            std::cout << "Applying BC... " << std::flush;
        }
        bcManage (*systemMatrix, *rhs, *uFESpace->mesh(), uFESpace->dof(), bcHandler, uFESpace->feBd(), 1.0, currentTime);


        if (verbose)
        {
            std::cout << "Solving the system... " << std::endl;
        }
        *solution *= 0;

        // LinearSolver needs the monolithic matrix
        boost::shared_ptr<matrix_Type> systemMatrixNoBlock (new matrix_Type ( systemMatrix->matrixPtr() ) );
        linearSolver.setMatrix (*systemMatrix);
        linearSolver.solveSystem (*rhs, *solution, systemMatrixNoBlock);

        // Updating the BDF scheme
        bdf.shiftRight ( *solution );

        // Exporting the solution
        exporter.postProcess ( currentTime );

        // Store previous solution, its time derivative and the RHS for the explicit turbulence model
        *prevSolution               = *solution;
        *prevRhs                    = *rhs;

        iterChrono.stop();
        if (verbose)
        {
            std::cout << "Iteration time: " << iterChrono.diff() << " s." << std::endl;
        }

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


