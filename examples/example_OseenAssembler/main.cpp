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
    @date 14-03-2011
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


#include <life/lifealg/PreconditionerIfpack.hpp>
#include <life/lifealg/PreconditionerML.hpp>
#include <life/lifealg/SolverAztecOO.hpp>
#include <life/lifearray/MatrixEpetra.hpp>
#include <life/lifefilters/ExporterHDF5.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/BCManage.hpp>
#include <life/lifemesh/MeshPartitioner.hpp>
#include <life/lifemesh/RegionMesh3DStructured.hpp>
#include <life/lifemesh/RegionMesh3D.hpp>
#include <life/lifesolver/OseenAssembler.hpp>
#include <life/lifemesh/MeshData.hpp>

#include "EthierSteinmanUnsteady.hpp"

using namespace LifeV;

namespace
{
static bool regIF = (PRECFactory::instance().registerProduct( "Ifpack", &createIfpack ));
static bool regML = (PRECFactory::instance().registerProduct( "ML", &createML ));

enum DiffusionType{ViscousStress,StiffStrain};
enum MeshType{RegularMesh,File};

typedef RegionMesh3D<LinearTetra> mesh_type;
typedef MatrixEpetra<Real> matrix_type;
typedef VectorEpetra vector_type;
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
            << " |             Date: 2010-03-14                  |" << std::endl
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

    // **** Stupid GetPot stuff ****
    GetPot command_line(argc,argv);
    const std::string dataFileName = command_line.follow("data", 2, "-f","--file");
    GetPot dataFile(dataFileName);
    // *****************************

    // Physical quantity
    const Real viscosity      = 0.01;
    const Real density        = 1.0;

    // Time discretization
    const Real initialtime    = 0.0;
    const Real endtime        = 1e-4;
    const Real timestep       = 1e-5;

    // Space discretization
    const MeshType meshSource = RegularMesh;
    const UInt numMeshElem    = 10;

    // Numerical scheme
    const DiffusionType diffusionType = ViscousStress;

    // EthierSteinman data
    EthierSteinmanUnsteady ethierSteinmanData;
    ethierSteinmanData.setA(1.0);
    ethierSteinmanData.setD(1.0);
    ethierSteinmanData.setViscosity(viscosity);
    ethierSteinmanData.setDensity(density);

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

    if (verbose) std::cout << "Partitioning the mesh ... " << std::flush;
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
    boost::shared_ptr<FESpace< mesh_type, MapEpetra > > uFESpace( new FESpace< mesh_type, MapEpetra >(meshPart,uOrder, 3, Comm));
    if (verbose) std::cout << "ok." << std::endl;

    if (verbose) std::cout << "Building the pressure FE space ... " << std::flush;
    boost::shared_ptr<FESpace< mesh_type, MapEpetra > > pFESpace( new FESpace< mesh_type, MapEpetra >(meshPart,pOrder, 1, Comm));
    if (verbose) std::cout << "ok." << std::endl;

    if (verbose) std::cout << "Total Velocity Dof = " << 3*uFESpace->dof().numTotalDof() << std::endl;
    if (verbose) std::cout << "Total Pressure Dof = " << pFESpace->dof().numTotalDof() << std::endl;

    // +-----------------------------------------------+
    // |               Matrix Assembly                 |
    // +-----------------------------------------------+
    if (verbose) std::cout << std::endl << "[Matrix Assembly]" << std::endl;

    if (verbose) std::cout << "Setting up assembler... " << std::flush;
    OseenAssembler<mesh_type,matrix_type,vector_type> oseenAssembler;
    oseenAssembler.setup(uFESpace,pFESpace);
    if (verbose) std::cout << "done" << std::endl;

    if (verbose) std::cout << "Defining the matrix... " << std::flush;
    MapEpetra solutionMap(uFESpace->map()+pFESpace->map());
    boost::shared_ptr<matrix_type> systemMatrix(new matrix_type( solutionMap ));
    *systemMatrix *=0.0;
    if (verbose) std::cout << "done" << std::endl;

    // Perform the assembly of the matrix
    switch(diffusionType)
    {
        case ViscousStress:
            if (verbose) std::cout << "Adding the viscous stress... " << std::flush;
            oseenAssembler.addViscousStress(systemMatrix,viscosity/density);
            if (verbose) std::cout << "done" << std::endl;
            break;
        case StiffStrain:
            if (verbose) std::cout << "Adding the stiff strain... " << std::flush;
            oseenAssembler.addStiffStrain(systemMatrix,viscosity/density);
            if (verbose) std::cout << "done" << std::endl;
            break;
        default:
            cerr << "[Error] Diffusion type unknown" << std::endl;
            exit(1);
            break;
    }

    if (verbose) std::cout << "Adding the convection... " << std::flush;
    vector_type beta(uFESpace->map(),Repeated);
    beta *= 0;
    beta += 1;
    oseenAssembler.addConvection(systemMatrix,beta);
    if (verbose) std::cout << "done" << std::endl;

    if (verbose) std::cout << "Adding the mass... " << std::flush;
    oseenAssembler.addMass(systemMatrix,1.0/timestep);
    if (verbose) std::cout << "done" << std::endl;

    if (verbose) std::cout << "Adding the gradient of the pressure... " << std::flush;
    oseenAssembler.addGradPressure(systemMatrix);
    if (verbose) std::cout << "done" << std::endl;

    if (verbose) std::cout << "Adding the divergence free constraint... " << std::flush;
    oseenAssembler.addDivergence(systemMatrix);
    if (verbose) std::cout << "done" << std::endl;

    if (verbose) std::cout << "Closing the matrix... " << std::flush;
    systemMatrix->globalAssemble();
    if (verbose) std::cout << "done" << std::endl;

    // +-----------------------------------------------+
    // |                 RHS Assembly                  |
    // +-----------------------------------------------+
    if (verbose) std::cout << std::endl << "[RHS Assembly]" << std::endl;
    vector_type rhs(uFESpace->map(),Repeated);
    rhs*=0.0;

    /*
#ifdef TEST_RHS
    vector_type fInterpolated(uFESpace->map(),Repeated);
    fInterpolated*=0.0;
    uFESpace->interpolate(fRhs,fInterpolated,0.0);
    adrAssembler.addMassRhs(rhs,fInterpolated);
    rhs.globalAssemble();
#endif
*/

    if (verbose) std::cout << " done ! " << std::endl;

    /*
    // +-----------------------------------------------+
    // |             Boundary conditions               |
    // +-----------------------------------------------+
    if (verbose) std::cout << " -- Building the BCHandler ... " << std::flush;
    BCHandler bchandler;
    BCFunctionBase BCu( exactSolution );
    bchandler.addBC("Dirichlet",1,Essential,Full,BCu,1);
    for (UInt i(2); i<=6; ++i)
    {
        bchandler.addBC("Dirichlet",i,Essential,Full,BCu,1);
    }
    if (verbose) std::cout << " done ! " << std::endl;

    if (verbose) std::cout << " -- Updating the BCs ... " << std::flush;
    bchandler.bcUpdate(*uFESpace->mesh(),uFESpace->feBd(),uFESpace->dof());
    if (verbose) std::cout << " done ! " << std::endl;

    if (verbose) std::cout << " -- Applying the BCs ... " << std::flush;
    vector_type rhsBC(rhs,Unique);
    bcManage(*systemMatrix,rhsBC,*uFESpace->mesh(),uFESpace->dof(),bchandler,uFESpace->feBd(),1.0,0.0);
    rhs = rhsBC;
    if (verbose) std::cout << " done ! " << std::endl;

    */

    //************* SPY ***********
    systemMatrix->spy("matrix");
    //rhs.spy("vector");
    //*****************************

    /*

    // +-----------------------------------------------+
    // |            Solver initialization              |
    // +-----------------------------------------------+
    if (verbose) std::cout << "[Solver initialization]" << std::flush;
    SolverAztecOO linearSolver;

    if (verbose) std::cout << "Setting up the solver... " << std::flush;
    linearSolver.setDataFromGetPot(dataFile,"solver");
    linearSolver.setupPreconditioner(dataFile,"prec");
    if (verbose) std::cout << "done" << std::endl;

    if (verbose) std::cout << "Setting matrix in the solver... " << std::flush;
    linearSolver.setMatrix(*systemMatrix);
    if (verbose) std::cout << "done" << std::endl;

    linearSolver.setCommunicator(Comm);

    // +-----------------------------------------------+
    // |             Solving the problem               |
    // +-----------------------------------------------+
    // Definition of the solution
    if (verbose) std::cout << "Defining the solution vector... " << std::flush;
    vector_type solution(uFESpace->map(),Unique);
    solution*=0.0;
    if (verbose) std::cout << "done" << std::endl;

    // Solve the problem
    if (verbose) std::cout << "Solving the system... " << std::flush;
    linearSolver.solveSystem(rhsBC,solution,systemMatrix);
    if (verbose) std::cout << "done" << std::endl;


    // Exporter definition and use
    if (verbose) std::cout << "Defining the exporter... " << std::flush;
    ExporterHDF5<mesh_type> exporter ( dataFile, "OseenAssembler");
    exporter->setPostDir( "./" ); // This is a test to see if M_post_dir is working
    exporter->setMeshProcId( meshPart.meshPartition(), d->comm->MyPID() );
    if (verbose) std::cout << "done" << std::endl;

    if (verbose) std::cout << "Defining the exported quantities... " << std::flush;
    boost::shared_ptr<vector_type> solutionPtr (new vector_type(solution,Repeated));
    if (verbose) std::cout << "done" << std::endl;

    if (verbose) std::cout << "Updating the exporter... " << std::flush;
    exporter->addVariable( ExporterData::Vector, "velocity", solutionPtr,
                                       UInt(0), uFESpace.dof().numTotalDof() );
    exporter->addVariable( ExporterData::Scalar, "pressure", solutionPtr,
                                       UInt(3*uFESpace.dof().numTotalDof()),
                                       UInt(pFESpace.dof().numTotalDof()) );
    if (verbose) std::cout << "done" << std::endl;

    if (verbose) std::cout << "Exporting... " << std::flush;
    exporter.postProcess(0);
    if (verbose) std::cout << "done" << std::endl;
    */

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return( EXIT_SUCCESS );
}


