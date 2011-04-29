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
/**
   \file main.cpp
   \author Simone Rossi <simone.rossi@epfl.ch>
   \date 2011-04
 */

// LifeV definition files
#include <life/lifecore/LifeV.hpp>
#include <life/lifecore/LifeChrono.hpp>
#include <life/lifesolver/IncompressibleStructureSolver.hpp>
#include <life/lifemesh/MeshData.hpp>
#include <life/lifesolver/IncompressibleStructureData.hpp>
#include <life/lifearray/MapEpetra.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/TimeAdvanceBDFNavierStokes.hpp>
#include <life/lifemesh/MeshPartitioner.hpp>
#include <life/lifefilters/ExporterEnsight.hpp>
#include <life/lifefilters/ExporterHDF5.hpp>
#include <life/lifefilters/ExporterEmpty.hpp>
//#include <fstream> // To create an output for the flux

// Trilinos-MPI communication definitions
#include <Epetra_ConfigDefs.h>
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

// Object type definitions
typedef LifeV::RegionMesh3D<LifeV::LinearTetra>                 mesh_Type;
typedef LifeV::IncompressibleStructureSolver< mesh_Type >       fluid_Type;
typedef fluid_Type::vector_Type                                 vector_Type;
typedef boost::shared_ptr<vector_Type>                          vectorPtr_Type;   //Pointer

// +-----------------------------------------------+
// | Data and functions for the boundary conditions|
// +-----------------------------------------------+
const int BACK   = 1;
const int FRONT  = 2;
const int LEFT   = 3;
const int RIGHT  = 4;
const int BOTTOM = 5;
const int TOP    = 6;

//const int BACK   = 20;
//const int FRONT  = 20;
//const int LEFT   = 1;
//const int RIGHT  = 1;
//const int BOTTOM = 3;
//const int TOP    = 2;

LifeV::Real boundaryLoadBC(const LifeV::Real& /*t*/, const LifeV::Real& /*x*/, const LifeV::Real& /*y*/, const LifeV::Real& /*z*/, const LifeV::ID& i)
{
    switch (i)
    {
    case 2:
        return -496.0;
    default:
        return 0.0;
    }
}

LifeV::Real stressFreeBC( const LifeV::Real& /* t */,
                   const LifeV::Real& /* x */,
                   const LifeV::Real& /* y */,
                   const LifeV::Real& /* z */,
                   const LifeV::ID& /* i */ )
{
    return 0.0;
}
LifeV::Real noDisplacementBC( const LifeV::Real& /* t */,
                   const LifeV::Real& /* x */,
                   const LifeV::Real& /* y */,
                   const LifeV::Real& /* z */,
                   const LifeV::ID& /* i */ )
{
    return 0.0;
}

int main(int argc, char** argv)
{

    // +-----------------------------------------------+
    // |            Initialization of MPI              |
    // +-----------------------------------------------+
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
#endif

    boost::shared_ptr<Epetra_Comm>   comm;
#ifdef EPETRA_MPI
    comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#else
    comm.reset( new Epetra_SerialComm() );
#endif

    bool verbose(false);
    if (comm->MyPID() == 0)
    {
        verbose = true;
        std::cout
            << " +-----------------------------------------------+" << std::endl
            << " |  Incompressible Structure example for LifeV   |" << std::endl
            << " +-----------------------------------------------+" << std::endl
            << std::endl
            << " +-----------------------------------------------+" << std::endl
            << " |           Author: Simone Rossi                |" << std::endl
            << " |             Date: April, 2011                 |" << std::endl
            << " +-----------------------------------------------+" << std::endl
            << std::endl;

        std::cout << "[Initilization of MPI]" << std::endl;
#ifdef HAVE_MPI
        std::cout << "Using MPI (" << nproc << " proc.)" << std::endl;
#else
        std::cout << "Using serial version" << std::endl;
#endif
    }

    // Chronometer
    LifeV::LifeChrono globalChrono;
    LifeV::LifeChrono initChrono;
    LifeV::LifeChrono iterChrono;
    globalChrono.start();
    initChrono.start();

    // +-----------------------------------------------+
    // |               Loading the data                |
    // +-----------------------------------------------+
    if (verbose) std::cout << std::endl << "[Loading the data]" << std::endl;
    GetPot command_line(argc,argv);
    const std::string dataFileName = command_line.follow("data", 2, "-f","--file");
    GetPot dataFile(dataFileName);


    // +-----------------------------------------------+
    // |               Loading the mesh                |
    // +-----------------------------------------------+
    if (verbose) std::cout << std::endl << "[Loading the mesh]" << std::endl;
    LifeV::MeshData meshData;
    meshData.setup(dataFile, "structure/space_discretization");

    if (verbose) std::cout << "Mesh file: " << meshData.meshDir() << meshData.meshFile() << std::endl;
    boost::shared_ptr< LifeV::RegionMesh3D<LifeV::LinearTetra> > fullMeshPtr(new LifeV::RegionMesh3D<LifeV::LinearTetra>);
    LifeV::readMesh(*fullMeshPtr, meshData);
    // Split the mesh between processors
    LifeV::MeshPartitioner< LifeV::RegionMesh3D<LifeV::LinearTetra> >   meshPart(fullMeshPtr, comm);

    // +-----------------------------------------------+
    // |            Creating the FE spaces             |
    // +-----------------------------------------------+
    if (verbose) std::cout << std::endl << "[Creating the FE spaces]" << std::endl;
    std::string uOrder =  dataFile( "structure/space_discretization/disp_order",   "P2");
    std::string pOrder =  dataFile( "structure/space_discretization/press_order", "P1");

    if (verbose) std::cout << "FE for the displacements: " << uOrder << std::endl
                               << "FE for the pressure: " << pOrder << std::endl;

    if (verbose) std::cout << "Building the displacements FE space... " << std::flush;
    LifeV::FESpace< LifeV::RegionMesh3D<LifeV::LinearTetra>, LifeV::MapEpetra > uFESpace(meshPart, uOrder, 3, comm);
    if (verbose)
        std::cout << "ok." << std::endl;

    if (verbose) std::cout << "Building the pressure FE space... " << std::flush;
    LifeV::FESpace< LifeV::RegionMesh3D<LifeV::LinearTetra>, LifeV::MapEpetra > pFESpace(meshPart,pOrder,1,comm);
    if (verbose) std::cout << "ok." << std::endl;

    // Total degrees of freedom (elements of matrix)
    LifeV::UInt totalDispDof   = uFESpace.map().map(LifeV::Unique)->NumGlobalElements();
    LifeV::UInt totalPressDof = pFESpace.map().map(LifeV::Unique)->NumGlobalElements();

    if (verbose) std::cout << "Total Displacements Dof: " << totalDispDof << std::endl;
    if (verbose) std::cout << "Total Pressure Dof: " << totalPressDof << std::endl;

    // +-----------------------------------------------+
    // |             Boundary conditions               |
    // +-----------------------------------------------+
    if (verbose) std::cout<< std::endl << "[Boundary conditions]" << std::endl;

    LifeV::BCFunctionBase stressFree(stressFreeBC);
    LifeV::BCFunctionBase boundaryLoad(boundaryLoadBC);
    LifeV::BCFunctionBase noDisplacement(noDisplacementBC);



    LifeV::BCHandler bcH;
    // A boundary condition in every face

    bcH.addBC( "Top"   , TOP   , LifeV::Natural  , LifeV::Full     , boundaryLoad  , 3     );
//    bcH.addBC( "Top"   , TOP   , LifeV::Natural  , LifeV::Full     , stressFree    , 3     );
    bcH.addBC( "Left"  , LEFT  , LifeV::Natural  , LifeV::Full     , stressFree    , 3     );
    bcH.addBC( "Front" , FRONT , LifeV::Natural  , LifeV::Full     , stressFree    , 3     );
    bcH.addBC( "Right" , RIGHT , LifeV::Natural  , LifeV::Full     , stressFree    , 3     );
    bcH.addBC( "Back"  , BACK  , LifeV::Natural  , LifeV::Full     , stressFree    , 3     );
    bcH.addBC( "Bottom", BOTTOM, LifeV::Essential, LifeV::Full     , noDisplacement, 3     );

    // Get the number of Lagrange Multiplyers (LM) and set the offsets
    std::vector<LifeV::bcName_Type> fluxVector = bcH.findAllBCWithType( LifeV::Flux );
    LifeV::UInt numLM = static_cast<LifeV::UInt>( fluxVector.size() );

    LifeV::UInt offset = uFESpace.map().map(LifeV::Unique)->NumGlobalElements()
                         + pFESpace.map().map(LifeV::Unique)->NumGlobalElements();

    for ( LifeV::UInt i = 0; i < numLM; ++i )
        bcH.setOffset( fluxVector[i], offset + i );

    // +-----------------------------------------------+
    // |             Creating the problem              |
    // +-----------------------------------------------+
    if (verbose) std::cout<< std::endl << "[Creating the problem]" << std::endl;

    std::cout<< std::endl << "Sono qui A: create new data" << std::endl;

    boost::shared_ptr<LifeV::IncompressibleStructureData> incompressibleStructureData(new LifeV::IncompressibleStructureData());

    std::cout<< std::endl << "Sono qui B: setup file" << std::endl;

    incompressibleStructureData->setup( dataFile );
    incompressibleStructureData->showMe(cout);


    //if (verbose) std::cout << "Time discretization order " << incompressibleStructureData->dataTime()->orderBDF() << std::endl;


    MPI_Barrier(MPI_COMM_WORLD);
    // The problem (matrix and rhs) is packed in an object called fluid
    std::cout<< std::endl << "Sono qui C: create solver" << std::endl;

    LifeV::IncompressibleStructureSolver< LifeV::RegionMesh3D<LifeV::LinearTetra> > structure (incompressibleStructureData,
                                                                                               uFESpace,
                                                                                               pFESpace,
                                                                                               comm,
                                                                                               numLM);
    MPI_Barrier(MPI_COMM_WORLD);
    // Gets inputs from the data file
    std::cout<< std::endl << "Sono qui D: setup solver" << std::endl;

    structure.setUp(dataFile);

    // Assemble the matrices
    std::cout<< std::endl << "Sono qui E: build matrices" << std::endl;

    structure.buildSystem();


    std::cout<< std::endl << "Sono qui F: communication map" << std::endl;

    // Communication map
    LifeV::MapEpetra fullMap(structure.getMap());


    std::cout<< std::endl << "Sono qui G: MPI BARRIER" << std::endl;

    // Synchronization
    MPI_Barrier(MPI_COMM_WORLD);

    // +-----------------------------------------------+
    // |       Solve the problem                       |
    // +-----------------------------------------------+
    if (verbose) std::cout<< std::endl << "[Initialization of the simulation]" << std::endl;


//    vector_Type beta( fullMap );
    vector_Type rhs ( fullMap );

    MPI_Barrier(MPI_COMM_WORLD);

    std::string   RHS1="RHS1";
            rhs.spy(RHS1);

//    beta *= 0.;
    rhs *= 0.0;
    //rhs  -= 0.003;


    std::string   RHS="RHS";
    rhs.spy(RHS);
    //myfile.open ("example.txt");
    //myfile << "ciao";
    //myfile.close();

    std::cout<< std::endl << "Sono qui 1: update system" << std::endl;
    structure.updateSystem(rhs);

    std::cout<< std::endl << "Sono qui 2: iterate" << std::endl;
    structure.iterate(bcH);


    // +-----------------------------------------------+
    // |             Export solution                   |
    // +-----------------------------------------------+

    boost::shared_ptr< LifeV::ExporterHDF5<LifeV::RegionMesh3D<LifeV::LinearTetra> > > exporter;

    vectorPtr_Type dispAndPressure;

    std::string const exporterType =  dataFile( "exporter/type", "ensight");

    exporter.reset( new LifeV::ExporterHDF5<LifeV::RegionMesh3D<LifeV::LinearTetra> > ( dataFile, "example_incompressibleStructure" ) );
    exporter->setPostDir( "./" ); // This is a test to see if M_post_dir is working
    exporter->setMeshProcId( meshPart.meshPartition(), comm->MyPID() );

    dispAndPressure.reset( new vector_Type(*structure.solution(), exporter->mapType() ) );

    exporter->addVariable( LifeV::ExporterData::Vector, "displacement", dispAndPressure,
                           LifeV::UInt(0), uFESpace.dof().numTotalDof() );

    exporter->addVariable( LifeV::ExporterData::Scalar, "pressure", dispAndPressure,
                           LifeV::UInt(3*uFESpace.dof().numTotalDof()),
                           LifeV::UInt(pFESpace.dof().numTotalDof()) );
    exporter->postProcess( 0 );

    initChrono.stop();
    if (verbose) std::cout << "Initialization time:  " << initChrono.diff() << " s." << std::endl;

    // +-----------------------------------------------+
    // |             End of the simulation             |
    // +-----------------------------------------------+


    globalChrono.stop();
    if (verbose) std::cout << "Total simulation time:  " << globalChrono.diff() << " s." << std::endl;

    exporter->closeFile();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    std::cout << std::endl << " The End of the Simulation " << std::endl;
    return 0;
}

