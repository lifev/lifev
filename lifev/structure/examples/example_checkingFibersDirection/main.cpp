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

   This exampleis meant to help people checking if their analytical definition
   of the collagen fibers families in the computational domain is correct.

   Successfully used for cube, cylinder, torus, idealized aneurysm

   \author Paolo Tricerri <paolo.tricerri@epfl.ch>
   \date 1861-03-17
 */

#ifdef TWODIM
#error test_structure cannot be compiled in 2D
#endif

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

#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/PartitionIO.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <iostream>
#include <fstream>

#include "ud_functions.hpp"


using namespace LifeV;

std::set<UInt> parseList ( const std::string& list )
{
    std::string stringList = list;
    std::set<UInt> setList;
    if ( list == "" )
    {
        return setList;
    }
    size_t commaPos = 0;
    while ( commaPos != std::string::npos )
    {
        commaPos = stringList.find ( "," );
        setList.insert ( atoi ( stringList.substr ( 0, commaPos ).c_str() ) );
        stringList = stringList.substr ( commaPos + 1 );
    }
    setList.insert ( atoi ( stringList.c_str() ) );
    return setList;
}


class Structure
{
public:

    // Public typedefs
    typedef RegionMesh<LinearTetra>                                     mesh_Type;
    typedef VectorEpetra                                                vector_Type;
    typedef std::shared_ptr<vector_Type>                              vectorPtr_Type;
    typedef FESpace< RegionMesh<LinearTetra>, MapEpetra >               solidFESpace_Type;
    typedef std::shared_ptr<solidFESpace_Type>                        solidFESpacePtr_Type;

    // typedefs for fibers
    // std function for fiber direction
    typedef std::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > fiberFunction_Type;
    typedef std::shared_ptr<fiberFunction_Type> fiberFunctionPtr_Type;

    typedef std::vector<fiberFunctionPtr_Type>                          vectorFiberFunction_Type;
    typedef std::shared_ptr<vectorFiberFunction_Type>                 vectorFiberFunctionPtr_Type;

    typedef std::vector<vectorPtr_Type>                                 listOfFiberDirections_Type;

    typedef LifeV::ExporterEmpty<mesh_Type >                      emptyExporter_Type;
    typedef std::shared_ptr<emptyExporter_Type>                 emptyExporterPtr_Type;

    typedef LifeV::ExporterEnsight<mesh_Type >                    ensightFilter_Type;
    typedef std::shared_ptr<ensightFilter_Type>                 ensightFilterPtr_Type;

#ifdef HAVE_HDF5
    typedef LifeV::ExporterHDF5<mesh_Type >                       hdf5Filter_Type;
    typedef std::shared_ptr<hdf5Filter_Type>                    hdf5FilterPtr_Type;
#endif

    /** @name Constructors, destructor
     */
    //@{
    Structure ( int                                   argc,
                char**                                argv,
                std::shared_ptr<Epetra_Comm>        structComm );

    ~Structure()
    {}
    //@}

    //@{
    void run()
    {
        run3d();
    }
    //@}

protected:

private:

    /**
     * run the 2D driven cylinder simulation
     */
    void run2d();

    /**
     * run the 3D driven cylinder simulation
     */
    void run3d();

private:
    struct Private;
    std::shared_ptr<Private> parameters;
};



struct Structure::Private
{
    Private() :
        rho (1), poisson (1), young (1), bulk (1), alpha (1), gamma (1)
    {}
    double rho, poisson, young, bulk, alpha, gamma;

    std::string data_file_name;

    std::shared_ptr<Epetra_Comm>     comm;
};



Structure::Structure ( int                                   argc,
                       char**                                argv,
                       std::shared_ptr<Epetra_Comm>        structComm) :
    parameters ( new Private() )
{
    GetPot command_line (argc, argv);
    string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile ( data_file_name );
    parameters->data_file_name = data_file_name;

    parameters->comm = structComm;
    int ntasks = parameters->comm->NumProc();

    if (!parameters->comm->MyPID() )
    {
        std::cout << "My PID = " << parameters->comm->MyPID() << " out of " << ntasks << " running." << std::endl;
    }
}



void
Structure::run2d()
{
    std::cout << "2D cylinder test case is not available yet\n";
}



void
Structure::run3d()
{
    bool verbose = (parameters->comm->MyPID() == 0);

    //! dataElasticStructure
    GetPot dataFile ( parameters->data_file_name.c_str() );

    //Loading a partitoned mesh or reading a new one
    const std::string partitioningMesh = dataFile ( "partitioningOffline/loadMesh", "no");

    //Creation of pointers
    std::shared_ptr<MeshPartitioner<mesh_Type> > meshPart;
    std::shared_ptr<mesh_Type> pointerToMesh;

    if ( ! (partitioningMesh.compare ("no") ) )
    {
        std::shared_ptr<mesh_Type > fullMeshPtr (new mesh_Type ( ( parameters->comm ) ) );
        //Creating a new mesh from scratch
        MeshData             meshData;
        meshData.setup (dataFile, "solid/space_discretization");
        readMesh (*fullMeshPtr, meshData);

        meshPart.reset ( new MeshPartitioner<mesh_Type> (fullMeshPtr, parameters->comm ) );

        pointerToMesh = meshPart->meshPartition();
    }
    else
    {
        //Creating a mesh object from a partitioned mesh
        const std::string partsFileName (dataFile ("partitioningOffline/hdf5_file_name", "NO_DEFAULT_VALUE.h5") );

        std::shared_ptr<Epetra_MpiComm> mpiComm =
            std::dynamic_pointer_cast<Epetra_MpiComm>(parameters->comm);
        PartitionIO<mesh_Type> partitionIO (partsFileName, mpiComm);


        partitionIO.read (pointerToMesh);

    }

    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");

    //Mainly used for BCs assembling (Neumann type)
    solidFESpacePtr_Type dFESpace ( new solidFESpace_Type (pointerToMesh, dOrder, 3, parameters->comm) );
    solidFESpacePtr_Type dScalarFESpace ( new solidFESpace_Type (pointerToMesh, dOrder, 1, parameters->comm) );

    UInt numberOfFibersFamilies =  dataFile ( "solid/model/fibers/numberFamilies" , 0 );
    ASSERT( numberOfFibersFamilies, "This test should be used to check at least one fiber family!! Check the data file, please.")

    // Setting the fibers
    vectorFiberFunctionPtr_Type pointerToVectorOfFamilies( new vectorFiberFunction_Type( ) );
    (*pointerToVectorOfFamilies).resize( numberOfFibersFamilies );

    listOfFiberDirections_Type fiberDirections;
    fiberDirections.resize( numberOfFibersFamilies );

    if( verbose )
        std::cout << "Size of the number of families: " << (*pointerToVectorOfFamilies).size() << std::endl;

    fibersDirectionList setOfFiberFunctions;
    setOfFiberFunctions.setupFiberDefinitions( numberOfFibersFamilies );

    // Setting the vector of fibers functions
    for( UInt k(1); k <= pointerToVectorOfFamilies->size( ); k++ )
    {
        // Setting up the name of the function to define the family
        std::string family="Family";
        // adding the number of the family
        std::string familyNumber;
        std::ostringstream number;
        number << ( k );
        familyNumber = number.str();

        // Name of the function to create
        std::string creationString = family + familyNumber;
        (*pointerToVectorOfFamilies)[ k-1 ].reset( new fiberFunction_Type() );
        (*pointerToVectorOfFamilies)[ k-1 ] = setOfFiberFunctions.fiberDefinition( creationString );

        fiberDirections[ k-1 ].reset( new vector_Type( dFESpace->map() ) );
    }

    vectorPtr_Type thetaSection;
    thetaSection.reset( new vector_Type( dScalarFESpace->map() ) );
    vectorPtr_Type thetaRotation;
    thetaRotation.reset( new vector_Type( dScalarFESpace->map() ) );
    vectorPtr_Type sphereIndicator;
    sphereIndicator.reset( new vector_Type( dScalarFESpace->map() ) );
    vectorPtr_Type positionCenterVector;
    positionCenterVector.reset( new vector_Type( dFESpace->map() ) );
    vectorPtr_Type localPositionVector;
    localPositionVector.reset( new vector_Type( dFESpace->map() ) );

    MPI_Barrier (MPI_COMM_WORLD);

    std::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporter;

    std::string const exporterType =  dataFile ( "exporter/type", "ensight");
    std::string const exportFileName = dataFile ( "exporter/nameFile", "structure");
#ifdef HAVE_HDF5
    if (exporterType.compare ("hdf5") == 0)
    {
        exporter.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, exportFileName ) );
    }
    else
#endif
    {
        if (exporterType.compare ("none") == 0)
        {
            exporter.reset ( new ExporterEmpty<RegionMesh<LinearTetra> > ( dataFile, pointerToMesh, exportFileName, parameters->comm->MyPID() ) );
        }

        else
        {
            exporter.reset ( new ExporterEnsight<RegionMesh<LinearTetra> > ( dataFile, pointerToMesh, exportFileName, parameters->comm->MyPID() ) );
        }
    }

    exporter->setPostDir ( "./" );
    exporter->setMeshProcId ( pointerToMesh, parameters->comm->MyPID() );

    // Adding the fibers vectors
    // Setting the vector of fibers functions
    for( UInt k(1); k <= pointerToVectorOfFamilies->size( ); k++ )
      {
	// Setting up the name of the function to define the family
	std::string family="Family-";
	// adding the number of the family
	std::string familyNumber;
	std::ostringstream number;
	number << ( k );
	familyNumber = number.str();

	// Name of the function to create
	std::string creationString = family + familyNumber;
	exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, creationString, dFESpace, fiberDirections[ k-1 ], UInt (0) );
      }

    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "sphereIndicator", dScalarFESpace, sphereIndicator, UInt (0) );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "thetaSection", dScalarFESpace, thetaSection, UInt (0) );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "thetaRotation", dScalarFESpace, thetaRotation, UInt (0) );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "positionCenter", dFESpace, positionCenterVector, UInt (0) );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "localPosition", dFESpace, localPositionVector, UInt (0) );

    dScalarFESpace->interpolate ( static_cast<solidFESpace_Type::function_Type>( thetaFunction ),
				  *thetaSection,
				  0.0 );

    dScalarFESpace->interpolate ( static_cast<solidFESpace_Type::function_Type>( thetaRotationFunction ),
				  *thetaRotation,
				  0.0 );

    dScalarFESpace->interpolate ( static_cast<solidFESpace_Type::function_Type>( sphereIndicatorFunction ),
				  *sphereIndicator,
				  0.0 );

    dFESpace->interpolate ( static_cast<solidFESpace_Type::function_Type>( positionCenterSpherical ),
				  *positionCenterVector,
				  0.0 );

    dFESpace->interpolate ( static_cast<solidFESpace_Type::function_Type>( localPositionSpherical ),
				  *localPositionVector,
				  0.0 );

    // For loop on the families of fibers
    for( UInt k(1); k <= pointerToVectorOfFamilies->size( ); k++ )
      {
	// Interpolating the i-th fiber family
	dFESpace->interpolate ( static_cast<solidFESpace_Type::function_Type>( *(*pointerToVectorOfFamilies)[k - 1] ),
				*fiberDirections[ k-1 ],
				0.0 );
      }

    exporter->postProcess ( 0.0 );
    cout.precision(16);

    MPI_Barrier (MPI_COMM_WORLD);
}


int
main ( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    std::shared_ptr<Epetra_MpiComm> Comm (new Epetra_MpiComm ( MPI_COMM_WORLD ) );
    if ( Comm->MyPID() == 0 )
    {
        cout << "% using MPI" << endl;
    }
#else
    std::shared_ptr<Epetra_SerialComm> Comm ( new Epetra_SerialComm() );
    cout << "% using serial Version" << endl;
#endif

    Structure structure ( argc, argv, Comm );
    structure.run();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return EXIT_SUCCESS ;
}
