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
   \author Paolo Tricerri <paolo.tricerri@epfl.ch>

   This example shows how to select portion of the computational domain
   and change their flag using the MeshUtility features.
   Now the elements for which the flag is changed to 2 are the ones
   belonging to the intersection of the computational domain and the sphere
   defined by its center and

   I have always executed the example in serial since the regionmesh object that
   I pass to the elementlist function is the unpartitioned object.

   \date 2005-04-16
 */
#undef HAVE_HDF5
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
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>


//Include fils which were in the structure.cpp file
#include <lifev/core/array/MapEpetra.hpp>

#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/MeshWriter.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>

#include <lifev/structure/solver/isotropic/ExponentialMaterialNonLinear.hpp>


#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <iostream>



using namespace LifeV;

int returnValue = EXIT_FAILURE;


//Class to investigate the mesh
template < typename MeshEntityType,
         typename ComparisonPolicyType = std::function < bool (
         Real const&,
         Real const&) > >
class SphereInterrogator
{
public:
    typedef MeshEntityType       meshEntity_Type;
    typedef ComparisonPolicyType comparisonPolicy_Type;

    SphereInterrogator ( Vector3D const& center,
                         Real radius,
                         comparisonPolicy_Type const& policy = std::less<Real>() )
        : M_center ( center ),
          M_radius ( radius ),
          M_policy ( policy ) {}

    void operator() ( meshEntity_Type& entity ) const
    {
        // compute the barycenter of the entity
        // (this should be a method of the object)
        Vector3D barycenter;
        for ( UInt k = 0; k < meshEntity_Type::S_numPoints; k++ )
        {
            barycenter += entity.point ( k ).coordinates();
        }
        barycenter /= meshEntity_Type::S_numPoints;

        UInt newMarker (2);
        // check if the distance between the barycenter and the center of the circle
        // satisfies the policy (default: distance less than the radius)
        if ( M_policy ( ( barycenter - M_center ).norm(), M_radius ) )
        {
            entity.setMarkerID ( newMarker );
        }
    }

private:
    const Vector3D M_center;
    const Real M_radius;
    const comparisonPolicy_Type M_policy;

}; // SphereInterrogator

//This class is to count how many volumes are in the sphere. This is for debug purposes
//It cannot be used the same functor as SphereInspector since they perform different operation
template < typename MeshEntityType,
         typename ComparisonPolicyType = std::function < bool (
         Real const&,
         Real const&) > >
class SphereCounter
{
public:
    typedef MeshEntityType       meshEntity_Type;
    typedef ComparisonPolicyType comparisonPolicy_Type;

    SphereCounter ( Vector3D const& center,
                    Real radius,
                    comparisonPolicy_Type const& policy = std::less<Real>() )
        : M_center ( center ),
          M_radius ( radius ),
          M_policy ( policy ) {}

    bool operator() ( const meshEntity_Type& entity ) const
    {
        // compute the barycenter of the entity
        // (this should be a method of the object)
        Vector3D barycenter;
        for ( UInt k = 0; k < meshEntity_Type::S_numPoints; k++ )
        {
            barycenter += entity.point ( k ).coordinates();
        }
        barycenter /= meshEntity_Type::S_numPoints;

        // check if the distance between the barycenter and the center of the circle
        // satisfies the policy (default: distance less than the radius)
        return M_policy ( ( barycenter - M_center ).norm(), M_radius );
    }

private:
    const Vector3D M_center;
    const Real M_radius;
    const comparisonPolicy_Type M_policy;

}; // SphereCounter


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

    typedef RegionMesh<LifeV::LinearTetra >                       mesh_Type;

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
    typedef StructuralOperator<mesh_Type >::vector_Type                 vector_Type;
    typedef boost::shared_ptr<vector_Type>                              vectorPtr_Type;
    typedef FESpace< mesh_Type, MapEpetra >                             solidFESpace_Type;
    typedef boost::shared_ptr<solidFESpace_Type>                        solidFESpacePtr_Type;
    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 >       solidETFESpace_Type;
    typedef boost::shared_ptr<solidETFESpace_Type>                      solidETFESpacePtr_Type;

    typedef LifeV::RegionMesh<LinearTetra>                              mesh_Type;

    // Filters
    typedef LifeV::Exporter<mesh_Type  >                       filter_Type;
    typedef boost::shared_ptr< LifeV::Exporter<mesh_Type  > >           filterPtr_Type;

    typedef LifeV::ExporterEmpty<mesh_Type >                            emptyFilter_Type;
    typedef boost::shared_ptr<emptyFilter_Type>                         emptyFilterPtr_Type;
    typedef LifeV::ExporterEnsight<mesh_Type >                          ensightFilter_Type;
    typedef boost::shared_ptr<ensightFilter_Type>                       ensightFilterPtr_Type;

#ifdef HAVE_HDF5
    typedef LifeV::ExporterHDF5<mesh_Type >                             hdf5Filter_Type;
    typedef boost::shared_ptr<hdf5Filter_Type>                          hdf5FilterPtr_Type;
#endif

    bool verbose = (parameters->comm->MyPID() == 0);

    //! dataElasticStructure
    GetPot dataFile ( parameters->data_file_name.c_str() );
    boost::shared_ptr<StructuralConstitutiveLawData> dataStructure (new StructuralConstitutiveLawData( ) );
    dataStructure->setup (dataFile);

    MeshData             meshData;
    meshData.setup (dataFile, "solid/space_discretization");


    std::shared_ptr<RegionMesh<LinearTetra> > fullMeshPtr (new RegionMesh<LinearTetra> ( ( parameters->comm ) ) );
    std::shared_ptr<RegionMesh<LinearTetra> > localMeshPtr (new RegionMesh<LinearTetra> ( ( parameters->comm ) ) );

    readMesh (*fullMeshPtr, meshData);

    //fullMeshPtr->showMe( );

    if (verbose)
    {
        std::cout << std::endl;
    }

    //Geometrical Infos on the sphere
    Vector3D center (0.138, 0.0, 0.138);
    Real     radius (0.1);

    //Count how many volumes are in the sphere
    //Create the Predicate
    SphereCounter<mesh_Type::element_Type> countVolumesFunctor ( center, radius );
    //Count phase
    UInt numExtractedVolumes = fullMeshPtr->elementList().countAccordingToPredicate ( countVolumesFunctor );

    std::cout << " The Number of volumes inside the sphere is: " << numExtractedVolumes << std::endl;

    //Extract the list of elements inside the sphere.
    //Creation of the class to investigate the mesh
    SphereInterrogator<mesh_Type::element_Type> setSphereInterrogator ( center, radius );

    //This method changes the markerID for the volumes which are inside the defined sphere
    fullMeshPtr->elementList().changeAccordingToFunctor ( setSphereInterrogator );

    //Exporting the modified mesh
    std::string const nameExportMesh =  dataFile ( "exporter/modifiedMesh", "damagedMesh");
    MeshWriter::writeMeshMedit<RegionMesh<LinearTetra> > ( nameExportMesh , *fullMeshPtr );


    // Exporting the mesh to check region with changed flag
    MeshPartitioner< mesh_Type > meshPart ( fullMeshPtr, parameters->comm );
    localMeshPtr = meshPart.meshPartition();

    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");
    solidFESpacePtr_Type dFESpace ( new solidFESpace_Type (localMeshPtr, dOrder, 3, parameters->comm) );
    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (meshPart, & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), parameters->comm) );

    //! 1. Constructor of the class to compute the tensions
    StructuralOperator<RegionMesh<LinearTetra> >  solid;

    //! face BChandler object
    boost::shared_ptr<BCHandler> BCh ( new BCHandler() );

    //! 2. Its setup
    solid.setup (dataStructure,
                 dFESpace,
                 dETFESpace,
                 BCh,
                 parameters->comm);

    //! 6. Post-processing setting
    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporter;

    std::string const exporterType =  dataFile ( "exporter/type", "hdf5");
    std::string const nameExporter =  dataFile ( "exporter/name", "");

#ifdef HAVE_HDF5
    if (exporterType.compare ("hdf5") == 0)
    {
        exporter.reset ( new hdf5Filter_Type ( dataFile, nameExporter ) );
    }
    else
#endif
    {
        if (exporterType.compare ("none") == 0)
        {
            exporter.reset ( new emptyFilter_Type ( dataFile, meshPart.meshPartition(), nameExporter, parameters->comm->MyPID() ) ) ;
        }

        else
        {
            exporter.reset ( new ensightFilter_Type ( dataFile, meshPart.meshPartition(), nameExporter, parameters->comm->MyPID() ) ) ;
        }
    }

    exporter->setMeshProcId (dFESpace->mesh(), dFESpace->map().comm().MyPID() );

    vectorPtr_Type meshColors ( new vector_Type (solid.displacement(),  LifeV::Unique ) );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "colors", dFESpace, meshColors, UInt (0) );

    exporter->postProcess ( 0.0 );
    MPI_Barrier (MPI_COMM_WORLD);

    //color the mesh according to the marker of the volume
    solid.colorMesh ( *meshColors );

    exporter->postProcess ( 1000.0 );

    //Closing files
    exporter->closeFile();

    if (verbose )
    {
        std::cout << "finished" << std::endl;
    }

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
    return EXIT_SUCCESS; //Check has to be perfomed later.
}
