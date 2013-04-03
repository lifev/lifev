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


   Attention: At the moment the restart works only if the solution is saved at
   each time step
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

#include <iostream>


using namespace LifeV;

int returnValue = EXIT_FAILURE;


//Class to investigate the mesh
template < typename MeshEntityType,
         typename ComparisonPolicyType = boost::function2 < bool,
         Real const&,
         Real const& > >
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
         typename ComparisonPolicyType = boost::function2 < bool,
         Real const&,
         Real const& > >
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
                boost::shared_ptr<Epetra_Comm>        structComm );

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
    boost::shared_ptr<Private> parameters;

};



struct Structure::Private
{
    Private() :
        rho (1), poisson (1), young (1), bulk (1), alpha (1), gamma (1)
    {}
    double rho, poisson, young, bulk, alpha, gamma;

    std::string data_file_name;

    boost::shared_ptr<Epetra_Comm>     comm;

};



Structure::Structure ( int                                   argc,
                       char**                                argv,
                       boost::shared_ptr<Epetra_Comm>        structComm) :
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

    MeshData             meshData;
    meshData.setup (dataFile, "solid/space_discretization");

    boost::shared_ptr<RegionMesh<LinearTetra> > fullMeshPtr (new RegionMesh<LinearTetra> ( ( parameters->comm ) ) );
    readMesh (*fullMeshPtr, meshData);

    fullMeshPtr->showMe( );

    if (verbose)
    {
        std::cout << std::endl;
    }

    //Geometrical Infos on the sphere
    Vector3D center (0.0, 0.0335, 1.0);
    Real     radius (0.05);

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

    MPI_Barrier (MPI_COMM_WORLD);
}






int
main ( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_MpiComm> Comm (new Epetra_MpiComm ( MPI_COMM_WORLD ) );
    if ( Comm->MyPID() == 0 )
    {
        cout << "% using MPI" << endl;
    }
#else
    boost::shared_ptr<Epetra_SerialComm> Comm ( new Epetra_SerialComm() );
    cout << "% using serial Version" << endl;
#endif

    Structure structure ( argc, argv, Comm );
    structure.run();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return EXIT_SUCCESS; //Check has to be perfomed later.
}
