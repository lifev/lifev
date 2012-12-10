/*
 * run_MeshExtractor.cpp
 *
 *      Author: uvilla
 *
 *      Extract a 2d dimensional mesh from the boundary of a 3d mesh.
 */

#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/util/LifeChrono.hpp>


#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>

#include "MeshExtractor.hpp"

int run(GetPot & dataFile, bool /*verbose*/, boost::shared_ptr<Epetra_Comm>& comm)
{
    using namespace LifeV;
    typedef RegionMesh<LinearTetra> mesh_Type;
    typedef RegionMesh<LinearTriangle> mesh2d_Type;

    // LifeChrono chrono;

    std::string mesh_section("mesh");
    //FIXME At the moment we can extract only one marker at a time.
    //Only the first marker in boundaryMarkerListToExtract will be extracted
    std::list<UInt> boundaryMarkerListToExtract, otherBoundaryMarkerList;
    parseList( dataFile( (mesh_section+"/interfaceList").c_str(), ""), boundaryMarkerListToExtract );
    parseList( dataFile( (mesh_section+"/boundaryList").c_str(), ""), otherBoundaryMarkerList );


    //=======================================================================//
    // Mesh Stuff                                                            //
    //=======================================================================//
    // Read the 3d mesh
    boost::shared_ptr<mesh_Type> mesh;
    boost::shared_ptr< MeshPartitioner<mesh_Type> > meshPart;
    MeshData meshData(dataFile, mesh_section);
    mesh.reset( new mesh_Type( comm ) );
    readMesh(*mesh, meshData);

    //Extract the 2d mesh from the boundary with a given marker
    boost::shared_ptr<mesh2d_Type> mesh2d;
    boost::shared_ptr< MeshPartitioner<mesh2d_Type> > mesh2dPart;
    mesh2d.reset(extractBoundaryMesh(*mesh, *boundaryMarkerListToExtract.begin(), otherBoundaryMarkerList));
    mesh2d->showMe(false, std::cout);

    meshPart.reset(new MeshPartitioner<mesh_Type>(mesh, comm));
    mesh2dPart.reset(new MeshPartitioner<mesh2d_Type>(mesh2d, comm)); //Fails in debug mode!!!

    //========================================================================//
    // Create a dummy feSpace on the boundary mesh                            //
    //========================================================================//
    boost::shared_ptr< FESpace<mesh2d_Type, MapEpetra> > feSpace;
    feSpace.reset(new FESpace<mesh2d_Type, MapEpetra>(mesh2dPart->meshPartition(), "P1", 1, comm));

    //========================================================================//
    // post processing setup                                                  //
    //========================================================================//
    boost::shared_ptr<Exporter<mesh2d_Type> > exporter;
    std::string const exporterType =  dataFile( "exporter/type", "hdf5");

#ifdef HAVE_HDF5
    if (exporterType.compare("hdf5") == 0)
    {
        exporter.reset( new ExporterHDF5<mesh2d_Type > ( dataFile, "2dmesh" ) );
    }
    else
#endif
    {
        if (exporterType.compare("none") == 0)
        {
            exporter.reset( new ExporterEmpty<mesh2d_Type > ( dataFile, mesh2d, "2dmesh", comm->MyPID()) );
        }
        else
        {
            exporter.reset( new ExporterEnsight<mesh2d_Type > ( dataFile, mesh2d, "2dmesh", comm->MyPID()) );
        }
    }

    exporter->setPostDir( "./" );
    exporter->setMeshProcId( mesh2d, comm->MyPID() );


    //====================================================================//
    // Show the extracted mesh, color element according to the PID        //
    //====================================================================//
    boost::shared_ptr<VectorEpetra> u(new VectorEpetra(feSpace->map(), exporter->mapType()));
    exporter->addVariable(ExporterData<mesh2d_Type >::ScalarField, "u", feSpace, u, UInt(0));
    *u = comm->MyPID();
    exporter->postProcess(0);



    return 0;
}
