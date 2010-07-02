//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
  @file
  @brief This file provides the class  Hdf5exporter for post-processing with hdf5

  @author Simone Deparis <simone.deparis@epfl.ch>
  @author Radu Popescu <radu.popescu@epfl.ch>
  @date 11-11-2008
 */

#ifndef HDF5EXPORTER_H
#define HDF5EXPORTER_H 1

#include <vector>

#include <lifeconfig.h>

#ifndef HAVE_HDF5
#warning warning you should reconfigure with --with-hdf5=... flag

#else
#include <life/lifecore/util_string.hpp>
#include <life/lifefilters/exporter.hpp>
#include <EpetraExt_HDF5.h>
#include <Epetra_MultiVector.h>
#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>
#include <EpetraExt_DistArray.h>
#include <Epetra_IntVector.h>
#include <life/lifealg/EpetraMap.hpp>
#include <life/lifefem/refFE.hpp>
#include <Epetra_Comm.h>
#include <sstream>
#include <string>
#include <map>

namespace LifeV {

//! Hdf5 data exporter, implementation of Exporter
/*!
  @author Simone Deparis <simone.deparis@epfl.ch>
  @author Radu Popescu <radu.popescu@epfl.ch>

  Usage: two steps
  <ol>
    <li> first: add the variables using addVariable
    <li> second: call postProcess( time );
  </ol>

  Guidelines when operating with partition information:
  <ol>
    <li> Use different objects for import and export
    <li> When exporting graph and partition(s) create a Hdf5exporter to be used for output.
         Call methods: addPartitionGraph, addMeshPartitionAll, or addMyMeshPartition
         Call postProcess to write data
         Call CloseFile
    <li> When importing graph and/or partition create a Hdf5exporter to be used for input,
         Call loadGraph and partitionMesh::update, call loadMyPartition
         Call CloseFile
  </ol>
*/
template<typename Mesh>
class Hdf5exporter : public Exporter<Mesh> {

public:

    typedef Exporter<Mesh> super;
    typedef typename super::mesh_ptrtype mesh_ptrtype;
    typedef typename super::vector_rawtype vector_type;
    typedef typename super::vector_ptrtype vector_ptrtype;

    typedef EpetraExt::HDF5 hdf5_type;
    typedef boost::shared_ptr<hdf5_type> hdf5_ptrtype;
    typedef std::vector<std::vector<int> > graph_type;
    typedef boost::shared_ptr<graph_type> graph_ptrtype;
    typedef boost::shared_ptr<std::vector<mesh_ptrtype> > serial_mesh_ptrtype;


    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor for Hdf5exporter
    Hdf5exporter();

    //! Constructor for Hdf5exporter
    /*!
       @param dfile the GetPot data file where you must provide an [exporter] section with:
          "start"     (start index for sections in the hdf5 data structure 0 for 000, 1 for 001 etc.),
          "save"      (how many time steps per postprocessing)
          "multimesh" ( = true if the mesh has to be saved at each post-processing step)
       @param mesh the mesh
       @param the prefix for the case file (ex. "test" for test.case)
       @param the procId determines de CPU id. if negative, it ussemes there is only one processor
    */
    Hdf5exporter(const GetPot& dfile, mesh_ptrtype mesh, const std::string& prefix, const int& procId);

    //! Constructor for Hdf5exporter without prefix and procID
    /*!
       @param dfile the GetPot data file where you must provide an [exporter] section with:
          "start"     (start index for sections in the hdf5 data structure 0 for 000, 1 for 001 etc.),
          "save"      (how many time steps per postprocessing)
          "multimesh" ( = true if the mesh has to be saved at each post-processing step)
       @param mesh the mesh
    */
    Hdf5exporter(const GetPot& dfile, const std::string& prefix);

    //! Destructor for Hdf5exporter
    ~Hdf5exporter() {}

    //@}


    //! @name Methods
    //@{
    //! Add the partition graph to the post processing data file
    /*!
      Add the partition graph to the post processing data file.
      \param graph - shared_ptr<vector<vector<int> > > - shared pointer to the partition graph data structure
      (as returned by partitionMesh::graph() )
      \param comm - Epetra_Comm* - raw pointer to the Epetra communicator to be used
    */
    void addPartitionGraph(const graph_ptrtype& graph, Epetra_Comm* comm)
    {M_graph = graph; M_comm = comm;}

    //! Add all of the mesh partitions to the post processing data file (serial operation)
    /*!
      Add all of the mesh partitions to the post processing data file.
      \param meshPointer - shared_ptr<vector<shared_ptr<Mesh> > > - shared pointer to the vector storing
      pointers to the mesh partitions (as returned by partitionMesh::meshAllPartitions() )
      \param comm - Epetra_Comm* - raw pointer to the Epetra communicator to be used
    */
    void addMeshPartitionAll(const serial_mesh_ptrtype& meshPointer, Epetra_Comm* comm)
    {M_serialMesh = meshPointer; M_parallelMesh.reset(); M_comm = comm;}

    //! Add to HDF5 file the mesh partition that belongs to the current process (parallel operation)
    /*!
      After the mesh partition is loaded from the HDF5 file, the simulation is run and Hdf5exporter::postProcess()
      is called, the original contents of the HDF5 will be lost. To keep mesh partition, call
      Hdf5exporter::addMyMeshPartition() before calling Hdf5exporter::postProcess();
      \param meshPointer - shared_ptr<Mesh> - shared pointer to a mesh partition (as returned by
      partitionMesh::mesh() )
      \param comm - Epetra_Comm* - raw pointer to the Epetra communicator to be used
    */
    void addMyMeshPartition(const mesh_ptrtype& meshPointer, Epetra_Comm* comm)
    {/*M_parallelMesh = meshPointer; M_serialMesh.reset(); M_comm = comm;*/}

    //! Load the partitioned graph from a HDF5 file into a partitionMesh object
    /*!
      \param graph - shared_ptr<vector<vector<int> > > - a shared pointer to the graph data structure in the
      partitionMesh object (as returned by partitionMesh::graph() )
      \param comm - Epetra_Comm* - a raw pointer to the Epetra communicator to be used
    */
    void loadGraph(graph_ptrtype graph, Epetra_Comm* comm);

    //! Load a mesh partition according to the MPI PID
    /*!
      This method is to be used in parallel. The mesh partition corresponding to the
      MPI process id is loaded into the meshPartition object.
      Could be used with any object that contains a RegionMesh3D member but it is
      primarily used to load a mesh partition into the partitionMesh object.
      \param meshPartition - shared_ptr<Mesh> - shared pointer to mesh partition object
      \param comm -Epetra_Comm* - raw pointer to the Epetra communicator to be used
    */
    void loadMyPartition(mesh_ptrtype meshPartition, Epetra_Comm* comm);

    //! Post-process the variables added to the list
    /*!
        @param time the solver time
    */
    void postProcess(const Real& time);

    //! Import data from previous simulations at a certain time
    /*!
       @param Time the time of the data to be imported
       @return number of iteration corresponding at the time step
     */
    UInt importFromTime( const Real& Time );

    //! Import data from previous simulations at a certain time
    /*!
       @param Time the time of the data to be imported
       @return the simulation time corresponding to the iteration
     */
    double importFromIter( const UInt& );

    //! Import data from previous simulations
    /*!
       @param time the solver time
    */
    void import(const Real& /*Tstart*/, const Real& /*dt*/); // dt is used to rebuild the history up to now

    //! Import data from previous simulations
    /*!
       @param time the solver time
    */
    void import(const Real& /*Tstart*/);

    //! Close the Hdf5 file
    /*!
         Close the HDF5 file.
     */
    void CloseFile() {M_HDF5->Close();}

    //@}


    //! @name Get Methods
    //@{

    //! returns the type of the map to use for the EpetraVector
    EpetraMapType mapType() const;

    void rd_var( ExporterData& dvar);

    //@}

private:

    //! @name Private Methods
    //@{

    //! Define the shape of the elements
    void defineShape();
    //! write empty xdmf file
    void M_wr_initXdmf();
    //! append to xdmf file
    void M_wr_Xdmf(const Real& time);
    //! save position and write closing lines
    void M_wr_closeLinesXdmf();
    //! remove closing lines
    void M_wr_removeCloseLinesXdmf();

    void M_wr_topology  ( std::ofstream& xdmf );
    void M_wr_geometry  ( std::ofstream& xdmf );
    void M_wr_attributes( std::ofstream& xdmf );
    void M_wr_scalar_datastructure  ( std::ofstream& xdmf, const ExporterData& dvar );
    void M_wr_vector_datastructure  ( std::ofstream& xdmf, const ExporterData& dvar );

    void M_wr_var(const ExporterData& dvar);
    void M_wr_scalar(const ExporterData& dvar);
    void M_wr_vector(const ExporterData& dvar);

    void M_wr_geo();

    void M_rd_scalar( ExporterData& dvar);
    void M_rd_vector( ExporterData& dvar);

    // The following private methods are for writing the partitioned graph
    // and mesh to the output file
    void writeGraph();
    void writePartition(mesh_ptrtype partition, std::string& suffix);
    void writeParallelMesh();
    void writeSerialMesh();

    //@}

    hdf5_ptrtype      M_HDF5;
    std::ofstream     M_xdmf;

    const std::string M_closingLines;
    std::streampos    M_closingLinesPosition;
    std::string       M_outputFileName;

    serial_mesh_ptrtype   M_serialMesh;
    mesh_ptrtype          M_parallelMesh;
    graph_ptrtype         M_graph;
    Epetra_Comm*          M_comm;

};



// ===================================================
// Constructors
// ===================================================
template<typename Mesh>
Hdf5exporter<Mesh>::Hdf5exporter():
    super               (),
    M_HDF5              (),
    M_closingLines      ( "\n    </Grid>\n\n  </Domain>\n</Xdmf>\n"),
    M_outputFileName    ( "noninitialisedFileName" )
{
}

template<typename Mesh>
Hdf5exporter<Mesh>::Hdf5exporter(const GetPot& dfile, mesh_ptrtype mesh, const std::string& prefix,
                                 const int& procId) :
    super               ( dfile, prefix ),
    M_HDF5              (),
    M_closingLines      ( "\n    </Grid>\n\n  </Domain>\n</Xdmf>\n"),
    M_outputFileName    ( "noninitialisedFileName" )
{
    setMeshProcId( mesh, procId );
}

template<typename Mesh>
Hdf5exporter<Mesh>::Hdf5exporter(const GetPot& dfile, const std::string& prefix):
    super               ( dfile, prefix ),
    M_HDF5              (),
    M_closingLines      ( "\n    </Grid>\n\n  </Domain>\n</Xdmf>\n"),
    M_outputFileName    ( "noninitialisedFileName" )
{
}

// ===================================================
// Methods
// ===================================================
template<typename Mesh>
void Hdf5exporter<Mesh>::loadGraph(graph_ptrtype graph, Epetra_Comm *comm)
{
    if (M_HDF5.get() == 0)
    {
        M_HDF5.reset(new hdf5_type(*comm));
    }
    if (! M_HDF5->IsOpen())
    {
        M_HDF5->Open(this->M_post_dir + this->M_prefix + ".h5", H5F_ACC_RDONLY);
    }

    int nPartitions;

    M_HDF5->Read("Graph", "number_partitions", nPartitions);

    std::vector<int> partitionSizes(nPartitions);
    M_HDF5->Read("Graph", "partition_sizes", H5T_NATIVE_INT, nPartitions,
                 &partitionSizes[0]);

    graph->resize(0);
    graph->reserve(nPartitions);

    std::vector<int> partBuffer;
    std::stringstream index;

    for (UInt i = 0; i < nPartitions; ++i)
    {
        partBuffer.resize(partitionSizes[i]);
        index << i;
        M_HDF5->Read("Graph", "partition_graph_" + index.str(),
                     H5T_NATIVE_INT, partitionSizes[i],
                     &partBuffer[0]);
        graph->push_back(partBuffer);
        index.str(std::string());
        index.clear();
    }
}

template<typename Mesh>
void Hdf5exporter<Mesh>::loadMyPartition(mesh_ptrtype meshPartition, Epetra_Comm* comm)
{
    UInt elementNodes, faceNodes;
    switch (Mesh::ElementShape::Shape)
    {
    case HEXA:
        elementNodes = 8;
        faceNodes    = 4;
        break;
    case TETRA:
        elementNodes = 4;
        faceNodes    = 3;
    }

    if (M_HDF5.get() == 0)
    {
        M_HDF5.reset(new hdf5_type(*comm));
    }
    if (! M_HDF5->IsOpen())
    {
        M_HDF5->Open(this->M_post_dir + this->M_prefix + ".h5", H5F_ACC_RDONLY);
    }

    std::stringstream index;
    index << comm->MyPID();
    std::string suffix = "." + index.str();

    // Read counters from file and set them in the mesh

    // Points
    int numPoints;
    M_HDF5->Read("Mesh", "Counters.NumPoints" + suffix, numPoints);
    meshPartition->setMaxNumPoints(numPoints, true);
    int numBPoints;
    M_HDF5->Read("Mesh", "Counters.NumBPoints" + suffix, numBPoints);
    meshPartition->setNumBPoints(numBPoints);

    // Vertices
    int numVertices;
    M_HDF5->Read("Mesh", "Counters.NumVertices" + suffix, numVertices);
    meshPartition->setNumVertices(numVertices);
    int numBVertices;
    M_HDF5->Read("Mesh", "Counters.NumBVertices" + suffix, numBVertices);
    meshPartition->setNumBVertices(numBVertices);
    int numGlobalVertices;
    M_HDF5->Read("Mesh", "Counters.NumGlobalVertices" + suffix, numGlobalVertices);
    meshPartition->setNumGlobalVertices(numGlobalVertices);

    // Edges
    int numEdges;
    M_HDF5->Read("Mesh", "Counters.NumEdges" + suffix, numEdges);
    meshPartition->setNumEdges(numEdges);
    meshPartition->setMaxNumEdges(numEdges);
    int numBEdges;
    M_HDF5->Read("Mesh", "Counters.NumBEdges" + suffix, numBEdges);
    meshPartition->setNumBEdges(numBEdges);
    int numGlobalEdges;
    M_HDF5->Read("Mesh", "Counters.NumGlobalEdges" + suffix, numGlobalEdges);
    meshPartition->setMaxNumGlobalEdges(numGlobalEdges);

    // Faces
    int numFaces;
    M_HDF5->Read("Mesh", "Counters.NumFaces" + suffix, numFaces);
    meshPartition->setNumFaces(numFaces);
    meshPartition->setMaxNumFaces(numFaces);
    int numBFaces;
    M_HDF5->Read("Mesh", "Counters.NumBFaces" + suffix, numBFaces);
    meshPartition->setNumBFaces(numBFaces);
    int numGlobalFaces;
    M_HDF5->Read("Mesh", "Counters.NumGlobalFaces" + suffix, numGlobalFaces);
    meshPartition->setMaxNumGlobalFaces(numGlobalFaces);

    // Volumes
    int numVolumes;
    M_HDF5->Read("Mesh", "Counters.NumVolumes" + suffix, numVolumes);
    meshPartition->setMaxNumVolumes(numVolumes, true);
    int numGlobalVolumes;
    M_HDF5->Read("Mesh", "Counters.NumGlobalVolumes" + suffix, numGlobalVolumes);
    meshPartition->setMaxNumGlobalVolumes(numGlobalVolumes);

    // Read the list of points
    std::vector<double> tmpVectorDouble(numPoints);
    std::vector<std::vector<double> > pointCoordinates(3, tmpVectorDouble);

    std::vector<int> pointMarkers(numPoints);
    std::vector<int> pointBoundaryFlags(numPoints);
    std::vector<int> pointGlobalId(numPoints);

    M_HDF5->Read("Mesh", "Points.x" + suffix, H5T_NATIVE_DOUBLE, numPoints, &pointCoordinates[0][0]);
    M_HDF5->Read("Mesh", "Points.y" + suffix, H5T_NATIVE_DOUBLE, numPoints, &pointCoordinates[1][0]);
    M_HDF5->Read("Mesh", "Points.z" + suffix, H5T_NATIVE_DOUBLE, numPoints, &pointCoordinates[2][0]);
    M_HDF5->Read("Mesh", "Points.f" + suffix, H5T_NATIVE_INT, numPoints, &pointMarkers[0]);
    M_HDF5->Read("Mesh", "Points.BoundaryFlag" + suffix, H5T_NATIVE_INT, numPoints, &pointBoundaryFlags[0]);
    M_HDF5->Read("Mesh", "Points.GlobalId" + suffix, H5T_NATIVE_INT, numPoints, &pointGlobalId[0]);

    meshPartition->pointList.reserve(numPoints);
    meshPartition->_bPoints.reserve(meshPartition->numBPoints());

    typename Mesh::PointType *pp = 0;

    for (UInt j = 0; j < numPoints; ++j)
    {
        pp = &(meshPartition->addPoint(bool(pointBoundaryFlags[j])));
        pp->setMarker(pointMarkers[j]);
        pp->x() = pointCoordinates[0][j];
        pp->y() = pointCoordinates[1][j];
        pp->z() = pointCoordinates[2][j];
        pp->setLocalId(j + 1);
        pp->setId(pointGlobalId[j]);

        meshPartition->localToGlobalNode().insert(std::make_pair(j + 1, pointGlobalId[j]));
        meshPartition->globalToLocalNode().insert(std::make_pair(pointGlobalId[j], j + 1));
    }

    pointCoordinates.clear();
    pointMarkers.clear();
    pointBoundaryFlags.clear();
    pointGlobalId.clear();

    // Read the list of edges
    std::vector<int> tmpVectorInt(numEdges);
    std::vector<std::vector<int> > edgePoints(2, tmpVectorInt);

    std::vector<int> edgeMarkers(numEdges);
    std::vector<int> edgeGlobalId(numEdges);
    std::vector<int> edgeBoundaryFlags(numEdges);

    M_HDF5->Read("Mesh", "Edges.p1" + suffix, H5T_NATIVE_INT, numEdges, &edgePoints[0][0]);
    M_HDF5->Read("Mesh", "Edges.p2" + suffix, H5T_NATIVE_INT, numEdges, &edgePoints[1][0]);
    M_HDF5->Read("Mesh", "Edges.f" + suffix, H5T_NATIVE_INT, numEdges, &edgeMarkers[0]);
    M_HDF5->Read("Mesh", "Edges.GlobalId" + suffix, H5T_NATIVE_INT, numEdges, &edgeGlobalId[0]);
    M_HDF5->Read("Mesh", "Edges.BoundaryFlag" + suffix, H5T_NATIVE_INT, numEdges, &edgeBoundaryFlags[0]);

    meshPartition->edgeList.reserve(numEdges);

    typename Mesh::EdgeType *pe;

    for (UInt j = 0; j < numEdges; ++j)
    {
        pe = &(meshPartition->addEdge(edgeBoundaryFlags[j]));
        pe->setLocalId(j + 1);
        pe->setId(edgeGlobalId[j]);
        pe->setPoint(1, meshPartition->point(edgePoints[0][j]));
        pe->setPoint(2, meshPartition->point(edgePoints[1][j]));
        pe->setMarker(edgeMarkers[j]);
    }

    edgePoints.clear();
    edgeMarkers.clear();
    edgeGlobalId.clear();
    edgeBoundaryFlags.clear();

    // Read the list of faces
    tmpVectorInt.resize(numFaces);
    std::vector<std::vector<int> > facePoints(faceNodes, tmpVectorInt);

    std::vector<int> faceMarkers(numFaces);
    std::vector<int> faceGlobalId(numFaces);
    std::vector<int> faceBoundaryFlags(numFaces);

    std::vector<std::vector<int> > faceNeighbourId(2, tmpVectorInt);
    std::vector<std::vector<int> > faceNeighbourPos(2, tmpVectorInt);

    std::stringstream idx;
    for (UInt k = 0; k < faceNodes; ++k)
    {
        idx << k + 1;
        M_HDF5->Read("Mesh", "Faces.p" + idx.str() + suffix, H5T_NATIVE_INT, numFaces, &facePoints[k][0]);
        idx.str(std::string());
        idx.clear();
    }
    M_HDF5->Read("Mesh", "Faces.f" + suffix, H5T_NATIVE_INT, numFaces, &faceMarkers[0]);
    M_HDF5->Read("Mesh", "Faces.GlobalId" + suffix, H5T_NATIVE_INT, numFaces, &faceGlobalId[0]);
    M_HDF5->Read("Mesh", "Faces.BoundaryFlag" + suffix, H5T_NATIVE_INT, numFaces, &faceBoundaryFlags[0]);

    M_HDF5->Read("Mesh", "Faces.NeighbourId1" + suffix, H5T_NATIVE_INT, numFaces, &faceNeighbourId[0][0]);
    M_HDF5->Read("Mesh", "Faces.NeighbourId2" + suffix, H5T_NATIVE_INT, numFaces, &faceNeighbourId[1][0]);
    M_HDF5->Read("Mesh", "Faces.NeighbourPos1" + suffix, H5T_NATIVE_INT, numFaces, &faceNeighbourPos[0][0]);
    M_HDF5->Read("Mesh", "Faces.NeighbourPos2" + suffix, H5T_NATIVE_INT, numFaces, &faceNeighbourPos[1][0]);


    typename Mesh::FaceType *pf = 0;

    meshPartition->faceList.reserve(numFaces);

    for (UInt j = 0; j < numFaces; ++j)
    {
        pf = &(meshPartition->addFace(faceBoundaryFlags[j]));
        pf->setLocalId(j + 1);
        pf->setId(faceGlobalId[j]);

        pf->ad_first() = faceNeighbourId[0][j];
        pf->ad_second() = faceNeighbourId[1][j];
        pf->pos_first() = faceNeighbourPos[0][j];
        pf->pos_second() = faceNeighbourPos[1][j];

        pf->setMarker(faceMarkers[j]);
        for (UInt k = 0; k < faceNodes; ++k)
        {
            pf->setPoint(k + 1, meshPartition->point(facePoints[k][j]));
        }
    }

    meshPartition->setLinkSwitch("HAS_ALL_FACES");
    meshPartition->setLinkSwitch("FACES_HAVE_ADIACENCY");

    facePoints.clear();
    faceMarkers.clear();
    faceGlobalId.clear();
    faceBoundaryFlags.clear();
    faceNeighbourId.clear();
    faceNeighbourPos.clear();

    // Read the list of volumes
    tmpVectorInt.resize(numVolumes);
    std::vector<std::vector<int> > volumePoints(elementNodes, tmpVectorInt);

    std::vector<int> volumeMarkers(numVolumes);
    std::vector<int> volumeGlobalId(numVolumes);

    for (UInt k = 0; k < elementNodes; ++k)
    {
        idx << k + 1;
        M_HDF5->Read("Mesh", "Volumes.p" + idx.str() + suffix, H5T_NATIVE_INT, numVolumes, &volumePoints[k][0]);
        idx.str(std::string());
        idx.clear();
    }
    M_HDF5->Read("Mesh", "Volumes.f" + suffix, H5T_NATIVE_INT, numVolumes, &volumeMarkers[0]);
    M_HDF5->Read("Mesh", "Volumes.GlobalId" + suffix, H5T_NATIVE_INT, numVolumes, &volumeGlobalId[0]);

    meshPartition->volumeList.reserve(numVolumes);

    typename Mesh::VolumeType *pv = 0;

    for (UInt j = 0; j < numVolumes; ++j)
    {
        pv = &(meshPartition->addVolume());
        pv->setId(volumeGlobalId[j]);
        pv->setLocalId(j + 1);
        for (UInt k = 0; k < elementNodes; ++k)
        {
            pv->setPoint(k + 1, meshPartition->point(volumePoints[k][j]));
        }
        pv->setMarker(volumeMarkers[j]);
    }

    volumePoints.clear();
    volumeMarkers.clear();
    volumeGlobalId.clear();

    meshPartition->updateElementEdges(false, false);
    meshPartition->updateElementFaces(false, false);

}

template<typename Mesh>
void Hdf5exporter<Mesh>::postProcess(const Real& time)
{
    if ( M_HDF5.get() == 0)
    {
        if (this->M_listData.size() != 0)
        {
            M_HDF5.reset(new hdf5_type(this->M_listData.begin()->storedArray()->Comm()));
        }
        else
        {
            M_HDF5.reset(new hdf5_type(*M_comm));
        }
        M_outputFileName=this->M_prefix+".h5";
        M_HDF5->Create(this->M_post_dir+M_outputFileName);

        // write empty xdmf file
        M_wr_initXdmf();

        if (!this->M_multimesh)
        {
            if (this->M_listData.size() != 0)
            {
                M_wr_geo(); // see also M_wr_geometry
            }

            if (M_graph.get() != 0)
            {
                if (M_comm->MyPID() == 0)
                {
                    // Write the partition graph to file
                    writeGraph();
                }
            }
            if (M_serialMesh.get() != 0)
            {
                // Write all the mesh partitions to file
                writeSerialMesh();
            }
            if (M_parallelMesh.get() != 0)
            {
                // Write mesh partition that belongs to you (parallel op)
                writeParallelMesh();
            }

            M_HDF5->Flush();
        }
    }

    typedef std::list< ExporterData >::const_iterator Iterator;

    this->computePostfix();

    if ( this->M_postfix != "*****" )
    {
        if (!this->M_procId) std::cout << "  x-  HDF5 post-processing ...        " << std::flush;
        Chrono chrono;
        chrono.start();
        for (Iterator i=this->M_listData.begin(); i != this->M_listData.end(); ++i)
        {
            M_wr_var(*i);
        }
        // pushing time
        this->M_timeSteps.push_back(time);

        M_wr_Xdmf(time);

        if (this->M_multimesh) {
            M_wr_geo(); // see also M_wr_geometry
        }

        chrono.stop();

        // Write to file without closing the file
        M_HDF5->Flush();

        if (!this->M_procId) std::cout << "         done in " << chrono.diff() << " s." << std::endl;
    }
}

template<typename Mesh>
UInt Hdf5exporter<Mesh>::importFromTime( const Real& Time )
{
    // Container for the time and the postfix
    std::pair< Real, int > SelectedTimeAndPostfix;
    if ( !this->M_procId )
    {
        // Open the xmf file
        std::ifstream xmfFile;
        xmfFile.open( ( this->M_post_dir + this->M_prefix + ".xmf" ).c_str(), std::ios::in );

        // Vector of TimeStep
        std::vector< std::pair< Real, int > > TimeAndPostfix;
        if ( xmfFile.is_open() )
        {
            // Define some variables
            std::string line;
            std::vector<std::string> stringsVector;

            // Read one-by-one all the lines of the file
            while ( !xmfFile.eof() )
            {
                std::getline( xmfFile, line, '\n' );

                // If the line begin with "<!-- Time " it is the beginning of a new block
                if ( !line.compare( 0, 10, "<!-- Time " ) )
                {
                    boost::split( stringsVector, line, boost::is_any_of( " " ) );
                    TimeAndPostfix.push_back( make_pair( string2number( stringsVector[2] ), string2number( stringsVector[4] ) ) );
                }
            }
        }
        xmfFile.close();

        // Find the closest time step
        SelectedTimeAndPostfix = TimeAndPostfix.front();
        for ( std::vector< std::pair< Real, int > >::const_iterator i = TimeAndPostfix.begin(); i < TimeAndPostfix.end() ; ++i )
            if ( std::abs( SelectedTimeAndPostfix.first - Time ) >= std::abs( (*i).first - Time ) )
                SelectedTimeAndPostfix = *i;

        //std::cout << "  x-  HDF5 import from time " << SelectedTimeAndPostfix.first << " iteration " << SelectedTimeAndPostfix.second << std::endl;
    }

    this->M_listData.begin()->storedArray()->Comm().Broadcast( &SelectedTimeAndPostfix.second, 1, 0 );
    this->M_count = SelectedTimeAndPostfix.second;
    this->computePostfix();

    // Importing
    if ( !this->M_procId )
        std::cout << "  x-  HDF5 importing ...                       "<< std::flush;

    Chrono chrono;
    chrono.start();
    for ( std::list< ExporterData >::iterator i=this->M_listData.begin(); i != this->M_listData.end(); ++i )
        this->rd_var(*i);

    chrono.stop();
    if ( !this->M_procId )
        std::cout << "done in " << chrono.diff() << " s. (Time " << SelectedTimeAndPostfix.first
                                                 << ", Iteration " << SelectedTimeAndPostfix.second << " )" << std::endl;

    return static_cast <UInt> ( SelectedTimeAndPostfix.second );
}



template<typename Mesh>
double Hdf5exporter<Mesh>::importFromIter( const UInt& iter )
{
    // Container for the time and the postfix
    std::pair< Real, int > SelectedTimeAndPostfix;
    if ( !this->M_procId )
    {
        // Open the xmf file
        std::ifstream xmfFile;
        xmfFile.open( ( this->M_post_dir + this->M_prefix + ".xmf" ).c_str(), std::ios::in );

        // Vector of TimeStep
        std::vector< std::pair< Real, int > > TimeAndPostfix;

        if ( xmfFile.is_open() )
        {
            // Define some variables
            std::string line;
            std::vector<std::string> stringsVector;

            // Read one-by-one all the lines of the file
            while ( !xmfFile.eof() )
            {
                std::getline( xmfFile, line, '\n' );

                // If the line begin with "<!-- Time " it is the beginning of a new block
                if ( !line.compare( 0, 10, "<!-- Time " ) )
                {
                    boost::split( stringsVector, line, boost::is_any_of( " " ) );
                    TimeAndPostfix.push_back( make_pair( string2number( stringsVector[2] ), string2number( stringsVector[4] ) ) );
                }
            }
        }

        xmfFile.close();

        // Find the closest time step
        SelectedTimeAndPostfix = TimeAndPostfix.front();
        bool found             = false;

        for ( std::vector< std::pair< Real, int > >::const_iterator i = TimeAndPostfix.begin(); i < TimeAndPostfix.end()  ; ++i )
        {
            if ( i->second == iter )
            {
                SelectedTimeAndPostfix = *i;
                found = true;
                break;
            }
        }

        ASSERT(found, "Selected iteration not found");

    }



    //std::cout << "  x-  HDF5 import from time " << SelectedTimeAndPostfix.first << " iteration " << SelectedTimeAndPostfix.second << std::endl;
    this->M_listData.begin()->storedArray()->Comm().Broadcast( &SelectedTimeAndPostfix.second, 1, 0 );
    this->M_count = SelectedTimeAndPostfix.second;

    std::ostringstream index;
    index.fill('0');

    index << std::setw(5) << this->M_count;
    this->M_postfix = "." + index.str();

    //    this->computePostfix();

    // Importing
    if ( !this->M_procId )
        std::cout << "  x-  HDF5 importing iteration "
                  << index.str()
                  << " at time " << SelectedTimeAndPostfix.first
                  << " ... " << std::flush;

    Chrono chrono;
    chrono.start();
    for ( std::list< ExporterData >::iterator i = this->M_listData.begin(); i != this->M_listData.end(); ++i )
    {
        this->rd_var(*i);
    }
    chrono.stop();

    if ( !this->M_procId )
        std::cout << "done in " << chrono.diff() << " s. (Time " << SelectedTimeAndPostfix.first
                  << ", Iteration " << SelectedTimeAndPostfix.second << " )" << std::endl;

    return SelectedTimeAndPostfix.first;
}


template<typename Mesh>
void Hdf5exporter<Mesh>::import(const Real& Tstart, const Real& dt)
{
    // dt is used to rebuild the history up to now
    Real time(Tstart - this->M_count*dt);

    for ( UInt count(0); count < this->M_count; ++count )
    {
        this->M_timeSteps.push_back(time);
        time += dt;
    }

    time += dt;

    import(time);
}

template<typename Mesh>
void Hdf5exporter<Mesh>::import(const Real& time)
{
    if ( M_HDF5.get() == 0)
    {
        M_HDF5.reset(new hdf5_type(this->M_listData.begin()->storedArray()->Comm()));
        M_HDF5->Open(this->M_post_dir+this->M_prefix+".h5"); //!! Simone
    }

    this->M_timeSteps.push_back(time);

    this->computePostfix();

    assert( this->M_postfix != "*****" );

    if (!this->M_procId) std::cout << "  x-  HDF5 importing ..."<< std::endl;

    Chrono chrono;
    chrono.start();
    for (std::list< ExporterData >::iterator i=this->M_listData.begin(); i != this->M_listData.end(); ++i)
    {
        this->rd_var(*i); ///!!! Simone
    }
    chrono.stop();
    if (!this->M_procId) std::cout << "      done in " << chrono.diff() << " s." << std::endl;
}

// ===================================================
// Get Methods
// ===================================================
template<typename Mesh>
EpetraMapType Hdf5exporter<Mesh>::mapType() const
{
    return Unique;
}

// ===================================================
// Private Methods
// ===================================================
template<typename Mesh>
void Hdf5exporter<Mesh>::defineShape()
{
}

template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_var(const ExporterData& dvar)
{

    switch( dvar.type() )
    {
    case ExporterData::Scalar:
        M_wr_scalar(dvar);
        break;
    case ExporterData::Vector:
        M_wr_vector(dvar);
        break;
    }
}

template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_scalar(const ExporterData& dvar)
{
    /* Examples:
    M_HDF5->Write("map-" + toString(Comm.NumProc()), Map);
    M_HDF5->Write("matrix", Matrix);
    M_HDF5->Write("LHS", LHS);
    M_HDF5->Write("RHS", RHS);
    */

    UInt size  = dvar.size();
    UInt start = dvar.start();

    EpetraMap subMap(dvar.storedArray()->BlockMap(), start, size);
    vector_type subVar(subMap);
    subVar.subset(*dvar.storedArray(),start);

    std::string varname (dvar.variableName()+ this->M_postfix); // see also in M_wr_attributes
    bool writeTranspose (true);
    M_HDF5->Write(varname, subVar.getEpetraVector(), writeTranspose );
}

template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_vector(const ExporterData& dvar)
{

    UInt size  = dvar.size();
    UInt start = dvar.start();

    using namespace boost;

    // solution array has to be reordered and stored in a Multivector.
    // Using auxiliary arrays:
    //shared_array<double*>                   ArrayOfPointers(new double*[nDimensions]);
    double **                                 ArrayOfPointers(new double*[nDimensions]);
    shared_array< shared_ptr<vector_type> > ArrayOfVectors (new shared_ptr<vector_type>[nDimensions]);

    int MyLDA;


    // Building subsets (new vectors) of the original data and than taking a view of them to
    // build a multivector.
    // Note: the contents of ArrayOfPointers[0,1,2] must not be deleted explicitly, since their
    // content belongs to ArrayOfVectors[0,1,2].
    // ArrayOfVectors[0,1,2] are deleted when ArrayOfVectors is destroyed

    for (UInt d ( 0 ); d < nDimensions; ++d)
    {
        EpetraMap subMap(dvar.storedArray()->BlockMap(), start+d*size, size);
        ArrayOfVectors[d].reset(new  vector_type(subMap));
        ArrayOfVectors[d]->subset(*dvar.storedArray(),start+d*size);

        ArrayOfVectors[d]->getEpetraVector().ExtractView(&ArrayOfPointers[d], &MyLDA);
    }

    EpetraMap subMap(dvar.storedArray()->BlockMap(), start, size);
    Epetra_MultiVector multiVector(View, *subMap.getMap(Unique), ArrayOfPointers, nDimensions);


    bool writeTranspose (true);
    std::string varname (dvar.variableName() + this->M_postfix); // see also in M_wr_attributes
    M_HDF5->Write(varname, multiVector, writeTranspose);

    delete[] ArrayOfPointers;
}

template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_geo()
{


    /*
     2 variables:
        &DataFile;:/Connections/Values
        &DataFile;:/PointsX/Values
        &DataFile;:/PointsY/Values
        &DataFile;:/PointsZ/Values
     note:
      if (this->M_multimesh)
        &DataFile;:/Points+ this->M_postfix + X/Values
      see also M_wr_geometry
    */

    /* Steps:
       generate local  connections and local coordinates (not repeated)
       write out
    */

  // We need a map ,but it's not always possible to use that from the variables
  // (if we write out a P0 variable)
  // We build a map for the connections based on the element numbers and for the points we fake a P1 map

    ASSERT (this->M_listData.size() > 0 , "hdf5exporter: ListData is empty");

    // Connections
    // Need to use elements not dofs for this map. Recover local element lists

    std::vector<int> elementList;
    elementList.reserve(this->M_mesh->numElements()*Mesh::ElementShape::numPoints);
    for (ID i=1; i <= this->M_mesh->numElements(); ++i)
    {
        typename Mesh::ElementType const& element (this->M_mesh->element(i));
        UInt lid=(i-1)*Mesh::ElementShape::numPoints;
        for (ID j=1; j<= Mesh::ElementShape::numPoints; ++j, ++lid)
        {
            elementList[lid] = (element.id()-1)*Mesh::ElementShape::numPoints+(j-1);
        }
    }

    Epetra_Map connectionsMap(this->M_mesh->numGlobalElements()*Mesh::ElementShape::numPoints,
			      this->M_mesh->numElements()*Mesh::ElementShape::numPoints,
			      &elementList[0],
			      0, this->M_listData.begin()->storedArray()->Comm());

    Epetra_IntVector connections(connectionsMap);
    for (ID i=1; i <= this->M_mesh->numElements(); ++i)
    {
        typename Mesh::ElementType const& element (this->M_mesh->element(i));
	UInt lid=(i-1)*Mesh::ElementShape::numPoints;
        for (ID j=1; j<= Mesh::ElementShape::numPoints; ++j, ++lid)
        {
	  connections[lid] = element.point(j).id();
        }
    }

    this->M_listData.begin()->storedArray()->Comm().Barrier();

    // this offset is needed by hdf5 since it starts numbering from 0
    //int const hdf5Offset(this->M_listData.begin()->storedArray()->BlockMap().IndexBase());
    int const hdf5Offset(0);

    // Points

    // Build a map for linear elements, even though the origianl FE might be P0
    // This gives the right map for the coordinate arrays

    EpetraMap subMap;
    switch ( Mesh::ElementShape::Shape )
    {
    case TETRA:
      {
        const RefFE & refFEP1 = feTetraP1;
        EpetraMap tmpMapP1(refFEP1, *this->M_mesh,
		       const_cast<Epetra_Comm&>(this->M_listData.begin()->storedArray()->Comm()));
        subMap = tmpMapP1;
        break;
      }
    case HEXA:
      {
        const RefFE & refFEQ1 = feHexaQ1;
        EpetraMap tmpMapQ1(refFEQ1, *this->M_mesh,
		       const_cast<Epetra_Comm&>(this->M_listData.begin()->storedArray()->Comm()));
        subMap = tmpMapQ1;
        break;
      }
    case LINE:
      {
        const RefFE & refFEP11D = feSegP1;
        EpetraMap tmpMapQ11D(refFEP11D, *this->M_mesh,
		       const_cast<Epetra_Comm&>(this->M_listData.begin()->storedArray()->Comm()));
        subMap = tmpMapQ11D;
        break;
      }
    default:
        ERROR_MSG( "FE not allowed in HDF5 Exporter" );

    }
    //    EpetraMap subMap(refFE, *this->M_mesh,
    //		     const_cast<Epetra_Comm&>(this->M_listData.begin()->storedArray()->Comm()));

    EpetraVector pointsX(subMap);
    EpetraVector pointsY(subMap);
    EpetraVector pointsZ(subMap);

    int gid;
    for(ID i=1; i <= this->M_mesh->numVertices(); ++i)
    {
        typename Mesh::PointType const& point (this->M_mesh->pointList(i));
        gid = point.id() - hdf5Offset;

        bool insertedX(true);
        bool insertedY(true);
        bool insertedZ(true);

        insertedX = insertedX && pointsX.checkAndSet(gid, point.x());
        insertedY = insertedY && pointsY.checkAndSet(gid, point.y());
        insertedZ = insertedZ && pointsZ.checkAndSet(gid, point.z());
    }

    // Now we are ready to export the vectors to the hdf5 file

    std::string pointsXVarname("PointsX");
    std::string pointsYVarname("PointsY");
    std::string pointsZVarname("PointsZ");
    std::string connectionsVarname("Connections");

    if (this->M_multimesh)
    {
        connectionsVarname  += this->M_postfix; // see also in M_wr_topology
        pointsXVarname      += this->M_postfix; // see also in M_wr_geometry
        pointsYVarname      += this->M_postfix; // see also in M_wr_geometry
        pointsZVarname      += this->M_postfix; // see also in M_wr_geometry
    }

    M_HDF5->Write(connectionsVarname, connections);
    // bool writeTranspose (true);
    M_HDF5->Write(pointsXVarname, pointsX.getEpetraVector(), true);
    M_HDF5->Write(pointsYVarname, pointsY.getEpetraVector(), true);
    M_HDF5->Write(pointsZVarname, pointsZ.getEpetraVector(), true);
}


// write empty xdmf file
template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_initXdmf()
{
    if (this->M_procId == 0)
    {
        M_xdmf.open( (this->M_post_dir+this->M_prefix+".xmf").c_str(), std::ios_base::out );

        M_xdmf <<
            "<?xml version=\"1.0\" ?>\n" <<
            "<!DOCTYPE Xdmf SYSTEM \"" << this->M_prefix << ".xdmf\" [\n" <<
            "<!ENTITY DataFile \"" << this->M_prefix << ".h5\">\n" <<
            "]>\n" <<
            "<!-- " << this->M_prefix << ".h5 is generated by LifeV -->\n" <<
            "<Xdmf>\n" <<
            "  <Domain Name=\"" << this->M_prefix << "\">\n" <<
            "    <Grid Name=\"" << this->M_prefix << "Grid\" GridType=\"Collection\" CollectionType=\"Temporal\">\n" <<
            "\n";

        M_wr_closeLinesXdmf();
    }
}

// save position and write closing lines
template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_closeLinesXdmf()
{
    // save position
    M_closingLinesPosition = M_xdmf.tellp();

    // write closing lines
    M_xdmf << M_closingLines;
    M_xdmf.flush();

}

// remove closing lines
template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_removeCloseLinesXdmf()
{
    M_xdmf.seekp(M_closingLinesPosition);
}


template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_Xdmf(const Real& time)
{
    /*
      strategy: write the topology,
      <Topology
         Type="Tetrahedron"
         NumberOfElements="183"
         BaseOffset="1">
         <DataStructure Format="HDF"
                        Dimensions="602  4" ???
                        DataType="Int"
                        Precision="8">
             &DataFile;:/Connections/Values
         </DataStructure>
      </Topology>
      and then the geometry
      <Geometry Type="X_Y_Z">
         <DataStructure Format="HDF"
                        Dimensions="183"
                        DataType="Float"
                        Precision="8">
             &DataFile;:/PointsX/Values
         </DataStructure>
         <DataStructure Format="HDF"
                        Dimensions="183"
                        DataType="Float"
                        Precision="8">
             &DataFile;:/PointsY/Values
         </DataStructure>
         <DataStructure Format="HDF"
                        Dimensions="183"
                        DataType="Float"
                        Precision="8">
             &DataFile;:/PointsY/Values
         </DataStructure>
      </Geometry>


      In this aim we create
      two Epetra_IntVector, one with a repeated map and the second with a unique map
      two Epetra_MultiVector, one with a repeated map and the second with a unique map

      The Int vectors are for the connections, the double ones for the Vertices
    */

    if (this->M_procId == 0)
    {
        M_wr_removeCloseLinesXdmf();

        // write grid with time, topology, geometry and attributes
        // NOTE: The first line (<!-- Time t Iteration i -->) is used in function importFromTime.
        //       Check compatibility after any change on it!
        M_xdmf <<
            "<!-- Time " << time << " Iteration " << this->M_postfix.substr(1,5) << " -->\n" <<
            "    <Grid Name=\"Mesh " << time << "\">\n" <<
            "      <Time TimeType=\"Single\" Value=\"" << time << "\" />\n";
        M_wr_topology(M_xdmf);
        M_wr_geometry(M_xdmf);
        M_wr_attributes(M_xdmf);

        M_xdmf << "\n"
            "    </Grid>\n\n";



        // write closing lines
        M_wr_closeLinesXdmf();
    }
}

template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_topology  ( std::ofstream& xdmf )
{
    std::string FEstring;

    switch ( Mesh::ElementShape::Shape )
    {
    case TETRA:
        FEstring = "Tetrahedron";
        break;
    case HEXA:
        FEstring = "Hexahedron";
        break;
    case LINE:
        FEstring = "Polyline";
        break;
    default:
        ERROR_MSG( "FE not allowed in HDF5 Exporter" );
    }

    xdmf <<
        "      <Topology\n" <<
        "         Type=\"" << FEstring <<"\"\n" <<
        "         NumberOfElements=\"" << this->M_mesh->numGlobalElements() << "\"\n" <<
        "         BaseOffset=\"1\">\n" <<
        "         <DataStructure Format=\"HDF\"\n" <<
        "                        Dimensions=\""<< this->M_mesh->numGlobalElements() << " " << this->M_mesh->numLocalVertices() << "\"\n" <<
        "                        DataType=\"Int\"\n" <<
        "                        Precision=\"8\">\n" <<
        "             " << M_outputFileName << ":/Connections/Values\n" <<
        "         </DataStructure>\n" <<
        "      </Topology>\n";
}

template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_geometry  ( std::ofstream& xdmf )
{

  std::string postfix_string;

  // see also in postProcess
  if (this->M_multimesh)
    postfix_string = this->M_postfix;
  else
    postfix_string = "";


    xdmf <<
        "      <Geometry Type=\"X_Y_Z\">\n" <<
        "         <DataStructure Format=\"HDF\"\n" <<
        "                        Dimensions=\"" << this->M_mesh->numGlobalVertices() << "\"\n" <<
        "                        DataType=\"Float\"\n" <<
        "                        Precision=\"8\">\n" <<
      "             " << M_outputFileName << ":/" << "PointsX" << postfix_string << "/Values\n" <<
        "         </DataStructure>\n" <<
        "         <DataStructure Format=\"HDF\"\n" <<
        "                        Dimensions=\"" << this->M_mesh->numGlobalVertices() << "\"\n" <<
        "                        DataType=\"Float\"\n" <<
        "                        Precision=\"8\">\n" <<
        "             " << M_outputFileName << ":/" << "PointsY" << postfix_string << "/Values\n" <<
        "         </DataStructure>\n" <<
        "         <DataStructure Format=\"HDF\"\n" <<
        "                        Dimensions=\"" << this->M_mesh->numGlobalVertices() << "\"\n" <<
        "                        DataType=\"Float\"\n" <<
        "                        Precision=\"8\">\n" <<
        "             " << M_outputFileName << ":/" << "PointsZ" << postfix_string << "/Values\n" <<
        "         </DataStructure>\n" <<
        "      </Geometry>\n" <<
        "\n";
}

template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_attributes  ( std::ofstream& xdmf )
{

    // Loop on the variables to output
    for (std::list< ExporterData >::const_iterator i=this->M_listData.begin(); i != this->M_listData.end(); ++i)
    {
        xdmf <<
            "\n      <Attribute\n" <<
            "         Type=\"" << i->typeName() << "\"\n" <<
            "         Center=\"" << i->whereName() << "\"\n" <<
            "         Name=\"" << i->variableName()<<"\">\n";

        switch( i->type() )
        {
        case ExporterData::Scalar:
            M_wr_scalar_datastructure(xdmf, *i);
            break;
        case ExporterData::Vector:
            M_wr_scalar_datastructure(xdmf, *i);
            //M_wr_vector_datastructure(xdmf, *i);
            break;
        }

        xdmf <<
            "      </Attribute>\n";
    }
}

template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_scalar_datastructure  ( std::ofstream& xdmf, const ExporterData& dvar )
{

    Int globalUnknowns (0);
    switch ( dvar.where() )
    {
    case ExporterData::Node:
        globalUnknowns = this->M_mesh->numGlobalVertices();
        break;
    case ExporterData::Cell:
        globalUnknowns = this->M_mesh->numGlobalElements();
        break;
    }

    // First: hyperslab definition, then description of the data
    xdmf <<

        "         <DataStructure ItemType=\"HyperSlab\"\n" <<
        "                        Dimensions=\"" << globalUnknowns << " " << dvar.typeDim() << "\"\n" <<
        "                        Type=\"HyperSlab\">\n" <<
        "           <DataStructure  Dimensions=\"3 2\"\n" <<
        "                           Format=\"XML\">\n" <<
        "               0    0\n" <<
        "               1    1\n" <<
        "               " << globalUnknowns << " " << dvar.typeDim() << "\n" <<
        "           </DataStructure>\n" <<

        "           <DataStructure  Format=\"HDF\"\n" <<
        "                           Dimensions=\"" << dvar.size() << " " << dvar.typeDim() << "\"\n" <<
        "                           DataType=\"Float\"\n" <<
        "                           Precision=\"8\">\n" <<
        "               " << M_outputFileName << ":/" << dvar.variableName() << this->M_postfix  <<"/Values\n" << // see also in M_wr_vector/scalar
        "           </DataStructure>\n" <<
        "         </DataStructure>\n";

}

template <typename Mesh>
void Hdf5exporter<Mesh>::M_wr_vector_datastructure  ( std::ofstream& xdmf, const ExporterData& dvar )
{


    string coord[3]={"X","Y","Z"}; // see also wr_vector

    xdmf <<
                "         <DataStructure ItemType=\"Function\"\n" <<
                "                        Dimensions=\"" << this->M_mesh->numGlobalVertices() << " " << dvar.typeDim() << "\"\n" <<
                "                        Function=\"JOIN($0 , $1, $2)\">\n";

            for(int i(0); i < dvar.typeDim(); ++i)
                {
                    xdmf <<
                "           <DataStructure  Format=\"HDF\"\n" <<
                "                           Dimensions=\"" << this->M_mesh->numGlobalVertices() << " 1\"\n" <<
                "                           DataType=\"Float\"\n" <<
                "                           Precision=\"8\">\n" <<
                "               " << M_outputFileName << ":/" << dvar.variableName()<< coord[i] << this->M_postfix  <<"/Values\n" << // see also in M_wr_vector/scalar
                "           </DataStructure>\n";
                }

    xdmf <<
                "         </DataStructure>\n";


}

template <typename Mesh>
void Hdf5exporter<Mesh>::rd_var(ExporterData& dvar)
{
    if ( M_HDF5.get() == 0)
    {
        M_HDF5.reset(new hdf5_type(dvar.storedArray()->BlockMap().Comm()));
        M_HDF5->Open(this->M_post_dir+this->M_prefix+".h5"); //!! Simone
    }
    super::rd_var(dvar);
}

template <typename Mesh>
void Hdf5exporter<Mesh>::M_rd_scalar(ExporterData& dvar)
{

    UInt size  = dvar.size();
    UInt start = dvar.start();

    EpetraMap subMap(dvar.storedArray()->BlockMap(), start, size);
    Epetra_MultiVector* subVar(0);

    std::string varname (dvar.variableName()); // see also in M_wr_attributes
    if(this->M_postfix!="")
    {
        varname += this->M_postfix;
    }
    bool readTranspose (true);
    M_HDF5->Read(varname, *subMap.getMap(this->mapType()), subVar, readTranspose);

    dvar.storedArray()->subset(*subVar, subMap, 0, start );

    delete subVar;

}

template <typename Mesh>
void Hdf5exporter<Mesh>::M_rd_vector( ExporterData& dvar)
{
    UInt size  = dvar.size();
    UInt start = dvar.start();

    using namespace boost;

    // solution array has first to be read has Multivector.

    // first read the multivector:
    EpetraMap subMap(dvar.storedArray()->BlockMap(), start, size);
    Epetra_MultiVector* subVar(0);

    bool readTranspose (true);
    std::string varname (dvar.variableName()); // see also in M_wr_attributes

    if(this->M_postfix!="")
    {
        varname += this->M_postfix;
    }

    M_HDF5->Read(varname, *subMap.getMap(this->mapType()), subVar, readTranspose);


    // then put back value in our EpetraVector

    for (UInt d ( 0 ); d < nDimensions; ++d)
    {
        dvar.storedArray()->subset(*subVar, subMap,  0, start+d*size, d );
    }

    delete subVar;
}

template <typename Mesh>
void Hdf5exporter<Mesh>::writeGraph()
{
    std::vector<int> partitionSizes;
    int size, maxSize = 0;

    // Calculate the maximum size of the partitions and store partition
    // sizes
    partitionSizes.reserve(M_graph->size());
    for (std::vector<std::vector<int> >::iterator it = M_graph->begin();
         it != M_graph->end(); ++it)
    {
        size = it->size();
        if (size > maxSize)
        {
            maxSize = size;
        }
        partitionSizes.push_back(size);
    }

    M_HDF5->Write("Graph", "partition_sizes", H5T_NATIVE_INT,
                  M_graph->size(), &partitionSizes[0]);

    // Write partition array size
    M_HDF5->Write("Graph", "number_partitions", int(M_graph->size()));

    // Write partition array
    std::stringstream index;
    for (UInt i = 0; i < M_graph->size(); ++i)
    {
        index << i;
        M_HDF5->Write("Graph", "partition_graph_" + index.str(),
                      H5T_NATIVE_INT, partitionSizes[i],
                      &(*M_graph)[i][0]);
        index.str(std::string());
        index.clear();
    }
}

template <typename Mesh>
void Hdf5exporter<Mesh>::writePartition(mesh_ptrtype mesh, std::string& suffix)
{
    UInt elementNodes, faceNodes;
    switch (Mesh::ElementShape::Shape)
    {
    case HEXA:
        elementNodes = 8;
        faceNodes    = 4;
        break;
    case TETRA:
        elementNodes = 4;
        faceNodes    = 3;
    }

    M_HDF5->Write("Mesh", "Counters.NumPoints" + suffix, (int) mesh->numPoints());
    M_HDF5->Write("Mesh", "Counters.NumBPoints" + suffix, (int) mesh->numBPoints());

    M_HDF5->Write("Mesh", "Counters.NumVertices" + suffix, (int) mesh->numVertices());
    M_HDF5->Write("Mesh", "Counters.NumBVertices" + suffix, (int) mesh->numBVertices());
    M_HDF5->Write("Mesh", "Counters.NumGlobalVertices" + suffix, (int) mesh->numGlobalVertices());

    M_HDF5->Write("Mesh", "Counters.NumEdges" + suffix, (int) mesh->numEdges());
    M_HDF5->Write("Mesh", "Counters.NumBEdges" + suffix, (int) mesh->numBEdges());
    M_HDF5->Write("Mesh", "Counters.NumGlobalEdges" + suffix, (int) mesh->numGlobalEdges());

    M_HDF5->Write("Mesh", "Counters.NumFaces" + suffix, (int) mesh->numFaces());
    M_HDF5->Write("Mesh", "Counters.NumBFaces" + suffix, (int) mesh->numBFaces());
    M_HDF5->Write("Mesh", "Counters.NumGlobalFaces" + suffix, (int) mesh->numGlobalFaces());

    M_HDF5->Write("Mesh", "Counters.NumVolumes" + suffix, (int) mesh->numVolumes());
    M_HDF5->Write("Mesh", "Counters.NumGlobalVolumes" + suffix, (int) mesh->numGlobalVolumes());

    int numPoints = mesh->numPoints();
    std::vector<double> tmpVectorDouble(numPoints);
    std::vector<std::vector<double> > pointCoordinates(3, tmpVectorDouble);

    std::vector<int> pointMarkers(numPoints);
    std::vector<int> pointGlobalId(numPoints);
    std::vector<int> pointBoundaryFlags(numPoints);

    std::map<int, int>::iterator it = mesh->localToGlobalNode().begin();
    for (UInt j = 0; it != mesh->localToGlobalNode().end(); ++it, ++j)
    {
        pointCoordinates[0][j] = mesh->pointList[j].x();
        pointCoordinates[1][j] = mesh->pointList[j].y();
        pointCoordinates[2][j] = mesh->pointList[j].z();
        pointMarkers[j] = mesh->pointList[j].marker();
        pointGlobalId[j] = it->second;
        if (mesh->isBoundaryPoint(j+1))
        {
            pointBoundaryFlags[j] = 1;
        }
    }

    M_HDF5->Write("Mesh", "Points.x" + suffix, H5T_NATIVE_DOUBLE, numPoints, &pointCoordinates[0][0]);
    M_HDF5->Write("Mesh", "Points.y" + suffix, H5T_NATIVE_DOUBLE, numPoints, &pointCoordinates[1][0]);
    M_HDF5->Write("Mesh", "Points.z" + suffix, H5T_NATIVE_DOUBLE, numPoints, &pointCoordinates[2][0]);
    M_HDF5->Write("Mesh", "Points.f" + suffix, H5T_NATIVE_INT, numPoints, &pointMarkers[0]);
    M_HDF5->Write("Mesh", "Points.GlobalId" + suffix, H5T_NATIVE_INT, numPoints, &pointGlobalId[0]);
    M_HDF5->Write("Mesh", "Points.BoundaryFlag" + suffix, H5T_NATIVE_INT, numPoints, &pointBoundaryFlags[0]);

    pointCoordinates.clear();
    pointMarkers.clear();
    pointBoundaryFlags.clear();
    pointGlobalId.clear();

    int numEdges = mesh->numEdges();
    std::vector<int> tmpVectorInt(numEdges);
    std::vector<std::vector<int> > edgePoints(2,tmpVectorInt);

    std::vector<int> edgeMarkers(numEdges);
    std::vector<int> edgeGlobalId(numEdges);
    std::vector<int> edgeBoundaryFlags(numEdges);

    for (UInt j = 0; j < numEdges; ++j)
    {
        edgePoints[0][j] = mesh->edgeList[j].point(1).localId();
        edgePoints[1][j] = mesh->edgeList[j].point(2).localId();
        edgeMarkers[j] = mesh->edgeList[j].marker();
        edgeGlobalId[j] = mesh->edgeList[j].id();

        if (mesh->isBoundaryEdge(j+1))
        {
            edgeBoundaryFlags[j] = 1;
        }
    }

    M_HDF5->Write("Mesh", "Edges.p1" + suffix, H5T_NATIVE_INT, numEdges, &edgePoints[0][0]);
    M_HDF5->Write("Mesh", "Edges.p2" + suffix, H5T_NATIVE_INT, numEdges, &edgePoints[1][0]);
    M_HDF5->Write("Mesh", "Edges.f" + suffix, H5T_NATIVE_INT, numEdges, &edgeMarkers[0]);
    M_HDF5->Write("Mesh", "Edges.GlobalId" + suffix, H5T_NATIVE_INT, numEdges, &edgeGlobalId[0]);
    M_HDF5->Write("Mesh", "Edges.BoundaryFlag" + suffix, H5T_NATIVE_INT, numEdges, &edgeBoundaryFlags[0]);

    edgePoints.clear();
    edgeMarkers.clear();
    edgeGlobalId.clear();
    edgeBoundaryFlags.clear();

    int numFaces = mesh->numFaces();
    tmpVectorInt.resize(numFaces);
    std::vector<std::vector<int> > facePoints(faceNodes, tmpVectorInt);

    std::vector<int> faceMarkers(numFaces);
    std::vector<int> faceGlobalId(numFaces);
    std::vector<int> faceBoundaryFlags(numFaces);

    std::vector<std::vector<int> > faceNeighbourId(2, tmpVectorInt);
    std::vector<std::vector<int> > faceNeighbourPos(2, tmpVectorInt);

    for (UInt j = 0; j < numFaces; ++j)
    {
        for (UInt k = 0; k < faceNodes; ++k)
        {
            facePoints[k][j] = mesh->faceList[j].point(k + 1).localId();
        }
        faceMarkers[j] = mesh->faceList[j].marker();
        faceGlobalId[j] = mesh->faceList[j].id();

        faceNeighbourId[0][j] = mesh->faceList[j].ad_first();
        faceNeighbourId[1][j] = mesh->faceList[j].ad_second();
        faceNeighbourPos[0][j] = mesh->faceList[j].pos_first();
        faceNeighbourPos[1][j] = mesh->faceList[j].pos_second();

        if (mesh->isBoundaryFace(j+1))
        {
            faceBoundaryFlags[j] = 1;
        }
    }

    std::stringstream idx;
    for (UInt k = 0; k < faceNodes; ++k)
    {
        idx << k + 1;
        M_HDF5->Write("Mesh", "Faces.p" + idx.str() + suffix, H5T_NATIVE_INT, numFaces, &facePoints[k][0]);
        idx.str(std::string());
        idx.clear();
    }
    M_HDF5->Write("Mesh", "Faces.f" + suffix, H5T_NATIVE_INT, numFaces, &faceMarkers[0]);
    M_HDF5->Write("Mesh", "Faces.GlobalId" + suffix, H5T_NATIVE_INT, numFaces, &faceGlobalId[0]);
    M_HDF5->Write("Mesh", "Faces.BoundaryFlag" + suffix, H5T_NATIVE_INT, numFaces, &faceBoundaryFlags[0]);

    M_HDF5->Write("Mesh", "Faces.NeighbourId1" + suffix, H5T_NATIVE_INT, numFaces, &faceNeighbourId[0][0]);
    M_HDF5->Write("Mesh", "Faces.NeighbourId2" + suffix, H5T_NATIVE_INT, numFaces, &faceNeighbourId[1][0]);
    M_HDF5->Write("Mesh", "Faces.NeighbourPos1" + suffix, H5T_NATIVE_INT, numFaces, &faceNeighbourPos[0][0]);
    M_HDF5->Write("Mesh", "Faces.NeighbourPos2" + suffix, H5T_NATIVE_INT, numFaces, &faceNeighbourPos[1][0]);

    facePoints.clear();
    faceMarkers.clear();
    faceGlobalId.clear();
    faceBoundaryFlags.clear();
    faceNeighbourId.clear();
    faceNeighbourPos.clear();

    int numVolumes = mesh->numVolumes();
    tmpVectorInt.resize(numVolumes);
    std::vector<std::vector<int> > volumePoints(elementNodes, tmpVectorInt);

    std::vector<int> volumeMarkers(numVolumes);
    std::vector<int> volumeGlobalId(numVolumes);

    for (UInt j = 0; j < numVolumes; ++j)
    {
        for (UInt k = 0; k < elementNodes; ++k)
        {
            volumePoints[k][j] = mesh->volumeList[j].point(k + 1).localId();
        }
        volumeMarkers[j] = mesh->volumeList[j].marker();
        volumeGlobalId[j] = mesh->volumeList[j].id();
    }

    for (UInt k = 0; k < elementNodes; ++k)
    {
        idx << k + 1;
        M_HDF5->Write("Mesh", "Volumes.p" + idx.str() + suffix, H5T_NATIVE_INT, numVolumes, &volumePoints[k][0]);
        idx.str(std::string());
        idx.clear();
    }
    M_HDF5->Write("Mesh", "Volumes.f" + suffix, H5T_NATIVE_INT, numVolumes, &volumeMarkers[0]);
    M_HDF5->Write("Mesh", "Volumes.GlobalId" + suffix, H5T_NATIVE_INT, numVolumes, &volumeGlobalId[0]);

    volumePoints.clear();
    volumeMarkers.clear();
    volumeGlobalId.clear();
}

template <typename Mesh>
void Hdf5exporter<Mesh>::writeSerialMesh()
{
    std::stringstream index;
    std::string suffix;

    for (UInt i = 0; i < M_serialMesh->size(); ++i)
    {
        index << i;
        suffix = "." + index.str();

        writePartition((*M_serialMesh)[i], suffix);

        index.str(std::string());
        index.clear();
    }
}

template <typename Mesh>
void Hdf5exporter<Mesh>::writeParallelMesh()
{
    std::stringstream index;
    std::string suffix;

    index << M_comm->MyPID();
    suffix = "." + index.str();

    writePartition(M_parallelMesh, suffix);
}

}
#endif

#endif
