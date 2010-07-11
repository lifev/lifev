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
  @brief Class derived from Hdf5exporter to provide I/O for the mesh partitions (RegionMesh3D only)

  @author Radu Popescu <radu.popescu@epfl.ch>
  @date 9-07-2010
 */

#ifndef HDF5FILTER3DMESH_H
#define HDF5FILTER3DMESH_H 1

#include <life/lifefilters/hdf5exporter.hpp>

namespace LifeV {

//! Class derived from Hdf5exporter to provide I/O for the mesh partitions (RegionMesh3D only)
/*!
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
template <typename Mesh>
class HDF5Filter3DMesh : public Hdf5exporter<Mesh>
{

public:

    typedef Hdf5exporter<Mesh> base;
    typedef typename base::mesh_ptrtype mesh_ptrtype;
    typedef typename base::vector_rawtype vector_type;
    typedef typename base::vector_ptrtype vector_ptrtype;

    typedef EpetraExt::HDF5 hdf5_type;
    typedef boost::shared_ptr<hdf5_type> hdf5_ptrtype;
    typedef std::vector<std::vector<int> > graph_type;
    typedef boost::shared_ptr<graph_type> graph_ptrtype;
    typedef boost::shared_ptr<std::vector<mesh_ptrtype> > serial_mesh_ptrtype;

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor for HDF5Filter3DMesh
    HDF5Filter3DMesh() {}

    //! Constructor for HDF5Filter3DMesh
    /*!
       @param dfile the GetPot data file where you must provide an [exporter] section with:
          "start"     (start index for sections in the hdf5 data structure 0 for 000, 1 for 001 etc.),
          "save"      (how many time steps per postprocessing)
          "multimesh" ( = true if the mesh has to be saved at each post-processing step)
       @param mesh the mesh
       @param the prefix for the case file (ex. "test" for test.case)
       @param the procId determines de CPU id. if negative, it ussemes there is only one processor
    */
    HDF5Filter3DMesh(const GetPot& dfile, mesh_ptrtype mesh, const std::string& prefix, const int& procId);

    //! Constructor for HDF5Filter3DMesh without prefix and procID
    /*!
       @param dfile the GetPot data file where you must provide an [exporter] section with:
          "start"     (start index for sections in the hdf5 data structure 0 for 000, 1 for 001 etc.),
          "save"      (how many time steps per postprocessing)
          "multimesh" ( = true if the mesh has to be saved at each post-processing step)
       @param mesh the mesh
    */
    HDF5Filter3DMesh(const GetPot& dfile, const std::string& prefix);

    //! Destructor for HDF5Filter3DMesh
    ~HDF5Filter3DMesh() {}

    //@}

    //! @name Methods
    //@{
    virtual void postProcess(const Real& time);

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
      After the mesh partition is loaded from the HDF5 file, the simulation is run and HDF5Filter3DMesh::postProcess()
      is called, the original contents of the HDF5 will be lost. To keep mesh partition, call
      HDF5Filter3DMesh::addMyMeshPartition() before calling HDF5Filter3DMesh::postProcess();
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

private:

    // The following private methods are for writing the partitioned graph
    // and mesh to the output file
    void writeGraph();
    void writePartition(mesh_ptrtype partition, std::string& suffix);
    void writeParallelMesh();
    void writeSerialMesh();

    //@}

    serial_mesh_ptrtype   M_serialMesh;
    mesh_ptrtype          M_parallelMesh;
    graph_ptrtype         M_graph;
    Epetra_Comm*          M_comm;

};

// ===================================================
// Constructors
// ===================================================

template<typename Mesh>
HDF5Filter3DMesh<Mesh>::HDF5Filter3DMesh(const GetPot& dfile, mesh_ptrtype mesh, const std::string& prefix,
                                             const int& procId) :
    base                ( dfile, mesh, prefix, procId )
{
}

template<typename Mesh>
HDF5Filter3DMesh<Mesh>::HDF5Filter3DMesh(const GetPot& dfile, const std::string& prefix):
    base                ( dfile, prefix )
{
}

// ===================================================
// Methods
// ===================================================

template<typename Mesh>
void HDF5Filter3DMesh<Mesh>::loadGraph(graph_ptrtype graph, Epetra_Comm *comm)
{
    if (this->M_HDF5.get() == 0)
    {
        this->M_HDF5.reset(new hdf5_type(*comm));
    }
    if (! this->M_HDF5->IsOpen())
    {
        this->M_HDF5->Open(this->M_post_dir + this->M_prefix + ".h5", H5F_ACC_RDONLY);
    }

    int nPartitions;

    this->M_HDF5->Read("Graph", "number_partitions", nPartitions);

    std::vector<int> partitionSizes(nPartitions);
    this->M_HDF5->Read("Graph", "partition_sizes", H5T_NATIVE_INT, nPartitions,
                 &partitionSizes[0]);

    graph->resize(0);
    graph->reserve(nPartitions);

    std::vector<int> partBuffer;
    std::stringstream index;

    for (UInt i = 0; i < nPartitions; ++i)
    {
        partBuffer.resize(partitionSizes[i]);
        index << i;
        this->M_HDF5->Read("Graph", "partition_graph_" + index.str(),
                     H5T_NATIVE_INT, partitionSizes[i],
                     &partBuffer[0]);
        graph->push_back(partBuffer);
        index.str(std::string());
        index.clear();
    }
}

template<typename Mesh>
void HDF5Filter3DMesh<Mesh>::loadMyPartition(mesh_ptrtype meshPartition, Epetra_Comm* comm)
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

    if (this->M_HDF5.get() == 0)
    {
        this->M_HDF5.reset(new hdf5_type(*comm));
    }
    if (! this->M_HDF5->IsOpen())
    {
        this->M_HDF5->Open(this->M_post_dir + this->M_prefix + ".h5", H5F_ACC_RDONLY);
    }

    std::stringstream index;
    index << comm->MyPID();
    std::string suffix = "." + index.str();

    // Read counters from file and set them in the mesh

    // Points
    int numPoints;
    this->M_HDF5->Read("Mesh", "Counters.NumPoints" + suffix, numPoints);
    meshPartition->setMaxNumPoints(numPoints, true);
    int numBPoints;
    this->M_HDF5->Read("Mesh", "Counters.NumBPoints" + suffix, numBPoints);
    meshPartition->setNumBPoints(numBPoints);

    // Vertices
    int numVertices;
    this->M_HDF5->Read("Mesh", "Counters.NumVertices" + suffix, numVertices);
    meshPartition->setNumVertices(numVertices);
    int numBVertices;
    this->M_HDF5->Read("Mesh", "Counters.NumBVertices" + suffix, numBVertices);
    meshPartition->setNumBVertices(numBVertices);
    int numGlobalVertices;
    this->M_HDF5->Read("Mesh", "Counters.NumGlobalVertices" + suffix, numGlobalVertices);
    meshPartition->setNumGlobalVertices(numGlobalVertices);

    // Edges
    int numEdges;
    this->M_HDF5->Read("Mesh", "Counters.NumEdges" + suffix, numEdges);
    meshPartition->setNumEdges(numEdges);
    meshPartition->setMaxNumEdges(numEdges);
    int numBEdges;
    this->M_HDF5->Read("Mesh", "Counters.NumBEdges" + suffix, numBEdges);
    meshPartition->setNumBEdges(numBEdges);
    int numGlobalEdges;
    this->M_HDF5->Read("Mesh", "Counters.NumGlobalEdges" + suffix, numGlobalEdges);
    meshPartition->setMaxNumGlobalEdges(numGlobalEdges);

    // Faces
    int numFaces;
    this->M_HDF5->Read("Mesh", "Counters.NumFaces" + suffix, numFaces);
    meshPartition->setNumFaces(numFaces);
    meshPartition->setMaxNumFaces(numFaces);
    int numBFaces;
    this->M_HDF5->Read("Mesh", "Counters.NumBFaces" + suffix, numBFaces);
    meshPartition->setNumBFaces(numBFaces);
    int numGlobalFaces;
    this->M_HDF5->Read("Mesh", "Counters.NumGlobalFaces" + suffix, numGlobalFaces);
    meshPartition->setMaxNumGlobalFaces(numGlobalFaces);

    // Volumes
    int numVolumes;
    this->M_HDF5->Read("Mesh", "Counters.NumVolumes" + suffix, numVolumes);
    meshPartition->setMaxNumVolumes(numVolumes, true);
    int numGlobalVolumes;
    this->M_HDF5->Read("Mesh", "Counters.NumGlobalVolumes" + suffix, numGlobalVolumes);
    meshPartition->setMaxNumGlobalVolumes(numGlobalVolumes);

    // Read the list of points
    std::vector<double> tmpVectorDouble(numPoints);
    std::vector<std::vector<double> > pointCoordinates(3, tmpVectorDouble);

    std::vector<int> pointMarkers(numPoints);
    std::vector<int> pointBoundaryFlags(numPoints);
    std::vector<int> pointGlobalId(numPoints);

    this->M_HDF5->Read("Mesh", "Points.x" + suffix, H5T_NATIVE_DOUBLE, numPoints, &pointCoordinates[0][0]);
    this->M_HDF5->Read("Mesh", "Points.y" + suffix, H5T_NATIVE_DOUBLE, numPoints, &pointCoordinates[1][0]);
    this->M_HDF5->Read("Mesh", "Points.z" + suffix, H5T_NATIVE_DOUBLE, numPoints, &pointCoordinates[2][0]);
    this->M_HDF5->Read("Mesh", "Points.f" + suffix, H5T_NATIVE_INT, numPoints, &pointMarkers[0]);
    this->M_HDF5->Read("Mesh", "Points.BoundaryFlag" + suffix, H5T_NATIVE_INT, numPoints, &pointBoundaryFlags[0]);
    this->M_HDF5->Read("Mesh", "Points.GlobalId" + suffix, H5T_NATIVE_INT, numPoints, &pointGlobalId[0]);

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

    this->M_HDF5->Read("Mesh", "Edges.p1" + suffix, H5T_NATIVE_INT, numEdges, &edgePoints[0][0]);
    this->M_HDF5->Read("Mesh", "Edges.p2" + suffix, H5T_NATIVE_INT, numEdges, &edgePoints[1][0]);
    this->M_HDF5->Read("Mesh", "Edges.f" + suffix, H5T_NATIVE_INT, numEdges, &edgeMarkers[0]);
    this->M_HDF5->Read("Mesh", "Edges.GlobalId" + suffix, H5T_NATIVE_INT, numEdges, &edgeGlobalId[0]);
    this->M_HDF5->Read("Mesh", "Edges.BoundaryFlag" + suffix, H5T_NATIVE_INT, numEdges, &edgeBoundaryFlags[0]);

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
        this->M_HDF5->Read("Mesh", "Faces.p" + idx.str() + suffix, H5T_NATIVE_INT, numFaces, &facePoints[k][0]);
        idx.str(std::string());
        idx.clear();
    }
    this->M_HDF5->Read("Mesh", "Faces.f" + suffix, H5T_NATIVE_INT, numFaces, &faceMarkers[0]);
    this->M_HDF5->Read("Mesh", "Faces.GlobalId" + suffix, H5T_NATIVE_INT, numFaces, &faceGlobalId[0]);
    this->M_HDF5->Read("Mesh", "Faces.BoundaryFlag" + suffix, H5T_NATIVE_INT, numFaces, &faceBoundaryFlags[0]);

    this->M_HDF5->Read("Mesh", "Faces.NeighbourId1" + suffix, H5T_NATIVE_INT, numFaces, &faceNeighbourId[0][0]);
    this->M_HDF5->Read("Mesh", "Faces.NeighbourId2" + suffix, H5T_NATIVE_INT, numFaces, &faceNeighbourId[1][0]);
    this->M_HDF5->Read("Mesh", "Faces.NeighbourPos1" + suffix, H5T_NATIVE_INT, numFaces, &faceNeighbourPos[0][0]);
    this->M_HDF5->Read("Mesh", "Faces.NeighbourPos2" + suffix, H5T_NATIVE_INT, numFaces, &faceNeighbourPos[1][0]);


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
        this->M_HDF5->Read("Mesh", "Volumes.p" + idx.str() + suffix, H5T_NATIVE_INT, numVolumes, &volumePoints[k][0]);
        idx.str(std::string());
        idx.clear();
    }
    this->M_HDF5->Read("Mesh", "Volumes.f" + suffix, H5T_NATIVE_INT, numVolumes, &volumeMarkers[0]);
    this->M_HDF5->Read("Mesh", "Volumes.GlobalId" + suffix, H5T_NATIVE_INT, numVolumes, &volumeGlobalId[0]);

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
void HDF5Filter3DMesh<Mesh>::postProcess(const Real& time)
{
    if ( this->M_HDF5.get() == 0)
    {
        if (this->M_listData.size() != 0)
        {
            this->M_HDF5.reset(new hdf5_type(this->M_listData.begin()->storedArray()->Comm()));
        }
        else
        {
            this->M_HDF5.reset(new hdf5_type(*M_comm));
        }
        this->M_outputFileName=this->M_prefix+".h5";
        this->M_HDF5->Create(this->M_post_dir+this->M_outputFileName);

        // write empty xdmf file
        this->M_wr_initXdmf();

        if (!this->M_multimesh)
        {
            if (this->M_listData.size() != 0)
            {
                this->M_wr_geo(); // see also M_wr_geometry
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

            this->M_HDF5->Flush();
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
            this->M_wr_var(*i);
        }
        // pushing time
        this->M_timeSteps.push_back(time);

        this->M_wr_Xdmf(time);

        if (this->M_multimesh) {
            this->M_wr_geo(); // see also M_wr_geometry
        }

        chrono.stop();

        // Write to file without closing the file
        this->M_HDF5->Flush();

        if (!this->M_procId) std::cout << "         done in " << chrono.diff() << " s." << std::endl;
    }
}

template <typename Mesh>
void HDF5Filter3DMesh<Mesh>::writeGraph()
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

    this->M_HDF5->Write("Graph", "partition_sizes", H5T_NATIVE_INT,
                  M_graph->size(), &partitionSizes[0]);

    // Write partition array size
    this->M_HDF5->Write("Graph", "number_partitions", int(M_graph->size()));

    // Write partition array
    std::stringstream index;
    for (UInt i = 0; i < M_graph->size(); ++i)
    {
        index << i;
        this->M_HDF5->Write("Graph", "partition_graph_" + index.str(),
                      H5T_NATIVE_INT, partitionSizes[i],
                      &(*M_graph)[i][0]);
        index.str(std::string());
        index.clear();
    }
}

template <typename Mesh>
void HDF5Filter3DMesh<Mesh>::writePartition(mesh_ptrtype mesh, std::string& suffix)
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

    this->M_HDF5->Write("Mesh", "Counters.NumPoints" + suffix, (int) mesh->numPoints());
    this->M_HDF5->Write("Mesh", "Counters.NumBPoints" + suffix, (int) mesh->numBPoints());

    this->M_HDF5->Write("Mesh", "Counters.NumVertices" + suffix, (int) mesh->numVertices());
    this->M_HDF5->Write("Mesh", "Counters.NumBVertices" + suffix, (int) mesh->numBVertices());
    this->M_HDF5->Write("Mesh", "Counters.NumGlobalVertices" + suffix, (int) mesh->numGlobalVertices());

    this->M_HDF5->Write("Mesh", "Counters.NumEdges" + suffix, (int) mesh->numEdges());
    this->M_HDF5->Write("Mesh", "Counters.NumBEdges" + suffix, (int) mesh->numBEdges());
    this->M_HDF5->Write("Mesh", "Counters.NumGlobalEdges" + suffix, (int) mesh->numGlobalEdges());

    this->M_HDF5->Write("Mesh", "Counters.NumFaces" + suffix, (int) mesh->numFaces());
    this->M_HDF5->Write("Mesh", "Counters.NumBFaces" + suffix, (int) mesh->numBFaces());
    this->M_HDF5->Write("Mesh", "Counters.NumGlobalFaces" + suffix, (int) mesh->numGlobalFaces());

    this->M_HDF5->Write("Mesh", "Counters.NumVolumes" + suffix, (int) mesh->numVolumes());
    this->M_HDF5->Write("Mesh", "Counters.NumGlobalVolumes" + suffix, (int) mesh->numGlobalVolumes());

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

    this->M_HDF5->Write("Mesh", "Points.x" + suffix, H5T_NATIVE_DOUBLE, numPoints, &pointCoordinates[0][0]);
    this->M_HDF5->Write("Mesh", "Points.y" + suffix, H5T_NATIVE_DOUBLE, numPoints, &pointCoordinates[1][0]);
    this->M_HDF5->Write("Mesh", "Points.z" + suffix, H5T_NATIVE_DOUBLE, numPoints, &pointCoordinates[2][0]);
    this->M_HDF5->Write("Mesh", "Points.f" + suffix, H5T_NATIVE_INT, numPoints, &pointMarkers[0]);
    this->M_HDF5->Write("Mesh", "Points.GlobalId" + suffix, H5T_NATIVE_INT, numPoints, &pointGlobalId[0]);
    this->M_HDF5->Write("Mesh", "Points.BoundaryFlag" + suffix, H5T_NATIVE_INT, numPoints, &pointBoundaryFlags[0]);

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

    this->M_HDF5->Write("Mesh", "Edges.p1" + suffix, H5T_NATIVE_INT, numEdges, &edgePoints[0][0]);
    this->M_HDF5->Write("Mesh", "Edges.p2" + suffix, H5T_NATIVE_INT, numEdges, &edgePoints[1][0]);
    this->M_HDF5->Write("Mesh", "Edges.f" + suffix, H5T_NATIVE_INT, numEdges, &edgeMarkers[0]);
    this->M_HDF5->Write("Mesh", "Edges.GlobalId" + suffix, H5T_NATIVE_INT, numEdges, &edgeGlobalId[0]);
    this->M_HDF5->Write("Mesh", "Edges.BoundaryFlag" + suffix, H5T_NATIVE_INT, numEdges, &edgeBoundaryFlags[0]);

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
        this->M_HDF5->Write("Mesh", "Faces.p" + idx.str() + suffix, H5T_NATIVE_INT, numFaces, &facePoints[k][0]);
        idx.str(std::string());
        idx.clear();
    }
    this->M_HDF5->Write("Mesh", "Faces.f" + suffix, H5T_NATIVE_INT, numFaces, &faceMarkers[0]);
    this->M_HDF5->Write("Mesh", "Faces.GlobalId" + suffix, H5T_NATIVE_INT, numFaces, &faceGlobalId[0]);
    this->M_HDF5->Write("Mesh", "Faces.BoundaryFlag" + suffix, H5T_NATIVE_INT, numFaces, &faceBoundaryFlags[0]);

    this->M_HDF5->Write("Mesh", "Faces.NeighbourId1" + suffix, H5T_NATIVE_INT, numFaces, &faceNeighbourId[0][0]);
    this->M_HDF5->Write("Mesh", "Faces.NeighbourId2" + suffix, H5T_NATIVE_INT, numFaces, &faceNeighbourId[1][0]);
    this->M_HDF5->Write("Mesh", "Faces.NeighbourPos1" + suffix, H5T_NATIVE_INT, numFaces, &faceNeighbourPos[0][0]);
    this->M_HDF5->Write("Mesh", "Faces.NeighbourPos2" + suffix, H5T_NATIVE_INT, numFaces, &faceNeighbourPos[1][0]);

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
        this->M_HDF5->Write("Mesh", "Volumes.p" + idx.str() + suffix, H5T_NATIVE_INT, numVolumes, &volumePoints[k][0]);
        idx.str(std::string());
        idx.clear();
    }
    this->M_HDF5->Write("Mesh", "Volumes.f" + suffix, H5T_NATIVE_INT, numVolumes, &volumeMarkers[0]);
    this->M_HDF5->Write("Mesh", "Volumes.GlobalId" + suffix, H5T_NATIVE_INT, numVolumes, &volumeGlobalId[0]);

    volumePoints.clear();
    volumeMarkers.clear();
    volumeGlobalId.clear();
}

template <typename Mesh>
void HDF5Filter3DMesh<Mesh>::writeSerialMesh()
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
void HDF5Filter3DMesh<Mesh>::writeParallelMesh()
{
    std::stringstream index;
    std::string suffix;

    index << M_comm->MyPID();
    suffix = "." + index.str();

    writePartition(M_parallelMesh, suffix);
}

}

#endif // HDF5FILTER3DMESH
