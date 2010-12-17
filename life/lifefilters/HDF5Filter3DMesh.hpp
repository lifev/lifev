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
  @brief Class derived from Hdf5exporter to provide I/O for the mesh partitions (RegionMesh3D only)

  @date 9-07-2010
  @author Radu Popescu <radu.popescu@epfl.ch>

  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#ifndef HDF5FILTER3DMESH_H
#define HDF5FILTER3DMESH_H 1

#include <life/lifefilters/hdf5exporter.hpp>
#include <life/lifefem/dofInterface3Dto3D.hpp>

namespace LifeV
{

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
template <typename MeshType>
class HDF5Filter3DMesh : public Hdf5exporter<MeshType>
{
public:
    typedef MeshType mesh_Type;
    typedef Hdf5exporter<MeshType> base;
    typedef typename base::meshPtr_Type meshPtr_Type;
    typedef typename base::vectorRaw_Type vector_Type;
    typedef typename base::vectorPtr_Type vectorPtr_Type;

    typedef EpetraExt::HDF5 hdf5_Type;
    typedef boost::shared_ptr<hdf5_Type> hdf5Ptr_Type;
    typedef std::vector<std::vector<Int> > graph_Type;
    typedef boost::shared_ptr<graph_Type> graphPtr_Type;
    typedef boost::shared_ptr<std::vector<meshPtr_Type> > serialMeshPtr_Type;

    typedef DofInterface3Dto3D interface_Type;
    typedef boost::shared_ptr<interface_Type> interfacePtr_Type;
    typedef std::vector<interfacePtr_Type> interfaceVector_Type;
    // The vector contains pointers to each fluid partition's interface with
    // the solid.
    typedef boost::shared_ptr<interfaceVector_Type> interfaceVectorPtr_Type;


    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor for HDF5Filter3DMesh
    HDF5Filter3DMesh() {}

    //! Constructor for HDF5Filter3DMesh
    /*!
      @param dataFile the GetPot data file where you must provide an [exporter] section with:
      "start"     (start index for sections in the hdf5 data structure 0 for 000, 1 for 001 etc.),
      "save"      (how many time steps per postprocessing)
      "multimesh" ( = true if the mesh has to be saved at each post-processing step)
      @param mesh the mesh
      @param the prefix for the case file (ex. "test" for test.case)
      @param the procId determines de CPU id. if negative, it ussemes there is only one processor
    */
    HDF5Filter3DMesh(const GetPot& dataFile, meshPtr_Type mesh, const std::string& prefix, const Int& procId);

    //! Constructor for HDF5Filter3DMesh without prefix and procID
    /*!
      @param dataFile the GetPot data file where you must provide an [exporter] section with:
      "start"     (start index for sections in the hdf5 data structure 0 for 000, 1 for 001 etc.),
      "save"      (how many time steps per postprocessing)
      "multimesh" ( = true if the mesh has to be saved at each post-processing step)
      @param mesh the mesh
    */
    HDF5Filter3DMesh(const GetPot& dataFile, const std::string& prefix);

    //! Destructor for HDF5Filter3DMesh
    virtual ~HDF5Filter3DMesh() {}

    //@}

    //! @name Pubic Methods
    //@{
    virtual void postProcess(const Real& time);

    //! Add the partition graph to the post processing data file
    /*!
      Add the partition graph to the post processing data file.
      \param graph - shared_ptr<vector<vector<Int> > > - shared pointer to the partition graph data structure
      (as returned by partitionMesh::graph() )
      \param comm - Epetra_Comm* - raw pointer to the Epetra communicator to be used
    */
    void addPartitionGraph(const graphPtr_Type& graph, boost::shared_ptr<Epetra_Comm>& comm)
    {M_graph = graph; M_comm = comm;}

    //! Add all of the mesh partitions to the post processing data file (serial operation)
    /*!
      Add all of the mesh partitions to the post processing data file.
      \param meshPointer - shared_ptr<vector<shared_ptr<MeshType> > > - shared pointer to the vector storing
      pointers to the mesh partitions (as returned by partitionMesh::meshAllPartitions() )
      \param comm - Epetra_Comm* - raw pointer to the Epetra communicator to be used
    */
    void addMeshPartitionAll(const serialMeshPtr_Type& meshPointer, boost::shared_ptr<Epetra_Comm>& comm)
    {M_serialMesh = meshPointer; M_parallelMesh.reset(); M_comm = comm;}

    //! Add to HDF5 file the mesh partition that belongs to the current process (parallel operation)
    /*!
      After the mesh partition is loaded from the HDF5 file, the simulation is run and
      HDF5Filter3DMesh::postProcess() is called, the original contents of the HDF5 will be lost.
      To keep mesh partition, call:
      HDF5Filter3DMesh::addMyMeshPartition() before calling HDF5Filter3DMesh::postProcess();
      \param meshPointer - shared_ptr<Mesh> - shared pointer to a mesh partition (as returned by
      partitionMesh::mesh() )
      \param comm - Epetra_Comm* - raw pointer to the Epetra communicator to be used
    */
    void addMyMeshPartition(const meshPtr_Type& meshPointer, boost::shared_ptr<Epetra_Comm>& comm)
    {/*M_parallelMesh = meshPointer; M_serialMesh.reset(); M_comm = comm;*/}

    //! Add a DOF interface for writing to file
    /*!
      Add a DOF interface to the member vector M_interfaceData, for writing to the HDF5
      file. Call once for each interface that is to be written.
    */
    void addDOFInterface(const interfaceVectorPtr_Type& interfaces,
                         const std::string& type,
                         const Int& firstInterfaceFlag,
                         const Int& secondInterfaceFlag,
                         const boost::shared_ptr<Epetra_Comm>& comm);

    // TODO: Write a replacement methods that returns the graph vector directly,
    //       like with getMeshPartition and getStoredInterface
    //! Load the partitioned graph from a HDF5 file into a partitionMesh object
    /*!
      \param graph - shared_ptr<vector<vector<int> > > - a shared pointer to the graph data structure in the
      partitionMesh object (as returned by partitionMesh::graph() )
      \param comm - Epetra_Comm* - a raw pointer to the Epetra communicator to be used
    */
    void loadGraph(graphPtr_Type graph, boost::shared_ptr<Epetra_Comm>& comm);

    //! Load a mesh partition according to the MPI PID
    /*!
      This method is to be used in parallel. The mesh partition corresponding to the
      MPI process id is loaded into the meshPartition object.
      Could be used with any object that contains a RegionMesh3D member but it is
      primarily used to load a mesh partition into the partitionMesh object.
      \param meshPartition - shared_ptr<Mesh> - shared pointer to mesh partition object
      \param comm -shared_ptr<Epetra_Comm> - shared pointer to the Epetra communicator to be used
    */
    void __attribute__((__deprecated__)) loadMyPartition(meshPtr_Type meshPartition,
                                                         boost::shared_ptr<Epetra_Comm>& comm);

    //! Get the number of stored DOF interfaces
    Int queryStoredInterfaceNumber();

    //! Get the types of the stored DOF interfaces
    std::vector<std::string>& queryStoredInterfaceTypes();

    //! Return a pointer to the mesh partition that corresponds to the current MPI rank
    boost::shared_ptr<MeshType>& getMeshPartition();

    //! Return a pointer to the k-th interface stored inside the file
    boost::shared_ptr< std::map<UInt, UInt> >& getStoredInterface(Int k) ;

    // When reading back partitions and interfaces, the comm member must be set explicitly.
    // This is intended for use with getMeshPartition and getStoredInterface
    //! Set the M_comm data member.
    void setComm( const boost::shared_ptr<Epetra_Comm>& comm ) {M_comm = comm;}
    //@}

private:

//! @name Private methods
//@{
    // The following private methods are for writing the partitioned graph
    // and mesh to the output file
    void writeGraph();
    void writePartition(meshPtr_Type partition, std::string& suffix);
    void writeParallelMesh();
    void writeSerialMesh();
    void writeInterfaces();
//@}

    serialMeshPtr_Type                     M_serialMesh;
    meshPtr_Type                           M_parallelMesh;
    graphPtr_Type                          M_graph;
    // Use a vector to store the data of all interfaces we wish to write to disk
    std::vector<interfaceVectorPtr_Type>  M_DOFInterfaces;
    std::vector<std::string>               M_interfaceTypes;
    std::vector<Int>                       M_firstInterfaceFlags;
    std::vector<Int>                       M_secondInterfaceFlags;
    boost::shared_ptr<Epetra_Comm>         M_comm;

};

// ===================================================
// Constructors
// ===================================================

template<typename MeshType>
HDF5Filter3DMesh<MeshType>::HDF5Filter3DMesh(const GetPot& dataFile, meshPtr_Type mesh,
                                             const std::string& prefix,
                                             const Int& procId) :
    base                ( dataFile, mesh, prefix, procId )
{
}

template<typename MeshType>
HDF5Filter3DMesh<MeshType>::HDF5Filter3DMesh(const GetPot& dataFile, const std::string& prefix):
    base                ( dataFile, prefix )
{
}

// ===================================================
// Pubic Methods
// ===================================================

template<typename MeshType>
void HDF5Filter3DMesh<MeshType>::addDOFInterface(const interfaceVectorPtr_Type& interfaces,
                                                 const std::string& type,
                                                 const Int& firstInterfaceFlag,
                                                 const Int& secondInterfaceFlag,
                                                 const boost::shared_ptr<Epetra_Comm>& comm)
{
    M_DOFInterfaces.push_back(interfaces);
    M_interfaceTypes.push_back(type);
    M_firstInterfaceFlags.push_back(firstInterfaceFlag);
    M_secondInterfaceFlags.push_back(secondInterfaceFlag);

    M_comm = comm;
}

template<typename MeshType>
void HDF5Filter3DMesh<MeshType>::loadGraph(graphPtr_Type graph, boost::shared_ptr<Epetra_Comm>& comm)
{
    if (this->M_HDF5.get() == 0)
    {
        this->M_HDF5.reset(new hdf5_Type(*comm));
    }
    if (! this->M_HDF5->IsOpen())
    {
        this->M_HDF5->Open(this->M_postDir + this->M_prefix + ".h5", H5F_ACC_RDONLY);
    }

    Int nPartitions;

    this->M_HDF5->Read("Graph", "number_partitions", nPartitions);

    std::vector<Int> partitionSizes(nPartitions);
    this->M_HDF5->Read("Graph", "partition_sizes", H5T_NATIVE_INT, nPartitions,
                       &partitionSizes[0]);

    graph->resize(0);
    graph->reserve(nPartitions);

    std::vector<Int> partBuffer;
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

// TODO: MARK AS DEPRECATED AND REMOVE -> FIX TEST PARTITION AND MAYBE THE FSI PARTITIONER
/*template<typename MeshType>
void HDF5Filter3DMesh<MeshType>::loadMyPartition(meshPtr_Type meshPartition,
                                                 boost::shared_ptr<Epetra_Comm>& comm)
{
    UInt elementNodes, faceNodes;
    switch (MeshType::ElementShape::Shape)
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
        this->M_HDF5.reset(new hdf5_Type(*comm));
    }
    if (! this->M_HDF5->IsOpen())
    {
        this->M_HDF5->Open(this->M_postDir + this->M_prefix + ".h5", H5F_ACC_RDONLY);
    }

    std::stringstream index;
    index << comm->MyPID();
    std::string suffix = "." + index.str();

    // Read counters from file and set them in the mesh

    // Points
    Int numPoints;
    this->M_HDF5->Read("Mesh", "Counters.NumPoints" + suffix, numPoints);
    meshPartition->setMaxNumPoints(numPoints, true);
    Int numBPoints;
    this->M_HDF5->Read("Mesh", "Counters.NumBPoints" + suffix, numBPoints);
    meshPartition->setNumBPoints(numBPoints);

    // Vertices
    Int numVertices;
    this->M_HDF5->Read("Mesh", "Counters.NumVertices" + suffix, numVertices);
    meshPartition->setNumVertices(numVertices);
    Int numBVertices;
    this->M_HDF5->Read("Mesh", "Counters.NumBVertices" + suffix, numBVertices);
    meshPartition->setNumBVertices(numBVertices);
    Int numGlobalVertices;
    this->M_HDF5->Read("Mesh", "Counters.NumGlobalVertices" + suffix, numGlobalVertices);
    meshPartition->setNumGlobalVertices(numGlobalVertices);

    // Edges
    Int numEdges;
    this->M_HDF5->Read("Mesh", "Counters.NumEdges" + suffix, numEdges);
    meshPartition->setNumEdges(numEdges);
    meshPartition->setMaxNumEdges(numEdges);
    Int numBEdges;
    this->M_HDF5->Read("Mesh", "Counters.NumBEdges" + suffix, numBEdges);
    meshPartition->setNumBEdges(numBEdges);
    Int numGlobalEdges;
    this->M_HDF5->Read("Mesh", "Counters.NumGlobalEdges" + suffix, numGlobalEdges);
    meshPartition->setMaxNumGlobalEdges(numGlobalEdges);

    // Faces
    Int numFaces;
    this->M_HDF5->Read("Mesh", "Counters.NumFaces" + suffix, numFaces);
    meshPartition->setNumFaces(numFaces);
    meshPartition->setMaxNumFaces(numFaces);
    Int numBFaces;
    this->M_HDF5->Read("Mesh", "Counters.NumBFaces" + suffix, numBFaces);
    meshPartition->setNumBFaces(numBFaces);
    Int numGlobalFaces;
    this->M_HDF5->Read("Mesh", "Counters.NumGlobalFaces" + suffix, numGlobalFaces);
    meshPartition->setMaxNumGlobalFaces(numGlobalFaces);

    // Volumes
    Int numVolumes;
    this->M_HDF5->Read("Mesh", "Counters.NumVolumes" + suffix, numVolumes);
    meshPartition->setMaxNumVolumes(numVolumes, true);
    Int numGlobalVolumes;
    this->M_HDF5->Read("Mesh", "Counters.NumGlobalVolumes" + suffix, numGlobalVolumes);
    meshPartition->setMaxNumGlobalVolumes(numGlobalVolumes);

    // Read the list of points
    std::vector<Real> tmpVectorDouble(numPoints);
    std::vector<std::vector<Real> > pointCoordinates(3, tmpVectorDouble);

    std::vector<Int> pointMarkers(numPoints);
    std::vector<Int> pointBoundaryFlags(numPoints);
    std::vector<Int> pointGlobalId(numPoints);

    this->M_HDF5->Read("Mesh", "Points.x" + suffix, H5T_NATIVE_DOUBLE, numPoints, &pointCoordinates[0][0]);
    this->M_HDF5->Read("Mesh", "Points.y" + suffix, H5T_NATIVE_DOUBLE, numPoints, &pointCoordinates[1][0]);
    this->M_HDF5->Read("Mesh", "Points.z" + suffix, H5T_NATIVE_DOUBLE, numPoints, &pointCoordinates[2][0]);
    this->M_HDF5->Read("Mesh", "Points.f" + suffix, H5T_NATIVE_INT, numPoints, &pointMarkers[0]);
    this->M_HDF5->Read("Mesh", "Points.BoundaryFlag" + suffix, H5T_NATIVE_INT, numPoints,
                       &pointBoundaryFlags[0]);
    this->M_HDF5->Read("Mesh", "Points.GlobalId" + suffix, H5T_NATIVE_INT, numPoints, &pointGlobalId[0]);

    meshPartition->pointList.reserve(numPoints);
    meshPartition->_bPoints.reserve(meshPartition->numBPoints());

    typename MeshType::PointType *pp = 0;

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
    std::vector<Int> tmpVectorInt(numEdges);
    std::vector<std::vector<Int> > edgePoints(2, tmpVectorInt);

    std::vector<Int> edgeMarkers(numEdges);
    std::vector<Int> edgeGlobalId(numEdges);
    std::vector<Int> edgeBoundaryFlags(numEdges);

    this->M_HDF5->Read("Mesh", "Edges.p1" + suffix, H5T_NATIVE_INT, numEdges, &edgePoints[0][0]);
    this->M_HDF5->Read("Mesh", "Edges.p2" + suffix, H5T_NATIVE_INT, numEdges, &edgePoints[1][0]);
    this->M_HDF5->Read("Mesh", "Edges.f" + suffix, H5T_NATIVE_INT, numEdges, &edgeMarkers[0]);
    this->M_HDF5->Read("Mesh", "Edges.GlobalId" + suffix, H5T_NATIVE_INT, numEdges, &edgeGlobalId[0]);
    this->M_HDF5->Read("Mesh", "Edges.BoundaryFlag" + suffix, H5T_NATIVE_INT, numEdges,
                       &edgeBoundaryFlags[0]);

    meshPartition->edgeList.reserve(numEdges);

    typename MeshType::EdgeType *pe;

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
    std::vector<std::vector<Int> > facePoints(faceNodes, tmpVectorInt);

    std::vector<Int> faceMarkers(numFaces);
    std::vector<Int> faceGlobalId(numFaces);
    std::vector<Int> faceBoundaryFlags(numFaces);

    std::vector<std::vector<Int> > faceNeighbourId(2, tmpVectorInt);
    std::vector<std::vector<Int> > faceNeighbourPos(2, tmpVectorInt);

    std::stringstream idx;
    for (UInt k = 0; k < faceNodes; ++k)
    {
        idx << k + 1;
        this->M_HDF5->Read("Mesh", "Faces.p" + idx.str() + suffix, H5T_NATIVE_INT, numFaces,
                           &facePoints[k][0]);
        idx.str(std::string());
        idx.clear();
    }
    this->M_HDF5->Read("Mesh", "Faces.f" + suffix, H5T_NATIVE_INT, numFaces, &faceMarkers[0]);
    this->M_HDF5->Read("Mesh", "Faces.GlobalId" + suffix, H5T_NATIVE_INT, numFaces, &faceGlobalId[0]);
    this->M_HDF5->Read("Mesh", "Faces.BoundaryFlag" + suffix, H5T_NATIVE_INT, numFaces,
                       &faceBoundaryFlags[0]);

    this->M_HDF5->Read("Mesh", "Faces.NeighbourId1" + suffix, H5T_NATIVE_INT, numFaces,
                       &faceNeighbourId[0][0]);
    this->M_HDF5->Read("Mesh", "Faces.NeighbourId2" + suffix, H5T_NATIVE_INT, numFaces,
                       &faceNeighbourId[1][0]);
    this->M_HDF5->Read("Mesh", "Faces.NeighbourPos1" + suffix, H5T_NATIVE_INT, numFaces,
                       &faceNeighbourPos[0][0]);
    this->M_HDF5->Read("Mesh", "Faces.NeighbourPos2" + suffix, H5T_NATIVE_INT, numFaces,
                       &faceNeighbourPos[1][0]);


    typename MeshType::FaceType *pf = 0;

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
    std::vector<std::vector<Int> > volumePoints(elementNodes, tmpVectorInt);

    std::vector<Int> volumeMarkers(numVolumes);
    std::vector<Int> volumeGlobalId(numVolumes);

    for (UInt k = 0; k < elementNodes; ++k)
    {
        idx << k + 1;
        this->M_HDF5->Read("Mesh", "Volumes.p" + idx.str() + suffix, H5T_NATIVE_INT, numVolumes,
                           &volumePoints[k][0]);
        idx.str(std::string());
        idx.clear();
    }
    this->M_HDF5->Read("Mesh", "Volumes.f" + suffix, H5T_NATIVE_INT, numVolumes, &volumeMarkers[0]);
    this->M_HDF5->Read("Mesh", "Volumes.GlobalId" + suffix, H5T_NATIVE_INT, numVolumes, &volumeGlobalId[0]);

    meshPartition->volumeList.reserve(numVolumes);

    typename MeshType::VolumeType *pv = 0;

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

}*/

template<typename MeshType>
void HDF5Filter3DMesh<MeshType>::postProcess(const Real& time)
{
    if ( this->M_HDF5.get() == 0)
    {
        if (this->M_listData.size() != 0)
        {
            this->M_HDF5.reset(new hdf5_Type(this->M_listData.begin()->storedArray()->Comm()));
        }
        else
        {
            this->M_HDF5.reset(new hdf5_Type(*M_comm));
        }
        this->M_outputFileName=this->M_prefix+".h5";
        this->M_HDF5->Create(this->M_postDir+this->M_outputFileName);

        // write empty xdmf file
        this->writeInitXdmf();

        if (!this->M_multimesh)
        {
            if (this->M_listData.size() != 0)
            {
                this->writeGeometry(); // see also writeGeometrymetry
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
            if (M_DOFInterfaces.size() != 0)
            {
                writeInterfaces();
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
            this->writeVariable(*i);
        }
        // pushing time
        this->M_timeSteps.push_back(time);

        this->writeXdmf(time);

        if (this->M_multimesh)
        {
            this->writeGeometry(); // see also writeGeometry
        }

        chrono.stop();

        // Write to file without closing the file
        this->M_HDF5->Flush();

        if (!this->M_procId) std::cout << "         done in " << chrono.diff() << " s." << std::endl;
    }
}

template <typename MeshType>
int HDF5Filter3DMesh<MeshType>::queryStoredInterfaceNumber()
{
    if (this->M_HDF5.get() == 0)
    {
        this->M_HDF5.reset(new hdf5_Type(*M_comm));
    }
    if (! this->M_HDF5->IsOpen())
    {
        this->M_HDF5->Open(this->M_postDir + this->M_prefix + ".h5", H5F_ACC_RDONLY);
    }

    Int storedInterfaceNumber;
    this->M_HDF5->Read("Interfaces", "Number", storedInterfaceNumber);

    return storedInterfaceNumber;
}

template <typename MeshType>
std::vector<std::string>& HDF5Filter3DMesh<MeshType>::queryStoredInterfaceTypes()
{
    if (this->M_HDF5.get() == 0)
    {
        this->M_HDF5.reset(new hdf5_Type(*M_comm));
    }
    if (! this->M_HDF5->IsOpen())
    {
        this->M_HDF5->Open(this->M_postDir + this->M_prefix + ".h5", H5F_ACC_RDONLY);
    }

    int storedInterfaceNumber;
    this->M_HDF5->Read("Interfaces", "Number", storedInterfaceNumber);

    std::vector<std::string> storedInterfaceTypes(storedInterfaceNumber);
    for (Int k = 0; k < storedInterfaceNumber; ++k)
    {
        std::stringstream idx;
        idx << k;
        this->M_HDF5->Read("Interfaces", "Type." + idx.str(), storedInterfaceTypes[k]);
    }

    return storedInterfaceTypes;
}

template <typename MeshType>
boost::shared_ptr<MeshType>& HDF5Filter3DMesh<MeshType>::getMeshPartition()
{
    boost::shared_ptr<MeshType> tempMesh(new MeshType);

    UInt elementNodes, faceNodes;
    switch (MeshType::ElementShape::Shape)
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
        this->M_HDF5.reset(new hdf5_Type(*M_comm));
    }
    if (! this->M_HDF5->IsOpen())
    {
        this->M_HDF5->Open(this->M_postDir + this->M_prefix + ".h5", H5F_ACC_RDONLY);
    }

    std::stringstream index;
    index << M_comm->MyPID();
    std::string suffix = "." + index.str();

    // Read counters from file and set them in the mesh

    // Points
    Int numPoints;
    this->M_HDF5->Read("Mesh", "Counters.NumPoints" + suffix, numPoints);
    tempMesh->setMaxNumPoints(numPoints, true);
    Int numBPoints;
    this->M_HDF5->Read("Mesh", "Counters.NumBPoints" + suffix, numBPoints);
    tempMesh->setNumBPoints(numBPoints);

    // Vertices
    Int numVertices;
    this->M_HDF5->Read("Mesh", "Counters.NumVertices" + suffix, numVertices);
    tempMesh->setNumVertices(numVertices);
    Int numBVertices;
    this->M_HDF5->Read("Mesh", "Counters.NumBVertices" + suffix, numBVertices);
    tempMesh->setNumBVertices(numBVertices);
    Int numGlobalVertices;
    this->M_HDF5->Read("Mesh", "Counters.NumGlobalVertices" + suffix, numGlobalVertices);
    tempMesh->setNumGlobalVertices(numGlobalVertices);

    // Edges
    Int numEdges;
    this->M_HDF5->Read("Mesh", "Counters.NumEdges" + suffix, numEdges);
    tempMesh->setNumEdges(numEdges);
    tempMesh->setMaxNumEdges(numEdges);
    Int numBEdges;
    this->M_HDF5->Read("Mesh", "Counters.NumBEdges" + suffix, numBEdges);
    tempMesh->setNumBEdges(numBEdges);
    Int numGlobalEdges;
    this->M_HDF5->Read("Mesh", "Counters.NumGlobalEdges" + suffix, numGlobalEdges);
    tempMesh->setMaxNumGlobalEdges(numGlobalEdges);

    // Faces
    Int numFaces;
    this->M_HDF5->Read("Mesh", "Counters.NumFaces" + suffix, numFaces);
    tempMesh->setNumFaces(numFaces);
    tempMesh->setMaxNumFaces(numFaces);
    Int numBFaces;
    this->M_HDF5->Read("Mesh", "Counters.NumBFaces" + suffix, numBFaces);
    tempMesh->setNumBFaces(numBFaces);
    Int numGlobalFaces;
    this->M_HDF5->Read("Mesh", "Counters.NumGlobalFaces" + suffix, numGlobalFaces);
    tempMesh->setMaxNumGlobalFaces(numGlobalFaces);

    // Volumes
    Int numVolumes;
    this->M_HDF5->Read("Mesh", "Counters.NumVolumes" + suffix, numVolumes);
    tempMesh->setMaxNumVolumes(numVolumes, true);
    Int numGlobalVolumes;
    this->M_HDF5->Read("Mesh", "Counters.NumGlobalVolumes" + suffix, numGlobalVolumes);
    tempMesh->setMaxNumGlobalVolumes(numGlobalVolumes);

    // Read the list of points
    std::vector<Real> tmpVectorDouble(numPoints);
    std::vector<std::vector<Real> > pointCoordinates(3, tmpVectorDouble);

    std::vector<Int> pointMarkers(numPoints);
    std::vector<Int> pointBoundaryFlags(numPoints);
    std::vector<Int> pointGlobalId(numPoints);

    this->M_HDF5->Read("Mesh", "Points.x" + suffix, H5T_NATIVE_DOUBLE, numPoints, &pointCoordinates[0][0]);
    this->M_HDF5->Read("Mesh", "Points.y" + suffix, H5T_NATIVE_DOUBLE, numPoints, &pointCoordinates[1][0]);
    this->M_HDF5->Read("Mesh", "Points.z" + suffix, H5T_NATIVE_DOUBLE, numPoints, &pointCoordinates[2][0]);
    this->M_HDF5->Read("Mesh", "Points.f" + suffix, H5T_NATIVE_INT, numPoints, &pointMarkers[0]);
    this->M_HDF5->Read("Mesh", "Points.BoundaryFlag" + suffix, H5T_NATIVE_INT, numPoints,
                       &pointBoundaryFlags[0]);
    this->M_HDF5->Read("Mesh", "Points.GlobalId" + suffix, H5T_NATIVE_INT, numPoints, &pointGlobalId[0]);

    tempMesh->pointList.reserve(numPoints);
    tempMesh->_bPoints.reserve(tempMesh->numBPoints());

    typename MeshType::PointType *pp = 0;

    for (UInt j = 0; j < numPoints; ++j)
    {
        pp = &(tempMesh->addPoint(bool(pointBoundaryFlags[j])));
        pp->setMarker(pointMarkers[j]);
        pp->x() = pointCoordinates[0][j];
        pp->y() = pointCoordinates[1][j];
        pp->z() = pointCoordinates[2][j];
        pp->setLocalId(j + 1);
        pp->setId(pointGlobalId[j]);

        tempMesh->localToGlobalNode().insert(std::make_pair(j + 1, pointGlobalId[j]));
        tempMesh->globalToLocalNode().insert(std::make_pair(pointGlobalId[j], j + 1));
    }

    pointCoordinates.clear();
    pointMarkers.clear();
    pointBoundaryFlags.clear();
    pointGlobalId.clear();

    // Read the list of edges
    std::vector<Int> tmpVectorInt(numEdges);
    std::vector<std::vector<Int> > edgePoints(2, tmpVectorInt);

    std::vector<Int> edgeMarkers(numEdges);
    std::vector<Int> edgeGlobalId(numEdges);
    std::vector<Int> edgeBoundaryFlags(numEdges);

    this->M_HDF5->Read("Mesh", "Edges.p1" + suffix, H5T_NATIVE_INT, numEdges, &edgePoints[0][0]);
    this->M_HDF5->Read("Mesh", "Edges.p2" + suffix, H5T_NATIVE_INT, numEdges, &edgePoints[1][0]);
    this->M_HDF5->Read("Mesh", "Edges.f" + suffix, H5T_NATIVE_INT, numEdges, &edgeMarkers[0]);
    this->M_HDF5->Read("Mesh", "Edges.GlobalId" + suffix, H5T_NATIVE_INT, numEdges, &edgeGlobalId[0]);
    this->M_HDF5->Read("Mesh", "Edges.BoundaryFlag" + suffix, H5T_NATIVE_INT, numEdges,
                       &edgeBoundaryFlags[0]);

    tempMesh->edgeList.reserve(numEdges);

    typename MeshType::EdgeType *pe;

    for (UInt j = 0; j < numEdges; ++j)
    {
        pe = &(tempMesh->addEdge(edgeBoundaryFlags[j]));
        pe->setLocalId(j + 1);
        pe->setId(edgeGlobalId[j]);
        pe->setPoint(1, tempMesh->point(edgePoints[0][j]));
        pe->setPoint(2, tempMesh->point(edgePoints[1][j]));
        pe->setMarker(edgeMarkers[j]);
    }

    edgePoints.clear();
    edgeMarkers.clear();
    edgeGlobalId.clear();
    edgeBoundaryFlags.clear();

    // Read the list of faces
    tmpVectorInt.resize(numFaces);
    std::vector<std::vector<Int> > facePoints(faceNodes, tmpVectorInt);

    std::vector<Int> faceMarkers(numFaces);
    std::vector<Int> faceGlobalId(numFaces);
    std::vector<Int> faceBoundaryFlags(numFaces);

    std::vector<std::vector<Int> > faceNeighbourId(2, tmpVectorInt);
    std::vector<std::vector<Int> > faceNeighbourPos(2, tmpVectorInt);

    std::stringstream idx;
    for (UInt k = 0; k < faceNodes; ++k)
    {
        idx << k + 1;
        this->M_HDF5->Read("Mesh", "Faces.p" + idx.str() + suffix, H5T_NATIVE_INT, numFaces,
                           &facePoints[k][0]);
        idx.str(std::string());
        idx.clear();
    }
    this->M_HDF5->Read("Mesh", "Faces.f" + suffix, H5T_NATIVE_INT, numFaces, &faceMarkers[0]);
    this->M_HDF5->Read("Mesh", "Faces.GlobalId" + suffix, H5T_NATIVE_INT, numFaces, &faceGlobalId[0]);
    this->M_HDF5->Read("Mesh", "Faces.BoundaryFlag" + suffix, H5T_NATIVE_INT, numFaces,
                       &faceBoundaryFlags[0]);

    this->M_HDF5->Read("Mesh", "Faces.NeighbourId1" + suffix, H5T_NATIVE_INT, numFaces,
                       &faceNeighbourId[0][0]);
    this->M_HDF5->Read("Mesh", "Faces.NeighbourId2" + suffix, H5T_NATIVE_INT, numFaces,
                       &faceNeighbourId[1][0]);
    this->M_HDF5->Read("Mesh", "Faces.NeighbourPos1" + suffix, H5T_NATIVE_INT, numFaces,
                       &faceNeighbourPos[0][0]);
    this->M_HDF5->Read("Mesh", "Faces.NeighbourPos2" + suffix, H5T_NATIVE_INT, numFaces,
                       &faceNeighbourPos[1][0]);


    typename MeshType::FaceType *pf = 0;

    tempMesh->faceList.reserve(numFaces);

    for (UInt j = 0; j < numFaces; ++j)
    {
        pf = &(tempMesh->addFace(faceBoundaryFlags[j]));
        pf->setLocalId(j + 1);
        pf->setId(faceGlobalId[j]);

        pf->ad_first() = faceNeighbourId[0][j];
        pf->ad_second() = faceNeighbourId[1][j];
        pf->pos_first() = faceNeighbourPos[0][j];
        pf->pos_second() = faceNeighbourPos[1][j];

        pf->setMarker(faceMarkers[j]);
        for (UInt k = 0; k < faceNodes; ++k)
        {
            pf->setPoint(k + 1, tempMesh->point(facePoints[k][j]));
        }
    }

    tempMesh->setLinkSwitch("HAS_ALL_FACES");
    tempMesh->setLinkSwitch("FACES_HAVE_ADIACENCY");

    facePoints.clear();
    faceMarkers.clear();
    faceGlobalId.clear();
    faceBoundaryFlags.clear();
    faceNeighbourId.clear();
    faceNeighbourPos.clear();

    // Read the list of volumes
    tmpVectorInt.resize(numVolumes);
    std::vector<std::vector<Int> > volumePoints(elementNodes, tmpVectorInt);

    std::vector<Int> volumeMarkers(numVolumes);
    std::vector<Int> volumeGlobalId(numVolumes);

    for (UInt k = 0; k < elementNodes; ++k)
    {
        idx << k + 1;
        this->M_HDF5->Read("Mesh", "Volumes.p" + idx.str() + suffix, H5T_NATIVE_INT, numVolumes,
                           &volumePoints[k][0]);
        idx.str(std::string());
        idx.clear();
    }
    this->M_HDF5->Read("Mesh", "Volumes.f" + suffix, H5T_NATIVE_INT, numVolumes, &volumeMarkers[0]);
    this->M_HDF5->Read("Mesh", "Volumes.GlobalId" + suffix, H5T_NATIVE_INT, numVolumes, &volumeGlobalId[0]);

    tempMesh->volumeList.reserve(numVolumes);

    typename MeshType::VolumeType *pv = 0;

    for (UInt j = 0; j < numVolumes; ++j)
    {
        pv = &(tempMesh->addVolume());
        pv->setId(volumeGlobalId[j]);
        pv->setLocalId(j + 1);
        for (UInt k = 0; k < elementNodes; ++k)
        {
            pv->setPoint(k + 1, tempMesh->point(volumePoints[k][j]));
        }
        pv->setMarker(volumeMarkers[j]);
    }

    volumePoints.clear();
    volumeMarkers.clear();
    volumeGlobalId.clear();

    tempMesh->updateElementEdges(false, false);
    tempMesh->updateElementFaces(false, false);

    return tempMesh;
}

template <typename MeshType>
boost::shared_ptr< std::map<UInt, UInt> >& HDF5Filter3DMesh<MeshType>::getStoredInterface(int k)
{
    if (this->M_HDF5.get() == 0)
    {
        this->M_HDF5.reset(new hdf5_Type(*M_comm));
    }
    if (! this->M_HDF5->IsOpen())
    {
        this->M_HDF5->Open(this->M_postDir + this->M_prefix + ".h5", H5F_ACC_RDONLY);
    }

    int myRank = M_comm->MyPID();

    std::stringstream idx;
    idx << k << "." << myRank;

    boost::shared_ptr<std::map<UInt, UInt> > interface(new std::map<UInt, UInt>);

    Int size;
    this->M_HDF5->Read("Interfaces", "Size." + idx.str(), size);

    std::vector<UInt> keyVector(size);
    std::vector<UInt> valueVector(size);

    this->M_HDF5->Read("Interfaces", "Key." + idx.str(), H5T_NATIVE_INT, size,
                       &keyVector[0]);
    this->M_HDF5->Read("Interfaces", "Value." + idx.str(), H5T_NATIVE_INT, size,
                       &valueVector[0]);

    for (UInt i = 0; i < size; ++i)
    {
        interface->insert(std::make_pair(keyVector[i], valueVector[i]));
    }

    return interface;
}

// ========================
// Private methods
// ========================

template <typename MeshType>
void HDF5Filter3DMesh<MeshType>::writeGraph()
{
    std::vector<Int> partitionSizes;
    Int size, maxSize = 0;

    // Calculate the maximum size of the partitions and store partition
    // sizes
    partitionSizes.reserve(M_graph->size());
    for (std::vector<std::vector<Int> >::iterator it = M_graph->begin();
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
    this->M_HDF5->Write("Graph", "number_partitions", static_cast<Int>((M_graph->size())));

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

template <typename MeshType>
void HDF5Filter3DMesh<MeshType>::writePartition(meshPtr_Type mesh, std::string& suffix)
{
    UInt elementNodes, faceNodes;
    switch (MeshType::ElementShape::Shape)
    {
    case HEXA:
        elementNodes = 8;
        faceNodes    = 4;
        break;
    case TETRA:
        elementNodes = 4;
        faceNodes    = 3;
    }

    this->M_HDF5->Write("Mesh", "Counters.NumPoints" + suffix, static_cast<Int>(mesh->numPoints()));
    this->M_HDF5->Write("Mesh", "Counters.NumBPoints" + suffix, static_cast<Int>(mesh->numBPoints()));

    this->M_HDF5->Write("Mesh", "Counters.NumVertices" + suffix, static_cast<Int>(mesh->numVertices()));
    this->M_HDF5->Write("Mesh", "Counters.NumBVertices" + suffix, static_cast<Int>(mesh->numBVertices()));
    this->M_HDF5->Write("Mesh", "Counters.NumGlobalVertices" + suffix,
                        static_cast<Int>(mesh->numGlobalVertices()));

    this->M_HDF5->Write("Mesh", "Counters.NumEdges" + suffix, static_cast<Int>(mesh->numEdges()));
    this->M_HDF5->Write("Mesh", "Counters.NumBEdges" + suffix, static_cast<Int>(mesh->numBEdges()));
    this->M_HDF5->Write("Mesh", "Counters.NumGlobalEdges" + suffix, static_cast<Int>(mesh->numGlobalEdges()));

    this->M_HDF5->Write("Mesh", "Counters.NumFaces" + suffix, static_cast<Int>(mesh->numFaces()));
    this->M_HDF5->Write("Mesh", "Counters.NumBFaces" + suffix, static_cast<Int>(mesh->numBFaces()));
    this->M_HDF5->Write("Mesh", "Counters.NumGlobalFaces" + suffix, static_cast<Int>(mesh->numGlobalFaces()));

    this->M_HDF5->Write("Mesh", "Counters.NumVolumes" + suffix, static_cast<Int>(mesh->numVolumes()));
    this->M_HDF5->Write("Mesh", "Counters.NumGlobalVolumes" + suffix,
                        static_cast<Int>(mesh->numGlobalVolumes()));

    Int numPoints = mesh->numPoints();
    std::vector<Real> tmpVectorDouble(numPoints);
    std::vector<std::vector<Real> > pointCoordinates(3, tmpVectorDouble);

    std::vector<Int> pointMarkers(numPoints);
    std::vector<Int> pointGlobalId(numPoints);
    std::vector<Int> pointBoundaryFlags(numPoints);

    std::map<Int, Int>::iterator it = mesh->localToGlobalNode().begin();
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
    this->M_HDF5->Write("Mesh", "Points.BoundaryFlag" + suffix, H5T_NATIVE_INT, numPoints,
                        &pointBoundaryFlags[0]);

    pointCoordinates.clear();
    pointMarkers.clear();
    pointBoundaryFlags.clear();
    pointGlobalId.clear();

    Int numEdges = mesh->numEdges();
    std::vector<Int> tmpVectorInt(numEdges);
    std::vector<std::vector<Int> > edgePoints(2,tmpVectorInt);

    std::vector<Int> edgeMarkers(numEdges);
    std::vector<Int> edgeGlobalId(numEdges);
    std::vector<Int> edgeBoundaryFlags(numEdges);

    for (Int j = 0; j < numEdges; ++j)
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
    this->M_HDF5->Write("Mesh", "Edges.BoundaryFlag" + suffix, H5T_NATIVE_INT, numEdges,
                        &edgeBoundaryFlags[0]);

    edgePoints.clear();
    edgeMarkers.clear();
    edgeGlobalId.clear();
    edgeBoundaryFlags.clear();

    Int numFaces = mesh->numFaces();
    tmpVectorInt.resize(numFaces);
    std::vector<std::vector<Int> > facePoints(faceNodes, tmpVectorInt);

    std::vector<Int> faceMarkers(numFaces);
    std::vector<Int> faceGlobalId(numFaces);
    std::vector<Int> faceBoundaryFlags(numFaces);

    std::vector<std::vector<Int> > faceNeighbourId(2, tmpVectorInt);
    std::vector<std::vector<Int> > faceNeighbourPos(2, tmpVectorInt);

    for (Int j = 0; j < numFaces; ++j)
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
        this->M_HDF5->Write("Mesh", "Faces.p" + idx.str() + suffix, H5T_NATIVE_INT, numFaces,
                            &facePoints[k][0]);
        idx.str(std::string());
        idx.clear();
    }
    this->M_HDF5->Write("Mesh", "Faces.f" + suffix, H5T_NATIVE_INT, numFaces, &faceMarkers[0]);
    this->M_HDF5->Write("Mesh", "Faces.GlobalId" + suffix, H5T_NATIVE_INT, numFaces, &faceGlobalId[0]);
    this->M_HDF5->Write("Mesh", "Faces.BoundaryFlag" + suffix, H5T_NATIVE_INT, numFaces,
                        &faceBoundaryFlags[0]);

    this->M_HDF5->Write("Mesh", "Faces.NeighbourId1" + suffix, H5T_NATIVE_INT, numFaces,
                        &faceNeighbourId[0][0]);
    this->M_HDF5->Write("Mesh", "Faces.NeighbourId2" + suffix, H5T_NATIVE_INT, numFaces,
                        &faceNeighbourId[1][0]);
    this->M_HDF5->Write("Mesh", "Faces.NeighbourPos1" + suffix, H5T_NATIVE_INT, numFaces,
                        &faceNeighbourPos[0][0]);
    this->M_HDF5->Write("Mesh", "Faces.NeighbourPos2" + suffix, H5T_NATIVE_INT, numFaces,
                        &faceNeighbourPos[1][0]);

    facePoints.clear();
    faceMarkers.clear();
    faceGlobalId.clear();
    faceBoundaryFlags.clear();
    faceNeighbourId.clear();
    faceNeighbourPos.clear();

    Int numVolumes = mesh->numVolumes();
    tmpVectorInt.resize(numVolumes);
    std::vector<std::vector<Int> > volumePoints(elementNodes, tmpVectorInt);

    std::vector<Int> volumeMarkers(numVolumes);
    std::vector<Int> volumeGlobalId(numVolumes);

    for (Int j = 0; j < numVolumes; ++j)
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
        this->M_HDF5->Write("Mesh", "Volumes.p" + idx.str() + suffix, H5T_NATIVE_INT, numVolumes,
                            &volumePoints[k][0]);
        idx.str(std::string());
        idx.clear();
    }
    this->M_HDF5->Write("Mesh", "Volumes.f" + suffix, H5T_NATIVE_INT, numVolumes, &volumeMarkers[0]);
    this->M_HDF5->Write("Mesh", "Volumes.GlobalId" + suffix, H5T_NATIVE_INT, numVolumes, &volumeGlobalId[0]);

    volumePoints.clear();
    volumeMarkers.clear();
    volumeGlobalId.clear();
}

template <typename MeshType>
void HDF5Filter3DMesh<MeshType>::writeSerialMesh()
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

template <typename MeshType>
void HDF5Filter3DMesh<MeshType>::writeParallelMesh()
{
    std::stringstream index;
    std::string suffix;

    index << M_comm->MyPID();
    suffix = "." + index.str();

    writePartition(M_parallelMesh, suffix);
}

template <typename MeshType>
void HDF5Filter3DMesh<MeshType>::writeInterfaces()
{
    Int interfaceNumber = M_DOFInterfaces.size();

    std::vector<Int> firstDOF;
    std::vector<Int> secondDOF;

    this->M_HDF5->Write("Interfaces", "Number", interfaceNumber);

    for (Int i = 0; i < interfaceNumber; ++i)
    {
        interfaceVector_Type& currentInterfaceSet = *(M_DOFInterfaces[i]);

        Int partitionNumber = currentInterfaceSet.size();
        std::stringstream idx;
        idx << i;
        this->M_HDF5->Write("Interfaces", "Partitions." + idx.str(), partitionNumber);

        Int flag1 = M_firstInterfaceFlags[i];
        Int flag2 = M_secondInterfaceFlags[i];
        std::string type = M_interfaceTypes[i];

        this->M_HDF5->Write("Interfaces", "Flag1." + idx.str(), flag1);
        this->M_HDF5->Write("Interfaces", "Flag2." + idx.str(), flag2);
        this->M_HDF5->Write("Interfaces", "Type." + idx.str(), type);

        for (Int j = 0; j < partitionNumber; ++j)
        {
            interface_Type& currentInterface = *(currentInterfaceSet[j]);

            const std::map<UInt, UInt>& locDofMap = currentInterface.localDofMap();

            Int size = locDofMap.size();

            idx.str(std::string());
            idx.clear();
            idx << i << "." << j;
            this->M_HDF5->Write("Interfaces", "Size." + idx.str(), size);

            firstDOF.clear();
            secondDOF.clear();
            firstDOF.resize(size);
            secondDOF.resize(size);

            Int k = 0;
            for (std::map<UInt, UInt>::const_iterator it = locDofMap.begin();
                 it != locDofMap.end(); ++it)
            {
                firstDOF[k] = it->first;
                secondDOF[k] = it->second;
                ++k;
            }
            this->M_HDF5->Write("Interfaces", "Key." + idx.str(), H5T_NATIVE_INT, size, &firstDOF[0]);
            this->M_HDF5->Write("Interfaces", "Value." + idx.str(), H5T_NATIVE_INT, size, &secondDOF[0]);
        }

    }
}

} // Namespace LifeV

#endif // HDF5FILTER3DMESH_H
