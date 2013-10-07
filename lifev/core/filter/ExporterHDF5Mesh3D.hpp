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
  @brief Class derived from ExporterHDF5 to provide I/O for the mesh partitions (RegionMesh only)

  @date 9-07-2010
  @author Radu Popescu <radu.popescu@epfl.ch>

  @maintainer Radu Popescu <radu.popescu@epfl.ch>

  This class has been deprecated and will be removed in future releases.
  Please use the PartitionIO class in its place.
*/

#ifndef EXPORTER_HDF5_MESH_3D_H
#define EXPORTER_HDF5_MESH_3D_H 1


#include <lifev/core/LifeV.hpp>

#ifdef HAVE_HDF5

#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/fem/DOFInterface3Dto3D.hpp>

namespace LifeV
{

//! Class derived from ExporterHDF5 to provide I/O for the mesh partitions (RegionMesh only)
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
  <li> When exporting graph and partition(s) create a ExporterHDF5 to be used for output.
  Call methods: addPartitionGraph, addMeshPartitionAll, or addMyMeshPartition
  Call postProcess to write data
  Call CloseFile
  <li> When importing graph and/or partition create a ExporterHDF5 to be used for input,
  Call loadGraph and partitionMesh::update, call loadMyPartition
  Call CloseFile
  </ol>
*/
template <typename MeshType>
class ExporterHDF5Mesh3D : public ExporterHDF5<MeshType>
{
public:
    typedef MeshType mesh_Type;
    typedef ExporterHDF5<MeshType> base;
    typedef typename base::meshPtr_Type meshPtr_Type;
    typedef typename base::vector_Type    vector_Type;
    typedef typename base::vectorPtr_Type vectorPtr_Type;
    typedef typename base::exporterData_Type exporterData_Type;

    typedef EpetraExt::HDF5 hdf5_Type;
    typedef boost::shared_ptr<hdf5_Type> hdf5Ptr_Type;
    typedef std::vector<std::vector<Int> > graph_Type;
    typedef boost::shared_ptr<graph_Type> graphPtr_Type;
    typedef boost::shared_ptr<std::vector<meshPtr_Type> > serialMeshPtr_Type;

    typedef DOFInterface3Dto3D interface_Type;
    typedef boost::shared_ptr<interface_Type> interfacePtr_Type;
    typedef std::vector<interfacePtr_Type> interfaceVector_Type;
    // The vector contains pointers to each fluid partition's interface with
    // the solid.
    typedef boost::shared_ptr<interfaceVector_Type> interfaceVectorPtr_Type;


    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor for ExporterHDF5Mesh3D
    ExporterHDF5Mesh3D() {}

    //! Constructor for ExporterHDF5Mesh3D
    /*!
      @param dataFile the GetPot data file where you must provide an [exporter] section with:
      "start"     (start index for sections in the hdf5 data structure 0 for 000, 1 for 001 etc.),
      "save"      (how many time steps per postprocessing)
      "multimesh" ( = true if the mesh has to be saved at each post-processing step)
      @param mesh the mesh
      @param the prefix for the case file (ex. "test" for test.case)
      @param the procId determines de CPU id. if negative, it ussemes there is only one processor
    */
    ExporterHDF5Mesh3D (const GetPot& dataFile, meshPtr_Type mesh, const std::string& prefix, const Int& procId);

    //! Constructor for ExporterHDF5Mesh3D without prefix and procID
    /*!
      @param dataFile the GetPot data file where you must provide an [exporter] section with:
      "start"     (start index for sections in the hdf5 data structure 0 for 000, 1 for 001 etc.),
      "save"      (how many time steps per postprocessing)
      "multimesh" ( = true if the mesh has to be saved at each post-processing step)
      @param mesh the mesh
    */
    ExporterHDF5Mesh3D (const GetPot& dataFile, const std::string& prefix);

    //! Destructor for ExporterHDF5Mesh3D
    virtual ~ExporterHDF5Mesh3D() {}

    //@}

    //! @name Pubic Methods
    //@{
    virtual void postProcess (const Real& time);

    //! Add the partition graph to the post processing data file
    /*!
      Add the partition graph to the post processing data file.
      \param graph - shared_ptr<vector<vector<Int> > > - shared pointer to the partition graph data structure
      (as returned by partitionMesh::graph() )
      \param comm - Epetra_Comm* - raw pointer to the Epetra communicator to be used
    */
    void addPartitionGraph (const graphPtr_Type& graph, boost::shared_ptr<Epetra_Comm>& comm)
    {
        M_graph = graph;
        M_comm = comm;
    }

    //! Add all of the mesh partitions to the post processing data file (serial operation)
    /*!
      Add all of the mesh partitions to the post processing data file.
      \param meshPointer - shared_ptr<vector<shared_ptr<MeshType> > > - shared pointer to the vector storing
      pointers to the mesh partitions (as returned by partitionMesh::meshAllPartitions() )
      \param comm - Epetra_Comm* - raw pointer to the Epetra communicator to be used
    */
    void addMeshPartitionAll (const serialMeshPtr_Type& meshPointer, boost::shared_ptr<Epetra_Comm>& comm)
    {
        M_serialMesh = meshPointer;
        M_parallelMesh.reset();
        M_comm = comm;
    }

    //! Add to HDF5 file the mesh partition that belongs to the current process (parallel operation)
    /*!
      After the mesh partition is loaded from the HDF5 file, the simulation is run and
      ExporterHDF5Mesh3D::postProcess() is called, the original contents of the HDF5 will be lost.
      To keep mesh partition, call:
      ExporterHDF5Mesh3D::addMyMeshPartition() before calling ExporterHDF5Mesh3D::postProcess();
      \param meshPointer - shared_ptr<Mesh> - shared pointer to a mesh partition (as returned by
      partitionMesh::mesh() )
      \param comm - Epetra_Comm* - raw pointer to the Epetra communicator to be used
    */
    void addMyMeshPartition ( const meshPtr_Type& /*meshPointer*/, boost::shared_ptr<Epetra_Comm>& /*comm*/ )
    {
        /*M_parallelMesh = meshPointer; M_serialMesh.reset(); M_comm = comm;*/
    }

    //! Add a DOF interface for writing to file
    /*!
      Add a DOF interface to the member vector M_interfaceData, for writing to the HDF5
      file. Call once for each interface that is to be written.
    */
    void addDOFInterface (const interfaceVectorPtr_Type& interfaces,
                          const std::string& type,
                          const Int& firstInterfaceFlag,
                          const Int& secondInterfaceFlag,
                          const boost::shared_ptr<Epetra_Comm>& comm);

    //! Get the number of stored DOF interfaces
    Int queryStoredInterfaceNumber();

    //! Get the types of the stored DOF interfaces
    std::vector<std::string> queryStoredInterfaceTypes();

    //! Return a pointer to the graph
    graphPtr_Type getGraph();

    //! Return a pointer to the mesh partition that corresponds to the current MPI rank
    meshPtr_Type  getMeshPartition();

    //! Return a pointer to the k-th interface stored inside the file
    boost::shared_ptr< std::map<UInt, UInt> > getStoredInterface (Int k) ;

    // When reading back partitions and interfaces, the comm member must be set explicitly.
    // This is intended for use with getMeshPartition and getStoredInterface
    //! Set the M_comm data member.
    void setComm ( const boost::shared_ptr<Epetra_Comm>& comm )
    {
        M_comm = comm;
    }
    //@}

private:

    //! @name Private methods
    //@{
    // The following private methods are for writing the partitioned graph
    // and mesh to the output file
    void writeGraph();
    void writePartition (meshPtr_Type partition, std::string& suffix);
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
ExporterHDF5Mesh3D<MeshType>::ExporterHDF5Mesh3D (const GetPot& dataFile, meshPtr_Type mesh,
                                                  const std::string& prefix,
                                                  const Int& procId) :
    base                ( dataFile, mesh, prefix, procId )
{
}

template<typename MeshType>
ExporterHDF5Mesh3D<MeshType>::ExporterHDF5Mesh3D (const GetPot& dataFile, const std::string& prefix) :
    base                ( dataFile, prefix )
{
}

// ===================================================
// Pubic Methods
// ===================================================

template<typename MeshType>
void ExporterHDF5Mesh3D<MeshType>::addDOFInterface (const interfaceVectorPtr_Type& interfaces,
                                                    const std::string& type,
                                                    const Int& firstInterfaceFlag,
                                                    const Int& secondInterfaceFlag,
                                                    const boost::shared_ptr<Epetra_Comm>& comm)
{
    M_DOFInterfaces.push_back (interfaces);
    M_interfaceTypes.push_back (type);
    M_firstInterfaceFlags.push_back (firstInterfaceFlag);
    M_secondInterfaceFlags.push_back (secondInterfaceFlag);

    M_comm = comm;
}

template<typename MeshType>
void ExporterHDF5Mesh3D<MeshType>::postProcess (const Real& time)
{
    if ( this->M_HDF5.get() == 0)
    {
        if (this->M_dataVector.size() != 0)
        {
            this->M_HDF5.reset (new hdf5_Type (this->M_dataVector.begin()->storedArrayPtr()->comm() ) );
        }
        else
        {
            this->M_HDF5.reset (new hdf5_Type (*M_comm) );
        }
        this->M_outputFileName = this->M_prefix + ".h5";
        this->M_HDF5->Create (this->M_postDir + this->M_outputFileName);

        // write empty xdmf file
        this->writeInitXdmf();

        if (!this->M_multimesh)
        {
            if (this->M_dataVector.size() != 0)
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

    // typedef typename std::list< exporterData_Type >::const_iterator Iterator;

    this->computePostfix();

    if ( this->M_postfix != "*****" )
    {
        if (!this->M_procId)
        {
            std::cout << "  X-  HDF5 post-processing ...        " << std::flush;
        }
        LifeChrono chrono;
        chrono.start();
        for (typename base::dataVectorIterator_Type i = this->M_dataVector.begin(); i != this->M_dataVector.end(); ++i)
        {
            this->writeVariable (*i);
        }
        // pushing time
        this->M_timeSteps.push_back (time);

        this->writeXdmf (time);

        if (this->M_multimesh)
        {
            this->writeGeometry(); // see also writeGeometry
        }

        chrono.stop();

        // Write to file without closing the file
        this->M_HDF5->Flush();

        if (!this->M_procId)
        {
            std::cout << "         done in " << chrono.diff() << " s." << std::endl;
        }
    }
}

template <typename MeshType>
int ExporterHDF5Mesh3D<MeshType>::queryStoredInterfaceNumber()
{
    if (this->M_HDF5.get() == 0)
    {
        this->M_HDF5.reset (new hdf5_Type (*M_comm) );
    }
    if (! this->M_HDF5->IsOpen() )
    {
        this->M_HDF5->Open (this->M_postDir + this->M_prefix + ".h5", H5F_ACC_RDONLY);
    }

    Int storedInterfaceNumber;
    this->M_HDF5->Read ("Interfaces", "Number", storedInterfaceNumber);

    return storedInterfaceNumber;
}

template <typename MeshType>
std::vector<std::string> ExporterHDF5Mesh3D<MeshType>::queryStoredInterfaceTypes()
{
    if (this->M_HDF5.get() == 0)
    {
        this->M_HDF5.reset (new hdf5_Type (*M_comm) );
    }
    if (! this->M_HDF5->IsOpen() )
    {
        this->M_HDF5->Open (this->M_postDir + this->M_prefix + ".h5", H5F_ACC_RDONLY);
    }

    int storedInterfaceNumber;
    this->M_HDF5->Read ("Interfaces", "Number", storedInterfaceNumber);

    std::vector<std::string> storedInterfaceTypes (storedInterfaceNumber);
    for (Int k = 0; k < storedInterfaceNumber; ++k)
    {
        std::stringstream idx;
        idx << k;
        this->M_HDF5->Read ("Interfaces", "Type." + idx.str(), storedInterfaceTypes[k]);
    }

    return storedInterfaceTypes;
}

template <typename MeshType>
typename ExporterHDF5Mesh3D<MeshType>::graphPtr_Type ExporterHDF5Mesh3D<MeshType>::getGraph()
{
    graphPtr_Type tempGraph (new graphPtr_Type::element_type);

    if (this->M_HDF5.get() == 0)
    {
        this->M_HDF5.reset (new hdf5_Type (*M_comm) );
    }
    if (! this->M_HDF5->IsOpen() )
    {
        this->M_HDF5->Open (this->M_postDir + this->M_prefix + ".h5", H5F_ACC_RDONLY);
    }

    Int nPartitions;

    this->M_HDF5->Read ("Graph", "number_partitions", nPartitions);

    std::vector<Int> partitionSizes (nPartitions);
    this->M_HDF5->Read ("Graph", "partition_sizes", H5T_NATIVE_INT, nPartitions,
                        &partitionSizes[0]);

    tempGraph->resize (0);
    tempGraph->reserve (nPartitions);

    std::vector<Int> partBuffer;
    std::stringstream index;

    for (Int i = 0; i < nPartitions; ++i)
    {
        partBuffer.resize (partitionSizes[i]);
        index << i;
        this->M_HDF5->Read ("Graph_partitions", "P" + index.str(),
                            H5T_NATIVE_INT, partitionSizes[i],
                            &partBuffer[0]);
        tempGraph->push_back (partBuffer);
        index.str (std::string() );
        index.clear();
    }

    return tempGraph;
}

template <typename MeshType>
typename ExporterHDF5Mesh3D<MeshType>::meshPtr_Type ExporterHDF5Mesh3D<MeshType>::getMeshPartition()
{
    meshPtr_Type tempMesh ( new MeshType ( M_comm ) );

    UInt elementNodes, faceNodes;
    switch (MeshType::elementShape_Type::S_shape)
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
        this->M_HDF5.reset (new hdf5_Type (*M_comm) );
    }
    if (! this->M_HDF5->IsOpen() )
    {
        this->M_HDF5->Open (this->M_postDir + this->M_prefix + ".h5", H5F_ACC_RDONLY);
    }

    std::stringstream index;
    index << M_comm->MyPID();
    std::string suffix = "." + index.str();

    // Read counters from file and set them in the mesh

    std::vector<Int> counters (13);
    this->M_HDF5->Read ("Mesh_counters", "P" + suffix, H5T_NATIVE_INT, 13, &counters[0]);
    // Points
    Int numPoints = counters[0];
    tempMesh->setMaxNumPoints (counters[0], true);
    tempMesh->setNumBPoints (counters[1]);

    // Vertices
    tempMesh->setNumVertices (counters[2]);
    tempMesh->setNumBVertices (counters[3]);
    tempMesh->setNumGlobalVertices (counters[4]);

    // Edges
    Int numEdges = counters[5];
    tempMesh->setNumEdges (counters[5]);
    tempMesh->setMaxNumEdges (counters[5]);
    tempMesh->setNumBEdges (counters[6]);
    tempMesh->setMaxNumGlobalEdges (counters[7]);

    // Faces
    Int numFaces = counters[8];
    tempMesh->setNumFaces (counters[8]);
    tempMesh->setMaxNumFaces (counters[8]);
    tempMesh->setNumBFaces (counters[9]);
    tempMesh->setMaxNumGlobalFaces (counters[10]);

    // Volumes
    Int numVolumes = counters[11];
    tempMesh->setMaxNumVolumes (counters[11], true);
    tempMesh->setMaxNumGlobalVolumes (counters[12]);

    // Read the list of points
    std::vector<Real> tmpVectorDouble (numPoints);
    std::vector<std::vector<Real> > pointCoordinates (3, tmpVectorDouble);

    std::vector<Int> pointMarkers (numPoints);
    std::vector<Int> pointFlags (numPoints);
    std::vector<Int> pointBoundaryFlags (numPoints);
    std::vector<Int> pointGlobalId (numPoints);

    this->M_HDF5->Read ("Mesh_points_coordinates_x", "P" + suffix, H5T_NATIVE_DOUBLE, numPoints, &pointCoordinates[0][0]);
    this->M_HDF5->Read ("Mesh_points_coordinates_y", "P" + suffix, H5T_NATIVE_DOUBLE, numPoints, &pointCoordinates[1][0]);
    this->M_HDF5->Read ("Mesh_points_coordinates_z", "P" + suffix, H5T_NATIVE_DOUBLE, numPoints, &pointCoordinates[2][0]);
    this->M_HDF5->Read ("Mesh_points_markers", "P" + suffix, H5T_NATIVE_INT, numPoints, &pointMarkers[0]);
    this->M_HDF5->Read ("Mesh_points_flags", "P" + suffix, H5T_NATIVE_INT, numPoints, &pointFlags[0]);
    this->M_HDF5->Read ("Mesh_points_boundary_flag", "P" + suffix, H5T_NATIVE_INT, numPoints, &pointBoundaryFlags[0]);
    this->M_HDF5->Read ("Mesh_points_gid", "P" + suffix, H5T_NATIVE_INT, numPoints, &pointGlobalId[0]);

    tempMesh->pointList.reserve (numPoints);
    tempMesh->_bPoints.reserve (tempMesh->numBPoints() );

    typename MeshType::point_Type* pp = 0;

    for (Int j = 0; j < numPoints; ++j)
    {
        pp = & (tempMesh->addPoint (false, false) );
        pp->setMarkerID (pointMarkers[j]);
        pp->setFlag (pointFlags[j]);
        pp->x() = pointCoordinates[0][j];
        pp->y() = pointCoordinates[1][j];
        pp->z() = pointCoordinates[2][j];
        pp->setId (pointGlobalId[j]);
    }

    pointCoordinates.clear();
    pointMarkers.clear();
    pointFlags.clear();
    pointBoundaryFlags.clear();
    pointGlobalId.clear();

    // Read the list of edges
    std::vector<Int> tmpVectorInt (numEdges);
    std::vector<std::vector<Int> > edgePoints (2, tmpVectorInt);

    std::vector<Int> edgeMarkers (numEdges);
    std::vector<Int> edgeFlags (numEdges);
    std::vector<Int> edgeGlobalId (numEdges);
    std::vector<Int> edgeBoundaryFlags (numEdges);

    this->M_HDF5->Read ("Mesh_edges_p1", "P" + suffix, H5T_NATIVE_INT, numEdges, &edgePoints[0][0]);
    this->M_HDF5->Read ("Mesh_edges_p2", "P" + suffix, H5T_NATIVE_INT, numEdges, &edgePoints[1][0]);
    this->M_HDF5->Read ("Mesh_edges_markers", "P" + suffix, H5T_NATIVE_INT, numEdges, &edgeMarkers[0]);
    this->M_HDF5->Read ("Mesh_edges_flags", "P" + suffix, H5T_NATIVE_INT, numEdges, &edgeFlags[0]);
    this->M_HDF5->Read ("Mesh_edges_gid", "P" + suffix, H5T_NATIVE_INT, numEdges, &edgeGlobalId[0]);
    this->M_HDF5->Read ("Mesh_edges_boundary_flag", "P" + suffix, H5T_NATIVE_INT, numEdges,
                        &edgeBoundaryFlags[0]);

    tempMesh->edgeList.reserve (numEdges);

    typename MeshType::edge_Type* pe;

    for (Int j = 0; j < numEdges; ++j)
    {
        pe = & (tempMesh->addEdge (edgeBoundaryFlags[j]) );
        pe->setId (edgeGlobalId[j]);
        pe->setPoint (0, tempMesh->point (edgePoints[0][j]) );
        pe->setPoint (1, tempMesh->point (edgePoints[1][j]) );
        pe->setMarkerID (edgeMarkers[j]);
        pe->setFlag (edgeFlags[j]);
    }

    edgePoints.clear();
    edgeMarkers.clear();
    edgeFlags.clear();
    edgeGlobalId.clear();
    edgeBoundaryFlags.clear();

    // Read the list of faces
    tmpVectorInt.resize (numFaces);
    std::vector<std::vector<Int> > facePoints (faceNodes, tmpVectorInt);

    std::vector<Int> faceMarkers (numFaces);
    std::vector<Int> faceFlags (numFaces);
    std::vector<Int> faceGlobalId (numFaces);
    std::vector<Int> faceBoundaryFlags (numFaces);

    std::vector<std::vector<Int> > faceNeighbourId (2, tmpVectorInt);
    std::vector<std::vector<Int> > faceNeighbourPos (2, tmpVectorInt);

    std::stringstream idx;
    for (UInt k = 0; k < faceNodes; ++k)
    {
        idx << k + 1;
        this->M_HDF5->Read ("Mesh_faces_points", "P" + idx.str() + suffix, H5T_NATIVE_INT, numFaces,
                            &facePoints[k][0]);
        idx.str (std::string() );
        idx.clear();
    }
    this->M_HDF5->Read ("Mesh_faces_markers", "P" + suffix, H5T_NATIVE_INT, numFaces, &faceMarkers[0]);
    this->M_HDF5->Read ("Mesh_faces_flags", "P" + suffix, H5T_NATIVE_INT, numFaces, &faceFlags[0]);
    this->M_HDF5->Read ("Mesh_faces_gid", "P" + suffix, H5T_NATIVE_INT, numFaces, &faceGlobalId[0]);
    this->M_HDF5->Read ("Mesh_faces_boundary_flags", "P" + suffix, H5T_NATIVE_INT, numFaces,
                        &faceBoundaryFlags[0]);

    this->M_HDF5->Read ("Mesh_faces_neighbour1", "P" + suffix, H5T_NATIVE_INT, numFaces,
                        &faceNeighbourId[0][0]);
    this->M_HDF5->Read ("Mesh_faces_neighbour2", "P" + suffix, H5T_NATIVE_INT, numFaces,
                        &faceNeighbourId[1][0]);
    this->M_HDF5->Read ("Mesh_faces_neighbour_pos1", "P" + suffix, H5T_NATIVE_INT, numFaces,
                        &faceNeighbourPos[0][0]);
    this->M_HDF5->Read ("Mesh_faces_neighbour_pos2", "P" + suffix, H5T_NATIVE_INT, numFaces,
                        &faceNeighbourPos[1][0]);


    typename MeshType::face_Type* pf = 0;

    tempMesh->faceList.reserve (numFaces);

    for (Int j = 0; j < numFaces; ++j)
    {
        pf = & (tempMesh->addFace (faceBoundaryFlags[j]) );
        pf->setId (faceGlobalId[j]);

        pf->firstAdjacentElementIdentity() = faceNeighbourId[0][j];
        pf->secondAdjacentElementIdentity() = faceNeighbourId[1][j];
        pf->firstAdjacentElementPosition() = faceNeighbourPos[0][j];
        pf->secondAdjacentElementPosition() = faceNeighbourPos[1][j];

        pf->setMarkerID (faceMarkers[j]);
        pf->setFlag (faceFlags[j]);
        for (UInt k = 0; k < faceNodes; ++k)
        {
            pf->setPoint (k, tempMesh->point (facePoints[k][j]) );
        }
    }

    tempMesh->setLinkSwitch ("HAS_ALL_FACETS");
    tempMesh->setLinkSwitch ("FACETS_HAVE_ADIACENCY");

    facePoints.clear();
    faceMarkers.clear();
    faceFlags.clear();
    faceGlobalId.clear();
    faceBoundaryFlags.clear();
    faceNeighbourId.clear();
    faceNeighbourPos.clear();

    // Read the list of volumes
    tmpVectorInt.resize (numVolumes);
    std::vector<std::vector<Int> > volumePoints (elementNodes, tmpVectorInt);

    std::vector<Int> volumeMarkers (numVolumes);
    std::vector<Int> volumeFlags (numVolumes);
    std::vector<Int> volumeGlobalId (numVolumes);

    for (UInt k = 0; k < elementNodes; ++k)
    {
        idx << k + 1;
        this->M_HDF5->Read ("Mesh_volumes_points", "P" + idx.str() + suffix, H5T_NATIVE_INT, numVolumes,
                            &volumePoints[k][0]);
        idx.str (std::string() );
        idx.clear();
    }
    this->M_HDF5->Read ("Mesh_volumes_markers", "P" + suffix, H5T_NATIVE_INT, numVolumes, &volumeMarkers[0]);
    this->M_HDF5->Read ("Mesh_volumes_flags", "P" + suffix, H5T_NATIVE_INT, numVolumes, &volumeFlags[0]);
    this->M_HDF5->Read ("Mesh_volumes_gid", "P" + suffix, H5T_NATIVE_INT, numVolumes, &volumeGlobalId[0]);

    tempMesh->volumeList.reserve (numVolumes);

    typename MeshType::volume_Type* pv = 0;

    for (Int j = 0; j < numVolumes; ++j)
    {
        pv = & (tempMesh->addVolume() );
        pv->setId (volumeGlobalId[j]);
        for (UInt k = 0; k < elementNodes; ++k)
        {
            pv->setPoint (k, tempMesh->point (volumePoints[k][j]) );
        }
        pv->setMarkerID (volumeMarkers[j]);
        pv->setFlag (volumeFlags[j]);
    }

    volumePoints.clear();
    volumeMarkers.clear();
    volumeFlags.clear();
    volumeGlobalId.clear();

    tempMesh->updateElementEdges (false, false);
    tempMesh->updateElementFaces (false, false);

    return tempMesh;
}

template <typename MeshType>
boost::shared_ptr< std::map<UInt, UInt> > ExporterHDF5Mesh3D<MeshType>::getStoredInterface (int k)
{
    if (this->M_HDF5.get() == 0)
    {
        this->M_HDF5.reset (new hdf5_Type (*M_comm) );
    }
    if (! this->M_HDF5->IsOpen() )
    {
        this->M_HDF5->Open (this->M_postDir + this->M_prefix + ".h5", H5F_ACC_RDONLY);
    }

    int myRank = M_comm->MyPID();

    std::stringstream idx;
    idx << k << "." << myRank;

    boost::shared_ptr<std::map<UInt, UInt> > interface (new std::map<UInt, UInt>);

    Int size;
    this->M_HDF5->Read ("Interfaces", "Size." + idx.str(), size);

    std::vector<UInt> keyVector (size);
    std::vector<UInt> valueVector (size);

    this->M_HDF5->Read ("Interfaces", "Key." + idx.str(), H5T_NATIVE_INT, size,
                        &keyVector[0]);
    this->M_HDF5->Read ("Interfaces", "Value." + idx.str(), H5T_NATIVE_INT, size,
                        &valueVector[0]);

    for (Int i = 0; i < size; ++i)
    {
        interface->insert (std::make_pair (keyVector[i], valueVector[i]) );
    }

    return interface;
}

// ========================
// Private methods
// ========================

template <typename MeshType>
void ExporterHDF5Mesh3D<MeshType>::writeGraph()
{
    std::vector<Int> partitionSizes;
    Int size, maxSize = 0;

    // Calculate the maximum size of the partitions and store partition
    // sizes
    partitionSizes.reserve (M_graph->size() );
    for (std::vector<std::vector<Int> >::iterator it = M_graph->begin();
            it != M_graph->end(); ++it)
    {
        size = it->size();
        if (size > maxSize)
        {
            maxSize = size;
        }
        partitionSizes.push_back (size);
    }

    this->M_HDF5->Write ("Graph", "partition_sizes", H5T_NATIVE_INT,
                         M_graph->size(), &partitionSizes[0]);

    // Write partition array size
    this->M_HDF5->Write ("Graph", "number_partitions", static_cast<Int> ( (M_graph->size() ) ) );

    // Write partition array
    std::stringstream index;
    for (UInt i = 0; i < M_graph->size(); ++i)
    {
        index << i;
        this->M_HDF5->Write ("Graph_partitions", "P" + index.str(),
                             H5T_NATIVE_INT, partitionSizes[i],
                             & (*M_graph) [i][0]);
        index.str (std::string() );
        index.clear();

        this->M_HDF5->Flush();

        std::cout << "Wrote graph for partition " << i << std::endl;
    }
}

template <typename MeshType>
void ExporterHDF5Mesh3D<MeshType>::writePartition (meshPtr_Type mesh, std::string& suffix)
{
    UInt elementNodes, faceNodes;
    switch (MeshType::elementShape_Type::S_shape)
    {
        case HEXA:
            elementNodes = 8;
            faceNodes    = 4;
            break;
        case TETRA:
            elementNodes = 4;
            faceNodes    = 3;
        default:
            ERROR_MSG("element type not supported");
    }

    std::vector<Int> counters;
    counters.push_back (static_cast<Int> (mesh->numPoints() ) );
    counters.push_back (static_cast<Int> (mesh->numBPoints() ) );
    counters.push_back (static_cast<Int> (mesh->numVertices() ) );
    counters.push_back (static_cast<Int> (mesh->numBVertices() ) );
    counters.push_back (static_cast<Int> (mesh->numGlobalVertices() ) );
    counters.push_back (static_cast<Int> (mesh->numEdges() ) );
    counters.push_back (static_cast<Int> (mesh->numBEdges() ) );
    counters.push_back (static_cast<Int> (mesh->numGlobalEdges() ) );
    counters.push_back (static_cast<Int> (mesh->numFaces() ) );
    counters.push_back (static_cast<Int> (mesh->numBFaces() ) );
    counters.push_back (static_cast<Int> (mesh->numGlobalFaces() ) );
    counters.push_back (static_cast<Int> (mesh->numVolumes() ) );
    counters.push_back (static_cast<Int> (mesh->numGlobalVolumes() ) );

    this->M_HDF5->Write ("Mesh_counters", "P" + suffix, H5T_NATIVE_INT, 13, &counters[0]);
    this->M_HDF5->Flush();

    UInt numPoints = mesh->numPoints();
    std::vector<Real> tmpVectorDouble (numPoints);
    std::vector<std::vector<Real> > pointCoordinates (3, tmpVectorDouble);

    std::vector<Int> pointMarkers (numPoints);
    std::vector<Int> pointFlags (numPoints);
    std::vector<Int> pointGlobalId (numPoints);
    std::vector<Int> pointBoundaryFlags (numPoints);

    for (UInt j = 0; j < numPoints; ++j)
    {
        pointCoordinates[0][j] = mesh->pointList[j].x();
        pointCoordinates[1][j] = mesh->pointList[j].y();
        pointCoordinates[2][j] = mesh->pointList[j].z();
        pointMarkers[j] = mesh->pointList[j].markerID();
        pointFlags[j] = mesh->pointList[j].flag();
        pointGlobalId[j] = mesh->point (j).id();
        if (mesh->isBoundaryPoint (j) )
        {
            pointBoundaryFlags[j] = 1;
        }
    }

    this->M_HDF5->Write ("Mesh_points_coordinates_x", "P" + suffix, H5T_NATIVE_DOUBLE, numPoints, &pointCoordinates[0][0]);
    this->M_HDF5->Write ("Mesh_points_coordinates_y", "P" + suffix, H5T_NATIVE_DOUBLE, numPoints, &pointCoordinates[1][0]);
    this->M_HDF5->Write ("Mesh_points_coordinates_z", "P" + suffix, H5T_NATIVE_DOUBLE, numPoints, &pointCoordinates[2][0]);
    this->M_HDF5->Write ("Mesh_points_markers", "P" + suffix, H5T_NATIVE_INT, numPoints, &pointMarkers[0]);
    this->M_HDF5->Write ("Mesh_points_flags", "P" + suffix, H5T_NATIVE_INT, numPoints, &pointFlags[0]);
    this->M_HDF5->Write ("Mesh_points_gid", "P" + suffix, H5T_NATIVE_INT, numPoints, &pointGlobalId[0]);
    this->M_HDF5->Write ("Mesh_points_boundary_flag", "P" + suffix, H5T_NATIVE_INT, numPoints,
                         &pointBoundaryFlags[0]);

    this->M_HDF5->Flush();

    pointCoordinates.clear();
    pointMarkers.clear();
    pointFlags.clear();
    pointBoundaryFlags.clear();
    pointGlobalId.clear();

    Int numEdges = mesh->numEdges();
    std::vector<Int> tmpVectorInt (numEdges);
    std::vector<std::vector<Int> > edgePoints (2, tmpVectorInt);

    std::vector<Int> edgeMarkers (numEdges);
    std::vector<Int> edgeFlags (numEdges);
    std::vector<Int> edgeGlobalId (numEdges);
    std::vector<Int> edgeBoundaryFlags (numEdges);

    for (Int j = 0; j < numEdges; ++j)
    {
        edgePoints[0][j] = mesh->edgeList[j].point (0).localId();
        edgePoints[1][j] = mesh->edgeList[j].point (1).localId();
        edgeMarkers[j] = mesh->edgeList[j].markerID();
        edgeFlags[j] = mesh->edgeList[j].flag();
        edgeGlobalId[j] = mesh->edgeList[j].id();

        if (mesh->isBoundaryEdge (j) )
        {
            edgeBoundaryFlags[j] = 1;
        }
    }

    this->M_HDF5->Write ("Mesh_edges_p1", "P" + suffix, H5T_NATIVE_INT, numEdges, &edgePoints[0][0]);
    this->M_HDF5->Write ("Mesh_edges_p2", "P" + suffix, H5T_NATIVE_INT, numEdges, &edgePoints[1][0]);
    this->M_HDF5->Write ("Mesh_edges_markers", "P" + suffix, H5T_NATIVE_INT, numEdges, &edgeMarkers[0]);
    this->M_HDF5->Write ("Mesh_edges_flags", "P" + suffix, H5T_NATIVE_INT, numEdges, &edgeFlags[0]);
    this->M_HDF5->Write ("Mesh_edges_gid", "P" + suffix, H5T_NATIVE_INT, numEdges, &edgeGlobalId[0]);
    this->M_HDF5->Write ("Mesh_edges_boundary_flag", "P" + suffix, H5T_NATIVE_INT, numEdges,
                         &edgeBoundaryFlags[0]);

    this->M_HDF5->Flush();

    edgePoints.clear();
    edgeMarkers.clear();
    edgeFlags.clear();
    edgeGlobalId.clear();
    edgeBoundaryFlags.clear();

    Int numFaces = mesh->numFaces();
    tmpVectorInt.resize (numFaces);
    std::vector<std::vector<Int> > facePoints (faceNodes, tmpVectorInt);

    std::vector<Int> faceMarkers (numFaces);
    std::vector<Int> faceFlags (numFaces);
    std::vector<Int> faceGlobalId (numFaces);
    std::vector<Int> faceBoundaryFlags (numFaces);

    std::vector<std::vector<Int> > faceNeighbourId (2, tmpVectorInt);
    std::vector<std::vector<Int> > faceNeighbourPos (2, tmpVectorInt);

    for (Int j = 0; j < numFaces; ++j)
    {
        for (UInt k = 0; k < faceNodes; ++k)
        {
            facePoints[k][j] = mesh->faceList[j].point (k).localId();
        }
        faceMarkers[j] = mesh->faceList[j].markerID();
        faceFlags[j] = mesh->faceList[j].flag();
        faceGlobalId[j] = mesh->faceList[j].id();

        faceNeighbourId[0][j] = mesh->faceList[j].firstAdjacentElementIdentity();
        faceNeighbourId[1][j] = mesh->faceList[j].secondAdjacentElementIdentity();
        faceNeighbourPos[0][j] = mesh->faceList[j].firstAdjacentElementPosition();
        faceNeighbourPos[1][j] = mesh->faceList[j].secondAdjacentElementPosition();

        if (mesh->isBoundaryFace (j) )
        {
            faceBoundaryFlags[j] = 1;
        }
    }

    std::stringstream idx;
    for (UInt k = 0; k < faceNodes; ++k)
    {
        idx << k + 1;
        this->M_HDF5->Write ("Mesh_faces_points", "P" + idx.str() + suffix, H5T_NATIVE_INT, numFaces,
                             &facePoints[k][0]);
        idx.str (std::string() );
        idx.clear();
    }
    this->M_HDF5->Write ("Mesh_faces_markers", "P" + suffix, H5T_NATIVE_INT, numFaces, &faceMarkers[0]);
    this->M_HDF5->Write ("Mesh_faces_flags", "P" + suffix, H5T_NATIVE_INT, numFaces, &faceFlags[0]);
    this->M_HDF5->Write ("Mesh_faces_gid", "P" + suffix, H5T_NATIVE_INT, numFaces, &faceGlobalId[0]);
    this->M_HDF5->Write ("Mesh_faces_boundary_flags", "P" + suffix, H5T_NATIVE_INT, numFaces,
                         &faceBoundaryFlags[0]);

    this->M_HDF5->Write ("Mesh_faces_neighbour1", "P" + suffix, H5T_NATIVE_INT, numFaces,
                         &faceNeighbourId[0][0]);
    this->M_HDF5->Write ("Mesh_faces_neighbour2", "P" + suffix, H5T_NATIVE_INT, numFaces,
                         &faceNeighbourId[1][0]);
    this->M_HDF5->Write ("Mesh_faces_neighbour_pos1", "P" + suffix, H5T_NATIVE_INT, numFaces,
                         &faceNeighbourPos[0][0]);
    this->M_HDF5->Write ("Mesh_faces_neighbour_pos2", "P" + suffix, H5T_NATIVE_INT, numFaces,
                         &faceNeighbourPos[1][0]);

    this->M_HDF5->Flush();

    facePoints.clear();
    faceMarkers.clear();
    faceFlags.clear();
    faceGlobalId.clear();
    faceBoundaryFlags.clear();
    faceNeighbourId.clear();
    faceNeighbourPos.clear();

    Int numVolumes = mesh->numVolumes();
    tmpVectorInt.resize (numVolumes);
    std::vector<std::vector<Int> > volumePoints (elementNodes, tmpVectorInt);

    std::vector<Int> volumeMarkers (numVolumes);
    std::vector<Int> volumeFlags (numVolumes);
    std::vector<Int> volumeGlobalId (numVolumes);

    for (Int j = 0; j < numVolumes; ++j)
    {
        for (UInt k = 0; k < elementNodes; ++k)
        {
            volumePoints[k][j] = mesh->volumeList[j].point (k).localId();
        }
        volumeMarkers[j] = mesh->volumeList[j].markerID();
        volumeFlags[j] = mesh->volumeList[j].flag();
        volumeGlobalId[j] = mesh->volumeList[j].id();
    }

    for (UInt k = 0; k < elementNodes; ++k)
    {
        idx << k + 1;
        this->M_HDF5->Write ("Mesh_volumes_points", "P" + idx.str() + suffix, H5T_NATIVE_INT, numVolumes,
                             &volumePoints[k][0]);
        idx.str (std::string() );
        idx.clear();
    }
    this->M_HDF5->Write ("Mesh_volumes_markers", "P" + suffix, H5T_NATIVE_INT, numVolumes, &volumeMarkers[0]);
    this->M_HDF5->Write ("Mesh_volumes_flags", "P" + suffix, H5T_NATIVE_INT, numVolumes, &volumeFlags[0]);
    this->M_HDF5->Write ("Mesh_volumes_gid", "P" + suffix, H5T_NATIVE_INT, numVolumes, &volumeGlobalId[0]);

    this->M_HDF5->Flush();

    volumePoints.clear();
    volumeMarkers.clear();
    volumeFlags.clear();
    volumeGlobalId.clear();
}

template <typename MeshType>
void ExporterHDF5Mesh3D<MeshType>::writeSerialMesh()
{
    std::stringstream index;
    std::string suffix;

    for (UInt i = 0; i < M_serialMesh->size(); ++i)
    {
        index << i;
        suffix = "." + index.str();

        writePartition ( (*M_serialMesh) [i], suffix);

        index.str (std::string() );
        index.clear();

        std::cout << "Wrote mesh partition " << i << std::endl;
    }
}

template <typename MeshType>
void ExporterHDF5Mesh3D<MeshType>::writeParallelMesh()
{
    std::stringstream index;
    std::string suffix;

    index << M_comm->MyPID();
    suffix = "." + index.str();

    writePartition (M_parallelMesh, suffix);
}

template <typename MeshType>
void ExporterHDF5Mesh3D<MeshType>::writeInterfaces()
{
    Int interfaceNumber = M_DOFInterfaces.size();

    std::vector<Int> firstDOF;
    std::vector<Int> secondDOF;

    this->M_HDF5->Write ("Interfaces", "Number", interfaceNumber);

    for (Int i = 0; i < interfaceNumber; ++i)
    {
        interfaceVector_Type& currentInterfaceSet = * (M_DOFInterfaces[i]);

        Int partitionNumber = currentInterfaceSet.size();
        std::stringstream idx;
        idx << i;
        this->M_HDF5->Write ("Interfaces", "Partitions." + idx.str(), partitionNumber);

        Int flag1 = M_firstInterfaceFlags[i];
        Int flag2 = M_secondInterfaceFlags[i];
        std::string type = M_interfaceTypes[i];

        this->M_HDF5->Write ("Interfaces", "Flag1." + idx.str(), flag1);
        this->M_HDF5->Write ("Interfaces", "Flag2." + idx.str(), flag2);
        this->M_HDF5->Write ("Interfaces", "Type." + idx.str(), type);

        for (Int j = 0; j < partitionNumber; ++j)
        {
            interface_Type& currentInterface = * (currentInterfaceSet[j]);

            const std::map<UInt, UInt>& locDofMap = currentInterface.localDofMap();

            Int size = locDofMap.size();

            idx.str (std::string() );
            idx.clear();
            idx << i << "." << j;
            this->M_HDF5->Write ("Interfaces", "Size." + idx.str(), size);

            firstDOF.clear();
            secondDOF.clear();
            firstDOF.resize (size);
            secondDOF.resize (size);

            Int k = 0;
            for (std::map<UInt, UInt>::const_iterator it = locDofMap.begin();
                    it != locDofMap.end(); ++it)
            {
                firstDOF[k] = it->first;
                secondDOF[k] = it->second;
                ++k;
            }
            this->M_HDF5->Write ("Interfaces", "Key." + idx.str(), H5T_NATIVE_INT, size, &firstDOF[0]);
            this->M_HDF5->Write ("Interfaces", "Value." + idx.str(), H5T_NATIVE_INT, size, &secondDOF[0]);
        }

    }
}

} // Namespace LifeV

#endif // HAVE_HDF5

#endif // EXPORTER_HDF5_MESH_3D_H
