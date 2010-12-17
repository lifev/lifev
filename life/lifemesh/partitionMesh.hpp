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
  @brief Class that handles mesh partitioning

  @date 08-12-2010
  @author Gilles Fourestey <gilles.fourestey@epfl.ch>

  @contributor Radu Popescu <radu.popescu@epfl.ch>
  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#ifndef PARTMESH_H
#define PARTMESH_H 1

#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/shared_ptr.hpp>
#include <parmetis.h>
#include <Epetra_MpiComm.h>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifecore/life.hpp>
#include <life/lifefem/dof.hpp>
#include <life/lifemesh/regionMesh3D.hpp>

namespace LifeV
{
/*!
  @brief Class that handles mesh partitioning
  @author Gilles Fourestey gilles.fourestey@epfl.ch
  @contributor Radu Popescu radu.popescu@epfl.ch

  This class implements the partitioning of a mesh using (par)Metis. Currently
  it is possible to do both online and offline partitioning. In online
  partitioning, each MPI rank stores only its own mesh partition. In offline
  partitioning (run as a single process, on a workstation), all the mesh
  partitions are stored and can be saved to disk, using the HDF5 filter,
  for later use during a parallel run.
*/

template<typename MeshType>
class partitionMesh
{
public:
    //@{
    // Make the template's type available to the outside
    typedef MeshType mesh_Type;
    typedef boost::shared_ptr<MeshType> meshPtr_Type;
    typedef std::vector<std::vector<Int> > graph_Type;
    typedef boost::shared_ptr<graph_Type> graphPtr_Type;
    typedef std::vector<meshPtr_Type> partMesh_Type;
    typedef boost::shared_ptr<partMesh_Type> partMeshPtr_Type;
    //@}
    //! \name Constructors & Destructors
    //@{
    //! Default empty constructor
    partitionMesh();
    //! Constructor
    /*!
      This is a non-empty constructor. It takes as parameters the
      unpartitioned mesh (reference), the Epetra_Comm object in
      use (reference) and pointers to the Epetra interface and
      repeated interface maps. The constructor initializes the
      data members and calls a private method
      partitionMesh::execute which handles the mesh partitioning.
      \param mesh - Mesh& - the unpartitioned mesh
      \param _comm - Epetra_Comm& - Epetra communicator object
      \param interfaceMap - Epetra_Map*
      \param interfaceMapRep - Epetra_Map*
    */
    partitionMesh(meshPtr_Type &mesh, boost::shared_ptr<Epetra_Comm> comm, Epetra_Map* interfaceMap = 0,
                  Epetra_Map* interfaceMapRep = 0);
    //! Empty destructor
    virtual ~partitionMesh() {}
    //@}

    //! \name Public Methods
    //@{
    //! To be used with the new constructor.
    /*!
      Loads the parameters of the partitioning process from the simulation data file.
      Allocates and initializes data members according to the number of partitions specified in
      the data file.
      \param numPartitions - UInt - the number of partitions desired, in the offline case
      \param _comm - Epetra_Comm& - reference of the Epetra communicator used
    */
    void setup(UInt numPartitions, boost::shared_ptr<Epetra_Comm> comm);
    //! Call update() method after loading the graph, to rebuild all data structures
    /*!
      This method is to be called after the partitioned graph is LOADED from a HDF5 file.
      Set M_nPartitions and rebuilds the M_graphVertexLocations data member, required to do
      mesh partitioning.
      !!! This method should be the first one to be called after loading the graph. !!!
    */
    void update();
    //!  Stores a pointer to the unpartitioned mesh in the M_originalMesh data member.
    /*!
      Stores a pointer to the unpartitioned mesh in the M_originalMesh data member.
      \param mesh - Mesh& - the unpartitioned mesh (passed by reference)
      \param interfaceMap - Epetra_Map* - pointer to the interface map (default value 0)
      \param interfaceMapRep - Epetra_Map* - pointer to the repeated interface map (default value 0)
    */
    void attachUnpartitionedMesh(meshPtr_Type &mesh, Epetra_Map* interfaceMap = 0,
                                 Epetra_Map* interfaceMapRep = 0);
    //! Releases the original unpartitioned mesh
    /*!
      Releases the unpartitioned mesh so that it can be deleted, freein A LOT of memory
      in some cases.
    */
    void releaseUnpartitionedMesh();
    //! Executes the ParMETIS graph partitioning
    /*!
      Executes the ParMETIS graph partitioning
    */
    void doPartitionGraph();
    //! Builds the partitioned mesh using the partitioned graph
    /*!
      Builds the partitioned mesh using the partitioned graph
    */
    void doPartitionMesh();

    // Next method should be renamed and become a regular method
    //! Return a pointer to the mesh partition with rank k
    const meshPtr_Type&      getPartition(Int k)    const {return (*M_meshPartitions)[k];}

    //! Prints information about the state (data) of the object
    void showMe(std::ostream& output = std::cout) const;
    //@}

    //! \name Get Methods
    //@{
    //! Return a reference to M_vertexDistribution
    const std::vector<Int>&  vertexDistribution()   const {return M_vertexDistribution;};
    //! Return a const pointer to M_meshPartitions[0] - for parallel
    const meshPtr_Type&      meshPartition()        const {return (*M_meshPartitions)[0];}
    meshPtr_Type&            meshPartition()              {return (*M_meshPartitions)[0];}
    //! Return a pointer to M_meshPartitions
    const partMeshPtr_Type&  meshPartitions()       const {return M_meshPartitions;}
    //! Return a pointer to M_graphVertexLocations
    const std::vector<Int>&  graphVertexLocations() const {return M_graphVertexLocations;}
    //! Return a pointer to M_elementDomains
    const graphPtr_Type&     elementDomains()       const {return M_elementDomains;}
    //! Return a reference to M_repeatedNodeVector
    const std::vector<Int>&  repeatedNodeVector()   const {return M_repeatedNodeVector[0];}
    //! Return a reference to M_repeatedEdgeVector
    const std::vector<Int>&  repeatedEdgeVector()   const {return M_repeatedEdgeVector[0];}
    //! Return a reference to M_repeatedFaceVector
    const std::vector<Int>&  repeatedFaceVector()   const {return M_repeatedFaceVector[0];}
    //! Return a reference to M_repeatedVolumeVector
    const std::vector<Int>&  repeatedVolumeVector() const {return M_repeatedVolumeVector[0];}

    // DEPRECATED METHODS
    const std::vector<Int>& __attribute__((__deprecated__)) vertexDist() const {return M_vertexDistribution;};
    const meshPtr_Type& __attribute__((__deprecated__)) mesh() const {return (*M_meshPartitions)[0];}
    const partMeshPtr_Type& __attribute__((__deprecated__)) meshAllPartitions() const {return M_meshPartitions;}
    const std::vector<Int>& __attribute__((__deprecated__)) part() const {return M_graphVertexLocations;}
    const graphPtr_Type& __attribute__((__deprecated__)) graph() const {return M_elementDomains;}
    const meshPtr_Type& __attribute__((__deprecated__)) mesh(Int k) const {return (*M_meshPartitions)[k];}
    //@}

private:
    // Private copy constructor and assignment operator. No implementation
    partitionMesh(const partitionMesh&);
    partitionMesh& operator=(const partitionMesh&);
    //! Private Methods
    //@{
    //! Execute mesh partitioning using the configured MPI processes (online partitioning)
    /*!
      Executed the mesh partitioning using the number of MPI processes as the number of partitions.
      Sets current mesh element parameters: M_elementNodes, M_elementEdges, M_elementFaces, M_faceNodes
      Updates: M_elementDomains (indirectly)
      Other data members are changed indirectly by calling other private methods.
    */
    void execute();
    //! Sets the element parameters according to the type of mesh element used.
    /*!
      Sets element parameters (nodes, faces, edges and number of nodes on each
      face according to the type of mesh element used (Mesh::ElementShape::Shape).
      Updates M_elementNodes, M_elementFaces, M_elementEdges, M_faceNodes.
    */
    void setElementParameters();
    //! Build the graph vertex distribution vector
    /*!
      Updates the member M_vertexDistribution according to the number of processors to be
      used by ParMETIS (the number of processes started for MPI
      \param numElements - UInt - number of elements in the mesh
    */
    void distributeElements(UInt numElements);
    //! Find faces on the boundaries between domains (FSI)
    /*!
      Identifies the element faces that are common to both the fluid and the solid
      meshes and creates a map between the faces that reside on the boundary between
      the two meshes and the processors used. Updates the members M_repeatedFace and
      M_isOnProc.
    */
    void findRepeatedFacesFSI();
    //! Partition the connectivity graph using ParMETIS
    /*!
      Partitions the connectivity graph using ParMETIS. The result is stored in the
      member M_graphVertexLocations: this is a vector of integer values; its size is the number of
      elements in the unpartitioned mesh. Each value represents the number of the
      partition to which the element was assigned. Also creates M_elementDomains, the
      vector of elements in each subdomain.
      Updates: M_graphVertexLocations, M_elementDomains
      \param numParts - unsigned int - number of partitions for the graph cutting process
    */
    void partitionConnectivityGraph(UInt numParts);
    //! Updates the map between elements and processors in FSI
    /*!
      Updates M_elementDomains during FSI modeling.
    */
    void matchFluidPartitionsFSI();
    //! Redistribute elements among processes
    /*!
      Redistributes elements among processes, when needed, after the connectivity
      graph partitioning phase. Updates M_elementDomains
    */
    void redistributeElements();
    //! Construct local mesh
    /*!
      Constructs the data structures for the local mesh partition.
      Updates M_localNodes, M_localEdges, M_localFaces, M_localVolumes,
      M_globalToLocalNode.
    */
    void constructLocalMesh();
    //! Construct nodes
    /*!
      Adds nodes to the partitioned mesh object. Updates M_nBoundaryPoints,
      M_meshPartitions.
    */
    void constructNodes();
    //! Construct volumes
    /*!
      Adds volumes to the partitioned mesh object. Updates M_globalToLocalVolume,
      M_meshPartitions.
    */
    void constructVolumes();
    //! Construct edges
    /*!
      Adds edges to the partitioned mesh object. Updates M_nBoundaryEdges,
      M_meshPartitions.
    */
    void constructEdges();
    //! Construct faces
    /*!
      Adds faces to the partitioned mesh object. Updates M_nBoundaryFaces,
      M_meshPartitions.
    */
    void constructFaces();
    //! Final setup of local mesh
    /*!
      Updates the partitioned mesh object data members after adding the mesh
      elements (nodes, edges, faces, volumes).
      Updates M_meshPartitions.
    */
    void finalSetup();
    //! Create repeated element map
    /*!
      Creates a map of the boundary elements (nodes, edges, faces, volumes).
      Updates M_repeatedNodeVector, M_repeatedEdgeVector, M_repeatedFaceVector,
      M_repeatedVolumeVector.
    */
    void createRepeatedMap();

    //@}
    //! Private Data Members
    //@{
    UInt                                 M_numPartitions;
    partMeshPtr_Type                     M_meshPartitions;
    std::vector<Int>                     M_vertexDistribution;
    std::vector<Int>                     M_adjacencyGraphKeys;
    std::vector<Int>                     M_adjacencyGraphValues;
    boost::shared_ptr<Epetra_Comm>       M_comm;
    UInt                                 M_me;

    std::vector<std::vector<Int> >       M_localNodes;
    std::vector<std::set<Int> >          M_localEdges;
    std::vector<std::set<Int> >          M_localFaces;
    std::vector<std::vector<Int> >       M_localVolumes;
    std::vector<std::vector<Int> >       M_repeatedNodeVector;
    std::vector<std::vector<Int> >       M_repeatedEdgeVector;
    std::vector<std::vector<Int> >       M_repeatedFaceVector;
    std::vector<std::vector<Int> >       M_repeatedVolumeVector;
    std::vector<std::map<Int, Int> >     M_globalToLocalNode;
    std::vector<std::map<Int, Int> >     M_globalToLocalVolume;
    std::vector<UInt>                    M_nBoundaryPoints;
    std::vector<UInt>                    M_nBoundaryEdges;
    std::vector<UInt>                    M_nBoundaryFaces;
    // The following are utility variables used throughout the partitioning
    // process
    meshPtr_Type                         M_originalMesh;
    Epetra_Map*                          M_interfaceMap;
    Epetra_Map*                          M_interfaceMapRep;
    //! Number of partitions handled. 1 for parallel (old way), != 1 for serial
    UInt                                 M_elementNodes;
    UInt                                 M_elementFaces;
    UInt                                 M_elementEdges;
    UInt                                 M_faceNodes;
    boost::shared_ptr<std::vector<Int> > M_repeatedFace;
    boost::shared_ptr<std::vector<Int> > M_isOnProc;
    std::vector<Int>                     M_graphVertexLocations;
    graphPtr_Type                        M_elementDomains;
    bool                                 M_serialMode; // how to tell if running serial partition mode
    //@}
}; // class partitionMesh

//
// IMPLEMENTATION
//


// =================================
// Constructors and destructor
// =================================

template<typename MeshType>
partitionMesh<MeshType>::partitionMesh()
{
}

template<typename MeshType>
partitionMesh<MeshType>::partitionMesh(meshPtr_Type &mesh, boost::shared_ptr<Epetra_Comm> comm,
                                   Epetra_Map* interfaceMap,
                                   Epetra_Map* interfaceMapRep):
    M_numPartitions (1),
    M_comm (comm),
    M_originalMesh (mesh),
    M_interfaceMap (interfaceMap),
    M_interfaceMapRep (interfaceMapRep),
    M_elementDomains (new graph_Type),
    M_serialMode (false)
{
    M_me = M_comm->MyPID();

    meshPtr_Type newMesh (new MeshType);
    M_meshPartitions.reset(new partMesh_Type(M_numPartitions, newMesh));
    newMesh.reset();

    M_localNodes.resize(1);
    M_localEdges.resize(1);
    M_localFaces.resize(1);
    M_localVolumes.resize(1);
    M_repeatedNodeVector.resize(1);
    M_repeatedEdgeVector.resize(1);
    M_repeatedFaceVector.resize(1);
    M_repeatedVolumeVector.resize(1);
    M_globalToLocalNode.resize(1);
    M_globalToLocalVolume.resize(1);
    M_nBoundaryPoints.resize(1);
    M_nBoundaryEdges.resize(1);
    M_nBoundaryFaces.resize(1);

    execute();
}

// =================================
// Public methods
// =================================

template<typename MeshType>
void partitionMesh<MeshType>::setup(UInt numPartitions, boost::shared_ptr<Epetra_Comm> comm)
{
    M_serialMode = true;
    M_comm = comm;
    M_me = M_comm->MyPID();
    setElementParameters();

    M_numPartitions = numPartitions;

    M_meshPartitions.reset(new partMesh_Type);
    meshPtr_Type newMesh;
    for (UInt i = 0; i < M_numPartitions; ++i)
    {
        newMesh.reset(new MeshType);
        M_meshPartitions->push_back(newMesh);
    }
    newMesh.reset();

    M_elementDomains.reset(new graph_Type);

    M_localNodes.resize(M_numPartitions);
    M_localEdges.resize(M_numPartitions);
    M_localFaces.resize(M_numPartitions);
    M_localVolumes.resize(M_numPartitions);
    M_repeatedNodeVector.resize(M_numPartitions);
    M_repeatedEdgeVector.resize(M_numPartitions);
    M_repeatedFaceVector.resize(M_numPartitions);
    M_repeatedVolumeVector.resize(M_numPartitions);
    M_globalToLocalNode.resize(M_numPartitions);
    M_globalToLocalVolume.resize(M_numPartitions);
    M_nBoundaryPoints.resize(M_numPartitions);
    M_nBoundaryEdges.resize(M_numPartitions);
    M_nBoundaryFaces.resize(M_numPartitions);
}

template<typename MeshType>
void partitionMesh<MeshType>::update()
{
    M_numPartitions = M_elementDomains->size();

    Int numElements = 0;

    for (UInt i = 0; i < M_numPartitions; ++i)
    {
        numElements += (*M_elementDomains)[i].size();
    }

    // Rebuild M_graphVertexLocations
    M_graphVertexLocations.resize(numElements);
    for (std::vector<std::vector<Int> >::iterator it1 = M_elementDomains->begin();
         it1 != M_elementDomains->end(); ++it1)
    {
        for (std::vector<Int>::iterator it2 = it1->begin();
             it2 != it1->end(); ++it2)
        {
            M_graphVertexLocations[*it2] = static_cast<Int>((it1 - M_elementDomains->begin()));
        }
    }
}

template<typename MeshType>
void partitionMesh<MeshType>::attachUnpartitionedMesh(meshPtr_Type &mesh,
                                                  Epetra_Map* interfaceMap,
                                                  Epetra_Map* interfaceMapRep)
{
    M_originalMesh = mesh;
    M_interfaceMap = interfaceMap;
    M_interfaceMapRep = interfaceMapRep;
}

template<typename MeshType>
void partitionMesh<MeshType>::releaseUnpartitionedMesh()
{
    M_originalMesh.reset();
    M_interfaceMap = 0;
    M_interfaceMapRep = 0;
}

template<typename MeshType>
void partitionMesh<MeshType>::doPartitionGraph()
{
    distributeElements(M_originalMesh->numElements());
    if (M_interfaceMap)
    {
        findRepeatedFacesFSI();
    }
    partitionConnectivityGraph(M_numPartitions);
    if (M_interfaceMap)
    {
        matchFluidPartitionsFSI();
    }
    redistributeElements();
}

template<typename MeshType>
void partitionMesh<MeshType>::doPartitionMesh()
{
    constructLocalMesh();
    constructNodes();
    constructVolumes();
    constructEdges();
    constructFaces();
    finalSetup();
    createRepeatedMap();
}

template<typename MeshType>
void partitionMesh<MeshType>::showMe(std::ostream& output) const
{
    output << "Number of partitions: " << M_numPartitions << std::endl;
    output << "Serial mode:" << M_serialMode << std::endl;
}

// =================================
// Private methods
// =================================

template<typename MeshType>
void partitionMesh<MeshType>::setElementParameters()
{
    switch (MeshType::ElementShape::Shape)
    {
    case HEXA:
        M_elementNodes = 8;
        M_elementFaces = 6;
        M_elementEdges = 12;
        M_faceNodes    = 4;
        break;
    case TETRA:
        M_elementNodes = 4;
        M_elementFaces = 4;
        M_elementEdges = 6;
        M_faceNodes    = 3;
        break;
    default:
        ERROR_MSG( "Face Shape not implemented in partitionMesh" );
    }
}

template<typename MeshType>
void partitionMesh<MeshType>::distributeElements(UInt numElements)
{
    // ParMETIS is able to work in parallel: how many processors does it have at hand?
    Int numProcessors = M_comm->NumProc();
    M_me              = M_comm->MyPID();

    // CAREFUL: ParMetis works on a graph abstraction.
    // A graph is built over the data structure to be split, each vertex being a mesh element
    // so hereby a "vertex" is actually a _graph_ vertex, i. e. a mesh element
    M_vertexDistribution.resize(numProcessors + 1);
    M_vertexDistribution[0] = 0;

    UInt k = numElements;

    // Evenly distributed graph vertices
    for (Int i = 0; i < numProcessors; ++i)
    {
        UInt l = k / (numProcessors - i);
        M_vertexDistribution[i + 1] = M_vertexDistribution[i] + l;
        k -= l;
    }
    ASSERT(k == 0, "At this point we should have 0 volumes left") ;
}

template<typename MeshType>
void partitionMesh<MeshType>::findRepeatedFacesFSI()
{
    std::vector<Int>                    myRepeatedFace; // used for the solid partitioning
    boost::shared_ptr<std::vector<Int> > myIsOnProc;     // used for the solid partitioning

    myIsOnProc.reset(new std::vector<Int>(M_originalMesh->numVolumes()));

    bool myFaceRep;
    bool myFace(false);
    short count;
    for (UInt h = 0; h < M_originalMesh->numVolumes(); ++h)
    {
        (*myIsOnProc)[h] = -1;
    }

    // This loop is throughout the whole unpartitioned mesh,
    // it is expensive and not scalable.
    // Bad, this part should be done offline

    for (UInt ie = 1; ie <= M_originalMesh->numVolumes(); ++ie)
    {
        for (UInt iface = 1; iface <= M_elementFaces; ++iface)
        {
            UInt face = M_originalMesh->localFaceId(ie, iface);
            UInt vol  = M_originalMesh->face(face).ad_first();
            if (vol == ie)
            {
                vol = M_originalMesh->face(face).ad_second();
            }
            if (vol != 0)
            {
                myFace = false;
                myFaceRep = false;
                count = 0;
                for (Int ipoint = 1; ipoint <= static_cast<Int>(M_faceNodes); ++ipoint) // vertex-based dofs
                {
                    myFaceRep = ((M_interfaceMap->LID(M_originalMesh->face(face).point(ipoint).id())
                                  /* first is fluid */ == -1) &&
                                 (M_interfaceMapRep->LID(M_originalMesh->face(face).point(ipoint).id())
                                  /* first is fluid */ != -1));
                    myFace = myFace ||
                        (M_interfaceMap->LID(M_originalMesh->face(face).point(ipoint).id()) != -1);
                    if (myFaceRep)
                    {
                        ++count;
                    }
                }
                if (count > 1)
                {
                    myRepeatedFace.push_back(1);
                }
                else
                {
                    myRepeatedFace.push_back(0);
                }
            }
            if (myFace)
            {
                (*myIsOnProc)[ie-1] = M_me;
            }
        }
    }

    M_repeatedFace.reset(new std::vector<Int> (myRepeatedFace.size()));
    M_isOnProc.reset(new std::vector<Int> (*myIsOnProc));

    // Lot of communication here!!
    MPI_Allreduce(&myRepeatedFace[0], &(*M_repeatedFace)[0], myRepeatedFace.size(),
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&(*myIsOnProc)[0], &(*M_isOnProc)[0], myIsOnProc->size(),
                  MPI_INT, MPI_MAX, MPI_COMM_WORLD);
}

template<typename MeshType>
void partitionMesh<MeshType>::partitionConnectivityGraph(UInt numParts)
{
    // This array's size is equal to the number of locally-stored vertices:
    // at the end of the partitioning process, "M_graphVertexLocations" will contain the partitioning array:
    // M_graphVertexLocations[m] = n; means that graph vertex m belongs to subdomain n
    M_graphVertexLocations.resize(M_vertexDistribution[M_me + 1] - M_vertexDistribution[M_me]);

    // Now each processor will take care of its own graph vertices (i. e. mesh elements).
    // Nothing guarantees about the neighbor elements distribution across the processors,
    // since as of now we just split the set of volumes based on IDs.
    // Here we building up the neighbor arrays.

    UInt localStart = M_vertexDistribution[M_me] + 1;
    UInt localEnd   = M_vertexDistribution[M_me + 1] + 1;

    // this vector contains the weights for the edges of the graph,
    // it is set to null if it is not used.
    std::vector<Int> edgeWeights;

    M_adjacencyGraphKeys.resize(0);
    M_adjacencyGraphKeys.push_back(0);

    UInt sum = 0;

    for (UInt ie = localStart; ie < localEnd; ++ie)
    {
        for (UInt iface = 1; iface <= M_elementFaces; ++iface)
        {
            // global ID of the iface-th face in element ie
            UInt face = M_originalMesh->localFaceId(ie, iface);
            // first adjacent element to face "face"
            UInt elem = M_originalMesh->face(face).ad_first();
            if (elem == ie)
            {
                elem = M_originalMesh->face(face).ad_second();
            }
            if (elem != 0)
            {
                // this is the list of adjacency
                // for each graph vertex, simply push back the ID of its neighbors
                M_adjacencyGraphValues.push_back(elem - 1);
                ++sum;
                if (M_interfaceMap) // if I'm partitioning the solid in FSI
                {
                    if ((*M_repeatedFace)[sum])
                    {
                        edgeWeights.push_back(0);
                    }
                    else
                    {
                        edgeWeights.push_back(10);
                    }
                }
            }
        }
        // this is the list of "keys" to access M_adjacencyGraphValues
        // graph element i has neighbors M_adjacencyGraphValues[ k ],
        // with M_adjacencyGraphKeys[i] <= k < M_adjacencyGraphKeys[i+1]
        M_adjacencyGraphKeys.push_back(sum);
    }

    // **************
    // parMetis part

    // this array is to be used for weighted nodes on the graph:
    // usually we will set it to NULL

    Int* weightVector = 0;

    Int weightFlag;
    if (M_interfaceMap)
    {
        weightFlag = 1;
    }
    else
    {
        weightFlag = 0;
    }

    Int ncon = 1;
    Int numflag = 0;

    Int cutEdges; // here will be stored the number of edges cut in the partitioning process

    // additional options
    std::vector<Int>  options(3,0);
    options[0] = 1; // means that additional options are actually passed
    options[1] = 3; // level of information to be returned during execution (see ParMETIS's defs.h file)
    options[2] = 1; // random number seed for the ParMETIS routine

    // fraction of vertex weight to be distributed to each subdomain.
    // here we want the subdomains to be of the same size
    std::vector<float> tpwgts(ncon * numParts, 1. / numParts);
    // imbalance tolerance for each vertex weight
    std::vector<float> ubvec(ncon, 1.05);

    boost::shared_ptr<Epetra_MpiComm> mpiComm = boost::dynamic_pointer_cast <Epetra_MpiComm> (M_comm);
    MPI_Comm MPIcomm = mpiComm->Comm();

    Int nprocs;
    MPI_Comm_size(MPIcomm, &nprocs);

    /*
      (from ParMETIS v 3.1 manual)
      This routine is used to compute a k-way partitioning of a graph
      on p processors using the multilevel k-way multi-constraint
      partitioning algorithm.
    */

    Int numberParts = (Int) numParts;

    Int* adjwgtPtr(0);
    if (edgeWeights.size() > 0)
    {
        adjwgtPtr = static_cast<Int*>(&edgeWeights[0]);
    }
    ParMETIS_V3_PartKway(static_cast<Int*>(&M_vertexDistribution[0]),
                         static_cast<Int*>(&M_adjacencyGraphKeys[0]),
                         static_cast<Int*>(&M_adjacencyGraphValues[0]),
                         weightVector, adjwgtPtr, &weightFlag, &numflag,
                         &ncon, &numberParts, &tpwgts[0], &ubvec[0],
                         &options[0], &cutEdges, &M_graphVertexLocations[0],
                         &MPIcomm);

    M_comm->Barrier();

    Int nProc;
    nProc = M_comm->NumProc();

    // this is a vector of subdomains: each component is
    // the list of vertices belonging to the specific subdomain
    (*M_elementDomains).resize(numParts);

    // cycling on locally stored vertices
    for (UInt ii = 0; ii < M_graphVertexLocations.size(); ++ii)
    {
        // here we are associating the vertex global ID to the subdomain ID
        (*M_elementDomains)[M_graphVertexLocations[ii]].push_back(ii + M_vertexDistribution[M_me]);
    }
}

template<typename MeshType>
void partitionMesh<MeshType>::matchFluidPartitionsFSI()
{
    boost::shared_ptr<Epetra_MpiComm> mpiComm = boost::dynamic_pointer_cast <Epetra_MpiComm> (M_comm);
    MPI_Comm MPIcomm = mpiComm->Comm();
    Int numProcesses;
    MPI_Comm_size(MPIcomm, &numProcesses);

    Int procOrder[numProcesses];
    std::vector<std::vector<UInt> > myMatchesForProc(numProcesses);
    std::vector<std::vector<UInt> > matchesForProc(numProcesses);
    bool orderingError[numProcesses];

    for (Int i=0; i<numProcesses ; ++i)
    {
        orderingError[i]=false;
        for (Int j=0; j<numProcesses ; ++j)
        {
            myMatchesForProc[i].push_back(0);
            matchesForProc[i].push_back(0);
        }
    }

    for (UInt kk=0; kk<M_graphVertexLocations.size(); ++kk)
    {
        if ((*M_isOnProc)[kk+M_vertexDistribution[M_me]]!=-1)
        {
            ++myMatchesForProc[M_graphVertexLocations[kk]][(*M_isOnProc)[kk+M_vertexDistribution[M_me]]];
        }
    }

    for (UInt j=0; (Int)j<numProcesses; ++j)
    {
        MPI_Allreduce(&myMatchesForProc[j][0], &matchesForProc[j][0], numProcesses,
                      MPI_INT, MPI_SUM, MPIcomm);
    }

    M_comm->Barrier();

    Int suitableProcess = -1;
    UInt max = 0;

    for (Int ii = 0; ii<numProcesses; ++ii)
    {
        if (matchesForProc[M_me][ii] > max)
        {
            suitableProcess = ii;
            max = matchesForProc[M_me][ii];
        }
    }

    ASSERT(suitableProcess != -1, "one partition is without interface nodes!");
    procOrder[M_me] = suitableProcess;

    M_comm->Barrier();

    std::vector<UInt> maxs(numProcesses);
    maxs[M_me] = max;
    for (Int j = 0; j < numProcesses ; ++j) // Allgather
    {
        MPI_Bcast(&maxs[j], 1, MPI_INT, j, MPIcomm); // perhaps generates errors
    }

    std::vector<std::pair<UInt, Int> > procIndex(numProcesses);
    for (Int k = 0; k < numProcesses; ++k)
    {
        procIndex[k] = std::make_pair( maxs[k], k);
    }

    std::sort(procIndex.begin(), procIndex.end() /*, &booleanCondition::reordering*/);

    for (Int l=0; l<numProcesses; ++l)
    {
        for (Int l=0; l<numProcesses; ++l)
        {
            for (Int j=0; j<numProcesses ; ++j) // Allgather
            {
                MPI_Bcast( &procOrder[j], 1, MPI_INT, j, MPIcomm); // perhaps generates errors
            }
        }
    }

    std::vector< std::vector<Int> > locProc2((*M_elementDomains));
    for (Int j = numProcesses; j > 0 ; --j)
    {
        if (orderingError[procOrder[procIndex[j - 1].second]] == false)
        {
            (*M_elementDomains)[procOrder[procIndex[j - 1].second]] = locProc2[procIndex[j - 1].second];
        }
        else
        {
            std::cout << "Ordering error when assigning the processor"
                      << M_me << " to the partition," << std::endl
                      << " parmetis did a bad job." << std::endl;
            for (Int i = numProcesses; i > 0; --i)
            {
                if (orderingError[procIndex[i - 1].second] == false) // means that i is the first
                                                                     // proc not assigned
                {
                    procOrder[procIndex[j - 1].second] = procIndex[i - 1].second;
                    (*M_elementDomains)[procIndex[i - 1].second] = locProc2[procIndex[j - 1].second];
                    break;
                }
            }
        }
        orderingError[procOrder[procIndex[j - 1].second]] = true;
    }
}

template<typename MeshType>
void partitionMesh<MeshType>::redistributeElements()
{
    boost::shared_ptr<Epetra_MpiComm> mpiComm = boost::dynamic_pointer_cast <Epetra_MpiComm> (M_comm);
    MPI_Comm MPIcomm = mpiComm->Comm();
    Int numProcesses;
    MPI_Comm_size(MPIcomm, &numProcesses);

    Int maxInt (1000);
    Int sendSize[numProcesses];
    Int receiveSize[numProcesses];
    // cycling on subdomains
    // TODO: Matteo please comment this part :)

    MPI_Status status;
    Int size;

    // MPI_Status  recv_status;
    // MPI_Request send_request;

    for (Int iproc = 0; iproc < numProcesses; ++iproc)
    {
        sendSize[iproc] = (*M_elementDomains)[iproc].size();
    }
    MPI_Alltoall(sendSize, 1, MPI_INT, receiveSize, 1, MPI_INT, MPIcomm);

    for (Int iproc = 0; iproc < numProcesses; ++iproc)
    {
        if (static_cast<Int>(M_me) != iproc)
        {
            size = sendSize[iproc];
            // workaround for huge data to be passed
            if (size > maxInt)
            {
                Int incr = 1 ;
                Int pos = 0;
                Int sizePart = size;

                // divide the whole data set into smaller packets
                while (sizePart > maxInt)
                {
                    incr += 1;
                    sizePart = size / incr;
                }

                MPI_Send(&incr, 1, MPI_INT, iproc, 20, MPIcomm);
                MPI_Send(&sizePart, 1, MPI_INT, iproc, 30, MPIcomm);

                for (Int kk = 0; kk < incr; ++kk)
                {
                    MPI_Send(&pos, 1, MPI_INT, iproc, 100+kk, MPIcomm);
                    MPI_Send(&(*M_elementDomains)[iproc][pos], sizePart, MPI_INT, iproc,
                             5000000+kk, MPIcomm);
                    pos = pos + sizePart;
                }

                Int resto = size % incr;

                MPI_Send(&resto, 1, MPI_INT, iproc, 80, MPIcomm);

                if (resto != 0)
                {
                    MPI_Send(&pos, 1, MPI_INT, iproc, 40, MPIcomm);
                    MPI_Send(&(*M_elementDomains)[iproc][pos], resto, MPI_INT, iproc, 50, MPIcomm);
                }
            }
            else
            {
                if (size != 0)
                {
                    MPI_Send(&(*M_elementDomains)[iproc][0], size, MPI_INT, iproc, 60, MPIcomm);
                }
            }
        }
        else
        {
            for (Int jproc = 0; jproc < numProcesses; ++jproc)
            {
                if (jproc != iproc)
                {
                    size = receiveSize[jproc];
                    std::vector<Int> stack(size, 0);

                    if (size > maxInt)
                    {
                        Int sizePart, pos, incr;

                        MPI_Recv(&incr, 1, MPI_INT, jproc, 20, MPIcomm, &status);
                        MPI_Recv(&sizePart, 1, MPI_INT, jproc, 30, MPIcomm, &status);

                        for (Int kk = 0; kk < incr; ++kk)
                        {
                            MPI_Recv(&pos, 1, MPI_INT, jproc, 100+kk, MPIcomm, &status);
                            MPI_Recv(&stack[pos], sizePart , MPI_INT, jproc, 5000000+kk, MPIcomm, &status);
                        }
                        Int resto = 0;
                        MPI_Recv(&resto, 1, MPI_INT, jproc, 80, MPIcomm, &status);

                        if (resto != 0)
                        {
                            MPI_Recv(&pos, 1, MPI_INT, jproc, 40, MPIcomm, &status);
                            MPI_Recv(&stack[pos],  resto, MPI_INT, jproc, 50, MPIcomm, &status);
                        }
                    }
                    else
                    {
                        if (size != 0)
                        {
                            MPI_Recv(&stack[0], size , MPI_INT, jproc, 60, MPIcomm, &status);
                        }
                    }
                    for (Int jj = 0; jj < size; ++jj)
                    {
                        (*M_elementDomains)[M_me].push_back(stack[jj]);
                    }
                }
            }
        }
    }
}

template<typename MeshType>
void partitionMesh<MeshType>::constructLocalMesh()
{
    if (!M_me)
    {
        std::cout << "Building local mesh ..." << std::endl;
    }

    std::map<Int, Int>::iterator  im;
    std::set<Int>::iterator       is;

    Int count = 1;
    UInt ielem;
    UInt inode;

    for (UInt i = 0; i < M_numPartitions; ++i)
    {
        count = 1;
        // cycle on local element's ID

        UInt me = M_serialMode ? i : M_me;

        for (UInt jj = 0; jj < (*M_elementDomains)[me].size(); ++jj)
        {
            ielem = (*M_elementDomains)[me][jj];
            M_localVolumes[i].push_back(ielem);

            // cycle on element's nodes
            for (UInt ii = 1; ii <= M_elementNodes; ++ii)
            {
                inode = M_originalMesh->volume(ielem + 1).point(ii).id();
                im    = M_globalToLocalNode[i].find(inode);

                // if the node is not yet present in the list of local nodes, then add it
                // CAREFUL: also local numbering starts from 1 in RegionMesh
                if (im == M_globalToLocalNode[i].end())
                {
                    M_globalToLocalNode[i].insert(std::make_pair(inode, count));
                    ++count;
                    // store here the global numbering of the node
                    M_localNodes[i].push_back(M_originalMesh->volume(ielem + 1).point(ii).id());
                }
            }

            // cycle on element's edges
            for (UInt ii = 1; ii <= M_elementEdges; ++ii)
            {
                // store here the global numbering of the edge
                M_localEdges[i].insert(M_originalMesh->localEdgeId(ielem + 1, ii));
            }

            // cycle on element's faces
            for (UInt ii = 1; ii <= M_elementFaces; ++ii)
            {
                // store here the global numbering of the face
                M_localFaces[i].insert(M_originalMesh->localFaceId(ielem + 1, ii));
            }
        }
    }
}

template<typename MeshType>
void partitionMesh<MeshType>::constructNodes()
{
    UInt inode;
    for (UInt i = 0; i < M_numPartitions; ++i)
    {
        std::vector<Int>::iterator it;

        M_nBoundaryPoints[i] = 0;
        (*M_meshPartitions)[i]->pointList.reserve(M_localNodes[i].size());
        // guessing how many boundary points on this processor.
        (*M_meshPartitions)[i]->_bPoints.reserve(M_originalMesh->numBPoints() * M_localNodes[i].size() /
                                                 M_originalMesh->numBPoints());

        inode = 1;

        typename MeshType::PointType *pp = 0;

        // loop in the list of local nodes:
        // in this loop inode is the local numbering of the points
        for (it = M_localNodes[i].begin(); it != M_localNodes[i].end(); ++it, ++inode)
        {
            typename MeshType::PointType point = 0;

            // create a boundary point in the local mesh, if needed
            bool boundary = M_originalMesh->isBoundaryPoint(*it);
            if (boundary)
            {
                ++M_nBoundaryPoints[i];
            }

            pp = &(*M_meshPartitions)[i]->addPoint(boundary);
            pp->setMarker(M_originalMesh->point(*it).marker());

            pp->x() = M_originalMesh->point(*it).x();
            pp->y() = M_originalMesh->point(*it).y();
            pp->z() = M_originalMesh->point(*it).z();

            UInt id = M_originalMesh->point(*it).id();

            pp->setId(id);
            pp->setLocalId(inode);

            (*M_meshPartitions)[i]->localToGlobalNode().insert(std::make_pair(inode, id));
            (*M_meshPartitions)[i]->globalToLocalNode().insert(std::make_pair(id, inode));
        }
    }
}

template<typename MeshType>
void partitionMesh<MeshType>::constructVolumes()
{
    Int count;
    for (UInt i = 0; i < M_numPartitions; ++i)
    {
        std::map<Int, Int>::iterator im;
        std::vector<Int>::iterator it;
        count = 1;
        UInt inode;

        typename MeshType::VolumeType * pv = 0;

        (*M_meshPartitions)[i]->volumeList.reserve(M_localVolumes[i].size());

        // loop in the list of local elements
        // CAREFUL! in this loop inode is the global numbering of the points
        // We insert the local numbering of the nodes in the local volume list
        for (it = M_localVolumes[i].begin(); it != M_localVolumes[i].end(); ++it, ++count)
        {
            pv = &((*M_meshPartitions)[i]->addVolume());
            // CAREFUL! in ParMETIS data structures, numbering starts from 0
            pv->setId (M_originalMesh->volume(*it + 1).id());
            pv->setLocalId(count);

            M_globalToLocalVolume[i].insert(std::make_pair(M_originalMesh->volume(*it + 1).id(), count));

            for (ID id = 1; id <= M_elementNodes; ++id)
            {
                inode = M_originalMesh->volume(*it + 1).point(id).id();
                // im is an iterator to a map element
                // im->first is the key (i. e. the global ID "inode")
                // im->second is the value (i. e. the local ID "count")
                im = M_globalToLocalNode[i].find(inode);
                pv->setPoint(id, (*M_meshPartitions)[i]->pointList( (*im).second ));
            }

            Int ibc = M_originalMesh->volume(*it + 1).marker();

            pv->setMarker(entityFlag_Type( ibc ));
        }
    }
}

template<typename MeshType>
void partitionMesh<MeshType>::constructEdges()
{
    Int count;
    for (UInt i = 0; i < M_numPartitions; ++i)
    {
        std::map<Int, Int>::iterator im;
        std::set<Int>::iterator is;

        typename MeshType::EdgeType * pe;
        UInt inode;
        count = 1;

        M_nBoundaryEdges[i] = 0;
        (*M_meshPartitions)[i]->edgeList.reserve(M_localEdges[i].size());

        // loop in the list of local edges
        for (is = M_localEdges[i].begin(); is != M_localEdges[i].end(); ++is, ++count)
        {
            // create a boundary edge in the local mesh, if needed
            bool boundary = (M_originalMesh->isBoundaryEdge(*is));
            if (boundary)
            {
                // create a boundary edge in the local mesh, if needed
                ++M_nBoundaryEdges[i];
            }

            pe = &(*M_meshPartitions)[i]->addEdge(boundary);

            pe->setId (M_originalMesh->edge(*is).id());
            pe->setLocalId(count);

            for (ID id = 1; id <= 2; ++id)
            {
                inode = M_originalMesh->edge(*is).point(id).id();
                // im is an iterator to a map element
                // im->first is the key (i. e. the global ID "inode")
                // im->second is the value (i. e. the local ID "count")
                im = M_globalToLocalNode[i].find(inode);
                pe->setPoint(id, (*M_meshPartitions)[i]->pointList((*im).second));
            }
            pe->setMarker(M_originalMesh->edge(*is).marker());
        }
    }
}

template<typename MeshType>
void partitionMesh<MeshType>::constructFaces()
{
    Int count;
    for (UInt i = 0; i < M_numPartitions; ++i)
    {
        std::map<Int, Int>::iterator im;
        std::set<Int>::iterator      is;

        typename MeshType::FaceType * pf = 0;

        UInt inode;
        count = 1;

        M_nBoundaryFaces[i] = 0;
        (*M_meshPartitions)[i]->faceList.reserve(M_localFaces[i].size());

        // loop in the list of local faces
        for (is = M_localFaces[i].begin(); is != M_localFaces[i].end(); ++is, ++count)
        {
            // create a boundary face in the local mesh, if needed
            bool boundary = (M_originalMesh->isBoundaryFace(*is));
            if (boundary)
            {
                ++M_nBoundaryFaces[i];
            }

            pf =  &(*M_meshPartitions)[i]->addFace(boundary);

            pf->setId (M_originalMesh->face(*is).id());
            pf->setLocalId(count);

            Int elem1 = M_originalMesh->face(*is).ad_first();
            Int elem2 = M_originalMesh->face(*is).ad_second();

            // find the mesh elements adjacent to the face
            im =  M_globalToLocalVolume[i].find(elem1);

            Int localElem1;

            if (im == M_globalToLocalVolume[i].end())
            {
                localElem1 = 0;
            }
            else
            {
                localElem1 = (*im).second;
            }

            im =  M_globalToLocalVolume[i].find(elem2);

            Int localElem2;
            if (im == M_globalToLocalVolume[i].end())
            {
                localElem2 = 0;
            }
            else
            {
                localElem2 = (*im).second;
            }

            // if this process does not own either of the adjacent elements
            // then the two adjacent elements and the respective face positions coincide in the local mesh
            // possible bug fixed: not only the two adjacent elements face, but also the face
            // positions should coincide.
            // otherwise it could happen that a pair(element, position) is associated to different faces.
            // This can lead to a wrong treatment of the dofPerFace (in 2D of the dofPerEdge, as occurred
            // with P2)

            if ((localElem1 == 0) && !boundary)
            {
                pf->ad_first() = localElem2;
                pf->pos_first() = M_originalMesh->face(*is).pos_second();
            }
            else
            {
                pf->ad_first() = localElem1;
                pf->pos_first() = M_originalMesh->face(*is).pos_first();
            }

            if ((localElem2 == 0) && !boundary)
            {
                pf->ad_second() = localElem1;
                pf->pos_second() = M_originalMesh->face(*is).pos_first();
            }
            else
            {
                pf->ad_second()  = localElem2;
                pf->pos_second() = M_originalMesh->face(*is).pos_second();
            }


            for (ID id = 1; id <= M_originalMesh->face(*is).numLocalVertices; ++id)
            {
                inode = M_originalMesh->face(*is).point(id).id();
                im = M_globalToLocalNode[i].find(inode);
                pf->setPoint(id, (*M_meshPartitions)[i]->pointList((*im).second));
            }

            pf->setMarker(M_originalMesh->face(*is).marker());

            (*M_meshPartitions)[i]->setLinkSwitch("HAS_ALL_FACES");
            (*M_meshPartitions)[i]->setLinkSwitch("FACES_HAVE_ADIACENCY");
        }
    }
}

template<typename MeshType>
void partitionMesh<MeshType>::finalSetup()
{
    for (UInt i = 0; i < M_numPartitions; ++i)
    {
        UInt nVolumes = M_localVolumes[i].size();
        UInt nNodes   = M_localNodes[i].size();
        UInt nEdges   = M_localEdges[i].size();
        UInt nFaces   = M_localFaces[i].size();

        (*M_meshPartitions)[i]->setMaxNumPoints (nNodes, true);
        (*M_meshPartitions)[i]->setMaxNumEdges  (nEdges, true);
        (*M_meshPartitions)[i]->setMaxNumFaces  (nFaces, true);
        (*M_meshPartitions)[i]->setMaxNumVolumes( nVolumes, true);

        (*M_meshPartitions)[i]->setMaxNumGlobalPoints (M_originalMesh->numPoints());
        (*M_meshPartitions)[i]->setNumGlobalVertices  (M_originalMesh->numPoints());
        (*M_meshPartitions)[i]->setMaxNumGlobalEdges  (M_originalMesh->numEdges());
        (*M_meshPartitions)[i]->setMaxNumGlobalFaces  (M_originalMesh->numFaces());

        (*M_meshPartitions)[i]->setMaxNumGlobalVolumes(M_originalMesh->numVolumes());
        (*M_meshPartitions)[i]->setNumBFaces    (M_nBoundaryFaces[i]);

        (*M_meshPartitions)[i]->setNumBPoints   (M_nBoundaryPoints[i]);
        (*M_meshPartitions)[i]->setNumBEdges    (M_nBoundaryEdges[i]);

        (*M_meshPartitions)[i]->setNumVertices (nNodes );
        (*M_meshPartitions)[i]->setNumBVertices(M_nBoundaryPoints[i]);

        (*M_meshPartitions)[i]->updateElementEdges();

        (*M_meshPartitions)[i]->updateElementFaces();

        if (M_serialMode)
        {
            std::cout << "Created local mesh number " << i + 1 << std::endl;
        }
        else
        {
            std::cout << "Rank " << M_me << " created local mesh." << std::endl;
        }
    }
}

template<typename MeshType>
void partitionMesh<MeshType>::createRepeatedMap()
{
    std::set<Int>    repeatedNodeList;
    std::set<Int>    repeatedEdgeList;
    std::set<Int>    repeatedFaceList;

    if (! M_me)
    {
        std::cout << "Building repeated map... " << std::endl;
    }

    for (UInt i = 0; i < M_numPartitions; ++i)
    {
        std::set<Int>::iterator is;

        UInt me = M_serialMode ? i : M_me;

        std::vector<Int> elementList = (*M_elementDomains)[me];

        UInt inode, ielem;

        // repeated element map creation

        // use sets to store each entity only once
        repeatedNodeList.clear();
        repeatedEdgeList.clear();
        repeatedFaceList.clear();

        for (UInt ii = 0; ii < elementList.size(); ++ii)
        {
            ielem = elementList[ii];
            M_repeatedVolumeVector[i].push_back(ielem + 1);
            for (UInt jj = 1; jj <= M_elementNodes; ++jj)
            {
                inode = M_originalMesh->volume(ielem + 1).point(jj).id();
                repeatedNodeList.insert(inode);
            }
            for (UInt jj = 1; jj <= M_elementEdges; ++jj)
            {
                UInt iedge = M_originalMesh->localEdgeId(ielem + 1, jj);
                repeatedEdgeList.insert((Int) iedge);
            }
            for (UInt jj = 1; jj <= M_elementFaces; ++jj)
            {
                UInt iface = M_originalMesh->localFaceId(ielem + 1, jj);
                repeatedFaceList.insert(iface);
            }
        }

        // repeated node map creation
        M_repeatedNodeVector[i].reserve(repeatedNodeList.size());

        for (is = repeatedNodeList.begin(); is != repeatedNodeList.end(); ++is)
        {
            M_repeatedNodeVector[i].push_back(*is);
        }

        // repeated edge list creation
        M_repeatedEdgeVector[i].reserve(repeatedEdgeList.size());

        for (is = repeatedEdgeList.begin(); is != repeatedEdgeList.end(); ++is)
        {
            M_repeatedEdgeVector[i].push_back(*is);
        }

        // repeated face list creation
        M_repeatedFaceVector[i].reserve(repeatedFaceList.size());

        for (is = repeatedFaceList.begin(); is != repeatedFaceList.end(); ++is)
        {
            M_repeatedFaceVector[i].push_back(*is);
        }

        if (M_serialMode)
        {
            std::cout <<  "Created repeated map number " << i + 1 << std::endl;
        }
        else
        {
            std::cout << "Rank " << M_me << " created repeated map." << std::endl;
        }
    }
}

template<typename MeshType>
void partitionMesh<MeshType>::execute()
{
    // Set element parameters (number of nodes, faces, edges and number of nodes
    // on each face according to the type of mesh element used.
    setElementParameters();

    // Build graph vertex distribution vector. Graph vertex represents one element
    // in the mesh.
    distributeElements(M_originalMesh->numElements());


    // In fluid-structure interaction:
    // *    If the solid mesh is not partitioned the following part won't be
    //      executed
    // *    If the solid mesh is partitioned:
    //      - The fluid is partitioned first
    //      - The solid mesh partition tries to follow the partition of the fluid
    //      This is achieved by specifying a weight to some edge of the graph.
    //      The interface between two processors is the set of the nodes that for
    //      at least one processor are on the repeated map and not on the unique map.
    //      That's why the constructor needs both the unique and repeated maps
    //      on the interface


    //////////////////// BEGIN OF SOLID PARTITION PART ////////////////////////
    if (M_interfaceMap)
    {
        findRepeatedFacesFSI();
    }
    //////////////////// END OF SOLID PARTITION PART ////////////////////////

    // Partition connectivity graph
    partitionConnectivityGraph(M_comm->NumProc());

    //////////////// BEGIN OF SOLID PARTITION PART ////////////////
    if (M_interfaceMap)
    {
        matchFluidPartitionsFSI();
    }
    ////////////////// END OF SOLID PARTITION PART /////////////////////


    // Redistribute elements to appropriate processors before building the
    // partitioned mesh.
    redistributeElements();


#ifdef DEBUG
    Debug(4000) << M_me << " has " << (*M_elementDomains)[M_me].size() << " elements.\n";
#endif

    // ***********************
    // local mesh construction
    // ***********************
    constructLocalMesh();

    // ******************
    // nodes construction
    // ******************
    constructNodes();

    // ********************
    // volumes construction
    // ********************
    constructVolumes();

    // ******************
    // edges construction
    // ******************
    constructEdges();

    // ******************
    // faces construction
    // ******************
    constructFaces();

    // ******************
    // final setup
    // ******************
    finalSetup();

    // *********************
    // repeated map creation
    // *********************
    createRepeatedMap();

    // ***************************************
    // release the original unpartitioned mesh
    // allowing it to be deleted
    // ***************************************
    releaseUnpartitionedMesh();
}

}
#endif // PARTMESH_H
