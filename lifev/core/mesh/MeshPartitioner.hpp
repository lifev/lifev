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
  @contributor Mauro Perego <mperego@fsu.edu>
  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#ifndef MESH_PARTITIONER_H
#define MESH_PARTITIONER_H 1

#include <fstream>
#include <sstream>


#include <parmetis.h>
#include <Epetra_MpiComm.h>


#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/Switch.hpp>
#include <lifev/core/mesh/MeshElementMarked.hpp>
#include <lifev/core/util/LifeDebug.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/array/GhostHandler.hpp>

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
class MeshPartitioner
{
public:
    //@{
    // Make the template's type available to the outside
    typedef MeshType mesh_Type;
    typedef boost::shared_ptr<MeshType> meshPtr_Type;
    typedef std::vector<Int> idList_Type;
    typedef std::vector<idList_Type> graph_Type;
    typedef boost::shared_ptr<graph_Type> graphPtr_Type;
    typedef std::vector<meshPtr_Type> partMesh_Type;
    typedef boost::shared_ptr<partMesh_Type> partMeshPtr_Type;
    //@}
    //! \name Constructors & Destructors
    //@{
    //! Default empty constructor
    MeshPartitioner();

    //! Constructor
    /*!
      This is a non-empty constructor. It takes as parameters the
      unpartitioned mesh (reference), the Epetra_Comm object in
      use (reference) and pointers to the Epetra interface and
      repeated interface maps. The constructor initializes the
      data members and calls a private method
      MeshPartitioner::execute which handles the mesh partitioning.
      \param mesh - Mesh& - the unpartitioned mesh
      \param _comm - Epetra_Comm& - Epetra communicator object
      \param interfaceMap - Epetra_Map*
      \param interfaceMapRep - Epetra_Map*
    */
    MeshPartitioner ( meshPtr_Type& mesh,
                      boost::shared_ptr<Epetra_Comm>& comm,
                      Epetra_Map* interfaceMap = 0,
                      Epetra_Map* interfaceMapRep = 0 );

    //! Empty destructor
    virtual ~MeshPartitioner() {}
    //@}

    //! \name Public Methods
    //@{
    //! Partition the mesh.
    /*!
      It takes as parameters the unpartitioned mesh (reference),
      the Epetra_Comm object in use (reference) and pointers to
      the Epetra interface and repeated interface maps.
      The method initializes the data members and calls a private method
      MeshPartitioner::execute which handles the mesh partitioning.
      @param mesh - Mesh& - the unpartitioned mesh
      @param _comm - Epetra_Comm& - Epetra communicator object
      @param interfaceMap - Epetra_Map*
      @param interfaceMapRep - Epetra_Map*
      @note This method is meant to be used with the empty constructor.
    */
    void doPartition ( meshPtr_Type& mesh,
                       boost::shared_ptr<Epetra_Comm>& comm,
                       Epetra_Map* interfaceMap = 0,
                       Epetra_Map* interfaceMapRep = 0 );


    //! To be used with the new constructor.
    /*!
      Loads the parameters of the partitioning process from the simulation data file.
      Allocates and initializes data members according to the number of partitions specified in
      the data file.
      \param numPartitions - UInt - the number of partitions desired, in the offline case
      \param _comm - Epetra_Comm& - reference of the Epetra communicator used
    */
    void setup (UInt numPartitions, boost::shared_ptr<Epetra_Comm> comm);

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
    void attachUnpartitionedMesh (meshPtr_Type& mesh, Epetra_Map* interfaceMap = 0,
                                  Epetra_Map* interfaceMapRep = 0);

    //! Releases the original unpartitioned mesh
    /*!
      Releases the unpartitioned mesh so that it can be deleted, freeing A LOT of memory
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

    //! Initialize M_entityPID
    void fillEntityPID();

    // Next method should be renamed and become a regular method
    //! Return a pointer to the mesh partition with rank k
    const meshPtr_Type& LIFEV_DEPRECATED ( getPartition (Int k) const )
    {
        return (*M_meshPartitions) [k];
    }

    //! Clean structures that are not needed after partitioning
    void cleanUp();

    //! Prints information about the state (data) of the object
    void showMe (std::ostream& output = std::cout) const;
    //@}

    //! \name Get Methods
    //@{
    //! Return a reference to M_vertexDistribution
    const std::vector<Int>&  vertexDistribution()   const
    {
        return M_vertexDistribution;
    }
    //! Return a const pointer to M_meshPartitions[0] - for parallel
    const meshPtr_Type&      meshPartition()        const
    {
        return (*M_meshPartitions) [0];
    }
    meshPtr_Type&            meshPartition()
    {
        return (*M_meshPartitions) [0];
    }
    //! Return a pointer to M_meshPartitions
    const partMeshPtr_Type&  meshPartitions()       const
    {
        return M_meshPartitions;
    }
    //! Return a pointer to M_graphVertexLocations
    const std::vector<Int>&  graphVertexLocations() const
    {
        return M_graphVertexLocations;
    }
    //! Return a pointer to M_elementDomains
    const graphPtr_Type&     elementDomains()       const
    {
        return M_elementDomains;
    }
    graphPtr_Type&           elementDomains()
    {
        return M_elementDomains;
    }
    //! Return a pointer to the communicator M_comm
    const boost::shared_ptr<Epetra_Comm>& comm()   const
    {
        return M_comm;
    }
    //@}

    //! @name Set methos
    //@{

    //! Set M_buildOverlappingPartitions
    void setPartitionOverlap ( UInt const overlap )
    {
        M_partitionOverlap = overlap;
    }

    //@}

private:
    // Private copy constructor and assignment operator. No implementation
    MeshPartitioner (const MeshPartitioner&);
    MeshPartitioner& operator= (const MeshPartitioner&);
    //! Private Methods
    //@{

    //! Initialize the parameters with default value.
    void init ();

    //! Execute mesh partitioning using the configured MPI processes (online partitioning)
    /*!
      Executed the mesh partitioning using the number of MPI processes as the number of partitions.
      Sets current mesh element parameters: M_elementVertices, M_elementRidges, M_elementFacets, M_facetVertices
      Updates: M_elementDomains (indirectly)
      Other data members are changed indirectly by calling other private methods.
    */
    void execute();
    //! Build the graph vertex distribution vector
    /*!
      Updates the member M_vertexDistribution according to the number of processors to be
      used by ParMETIS (the number of processes started for MPI
      \param numElements - UInt - number of elements in the mesh
    */
    void distributeElements (UInt numElements);
    //! Find faces on the boundaries between domains (FSI)
    /*!
      Identifies the element faces that are common to both the fluid and the solid
      meshes and creates a map between the faces that reside on the boundary between
      the two meshes and the processors used. Updates the members M_repeatedFacet and
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
    void partitionConnectivityGraph (UInt numParts);

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
      Updates M_localNodes, M_localRidges, M_localFacets, M_localElements,
      M_globalToLocalNode.
    */
    void constructLocalMesh();
    //! Construct vertices
    /*!
      Adds vertices to the partitioned mesh object. Updates M_nBoundaryPoints,
      M_meshPartitions.
    */
    void constructNodes();
    //! Construct volumes
    /*!
      Adds volumes to the partitioned mesh object. Updates M_globalToLocalElement,
      M_meshPartitions.
    */
    void constructElements();
    //! Construct ridges
    /*!
      Adds ridges to the partitioned mesh object. Updates M_nBoundaryRidges,
      M_meshPartitions.
    */
    void constructRidges();
    //! Construct facets
    /*!
      Adds facets to the partitioned mesh object. Updates M_nBoundaryFacets,
      M_meshPartitions.
    */
    void constructFacets();
    //! Final setup of local mesh
    /*!
      Updates the partitioned mesh object data members after adding the mesh
      elements (vertices, ridges, facets, volumes).
      Updates M_meshPartitions.
    */
    void finalSetup();

    //! Mark ghost entities
    /*!
      Mark all ghost entities that have been added for overlapping partitions
      with the EntityFlags::GHOST flag in order to properly build the dof
      map in DOF::GlobalElements().
    */
    void markGhostEntities();

    //@}
    //! Private Data Members
    //@{
    UInt                                 M_numPartitions;
    partMeshPtr_Type                     M_meshPartitions;
    std::vector<Int>                     M_vertexDistribution;
    std::vector<Int>                     M_adjacencyGraphKeys;
    std::vector<Int>                     M_adjacencyGraphValues;
    boost::shared_ptr<Epetra_Comm>       M_comm;
    Int                                  M_me;

    std::vector<std::vector<Int> >       M_localNodes;
    std::vector<std::set<Int> >          M_localRidges;
    std::vector<std::set<Int> >          M_localFacets;
    std::vector<std::vector<Int> >       M_localElements;
    std::vector<std::map<Int, Int> >     M_globalToLocalNode;
    std::vector<std::map<Int, Int> >     M_globalToLocalElement;
    std::vector<UInt>                    M_nBoundaryPoints;
    std::vector<UInt>                    M_nBoundaryRidges;
    std::vector<UInt>                    M_nBoundaryFacets;
    // The following are utility variables used throughout the partitioning
    // process
    meshPtr_Type                         M_originalMesh;
    Epetra_Map*                          M_interfaceMap;
    Epetra_Map*                          M_interfaceMapRep;
    //! Number of partitions handled. 1 for parallel (old way), != 1 for serial
    UInt                                 M_elementVertices;
    UInt                                 M_elementFacets;
    UInt                                 M_elementRidges;
    UInt                                 M_facetVertices;
    boost::shared_ptr<std::vector<Int> > M_repeatedFacet;
    boost::shared_ptr<std::vector<Int> > M_isOnProc;
    std::vector<Int>                     M_graphVertexLocations;
    graphPtr_Type                        M_elementDomains;
    bool                                 M_serialMode; // how to tell if running serial partition mode
    UInt                                 M_partitionOverlap;

    //! Store ownership for each entity, subdivided by entity type
    struct EntityPIDList
    {
        idList_Type elements;
        idList_Type facets;
        idList_Type ridges;
        idList_Type points;
    } M_entityPID;

    //@}
}; // class MeshPartitioner

//
// IMPLEMENTATION
//


// =================================
// Constructors and destructor
// =================================

template < typename MeshType >
MeshPartitioner < MeshType >::
MeshPartitioner()
{
    init ();
} // constructor

template < typename MeshType >
MeshPartitioner < MeshType >::
MeshPartitioner ( meshPtr_Type& mesh, boost::shared_ptr<Epetra_Comm>& comm,
                  Epetra_Map* interfaceMap, Epetra_Map* interfaceMapRep)
{
    init ();
    doPartition ( mesh, comm, interfaceMap, interfaceMapRep );
} // constructor

template < typename MeshType >
void
MeshPartitioner < MeshType >::
init ()
{
    M_numPartitions = 1;
    M_localNodes.resize ( M_numPartitions );
    M_localRidges.resize ( M_numPartitions );
    M_localFacets.resize ( M_numPartitions );
    M_localElements.resize ( M_numPartitions );
    M_globalToLocalNode.resize ( M_numPartitions );
    M_globalToLocalElement.resize ( M_numPartitions );
    M_nBoundaryPoints.resize ( M_numPartitions );
    M_nBoundaryRidges.resize ( M_numPartitions );
    M_nBoundaryFacets.resize ( M_numPartitions );
    M_elementDomains.reset ( new graph_Type );
    M_serialMode = false;
    M_partitionOverlap = 0;

    /*
      Sets element parameters (nodes, faces, ridges and number of nodes on each
      facet according to the type of mesh element used (Mesh::ElementShape::S_shape).
      Updates M_elementVertices, M_elementFaces, M_elementRidges, M_facetVertices.
    */
    M_elementVertices = MeshType::elementShape_Type::S_numVertices;
    M_elementFacets   = MeshType::elementShape_Type::S_numFacets;
    M_elementRidges   = MeshType::elementShape_Type::S_numRidges;
    M_facetVertices   = MeshType::facetShape_Type::S_numVertices;

} // init

template < typename MeshType >
void
MeshPartitioner < MeshType >::
doPartition ( meshPtr_Type& mesh, boost::shared_ptr<Epetra_Comm>& comm,
              Epetra_Map* interfaceMap, Epetra_Map* interfaceMapRep )
{
    M_comm = comm;
    M_originalMesh = mesh;
    M_interfaceMap = interfaceMap;
    M_interfaceMapRep = interfaceMapRep;

    M_me = M_comm->MyPID();

    meshPtr_Type newMesh ( new MeshType ( comm ) );
    newMesh->setIsPartitioned ( true );
    M_meshPartitions.reset ( new partMesh_Type ( M_numPartitions, newMesh ) );
    newMesh.reset();

    execute();

} // doPartiton

// =================================
// Public methods
// =================================

template<typename MeshType>
void MeshPartitioner<MeshType>::setup (UInt numPartitions, boost::shared_ptr<Epetra_Comm> comm)
{
    M_serialMode = true;
    M_comm = comm;
    M_me = M_comm->MyPID();

    M_numPartitions = numPartitions;

    M_meshPartitions.reset (new partMesh_Type);
    meshPtr_Type newMesh;
    for (UInt i = 0; i < M_numPartitions; ++i)
    {
        newMesh.reset ( new MeshType ( comm ) );
        newMesh->setIsPartitioned ( true );
        M_meshPartitions->push_back (newMesh);
    }
    newMesh.reset();

    M_elementDomains.reset (new graph_Type);

    M_localNodes.resize (M_numPartitions);
    M_localRidges.resize (M_numPartitions);
    M_localFacets.resize (M_numPartitions);
    M_localElements.resize (M_numPartitions);
    M_globalToLocalNode.resize (M_numPartitions);
    M_globalToLocalElement.resize (M_numPartitions);
    M_nBoundaryPoints.resize (M_numPartitions);
    M_nBoundaryRidges.resize (M_numPartitions);
    M_nBoundaryFacets.resize (M_numPartitions);
}

template<typename MeshType>
void MeshPartitioner<MeshType>::update()
{
    M_numPartitions = M_elementDomains->size();

    Int numElements = 0;

    for (UInt i = 0; i < M_numPartitions; ++i)
    {
        numElements += (*M_elementDomains) [i].size();
    }

    // Rebuild M_graphVertexLocations
    M_graphVertexLocations.resize (numElements);
    for (std::vector<std::vector<Int> >::iterator it1 = M_elementDomains->begin();
            it1 != M_elementDomains->end(); ++it1)
    {
        for (std::vector<Int>::iterator it2 = it1->begin();
                it2 != it1->end(); ++it2)
        {
            M_graphVertexLocations[*it2] = static_cast<Int> ( (it1 - M_elementDomains->begin() ) );
        }
    }
}

template<typename MeshType>
void MeshPartitioner<MeshType>::attachUnpartitionedMesh (meshPtr_Type& mesh,
                                                         Epetra_Map* interfaceMap,
                                                         Epetra_Map* interfaceMapRep)
{
    M_originalMesh = mesh;
    M_interfaceMap = interfaceMap;
    M_interfaceMapRep = interfaceMapRep;
}

template<typename MeshType>
void MeshPartitioner<MeshType>::releaseUnpartitionedMesh()
{
    M_originalMesh.reset();
    M_interfaceMap = 0;
    M_interfaceMapRep = 0;
}

template<typename MeshType>
void MeshPartitioner<MeshType>::doPartitionGraph()
{
    distributeElements (M_originalMesh->numElements() );
    if (M_interfaceMap)
    {
        findRepeatedFacesFSI();
    }
    partitionConnectivityGraph (M_numPartitions);
    if (M_interfaceMap)
    {
        matchFluidPartitionsFSI();
    }
}

template<typename MeshType>
void MeshPartitioner<MeshType>::doPartitionMesh()
{

    // ***********************
    // local mesh construction
    // ***********************
    constructLocalMesh();

    // ******************
    // nodes construction
    // ******************
    constructNodes();

    // ********************
    // element construction
    // ********************
    constructElements();

    // ******************
    // ridges construction
    // ******************
    constructRidges();

    // new faces can be built only after all local volumes are complete in order to get proper ghost faces data
    M_comm->Barrier();

    // ******************
    // faces construction
    // ******************
    constructFacets();

    // ******************
    // final setup
    // ******************
    finalSetup();

    markGhostEntities();
}

template<typename MeshType>
void MeshPartitioner<MeshType>::showMe (std::ostream& output) const
{
    output << "Number of partitions: " << M_numPartitions << std::endl;
    output << "Serial mode:" << M_serialMode << std::endl;
}

// =================================
// Private methods
// =================================

template<typename MeshType>
void MeshPartitioner<MeshType>::distributeElements (UInt numElements)
{
    // ParMETIS is able to work in parallel: how many processors does it have at hand?
    Int numProcessors = M_comm->NumProc();
    M_me              = M_comm->MyPID();

    // CAREFUL: ParMetis works on a graph abstraction.
    // A graph is built over the data structure to be split, each vertex being a mesh element
    // so hereby a "vertex" is actually a _graph_ vertex, i. e. a mesh element
    M_vertexDistribution.resize (numProcessors + 1);
    M_vertexDistribution[0] = 0;

    UInt k = numElements;

    // Evenly distributed graph vertices
    for (Int i = 0; i < numProcessors; ++i)
    {
        UInt l = k / (numProcessors - i);
        M_vertexDistribution[i + 1] = M_vertexDistribution[i] + l;
        k -= l;
    }
    ASSERT (k == 0, "At this point we should have 0 volumes left") ;
}

template<typename MeshType>
void MeshPartitioner<MeshType>::findRepeatedFacesFSI()
{
    std::vector<Int>                    myRepeatedFacet; // used for the solid partitioning
    boost::shared_ptr<std::vector<Int> > myIsOnProc;     // used for the solid partitioning

    myIsOnProc.reset (new std::vector<Int> (M_originalMesh->numElements() ) );

    bool myFacetRep;
    bool myFacet (false);
    short count;
    for (UInt h = 0; h < M_originalMesh->numElements(); ++h)
    {
        (*myIsOnProc) [h] = -1;
    }

    // This loop is throughout the whole unpartitioned mesh,
    // it is expensive and not scalable.
    // Bad, this part should be done offline

    for (UInt ie = 0; ie < M_originalMesh->numElements(); ++ie)
    {
        for (UInt ifacet = 0; ifacet < M_elementFacets; ++ifacet)
        {
            UInt facet = M_originalMesh->localFacetId (ie, ifacet);
            UInt vol  = M_originalMesh->facet (facet).firstAdjacentElementIdentity();
            if (vol == ie)
            {
                vol = M_originalMesh->facet (facet).secondAdjacentElementIdentity();
            }
            if (vol != NotAnId)
            {
                myFacet = false;
                myFacetRep = false;
                count = 0;
                for (Int ipoint = 0; ipoint < static_cast<Int> (M_facetVertices); ++ipoint) // vertex-based dofs
                {
                    myFacetRep = ( (M_interfaceMap->LID (M_originalMesh->facet (facet).point (ipoint).id() )
                                    /* first is fluid */ == -1) &&
                                   (M_interfaceMapRep->LID (M_originalMesh->facet (facet).point (ipoint).id() )
                                    /* first is fluid */ != -1) );
                    myFacet = myFacet ||
                              (M_interfaceMap->LID (M_originalMesh->facet (facet).point (ipoint).id() ) != -1);
                    if (myFacetRep)
                    {
                        ++count;
                    }
                }
                if (count > 1)
                {
                    myRepeatedFacet.push_back (1);
                }
                else
                {
                    myRepeatedFacet.push_back (0);
                }
            }
            if (myFacet)
            {
                (*myIsOnProc) [ie] = M_me;
            }
        }
    }

    M_repeatedFacet.reset (new std::vector<Int> (myRepeatedFacet.size() ) );
    M_isOnProc.reset (new std::vector<Int> (*myIsOnProc) );

    // Lot of communication here!!
    boost::shared_ptr<Epetra_MpiComm> mpiComm = boost::dynamic_pointer_cast <Epetra_MpiComm> ( M_comm );
    MPI_Allreduce ( &myRepeatedFacet[0], & (*M_repeatedFacet) [0], myRepeatedFacet.size(),
                    MPI_INT, MPI_SUM, mpiComm->Comm() );
    MPI_Allreduce ( & (*myIsOnProc) [0], & (*M_isOnProc) [0], myIsOnProc->size(),
                    MPI_INT, MPI_MAX, mpiComm->Comm() );
}

template<typename MeshType>
void MeshPartitioner<MeshType>::partitionConnectivityGraph (UInt numParts)
{
    // This array's size is equal to the number of locally-stored vertices:
    // at the end of the partitioning process, "M_graphVertexLocations" will contain the partitioning array:
    // M_graphVertexLocations[m] = n; means that graph vertex m belongs to subdomain n
    M_graphVertexLocations.resize ( M_vertexDistribution[M_comm->NumProc()] - M_vertexDistribution[0], M_comm->NumProc() );

    // Now each processor will take care of its own graph vertices (i. e. mesh elements).
    // Nothing guarantees about the neighbor elements distribution across the processors,
    // since as of now we just split the set of volumes based on IDs.
    // Here we building up the neighbor arrays.

    UInt localStart = M_vertexDistribution[M_me];
    UInt localEnd   = M_vertexDistribution[M_me + 1];

    // this vector contains the weights for the edges of the graph,
    // it is set to null if it is not used.
    std::vector<Int> graphEdgeWeights;

    M_adjacencyGraphKeys.resize (0);
    M_adjacencyGraphKeys.push_back (0);

    UInt sum = 0;

    for (UInt ie = localStart; ie < localEnd; ++ie)
    {
        for (UInt ifacet = 0; ifacet < M_elementFacets; ++ifacet)
        {
            // global ID of the ifacet-th facet in element ie
            UInt facet = M_originalMesh->localFacetId (ie, ifacet);
            // first adjacent element to face "facet"
            UInt elem = M_originalMesh->facet (facet).firstAdjacentElementIdentity();
            if (elem == ie)
            {
                elem = M_originalMesh->facet (facet).secondAdjacentElementIdentity();
            }
            if (elem != NotAnId)
            {
                // this is the list of adjacency
                // for each graph vertex, simply push back the ID of its neighbors
                M_adjacencyGraphValues.push_back (elem);
                ++sum;
                if (M_interfaceMap) // if I'm partitioning the solid in FSI
                {
                    if ( (*M_repeatedFacet) [sum])
                    {
                        graphEdgeWeights.push_back (0);
                    }
                    else
                    {
                        graphEdgeWeights.push_back (10);
                    }
                }
            }
        }
        // this is the list of "keys" to access M_adjacencyGraphValues
        // graph element i has neighbors M_adjacencyGraphValues[ k ],
        // with M_adjacencyGraphKeys[i] <= k < M_adjacencyGraphKeys[i+1]
        M_adjacencyGraphKeys.push_back (sum);
    }

    // **************
    // parMetis part

    // this array is to be used for weighted vertices on the graph:
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

    Int cutGraphEdges; // here will be stored the number of edges cut in the partitioning process

    // additional options
    std::vector<Int>  options (3, 0);
    options[0] = 1; // means that additional options are actually passed
    options[1] = 3; // level of information to be returned during execution (see ParMETIS's defs.h file)
    options[2] = 1; // random number seed for the ParMETIS routine

    // fraction of vertex weight to be distributed to each subdomain.
    // here we want the subdomains to be of the same size
    std::vector<float> tpwgts (ncon * numParts, 1. / numParts);
    // imbalance tolerance for each vertex weight
    std::vector<float> ubvec (ncon, 1.05);

    boost::shared_ptr<Epetra_MpiComm> mpiComm = boost::dynamic_pointer_cast <Epetra_MpiComm> (M_comm);
    MPI_Comm MPIcomm = mpiComm->Comm();

    Int nprocs;
    MPI_Comm_size (MPIcomm, &nprocs);

    /*
      (from ParMETIS v 3.1 manual)
      This routine is used to compute a k-way partitioning of a graph
      on p processors using the multilevel k-way multi-constraint
      partitioning algorithm.
    */

    Int numberParts = (Int) numParts;

    Int* adjwgtPtr (0);
    if (graphEdgeWeights.size() > 0)
    {
        adjwgtPtr = static_cast<Int*> (&graphEdgeWeights[0]);
    }
    ParMETIS_V3_PartKway (static_cast<Int*> (&M_vertexDistribution[0]),
                          static_cast<Int*> (&M_adjacencyGraphKeys[0]),
                          static_cast<Int*> (&M_adjacencyGraphValues[0]),
                          weightVector, adjwgtPtr, &weightFlag, &numflag,
                          &ncon, &numberParts, &tpwgts[0], &ubvec[0],
                          &options[0], &cutGraphEdges, &M_graphVertexLocations[localStart],
                          &MPIcomm);

    M_comm->Barrier();

    Int nProc = M_comm->NumProc();

    // distribute the resulting partitioning stored in M_graphVertexLocations to all processors
    for ( Int proc = 0; proc < nProc; proc++ )
    {
        UInt procStart  = M_vertexDistribution[ proc ];
        UInt procLength = M_vertexDistribution[ proc + 1 ] - M_vertexDistribution[ proc ];
        M_comm->Broadcast ( &M_graphVertexLocations[ procStart ], procLength, proc );
    }

    // this is a vector of subdomains: each component is
    // the list of vertices belonging to the specific subdomain
    (*M_elementDomains).resize (numParts);

    // cycling on locally stored vertices
    for (UInt ii = 0; ii < M_graphVertexLocations.size(); ++ii)
    {
        // here we are associating the vertex global ID to the subdomain ID
        (*M_elementDomains) [ M_graphVertexLocations[ ii ] ].push_back ( ii );
    }
}

template<typename MeshType>
void MeshPartitioner<MeshType>::matchFluidPartitionsFSI()
{
    boost::shared_ptr<Epetra_MpiComm> mpiComm = boost::dynamic_pointer_cast <Epetra_MpiComm> (M_comm);
    MPI_Comm MPIcomm = mpiComm->Comm();
    Int numProcesses;
    MPI_Comm_size (MPIcomm, &numProcesses);

    std::vector<Int> procOrder (numProcesses);
    std::vector<std::vector<UInt> > myMatchesForProc (numProcesses);
    std::vector<std::vector<UInt> > matchesForProc (numProcesses);
    std::vector<bool> orderingError (numProcesses);

    for (Int i = 0; i < numProcesses ; ++i)
    {
        orderingError[i] = false;
        for (Int j = 0; j < numProcesses ; ++j)
        {
            myMatchesForProc[i].push_back (0);
            matchesForProc[i].push_back (0);
        }
    }

    for (UInt kk = 0; kk < M_graphVertexLocations.size(); ++kk)
    {
        if ( (*M_isOnProc) [kk + M_vertexDistribution[M_me]] != -1)
        {
            ++myMatchesForProc[M_graphVertexLocations[kk]][ (*M_isOnProc) [kk + M_vertexDistribution[M_me]]];
        }
    }

    for (UInt j = 0; (Int) j < numProcesses; ++j)
    {
        MPI_Allreduce (&myMatchesForProc[j][0], &matchesForProc[j][0], numProcesses,
                       MPI_INT, MPI_SUM, MPIcomm);
    }

    M_comm->Barrier();

    Int suitableProcess = -1;
    UInt max = 0;

    for (Int ii = 0; ii < numProcesses; ++ii)
    {
        if (matchesForProc[M_me][ii] > max)
        {
            suitableProcess = ii;
            max = matchesForProc[M_me][ii];
        }
    }

    ASSERT (suitableProcess != -1, "one partition is without interface nodes!");
    procOrder[M_me] = suitableProcess;

    M_comm->Barrier();

    std::vector<UInt> maxs (numProcesses);
    maxs[M_me] = max;
    for (Int j = 0; j < numProcesses ; ++j) // Allgather
    {
        MPI_Bcast (&maxs[j], 1, MPI_INT, j, MPIcomm); // perhaps generates errors
    }

    std::vector<std::pair<UInt, Int> > procIndex (numProcesses);
    for (Int k = 0; k < numProcesses; ++k)
    {
        procIndex[k] = std::make_pair ( maxs[k], k);
    }

    std::sort (procIndex.begin(), procIndex.end() /*, &booleanCondition::reordering*/);

    for (Int l = 0; l < numProcesses; ++l)
    {
        for (Int l = 0; l < numProcesses; ++l)
        {
            for (Int j = 0; j < numProcesses ; ++j) // Allgather
            {
                MPI_Bcast ( &procOrder[j], 1, MPI_INT, j, MPIcomm); // perhaps generates errors
            }
        }
    }

    std::vector< std::vector<Int> > locProc2 ( (*M_elementDomains) );
    for (Int j = numProcesses - 1; j >= 0 ; --j)
    {
        if (orderingError[procOrder[procIndex[j].second]] == false)
        {
            (*M_elementDomains) [procOrder[procIndex[j].second]] = locProc2[procIndex[j].second];
        }
        else
        {
            std::cout << "Ordering error when assigning the processor"
                      << M_me << " to the partition," << std::endl
                      << " parmetis did a bad job." << std::endl;
            for (Int i = numProcesses - 1; i >= 0; --i)
            {
                if (orderingError[procIndex[i].second] == false) // means that i is the first
                    // proc not assigned
                {
                    procOrder[procIndex[j].second] = procIndex[i].second;
                    (*M_elementDomains) [procIndex[i].second] = locProc2[procIndex[j].second];
                    break;
                }
            }
        }
        orderingError[procOrder[procIndex[j].second]] = true;
    }
}

template<typename MeshType>
void MeshPartitioner<MeshType>::redistributeElements()
{
    boost::shared_ptr<Epetra_MpiComm> mpiComm = boost::dynamic_pointer_cast <Epetra_MpiComm> (M_comm);
    MPI_Comm MPIcomm = mpiComm->Comm();
    Int numProcesses;
    MPI_Comm_size (MPIcomm, &numProcesses);

    Int maxInt (1000);
    std::vector<Int> sendSize ( numProcesses );
    std::vector<Int> receiveSize ( numProcesses );
    // cycling on subdomains
    // TODO: Matteo please comment this part :)

    MPI_Status status;
    Int size;

    // MPI_Status  recv_status;
    // MPI_Request send_request;

    for (Int iproc = 0; iproc < numProcesses; ++iproc)
    {
        sendSize[iproc] = (*M_elementDomains) [iproc].size();
    }
    MPI_Alltoall ( &sendSize[ 0 ], 1, MPI_INT, &receiveSize[ 0 ], 1, MPI_INT, MPIcomm );

    for (Int iproc = 0; iproc < numProcesses; ++iproc)
    {
        if (static_cast<Int> (M_me) != iproc)
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

                MPI_Send (&incr, 1, MPI_INT, iproc, 20, MPIcomm);
                MPI_Send (&sizePart, 1, MPI_INT, iproc, 30, MPIcomm);

                for (Int kk = 0; kk < incr; ++kk)
                {
                    MPI_Send (&pos, 1, MPI_INT, iproc, 100 + kk, MPIcomm);
                    MPI_Send (& (*M_elementDomains) [iproc][pos], sizePart, MPI_INT, iproc,
                              5000000 + kk, MPIcomm);
                    pos = pos + sizePart;
                }

                Int resto = size % incr;

                MPI_Send (&resto, 1, MPI_INT, iproc, 80, MPIcomm);

                if (resto != 0)
                {
                    MPI_Send (&pos, 1, MPI_INT, iproc, 40, MPIcomm);
                    MPI_Send (& (*M_elementDomains) [iproc][pos], resto, MPI_INT, iproc, 50, MPIcomm);
                }
            }
            else
            {
                if (size != 0)
                {
                    MPI_Send (& (*M_elementDomains) [iproc][0], size, MPI_INT, iproc, 60, MPIcomm);
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
                    std::vector<Int> stack (size, 0);

                    if (size > maxInt)
                    {
                        Int sizePart, pos, incr;

                        MPI_Recv (&incr, 1, MPI_INT, jproc, 20, MPIcomm, &status);
                        MPI_Recv (&sizePart, 1, MPI_INT, jproc, 30, MPIcomm, &status);

                        for (Int kk = 0; kk < incr; ++kk)
                        {
                            MPI_Recv (&pos, 1, MPI_INT, jproc, 100 + kk, MPIcomm, &status);
                            MPI_Recv (&stack[pos], sizePart , MPI_INT, jproc, 5000000 + kk, MPIcomm, &status);
                        }
                        Int resto = 0;
                        MPI_Recv (&resto, 1, MPI_INT, jproc, 80, MPIcomm, &status);

                        if (resto != 0)
                        {
                            MPI_Recv (&pos, 1, MPI_INT, jproc, 40, MPIcomm, &status);
                            MPI_Recv (&stack[pos],  resto, MPI_INT, jproc, 50, MPIcomm, &status);
                        }
                    }
                    else
                    {
                        if (size != 0)
                        {
                            MPI_Recv (&stack[0], size , MPI_INT, jproc, 60, MPIcomm, &status);
                        }
                    }
                    for (Int jj = 0; jj < size; ++jj)
                    {
                        (*M_elementDomains) [M_me].push_back (stack[jj]);
                    }
                }
            }
        }
    }
}

template<typename MeshType>
void MeshPartitioner<MeshType>::constructLocalMesh()
{
    if (!M_me)
    {
        std::cout << "Building local mesh ..." << std::endl;
    }

    std::map<Int, Int>::iterator  im;

    Int count = 0;
    UInt ielem;
    UInt inode;

    for (UInt i = 0; i < M_numPartitions; ++i)
    {
        count = 0;
        // cycle on local element's ID

        UInt me = M_serialMode ? i : M_me;

        for (UInt jj = 0; jj < (*M_elementDomains) [me].size(); ++jj)
        {
            ielem = (*M_elementDomains) [me][jj];
            M_localElements[i].push_back (ielem);

            // cycle on element's vertices
            for (UInt ii = 0; ii < M_elementVertices; ++ii)
            {
                inode = M_originalMesh->element (ielem).point (ii).id();
                im    = M_globalToLocalNode[i].find (inode);

                // if the node is not yet present in the list of local vertices, then add it
                if (im == M_globalToLocalNode[i].end() )
                {
                    M_globalToLocalNode[i].insert (std::make_pair (inode, count) );
                    ++count;
                    // store here the global numbering of the node
                    M_localNodes[i].push_back (M_originalMesh->element (ielem).point (ii).id() );
                }
            }

            // cycle on element's ridges
            for (UInt ii = 0; ii < M_elementRidges; ++ii)
            {
                // store here the global numbering of the ridge
                M_localRidges[i].insert (M_originalMesh->localRidgeId (ielem, ii) );
            }

            // cycle on element's facets
            for (UInt ii = 0; ii < M_elementFacets; ++ii)
            {
                // store here the global numbering of the facet
                M_localFacets[i].insert (M_originalMesh->localFacetId (ielem, ii) );
            }
        }
    }
}

template<typename MeshType>
void MeshPartitioner<MeshType>::constructNodes()
{
    UInt inode;
    for (UInt i = 0; i < M_numPartitions; ++i)
    {
        std::vector<Int>::iterator it;

        M_nBoundaryPoints[i] = 0;
        (*M_meshPartitions) [i]->pointList.reserve (M_localNodes[i].size() );
        // guessing how many boundary points on this processor.
        (*M_meshPartitions) [i]->_bPoints.reserve (M_originalMesh->numBPoints() * M_localNodes[i].size() /
                                                   M_originalMesh->numBPoints() );

        inode = 0;

        typename MeshType::point_Type* pp = 0;

        // loop in the list of local vertices:
        // in this loop inode is the local numbering of the points
        for (it = M_localNodes[i].begin(); it != M_localNodes[i].end(); ++it, ++inode)
        {
            // create a boundary point in the local mesh, if needed
            bool boundary = M_originalMesh->isBoundaryPoint (*it);
            M_nBoundaryPoints[i] += boundary;

            pp = & (*M_meshPartitions) [i]->addPoint ( boundary, false );
            *pp = M_originalMesh->point ( *it );

            pp->setLocalId ( inode );
        }
    }
}


template<typename MeshType>
void MeshPartitioner<MeshType>::constructElements()
{
    Int count;
    for (UInt i = 0; i < M_numPartitions; ++i)
    {
        std::map<Int, Int>::iterator im;
        std::vector<Int>::iterator it;
        count = 0;
        UInt inode;

        typename MeshType::element_Type* pv = 0;

        (*M_meshPartitions) [i]->elementList().reserve (M_localElements[i].size() );

        // loop in the list of local elements
        // CAREFUL! in this loop inode is the global numbering of the points
        // We insert the local numbering of the vertices in the local volume list
        for (it = M_localElements[i].begin(); it != M_localElements[i].end(); ++it, ++count)
        {
            pv = & ( (*M_meshPartitions) [i]->addElement() );
            *pv = M_originalMesh->element (*it);
            pv->setLocalId (count);

            M_globalToLocalElement[i].insert (std::make_pair ( pv->id(), pv -> localId() ) );

            for (ID id = 0; id < M_elementVertices; ++id)
            {
                inode = M_originalMesh->element (*it).point (id).id();
                // im is an iterator to a map element
                // im->first is the key (i. e. the global ID "inode")
                // im->second is the value (i. e. the local ID "count")
                im = M_globalToLocalNode[i].find (inode);
                pv->setPoint (id, (*M_meshPartitions) [i]->point ( (*im).second ) );
            }
        }
    }
}

template<typename MeshType>
void MeshPartitioner<MeshType>::constructRidges()
{
    if (mesh_Type::S_geoDimensions == 2)
    {
        M_nBoundaryRidges = M_nBoundaryPoints;
    }
    else if (mesh_Type::S_geoDimensions == 3)
    {
        Int count;
        for (UInt i = 0; i < M_numPartitions; ++i)
        {
            std::map<Int, Int>::iterator im;
            std::set<Int>::iterator is;

            typename MeshType::ridge_Type* pe;
            UInt inode;
            count = 0;

            M_nBoundaryRidges[i] = 0;
            (*M_meshPartitions) [i]->ridgeList().reserve (M_localRidges[i].size() );

            // loop in the list of local ridges
            for (is = M_localRidges[i].begin(); is != M_localRidges[i].end(); ++is, ++count)
            {
                // create a boundary ridge in the local mesh, if needed
                bool boundary = (M_originalMesh->isBoundaryRidge (*is) );

                // create a boundary ridge in the local mesh, if needed
                M_nBoundaryRidges[i] += boundary;

                pe = & (*M_meshPartitions) [i]->addRidge (boundary);
                *pe = M_originalMesh->ridge ( *is );

                pe->setLocalId (count);

                for (ID id = 0; id < 2; ++id)
                {
                    inode = M_originalMesh->ridge (*is).point (id).id();
                    // im is an iterator to a map element
                    // im->first is the key (i. e. the global ID "inode")
                    // im->second is the value (i. e. the local ID "count")
                    im = M_globalToLocalNode[i].find (inode);
                    pe->setPoint (id, (*M_meshPartitions) [i]->pointList ( (*im).second) );
                }
            }
        }
    }
}


template<typename MeshType>
void MeshPartitioner<MeshType>::constructFacets()
{
    Int count;
    for (UInt i = 0; i < M_numPartitions; ++i)
    {
        std::map<Int, Int>::iterator im;
        std::set<Int>::iterator      is;

        typename MeshType::facet_Type* pf = 0;

        UInt inode;
        count = 0;

        M_nBoundaryFacets[i] = 0;
        (*M_meshPartitions) [i]->facetList().reserve (M_localFacets[i].size() );

        // loop in the list of local faces
        for (is = M_localFacets[i].begin(); is != M_localFacets[i].end(); ++is, ++count)
        {
            // create a boundary facet in the local mesh, if needed
            bool boundary = (M_originalMesh->isBoundaryFacet (*is) );

            M_nBoundaryFacets[i] += boundary;


            Int elem1 = M_originalMesh->facet (*is).firstAdjacentElementIdentity();
            Int elem2 = M_originalMesh->facet (*is).secondAdjacentElementIdentity();

            // find the mesh elements adjacent to the facet
            im =  M_globalToLocalElement[i].find (elem1);

            ID localElem1;

            if (im == M_globalToLocalElement[i].end() )
            {
                localElem1 = NotAnId;
            }
            else
            {
                localElem1 = (*im).second;
            }

            im =  M_globalToLocalElement[i].find (elem2);

            ID localElem2;
            if (im == M_globalToLocalElement[i].end() )
            {
                localElem2 = NotAnId;
            }
            else
            {
                localElem2 = (*im).second;
            }

            pf =  & (*M_meshPartitions) [i]->addFacet (boundary);
            *pf = M_originalMesh->facet ( *is );

            pf->setLocalId ( count );

            for (ID id = 0; id < M_originalMesh->facet (*is).S_numLocalVertices; ++id)
            {
                inode = pf->point (id).id();
                im = M_globalToLocalNode[i].find (inode);
                pf->setPoint (id, (*M_meshPartitions) [i]->pointList ( (*im).second) );
            }

            // true if we are on a subdomain border
            ID ghostElem = ( localElem1 == NotAnId ) ? elem1 : elem2;

            if ( !boundary && ( localElem1 == NotAnId || localElem2 == NotAnId ) )
            {
                // set the flag for faces on the subdomain border
                pf->setFlag ( EntityFlags::SUBDOMAIN_INTERFACE );
                // set the flag for all points on that face
                for ( UInt pointOnFacet = 0; pointOnFacet < MeshType::facet_Type::S_numLocalPoints; pointOnFacet++ )
                {
                    (*M_meshPartitions) [i]->point ( pf->point ( pointOnFacet ).localId() ).setFlag ( EntityFlags::SUBDOMAIN_INTERFACE );
                }

            }

            // if this process does not own either of the adjacent elements
            // then the two adjacent elements and the respective facet positions coincide in the local mesh
            // possible bug fixed: not only the two adjacent elements facet, but also the facet
            // positions should coincide.
            // otherwise it could happen that a pair(element, position) is associated to different faces.
            // This can lead to a wrong treatment of the dofPerFace (in 2D of the dofPerRidge, as occurred
            // with P2)

            ASSERT ( (localElem1 != NotAnId) || (localElem2 != NotAnId), "A hanging facet in mesh partitioner!");

            if ( localElem1 == NotAnId )
            {
                pf->firstAdjacentElementIdentity()  = localElem2;
                pf->firstAdjacentElementPosition()  = M_originalMesh->facet (*is).secondAdjacentElementPosition();
                pf->secondAdjacentElementIdentity() = ghostElem;
                pf->secondAdjacentElementPosition() = NotAnId;
                pf->reversePoints();
            }
            else if ( localElem2 == NotAnId )
            {
                pf->firstAdjacentElementIdentity()  = localElem1;
                pf->firstAdjacentElementPosition()  = M_originalMesh->facet (*is).firstAdjacentElementPosition();
                pf->secondAdjacentElementIdentity() = ghostElem;
                pf->secondAdjacentElementPosition() = NotAnId;
            }
            else
            {
                pf->firstAdjacentElementIdentity()  = localElem1;
                pf->firstAdjacentElementPosition()  = M_originalMesh->facet (*is).firstAdjacentElementPosition();
                pf->secondAdjacentElementIdentity() = localElem2;
                pf->secondAdjacentElementPosition() = M_originalMesh->facet (*is).secondAdjacentElementPosition();
            }

        }
        (*M_meshPartitions) [i]->setLinkSwitch ("HAS_ALL_FACETS");
        (*M_meshPartitions) [i]->setLinkSwitch ("FACETS_HAVE_ADIACENCY");
    }

}

template<typename MeshType>
void MeshPartitioner<MeshType>::finalSetup()
{
    for (UInt i = 0; i < M_numPartitions; ++i)
    {
        UInt nElements = M_localElements[i].size();
        UInt nNodes   = M_localNodes[i].size();
        UInt nRidges   = M_localRidges[i].size();
        UInt nFacets   = M_localFacets[i].size();

        (*M_meshPartitions) [i]->setMaxNumPoints (nNodes, true);
        (*M_meshPartitions) [i]->setMaxNumRidges  (nRidges, true);
        (*M_meshPartitions) [i]->setMaxNumFacets  (nFacets, true);
        (*M_meshPartitions) [i]->setMaxNumElements ( nElements, true);

        (*M_meshPartitions) [i]->setMaxNumGlobalPoints (M_originalMesh->numPoints() );
        (*M_meshPartitions) [i]->setNumGlobalVertices  (M_originalMesh->numVertices() );
        (*M_meshPartitions) [i]->setMaxNumGlobalRidges  (M_originalMesh->numRidges() );
        (*M_meshPartitions) [i]->setMaxNumGlobalFacets  (M_originalMesh->numFacets() );

        (*M_meshPartitions) [i]->setMaxNumGlobalElements (M_originalMesh->numElements() );
        (*M_meshPartitions) [i]->setNumBoundaryFacets    (M_nBoundaryFacets[i]);

        (*M_meshPartitions) [i]->setNumBPoints   (M_nBoundaryPoints[i]);
        (*M_meshPartitions) [i]->setNumBoundaryRidges    (M_nBoundaryRidges[i]);

        (*M_meshPartitions) [i]->setNumVertices (nNodes );
        (*M_meshPartitions) [i]->setNumBVertices (M_nBoundaryPoints[i]);

        if (MeshType::S_geoDimensions == 3)
        {
            (*M_meshPartitions) [i]->updateElementRidges();
        }

        (*M_meshPartitions) [i]->updateElementFacets();

#ifdef HAVE_LIFEV_DEBUG
        if (M_serialMode)
        {
            debugStream (4000) << "Created local mesh number " << i << "\n";
        }
        else
        {
            debugStream (4000) << "Rank " << M_me << " created local mesh.\n";
        }
#endif
    }
}

template<typename MeshType>
void MeshPartitioner<MeshType>::execute()
{

    // Build graph vertex distribution vector. Graph vertex represents one element
    // in the mesh.
    distributeElements (M_originalMesh->numElements() );


    // In fluid-structure interaction:
    // *    If the solid mesh is not partitioned the following part won't be
    //      executed
    // *    If the solid mesh is partitioned:
    //      - The fluid is partitioned first
    //      - The solid mesh partition tries to follow the partition of the fluid
    //      This is achieved by specifying a weight to some ridge of the graph.
    //      The interface between two processors is the set of the vertices that for
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
    partitionConnectivityGraph (M_comm->NumProc() );

    //////////////// BEGIN OF SOLID PARTITION PART ////////////////
    if (M_interfaceMap)
    {
        matchFluidPartitionsFSI();
    }
    ////////////////// END OF SOLID PARTITION PART /////////////////////

#ifdef HAVE_LIFEV_DEBUG
    debugStream (4000) << M_me << " has " << (*M_elementDomains) [M_me].size() << " elements.\n";
#endif

    fillEntityPID ();
    if ( M_partitionOverlap > 0 )
    {
        GhostHandler<mesh_Type> gh ( M_originalMesh, M_comm );
        gh.extendGraphFE ( M_elementDomains, M_entityPID.points, M_partitionOverlap );
    }

    doPartitionMesh();

    // ***************************************
    // release the original unpartitioned mesh
    // allowing it to be deleted
    // ***************************************
    releaseUnpartitionedMesh();

    // ***************************************
    // clear all internal structures that are
    // not needed anymore
    // ***************************************
    cleanUp();
}

template<typename MeshType>
void MeshPartitioner<MeshType>::fillEntityPID ()
{
    Int numParts = (M_numPartitions > 1) ? M_numPartitions : M_comm->NumProc();

    // initialize entity PIDs to 0
    M_entityPID.points.resize   ( M_originalMesh->numPoints(),   numParts-1 );
    M_entityPID.elements.resize ( M_originalMesh->numElements(), numParts-1 );
    M_entityPID.facets.resize   ( M_originalMesh->numFacets(),   numParts-1 );
    M_entityPID.ridges.resize   ( M_originalMesh->numRidges(),   numParts-1 );

    // check: parallel algorithm seems to be slower for this
    // p = 0 can be skipped since M_entityPID is already initialized at that value
    for ( Int p = 0; p < numParts; p++ )
    {
        for ( UInt e = 0; e < (*M_elementDomains) [ p ].size(); e++ )
        {
            // point block
            for ( UInt k = 0; k < mesh_Type::element_Type::S_numPoints; k++ )
            {
                const ID& pointID = M_originalMesh->element ( (*M_elementDomains) [ p ][ e ] ).point ( k ).id();
                // pointPID should be the maximum between the procs that own it
                M_entityPID.points[ pointID ] = std::min ( M_entityPID.points[ pointID ], p );
            }

            // elem block
            const ID& elemID = M_originalMesh->element ( (*M_elementDomains) [ p ][ e ] ).id();
            // at his stage each element belongs to a single partition, overlap is not yet done.
            M_entityPID.elements[ elemID ] = p;

            // facet block
            for ( UInt k = 0; k < mesh_Type::element_Type::S_numFacets; k++ )
            {
                const ID& facetID = M_originalMesh->facet ( M_originalMesh->localFacetId ( elemID, k ) ).id();
                // facetPID should be the maximum between the proc that own it
                M_entityPID.facets[ facetID ] = std::min ( M_entityPID.facets[ facetID ], p );
            }

            // ridge block
            for ( UInt k = 0; k < mesh_Type::element_Type::S_numRidges; k++ )
            {
                const ID& ridgeID = M_originalMesh->ridge ( M_originalMesh->localRidgeId ( elemID, k ) ).id();
                // ridgePID should be the maximum between the proc that own it
                M_entityPID.ridges[ ridgeID ] = std::min ( M_entityPID.ridges[ ridgeID ], p );
            }
        }
    }
}

template<typename MeshType>
void MeshPartitioner<MeshType>::markGhostEntities()
{
    // mark ghost entities by each partition as described in M_entityPID
    //@todo: does not work for offline partitioning!
    //M_entityPID or flags should be exported and read back to make it work
    for (UInt i = 0; i < M_numPartitions; ++i)
    {
        Int const procId = (M_numPartitions > 1) ? i : M_me;
        for ( UInt e = 0; e < (*M_meshPartitions) [ i ]->numElements(); e++ )
        {
            typename MeshType::element_Type& element = (*M_meshPartitions) [ i ]->element ( e );
            if ( M_entityPID.elements[ element.id() ] != static_cast<UInt> ( procId ) )
            {
                element.setFlag ( EntityFlags::GHOST );
            }
        }

        for ( UInt f = 0; f < (*M_meshPartitions) [ i ]->numFacets(); f++ )
        {
            typename MeshType::facet_Type& facet = (*M_meshPartitions) [ i ]->facet ( f );
            if ( M_entityPID.facets[ facet.id() ] != static_cast<UInt> ( procId ) )
            {
                facet.setFlag ( EntityFlags::GHOST );
            }
        }

        for ( UInt r = 0; r < (*M_meshPartitions) [ i ]->numRidges(); r++ )
        {
            typename MeshType::ridge_Type& ridge = (*M_meshPartitions) [ i ]->ridge ( r );
            if ( M_entityPID.ridges[ ridge.id() ] != static_cast<UInt> ( procId ) )
            {
                ridge.setFlag ( EntityFlags::GHOST );
            }
        }

        for ( UInt p = 0; p < (*M_meshPartitions) [ i ]->numPoints(); p++ )
        {
            typename MeshType::point_Type& point = (*M_meshPartitions) [ i ]->point ( p );
            if ( M_entityPID.points[ point.id() ] != static_cast<UInt> ( procId ) )
            {
                point.setFlag ( EntityFlags::GHOST );
            }
        }
    }
    clearVector ( M_entityPID.elements );
    clearVector ( M_entityPID.facets );
    clearVector ( M_entityPID.ridges );
    clearVector ( M_entityPID.points );
}

template<typename MeshType>
void MeshPartitioner<MeshType>::cleanUp()
{
    clearVector ( M_vertexDistribution );
    clearVector ( M_adjacencyGraphKeys );
    clearVector ( M_adjacencyGraphValues );
    clearVector ( M_localNodes );
    clearVector ( M_localRidges );
    clearVector ( M_localFacets );
    clearVector ( M_localElements );
    clearVector ( M_nBoundaryPoints );
    clearVector ( M_nBoundaryRidges );
    clearVector ( M_nBoundaryFacets );
    clearVector ( M_graphVertexLocations );
    clearVector ( M_globalToLocalNode );
    clearVector ( M_globalToLocalElement );
}

}
#endif // MESH_PARTITIONER_H
