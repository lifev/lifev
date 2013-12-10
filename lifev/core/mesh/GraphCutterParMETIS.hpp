//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010, 2011, 2012 EPFL, Politecnico di Milano, Emory University

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
    @brief Class that partitions the graph associated with a mesh.
           Uses the ParMETIS graph processing library

    @date 28-11-2012
    @author Radu Popescu <radu.popescu@epfl.ch>
 */

#ifndef GRAPH_PARTITION_TOOL_PARMETIS_H
#define GRAPH_PARTITION_TOOL_PARMETIS_H 1

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <algorithm>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/bimap.hpp>
#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>
#include <parmetis.h>


#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/GraphCutterBase.hpp>
#include <lifev/core/mesh/GraphUtil.hpp>

namespace LifeV
{

//! Class that partitions the graph associated with a mesh (ParMETIS version)
/*!
    @author Radu Popescu <radu.popescu@epfl.ch>

    This class uses the ParMETIS package to partition the graph associated
    with a mesh. This class builds the dual graph of the mesh, partitions
    it according to a set of parameters and the stores the partition
    in a table (vector of vectors).
    At the end of the partition process, each vector will contain the
    GID of the elements in a part.

    While this class can be used stand-alone, it is used automatically by the
    MeshPartitionTool class during the mesh partition process.
 */
template<typename MeshType>
class GraphCutterParMETIS : public GraphCutterBase<MeshType>
{
public:
    //! @name Public Types
    //@{
    typedef Teuchos::ParameterList                 pList_Type;
    typedef boost::shared_ptr<Epetra_Comm>         commPtr_Type;
    typedef MeshType                               mesh_Type;
    typedef boost::shared_ptr<mesh_Type>           meshPtr_Type;

    typedef std::vector<Int>                       idList_Type;
    typedef boost::shared_ptr<idList_Type>         idListPtr_Type;
    typedef std::vector<idListPtr_Type>            vertexPartition_Type;
    typedef boost::shared_ptr<vertexPartition_Type> vertexPartitionPtr_Type;

    typedef boost::bimap<UInt, UInt>               biMap_Type;
    typedef biMap_Type::value_type                 biMapValue_Type;
    //@}

    //! @name Constructor & Destructor
    //! Constructor taking the original mesh, the MPI comm and parameters
    /*!
        This constructor can be used to build the object and perform the
        graph partitioning in one shot.

        @param mesh The original mesh whose graph this object will partition
        @param comm The Epetra MPI comm object which contains the processes
                    which participate
        @param parameters The Teuchos parameter list which contains the
                          partitioning parameters
     */
    GraphCutterParMETIS (meshPtr_Type& mesh,
                         commPtr_Type& comm,
                         pList_Type& parameters);

    //! Destructor
    virtual ~GraphCutterParMETIS() {}
    //@}

    //! @name Public methods
    //@{
    //! Performs the graph partitioning
    virtual Int run();
    //@}

    //! @name Get Methods
    //@{
    //! Get a pointer to one of the partitions
    virtual const idListPtr_Type& getPart (const UInt i) const
    {
        return M_vertexPartition->at(i);
    }
    virtual idListPtr_Type& getPart (const UInt i)
    {
        return M_vertexPartition->at(i);
    }

    //! Return the number of parts
    virtual const UInt numParts() const
    {
        return M_vertexPartition->size();
    }

    //! Get the entire partitioned graph, wrapped in a smart pointer
    virtual const vertexPartitionPtr_Type getGraph() const
    {
        vertexPartitionPtr_Type graph (new vertexPartition_Type(M_vertexPartition->size()));

        for (UInt i = 0; i < M_vertexPartition->size(); ++i)
        {
            graph->at(i) = getPart (i);
        }

        return graph;
    }

    //@}

private:
    //! @name Private methods
    //@{
    //! Set values for all the parameters, with default values where needed
    virtual void setParameters (pList_Type& parameters);

    //! Perform a flat, non-hierarchical partition
    Int partitionFlat();

    //! Perform a hierarchical (2-level) partition
    Int partitionHierarchical();

    //! Perform a partitioning on a given subset of elements
    Int partitionSubGraph (const biMap_Type& vertexMap,
                           const Int numParts,
                           vertexPartitionPtr_Type& vertexPartition);
    //@}

    // Private copy constructor and assignment operator are disabled
    GraphCutterParMETIS (const GraphCutterParMETIS&);
    GraphCutterParMETIS& operator= (const GraphCutterParMETIS&);

    //! @name Private Methods
    //@{
    //@}

    commPtr_Type                                       M_comm;
    Int                                                M_myPID;
    Int                                                M_numProcessors;
    Int                                                M_numParts;
    bool                                               M_hierarchical;
    Int                                                M_topology;
    pList_Type                                         M_parameters;
    meshPtr_Type                                       M_mesh;
    UInt                                               M_elementVertices;
    UInt                                               M_elementFacets;
    UInt                                               M_elementRidges;
    UInt                                               M_facetVertices;
    vertexPartitionPtr_Type                            M_vertexPartition;
};

//
// IMPLEMENTATION
//

// =================================
// Constructors and Destructor
// =================================

template<typename MeshType>
GraphCutterParMETIS<MeshType>::GraphCutterParMETIS (meshPtr_Type& mesh,
                                                    commPtr_Type& comm,
                                                    pList_Type& parameters) :
    M_comm (comm),
    M_myPID (M_comm->MyPID() ),
    M_numProcessors (M_comm->NumProc() ),
    M_numParts (0),
    M_parameters(),
    M_mesh (mesh),
    M_vertexPartition(new vertexPartition_Type(0))
{
    setParameters (parameters);
}

template<typename MeshType>
void GraphCutterParMETIS<MeshType>::setParameters (pList_Type& parameters)
{
    // Here put some default values for the parameters and then import
    // the user supplied list, overwriting the corresponding parameters
    M_parameters.set ("num-parts", static_cast<Int> (M_comm->NumProc() ),
                      "The desired number of parts");
    M_parameters.set<bool> ("hierarchical", false);
    M_parameters.set ("topology", "1",
                      "The topology of the mesh partition process.");

    M_parameters.setParameters (parameters);

    M_hierarchical = M_parameters.get<bool> ("hierarchical");
    M_topology = boost::lexical_cast<Int> (
                     M_parameters.get<std::string> ("topology") );
    M_numParts = M_parameters.get<Int> ("num-parts");
}

template<typename MeshType>
Int GraphCutterParMETIS<MeshType>::run()
{
    // Initialization
    /*
     * Sets element parameters (nodes, faces, ridges and number of nodes on each
     * facet according to the type of mesh element used
     * (Mesh::ElementShape::S_shape).
     * Updates M_elementVertices, M_elementFaces, M_elementRidges,
     * M_facetVertices.
     */
    M_elementVertices = MeshType::elementShape_Type::S_numVertices;
    M_elementFacets   = MeshType::elementShape_Type::S_numFacets;
    M_elementRidges   = MeshType::elementShape_Type::S_numRidges;
    M_facetVertices   = MeshType::facetShape_Type::S_numVertices;

    Int error;
    if ( (! M_hierarchical) || (M_topology == 1) )
    {
        error = partitionFlat();
    }
    else
    {
        error = partitionHierarchical();
    }

    return error;
}

template<typename MeshType>
Int GraphCutterParMETIS<MeshType>::partitionFlat()
{
    // In this case we want to partition the entire graph
    Int numVertices = M_mesh->numElements();

    // We need to build a bidirectional map between local and global IDs
    // vertexMap.left is the local-to-global map and vertexMap.right is
    // the global-to-local map
    biMap_Type vertexMap;
    for (UInt i = 0; i < numVertices; ++i)
    {
        vertexMap.insert (biMapValue_Type (i, i) );
    }

    // Call the partitionSubGraph method on the vertexList that was
    // prepared
    partitionSubGraph (vertexMap, M_numParts, M_vertexPartition);

    return 0;
}

template<typename MeshType>
Int GraphCutterParMETIS<MeshType>::partitionHierarchical()
{
    if (M_numProcessors != 1)
    {
        if (! M_myPID)
        {
            std::cout << "Hierarchical partitioning can only be performed with "
                      << "one MPI process." << std::endl;
        }
        return 1;
    }
    if (M_numParts % M_topology != 0)
    {
        if (! M_myPID)
        {
            std::cout << "Hierarchical partitioning can only be performed if "
                      << "the number of mesh parts is a multiple of the"
                      << " topology parameter." << std::endl;
        }
        return 2;
    }

    // First step is to partition the graph into the number of subdomains
    Int numSubdomains = M_numParts / M_topology;
    Int numVertices = M_mesh->numElements();

    // The vector contains the global IDs of the vertices in the graph
    biMap_Type vertexMap;
    for (Int i = 0; i < numVertices; ++i)
    {
        vertexMap.insert(biMapValue_Type(i, i));
    }

    vertexPartitionPtr_Type tempVertexPartition(new vertexPartition_Type(numSubdomains));

    /*
     * After calling partitionSubGraph, tempVertexPartition will contain
     * numSubdomains vectors with the graph vertices of each subdomain
     */
    partitionSubGraph(vertexMap, numSubdomains, tempVertexPartition);

    /*
     * Step two is to partition each subdomain into the number of sub parts
     * denoted by the M_topology parameter
     */
    M_vertexPartition->resize (M_numParts);
    Int currentPart = 0;
    for (Int i = 0; i < numSubdomains; ++i)
    {
        biMap_Type subdomainVertexMap;
        for (Int k = 0; k < tempVertexPartition->at(i)->size(); ++k) {
            subdomainVertexMap.insert(biMapValue_Type(k, tempVertexPartition->at(i)->at(k)));
        }
        vertexPartitionPtr_Type subdomainParts(new vertexPartition_Type(M_numParts));
              partitionSubGraph(subdomainVertexMap, M_topology,
                                subdomainParts);
        for (Int j = 0; j < M_topology; ++j)
        {
            M_vertexPartition->at(currentPart++) = subdomainParts->at(j);
        }
    }

    return 0;
}

template<typename MeshType>
Int GraphCutterParMETIS<MeshType>::partitionSubGraph (
    const biMap_Type& vertexMap,
    const Int numParts,
    vertexPartitionPtr_Type& vertexPartition)
{
    // Distribute elements
    UInt k = vertexMap.size();

    // CAREFUL: ParMetis works on a graph abstraction.
    // A graph is built over the data structure to be split, each vertex being
    // a mesh element so hereby a "vertex" is actually a _graph_ vertex,
    // i. e. a mesh element
    std::vector<Int> vertexDistribution (M_numProcessors + 1);
    vertexDistribution[0] = 0;
    // Evenly distributed graph vertices
    for (Int i = 0; i < M_numProcessors; ++i)
    {
        UInt l = k / (M_numProcessors - i);
        vertexDistribution[i + 1] = vertexDistribution[i] + l;
        k -= l;
    }
    ASSERT (k == 0, "At this point we should have 0 volumes left") ;


    /*
     * Partition Graph
     * This array's size is equal to the number of locally-stored vertices:
     * at the end of the partitioning process, "M_graphVertexLocations" will
     * contain the partitioning array:
     * M_graphVertexLocations[m] = n; means that graph vertex m belongs to
     * subdomain n
     */
    std::vector<Int> graphVertexLocations (
        vertexDistribution[M_numProcessors]
        - vertexDistribution[0], M_numProcessors);

    /*
     * Now each processor will take care of its own graph vertices
     * (i. e. mesh elements).
     * Nothing guarantees about the neighbor elements distribution across
     * the processors, since as of now we just split the set of volumes based
     * on IDs.
     * Here we building up the neighbor arrays.
     */
    UInt localStart = vertexDistribution[M_myPID];
    UInt localEnd   = vertexDistribution[M_myPID + 1];

    // this vector contains the weights for the edges of the graph,
    // it is set to null if it is not used.
    std::vector<Int> graphEdgeWeights;
    std::vector<Int> adjacencyGraphKeys (1, 0);
    std::vector<Int> adjacencyGraphValues (0);

    UInt sum = 0;

    for (UInt lid = localStart; lid < localEnd; ++lid)
    {
        for (UInt ifacet = 0; ifacet < M_elementFacets; ++ifacet)
        {
            UInt gid = vertexMap.left.at (lid);
            // global ID of the ifacet-th facet in element ie
            UInt facet = M_mesh->localFacetId (gid, ifacet);
            // first adjacent element to face "facet"
            UInt elem = M_mesh->facet (facet).firstAdjacentElementIdentity();
            if (elem == gid)
            {
                elem = M_mesh->facet (facet).secondAdjacentElementIdentity();
            }
            biMap_Type::right_const_iterator it = vertexMap.right.find (elem);

            bool inSubGraph = (vertexMap.right.end() != it);
            if ( (elem != NotAnId) && (inSubGraph) )
            {
                // this is the list of adjacency
                // for each graph vertex, push back the ID of its neighbors
                //adjacencyGraphValues.push_back(it - vertexMap.right.begin());
                adjacencyGraphValues.push_back (it->second);
                ++sum;
            }
        }
        // this is the list of "keys" to access M_adjacencyGraphValues
        // graph element i has neighbors M_adjacencyGraphValues[ k ],
        // with M_adjacencyGraphKeys[i] <= k < M_adjacencyGraphKeys[i+1]
        adjacencyGraphKeys.push_back (sum);
    }

    // **************
    // parMetis part

    // this array is to be used for weighted vertices on the graph:
    // usually we will set it to NULL

    Int* weightVector = 0;

    Int weightFlag = 0;

    Int ncon = 1;
    Int numflag = 0;

    Int cutGraphEdges; // here will be stored the number of edges cut in the
    // partitioning process

    // additional options
    std::vector<Int>  options (3, 0);
    options[0] = 1; // means that additional options are actually passed
    options[1] = 3; // level of information to be returned during execution
    options[2] = 1; // random number seed for the ParMETIS routine

    // fraction of vertex weight to be distributed to each subdomain.
    // here we want the subdomains to be of the same size
    std::vector<float> tpwgts (ncon * numParts, 1. / numParts);
    // imbalance tolerance for each vertex weight
    std::vector<float> ubvec (ncon, 1.05);

    boost::shared_ptr<Epetra_MpiComm> mpiComm
        = boost::dynamic_pointer_cast <Epetra_MpiComm> (M_comm);
    MPI_Comm MPIcomm = mpiComm->Comm();

    /*
      (from ParMETIS v 3.1 manual)
      This routine is used to compute a k-way partitioning of a graph
      on p processors using the multilevel k-way multi-constraint
      partitioning algorithm.
    */

    Int nParts = numParts;
    Int* adjwgtPtr (0);
    if (graphEdgeWeights.size() > 0)
    {
        adjwgtPtr = static_cast<Int*> (&graphEdgeWeights[0]);
    }
    ParMETIS_V3_PartKway (static_cast<Int*> (&vertexDistribution[0]),
                          static_cast<Int*> (&adjacencyGraphKeys[0]),
                          static_cast<Int*> (&adjacencyGraphValues[0]),
                          weightVector, adjwgtPtr, &weightFlag, &numflag,
                          &ncon, &nParts, &tpwgts[0], &ubvec[0],
                          &options[0], &cutGraphEdges,
                          &graphVertexLocations[localStart],
                          &MPIcomm);

    // distribute the resulting partitioning stored in M_graphVertexLocations
    // to all processors
    if (M_numProcessors != 1)
    {
        M_comm->Barrier();
        for (Int proc = 0; proc < M_numProcessors; proc++)
        {
            UInt procStart  = vertexDistribution[proc];
            UInt procLength = vertexDistribution[proc + 1]
                              - vertexDistribution[proc];
            M_comm->Broadcast (&graphVertexLocations[procStart],
                               procLength, proc);
        }
    }

    // cycling on locally stored vertices
    vertexPartition->resize (numParts);
    for (UInt i = 0; i < numParts; ++i)
    {
        vertexPartition->at(i).reset(new idList_Type(0));
    }
    for (UInt ii = 0; ii < graphVertexLocations.size(); ++ii)
    {
        // here we are associating the vertex global ID to the subdomain ID
        vertexPartition->at(graphVertexLocations[ii])->push_back (vertexMap.left.at (ii) );
    }

    return 0;
}

} // Namespace LifeV

#endif // GRAPH_PARTITION_TOOL_PARMETIS_H
