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
    @brief Utilitary functions for graph partitioning

    @date 2013-12-03
    @author Radu Popescu <radu.popescu@epfl.ch>
 */

#ifndef GRAPH_UTIL_H
#define GRAPH_UTIL_H 1

#include <vector>
//#include <algorithm>

#include <boost/shared_ptr.hpp>
#include <boost/bimap.hpp>
#include <Epetra_MpiComm.h>
#include <Epetra_SerialComm.h>
#include <Epetra_BlockMap.h>
#include <Epetra_CrsGraph.h>

#include <lifev/core/LifeV.hpp>


namespace LifeV
{
namespace GraphUtil
{

// Public typedefs
typedef std::vector <LifeV::Int>                idList_Type;
typedef std::vector <idList_Type>               graphPartition_Type;
typedef boost::shared_ptr<graphPartition_Type>  graphPartitionPtr_Type;
typedef boost::bimap<LifeV::UInt, LifeV::UInt>  biMap_Type;
typedef biMap_Type::value_type                  biMapValue_Type;

//! Function that partitions a graph of a subset of elements in a mesh
/*!
    @author Radu Popescu <radu.popescu@epfl.ch>

    TODO: add full description
 */
template <typename MeshType>
void partitionGraph (const idList_Type& vertexList,
					 const MeshType& mesh,
                     const Int numParts,
                     graphPartitionPtr_Type& vertexPartition);
}
};

template <typename MeshType>
void LifeV::GraphUtil::partitionGraph (const idList_Type& vertexList,
					 	 	 	 	   const MeshType& mesh,
					 	 	 	 	   const Int numParts,
					 	 	 	 	   graphPartitionPtr_Type& vertexPartition)
{
	// Build the CRS graph corresponding to the elements in the vertex list
#ifdef EPETRA_MPI
	Epetra_MpiComm comm(MPI_COMM_SELF);
#else
	Epetra_SerialComm comm();
#endif

	// We build a bidirectional map of vertex Id, for local-to-global and
	// global-to-local lookups
	UInt numVertices = vertexList.size();
	biMap_Type vertexIdMap;
	for (UInt i = 0; i < numVertices; ++i) {
		vertexIdMap.insert(biMapValue_Type(i, vertexList[i]));
	}

	Epetra_BlockMap map(static_cast<int>(numVertices), 1, 0, comm);
	Epetra_CrsGraph graph(Copy, map, 0, false);

	UInt elementFacets = MeshType::elementShape_Type::S_numFacets;
    for (UInt lid = 0; lid < numVertices; ++lid)
    {
    	int numAdjacentVertices = 0;
    	std::vector<int> adjacentVertices(0);
    	adjacentVertices.reserve(elementFacets);
        UInt gid = vertexIdMap.left.at (lid);
        for (UInt ifacet = 0; ifacet < elementFacets; ++ifacet)
        {
            UInt facet = mesh.localFacetId(gid, ifacet);
            UInt elem = mesh.facet(facet).firstAdjacentElementIdentity();
            if (elem == gid)
            {
                elem = mesh.facet(facet).secondAdjacentElementIdentity();
            }
            biMap_Type::right_const_iterator it = vertexIdMap.right.find (elem);

            bool inSubGraph = (vertexIdMap.right.end() != it);
            if ( (elem != NotAnId) && (inSubGraph) )
            {
            	adjacentVertices.push_back(it->second);
            	++numAdjacentVertices;
            }
        }
        graph.InsertGlobalIndices(lid, numAdjacentVertices, &adjacentVertices[0]);
    }
    graph.FillComplete();
    graph.Print(std::cout);

	// Partition the CRS graph using Isorropia

	// Return a list of elements in each graph part
}

//template<typename MeshType>
//Int GraphCutterParMETIS<MeshType>::partitionSubGraph (
//    const biMap_Type& vertexMap,
//    const Int numParts,
//    vertexPartition_Type& vertexPartition)
//{
//    // Distribute elements
//    UInt k = vertexMap.size();
//
//    // CAREFUL: ParMetis works on a graph abstraction.
//    // A graph is built over the data structure to be split, each vertex being
//    // a mesh element so hereby a "vertex" is actually a _graph_ vertex,
//    // i. e. a mesh element
//    std::vector<Int> vertexDistribution (M_numProcessors + 1);
//    vertexDistribution[0] = 0;
//    // Evenly distributed graph vertices
//    for (Int i = 0; i < M_numProcessors; ++i)
//    {
//        UInt l = k / (M_numProcessors - i);
//        vertexDistribution[i + 1] = vertexDistribution[i] + l;
//        k -= l;
//    }
//    ASSERT (k == 0, "At this point we should have 0 volumes left") ;
//
//
//    /*
//     * Partition Graph
//     * This array's size is equal to the number of locally-stored vertices:
//     * at the end of the partitioning process, "M_graphVertexLocations" will
//     * contain the partitioning array:
//     * M_graphVertexLocations[m] = n; means that graph vertex m belongs to
//     * subdomain n
//     */
//    std::vector<Int> graphVertexLocations (
//        vertexDistribution[M_numProcessors]
//        - vertexDistribution[0], M_numProcessors);
//
//    /*
//     * Now each processor will take care of its own graph vertices
//     * (i. e. mesh elements).
//     * Nothing guarantees about the neighbor elements distribution across
//     * the processors, since as of now we just split the set of volumes based
//     * on IDs.
//     * Here we building up the neighbor arrays.
//     */
//    UInt localStart = vertexDistribution[M_myPID];
//    UInt localEnd   = vertexDistribution[M_myPID + 1];
//
//    // this vector contains the weights for the edges of the graph,
//    // it is set to null if it is not used.
//    std::vector<Int> graphEdgeWeights;
//    std::vector<Int> adjacencyGraphKeys (1, 0);
//    std::vector<Int> adjacencyGraphValues (0);
//
//    UInt sum = 0;
//
//    for (UInt lid = localStart; lid < localEnd; ++lid)
//    {
//        for (UInt ifacet = 0; ifacet < M_elementFacets; ++ifacet)
//        {
//            UInt gid = vertexMap.left.at (lid);
//            // global ID of the ifacet-th facet in element ie
//            UInt facet = M_mesh->localFacetId (gid, ifacet);
//            // first adjacent element to face "facet"
//            UInt elem = M_mesh->facet (facet).firstAdjacentElementIdentity();
//            if (elem == gid)
//            {
//                elem = M_mesh->facet (facet).secondAdjacentElementIdentity();
//            }
//            biMap_Type::right_const_iterator it = vertexMap.right.find (elem);
//
//            bool inSubGraph = (vertexMap.right.end() != it);
//            if ( (elem != NotAnId) && (inSubGraph) )
//            {
//                // this is the list of adjacency
//                // for each graph vertex, push back the ID of its neighbors
//                //adjacencyGraphValues.push_back(it - vertexMap.right.begin());
//                adjacencyGraphValues.push_back (it->second);
//                ++sum;
//            }
//        }
//        // this is the list of "keys" to access M_adjacencyGraphValues
//        // graph element i has neighbors M_adjacencyGraphValues[ k ],
//        // with M_adjacencyGraphKeys[i] <= k < M_adjacencyGraphKeys[i+1]
//        adjacencyGraphKeys.push_back (sum);
//    }
//
//    // **************
//    // parMetis part
//
//    // this array is to be used for weighted vertices on the graph:
//    // usually we will set it to NULL
//
//    Int* weightVector = 0;
//
//    Int weightFlag = 0;
//
//    Int ncon = 1;
//    Int numflag = 0;
//
//    Int cutGraphEdges; // here will be stored the number of edges cut in the
//    // partitioning process
//
//    // additional options
//    std::vector<Int>  options (3, 0);
//    options[0] = 1; // means that additional options are actually passed
//    options[1] = 3; // level of information to be returned during execution
//    options[2] = 1; // random number seed for the ParMETIS routine
//
//    // fraction of vertex weight to be distributed to each subdomain.
//    // here we want the subdomains to be of the same size
//    std::vector<float> tpwgts (ncon * numParts, 1. / numParts);
//    // imbalance tolerance for each vertex weight
//    std::vector<float> ubvec (ncon, 1.05);
//
//    boost::shared_ptr<Epetra_MpiComm> mpiComm
//        = boost::dynamic_pointer_cast <Epetra_MpiComm> (M_comm);
//    MPI_Comm MPIcomm = mpiComm->Comm();
//
//    /*
//      (from ParMETIS v 3.1 manual)
//      This routine is used to compute a k-way partitioning of a graph
//      on p processors using the multilevel k-way multi-constraint
//      partitioning algorithm.
//    */
//
//    Int nParts = numParts;
//    Int* adjwgtPtr (0);
//    if (graphEdgeWeights.size() > 0)
//    {
//        adjwgtPtr = static_cast<Int*> (&graphEdgeWeights[0]);
//    }
//    ParMETIS_V3_PartKway (static_cast<Int*> (&vertexDistribution[0]),
//                          static_cast<Int*> (&adjacencyGraphKeys[0]),
//                          static_cast<Int*> (&adjacencyGraphValues[0]),
//                          weightVector, adjwgtPtr, &weightFlag, &numflag,
//                          &ncon, &nParts, &tpwgts[0], &ubvec[0],
//                          &options[0], &cutGraphEdges,
//                          &graphVertexLocations[localStart],
//                          &MPIcomm);
//
//    // distribute the resulting partitioning stored in M_graphVertexLocations
//    // to all processors
//    if (M_numProcessors != 1)
//    {
//        M_comm->Barrier();
//        for (Int proc = 0; proc < M_numProcessors; proc++)
//        {
//            UInt procStart  = vertexDistribution[proc];
//            UInt procLength = vertexDistribution[proc + 1]
//                              - vertexDistribution[proc];
//            M_comm->Broadcast (&graphVertexLocations[procStart],
//                               procLength, proc);
//        }
//    }
//
//    // cycling on locally stored vertices
//    vertexPartition.resize (numParts);
//    for (UInt i = 0; i < numParts; ++i)
//    {
//        vertexPartition[i].reset (new std::vector<Int> (0) );
//    }
//    for (UInt ii = 0; ii < graphVertexLocations.size(); ++ii)
//    {
//        // here we are associating the vertex global ID to the subdomain ID
//        vertexPartition[graphVertexLocations[ii]]->push_back (vertexMap.left.at (ii) );
//    }
//
//    return 0;
//}

#endif // GRAPH_UTIL_H
