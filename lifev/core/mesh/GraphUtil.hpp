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

#include <parmetis.h>

#include <boost/shared_ptr.hpp>
#include <boost/bimap.hpp>
#include <Epetra_MpiComm.h>
#include <Epetra_SerialComm.h>
#include <Epetra_BlockMap.h>
#include <Epetra_CrsGraph.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <lifev/core/LifeV.hpp>

namespace LifeV
{
namespace GraphUtil
{

// Public typedefs
typedef boost::shared_ptr<Epetra_Comm>          commPtr_Type;
typedef std::vector <LifeV::Int>                idList_Type;
typedef boost::shared_ptr<idList_Type>          idListPtr_Type;
typedef std::vector <idListPtr_Type>            idTable_Type;
typedef boost::shared_ptr<idTable_Type>         idTablePtr_Type;
typedef std::set<Int> 					        idSet_Type;
typedef boost::shared_ptr<idSet_Type>	        idSetPtr_Type;
typedef std::vector<idSetPtr_Type>              idSetGroup_Type;
typedef boost::shared_ptr<idSetGroup_Type>      idSetGroupPtr_Type;
typedef boost::bimap<LifeV::UInt, LifeV::UInt>  biMap_Type;
typedef biMap_Type::value_type                  biMapValue_Type;

//! Function that partitions a graph of a subset of elements in a mesh
/*!
    @author Radu Popescu <radu.popescu@epfl.ch>

    TODO: add full description
 */
template <typename MeshType>
void partitionGraphParMETIS (const idListPtr_Type& vertexList,
						  	 const MeshType& mesh,
						  	 const Teuchos::ParameterList& params,
						  	 idTablePtr_Type& vertexPartition,
						  	 commPtr_Type& comm);
}
};

template <typename MeshType>
void LifeV::GraphUtil::partitionGraphParMETIS (const idListPtr_Type& vertexList,
											   const MeshType& mesh,
											   const Teuchos::ParameterList& params,
											   idTablePtr_Type& vertexPartition,
											   commPtr_Type& comm)
{
	Int numProc = comm->NumProc();
	Int myPID = comm->MyPID();

	// We build a bidirectional map of vertex Id, for local-to-global and
	// global-to-local lookups
	UInt numVertices = vertexList->size();
	biMap_Type vertexIdMap;
	for (UInt i = 0; i < numVertices; ++i) {
		vertexIdMap.insert(biMapValue_Type(i, vertexList->at(i)));
	}

    // A graph is built over the data structure to be split, each vertex being
    // a mesh element so hereby a "vertex" is actually a _graph_ vertex,
    // i. e. a mesh element
    std::vector<Int> vertexDistribution (numProc + 1);
    vertexDistribution[0] = 0;
    UInt k = numVertices;
    // Evenly distributed graph vertices
    for (Int i = 0; i < numProc; ++i)
    {
        UInt l = k / (numProc - i);
        vertexDistribution[i + 1] = vertexDistribution[i] + l;
        k -= l;
    }

    /*
     * Partition Graph
     * This array's size is equal to the number of locally-stored vertices:
     * at the end of the partitioning process, "M_graphVertexLocations" will
     * contain the partitioning array:
     * M_graphVertexLocations[m] = n; means that graph vertex m belongs to
     * subdomain n
     */
    std::vector<Int> vertexLocations (vertexDistribution[numProc]
									  - vertexDistribution[0], numProc);
    UInt localStart = vertexDistribution[myPID];
    UInt localEnd   = vertexDistribution[myPID + 1];


    // this vector contains the weights for the edges of the graph,
    // it is set to null if it is not used.
    std::vector<Int> graphEdgeWeights;
    std::vector<Int> adjacencyGraphKeys (1, 0);
    std::vector<Int> adjacencyGraphValues (0);

    UInt sum = 0;

    UInt elementFacets = MeshType::elementShape_Type::S_numFacets;

    for (UInt lid = localStart; lid < localEnd; ++lid)
    {
        for (UInt ifacet = 0; ifacet < elementFacets; ++ifacet)
        {
            UInt gid = vertexIdMap.left.at (lid);
            // global ID of the ifacet-th facet in element ie
            UInt facet = mesh.localFacetId (gid, ifacet);
            // first adjacent element to face "facet"
            UInt elem = mesh.facet (facet).firstAdjacentElementIdentity();
            if (elem == gid)
            {
                elem = mesh.facet (facet).secondAdjacentElementIdentity();
            }
            biMap_Type::right_const_iterator it = vertexIdMap.right.find (elem);

            bool inSubGraph = (vertexIdMap.right.end() != it);
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

    Int* weightVector = 0;
    Int weightFlag = 0;
    Int ncon = 1;
    Int numflag = 0;
    Int cutGraphEdges;

    // additional options
    Int options[3] = {1, 3, 1};

    // fraction of vertex weight to be distributed to each subdomain.
    // here we want the subdomains to be of the same size
    Int numParts = params.get<Int>("num-parts");
    std::vector<float> tpwgts (ncon * numParts, 1. / numParts);
    // imbalance tolerance for each vertex weight
    std::vector<float> ubvec (ncon, 1.05);

    boost::shared_ptr<Epetra_MpiComm> mpiComm
        = boost::dynamic_pointer_cast <Epetra_MpiComm> (comm);
    MPI_Comm MPIcomm = mpiComm->Comm();

    Int* adjwgtPtr (0);
    if (graphEdgeWeights.size() > 0)
    {
        adjwgtPtr = static_cast<Int*> (&graphEdgeWeights[0]);
    }
    ParMETIS_V3_PartKway (static_cast<Int*> (&vertexDistribution[0]),
                          static_cast<Int*> (&adjacencyGraphKeys[0]),
                          static_cast<Int*> (&adjacencyGraphValues[0]),
                          weightVector, adjwgtPtr, &weightFlag, &numflag,
                          &ncon, &numParts, &tpwgts[0], &ubvec[0],
                          &options[0], &cutGraphEdges,
                          &vertexLocations[localStart],
                          &MPIcomm);

    // distribute the resulting partitioning stored in M_graphVertexLocations
    // to all processors
    if (numProc != 1)
    {
        comm->Barrier();
        for (Int proc = 0; proc < numProc; proc++)
        {
            UInt procStart  = vertexDistribution[proc];
            UInt procLength = vertexDistribution[proc + 1]
                              - vertexDistribution[proc];
            comm->Broadcast (&vertexLocations[procStart], procLength, proc);
        }
    }

    idTablePtr_Type vertexIds(new idTable_Type(numParts));
    // cycling on locally stored vertices
    for (UInt i = 0; i < numParts; ++i)
    {
        vertexIds->at(i).reset(new idList_Type(0));
        vertexIds->at(i)->reserve(vertexLocations.size() / numParts);
    }
    for (UInt ii = 0; ii < vertexLocations.size(); ++ii)
    {
        // here we are associating the vertex global ID to the subdomain ID
        vertexIds->at(vertexLocations[ii])->push_back(vertexIdMap.left.at(ii));
    }

    vertexPartition = vertexIds;
}

#endif // GRAPH_UTIL_H
