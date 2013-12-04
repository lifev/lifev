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
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <Isorropia_EpetraPartitioner.hpp>

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
void partitionGraphIsorropia (const idList_Type& vertexList,
					 const MeshType& mesh,
                     const Teuchos::ParameterList& params,
                     graphPartitionPtr_Type& vertexPartition);
}
};

template <typename MeshType>
void LifeV::GraphUtil::partitionGraphIsorropia (const idList_Type& vertexList,
					 	 	 	 	   	   	    const MeshType& mesh,
					 	 	 	 	   	   	    const Teuchos::ParameterList& params,
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

	// Build a graph of the elements, store it as an Epetra_CrsGraph
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
        graph.InsertGlobalIndices(lid, numAdjacentVertices,
        						  &adjacentVertices[0]);
    }
    graph.FillComplete();

	// Partition the CRS graph using Isorropia
    Int numParts = params.get<Int>("num-parts");
    Teuchos::ParameterList partitionParams;
    partitionParams.set<std::string>("NUM PARTS",
    								 boost::lexical_cast<std::string>(numParts));
        Teuchos::RCP<const Epetra_CrsGraph> graphPtr =
    		Teuchos::rcp(&graph, false);
    Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner =
    		Teuchos::rcp(new Isorropia::Epetra::Partitioner(graphPtr,
    														partitionParams,
    														true));

	// Build a list of elements in each graph part
    graphPartitionPtr_Type vertexIds(new graphPartition_Type);
    vertexIds->resize(numParts);
    for (int i = 0; i < numParts; ++i) {
    	// For each subpart get list of local element ids from the partitioner
    	int currentSize = partitioner->numElemsInPart(i);
    	std::vector<int> currentElements(currentSize);
    	partitioner->elemsInPart(i, &currentElements[0], currentSize);

        // Must do a local-to-global ID conversion
    	idList_Type& currentIdList = vertexIds->at(i);
		currentIdList.resize(currentSize);
    	for (int j = 0; j < currentSize; ++j) {
    		currentIdList[j] = vertexIdMap.left.at(currentElements[j]);
//    		currentIdList[j] = currentElements[j];
    	}
    }

    // Return the list of elements in each part
    vertexPartition = vertexIds;
}

#endif // GRAPH_UTIL_H
