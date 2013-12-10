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
    M_vertexPartition()
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

    idListPtr_Type vertexList(new idList_Type(numVertices));
    for (Int i = 0; i < numVertices; ++i) {
    	vertexList->at(i) = i;
    }

    // Call the partitionSubGraph method on the vertexList that was
    // prepared
    Teuchos::ParameterList pList;
    pList.set<Int>("num-parts", M_numParts);
    GraphUtil::partitionGraphParMETIS(vertexList, *(M_mesh), pList,
    								  M_vertexPartition, M_comm);

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
    idListPtr_Type vertexList(new idList_Type(numVertices));
    for (Int i = 0; i < numVertices; ++i)
    {
        vertexList->at(i) = i;
    }

    vertexPartitionPtr_Type tempVertexPartition;

    /*
     * After calling partitionSubGraph, tempVertexPartition will contain
     * numSubdomains vectors with the graph vertices of each subdomain
     */
    Teuchos::ParameterList pList;
    pList.set<Int>("num-parts", numSubdomains);
    GraphUtil::partitionGraphParMETIS(vertexList, *(M_mesh), pList,
    								  tempVertexPartition, M_comm);

    /*
     * Step two is to partition each subdomain into the number of sub parts
     * denoted by the M_topology parameter
     */
    pList.set<Int>("num-parts", M_topology);
    M_vertexPartition->resize (M_numParts);
    Int currentPart = 0;
    for (Int i = 0; i < numSubdomains; ++i)
    {
    	const idList_Type& currentVertices = *(tempVertexPartition->at(i));
        idListPtr_Type subdomainVertexMap(new idList_Type(currentVertices.size()));
        for (Int k = 0; k < currentVertices.size(); ++k) {
            subdomainVertexMap->at(k) = currentVertices[k];
        }

        vertexPartitionPtr_Type subdomainParts;

		GraphUtil::partitionGraphParMETIS(subdomainVertexMap, *(M_mesh),
										  pList, subdomainParts, M_comm);

        for (Int j = 0; j < M_topology; ++j)
        {
            M_vertexPartition->at(currentPart++) = subdomainParts->at(j);
        }
    }

    return 0;
}

} // Namespace LifeV

#endif // GRAPH_PARTITION_TOOL_PARMETIS_H
