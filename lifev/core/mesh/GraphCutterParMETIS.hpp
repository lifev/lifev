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

#include <boost/shared_ptr.hpp>
#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>
#include <parmetis.h>


#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>

namespace LifeV {

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
class GraphCutterParMETIS
{
public:
    //! @name Public Types
    //@{
	typedef Teuchos::ParameterList                          pList_Type;
	typedef boost::shared_ptr<Epetra_Comm>                  commPtr_Type;
    typedef MeshType                                        mesh_Type;
    typedef boost::shared_ptr<mesh_Type>                    meshPtr_Type;
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
    GraphCutterParMETIS(meshPtr_Type& mesh,
                       commPtr_Type& comm,
                       pList_Type& parameters);

    //! Destructor
    ~GraphCutterParMETIS() {}
    //@}

    //! @name Public methods
    //@{
    //! Performs the graph partitioning
    void run();
    //@}

    //! @name Get Methods
    //@{
    //! Get a pointer to one of the partitions
    const std::vector<Int>& getPart(const UInt i) const
	{
    	return M_elementDomains[i];
	}
    //@}

private:
    //! @name Private methods
    //@{
    //! Set values for all the parameters, with default values where needed
    void setParameters(pList_Type& parameters);
    //@}

    // Private copy constructor and assignment operator are disabled
    GraphCutterParMETIS(const GraphCutterParMETIS&);
    GraphCutterParMETIS& operator=(const GraphCutterParMETIS&);

    //! @name Private Methods
    //@{
    //@}
 
    commPtr_Type             				   M_comm;
    Int                                        M_myPID;
    Int                                        M_numProcessors;
    Int                                        M_numParts;
    Int                                        M_numPartsPerProcessor;
    pList_Type                                 M_parameters;
    meshPtr_Type                               M_mesh;
    std::vector<Int>                     	   M_vertexDistribution;
    std::vector<Int>                           M_graphVertexLocations;
    std::vector<Int>                           M_adjacencyGraphKeys;
    std::vector<Int>                           M_adjacencyGraphValues;
    UInt                                 	   M_elementVertices;
    UInt                                       M_elementFacets;
    UInt                                       M_elementRidges;
    UInt                                       M_facetVertices;
    std::vector<std::vector<Int> >			   M_elementDomains;
};

//
// IMPLEMENTATION
//

// =================================
// Constructors and Destructor
// =================================

template<typename MeshType>
GraphCutterParMETIS<MeshType>::GraphCutterParMETIS(meshPtr_Type& mesh,
                                   	   	   	   commPtr_Type& comm,
											   pList_Type& parameters) :
    M_comm(comm),
    M_myPID(M_comm->MyPID()),
    M_numProcessors(M_comm->NumProc()),
    M_numParts(0),
    M_numPartsPerProcessor(0),
    M_parameters(),
    M_mesh(mesh)
{
    setParameters(parameters);
}

template<typename MeshType>
void GraphCutterParMETIS<MeshType>::setParameters(pList_Type& parameters)
{
    // Here put some default values for the parameters and then import
    // the user supplied list, overwriting the corresponding parameters
    M_parameters.set("num_parts", static_cast<Int>(M_comm->NumProc()),
    				 "The desired number of parts");
    M_parameters.set("topology", "1",
    				 "The topology of the mesh partition process.");

    M_parameters.setParameters(parameters);

    M_numParts = M_parameters.get<Int>("num_parts");
    M_numPartsPerProcessor = M_numParts / M_numProcessors;
}

template<typename MeshType>
void GraphCutterParMETIS<MeshType>::run()
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

	// Distribute elements
    // ParMETIS is able to work in parallel: how many processors does it have
    // at hand?
    Int numProcessors = M_comm->NumProc();
    M_myPID = M_comm->MyPID();

    // CAREFUL: ParMetis works on a graph abstraction.
    // A graph is built over the data structure to be split, each vertex being
    // a mesh element so hereby a "vertex" is actually a _graph_ vertex,
    // i. e. a mesh element
    M_vertexDistribution.resize(numProcessors + 1);
    M_vertexDistribution[0] = 0;

    UInt k = M_mesh->numElements();

    // Evenly distributed graph vertices
    for (Int i = 0; i < numProcessors; ++i)
    {
        UInt l = k / (numProcessors - i);
        M_vertexDistribution[i + 1] = M_vertexDistribution[i] + l;
        k -= l;
    }
    ASSERT(k == 0, "At this point we should have 0 volumes left") ;


    /*
	 * Partition Graph
     * This array's size is equal to the number of locally-stored vertices:
     * at the end of the partitioning process, "M_graphVertexLocations" will
     * contain the partitioning array:
     * M_graphVertexLocations[m] = n; means that graph vertex m belongs to
     * subdomain n
     */
    M_graphVertexLocations.resize(M_vertexDistribution[M_comm->NumProc()]
								  - M_vertexDistribution[0], M_comm->NumProc());

    /*
     * Now each processor will take care of its own graph vertices
     * (i. e. mesh elements).
     * Nothing guarantees about the neighbor elements distribution across
     * the processors, since as of now we just split the set of volumes based
     * on IDs.
     * Here we building up the neighbor arrays.
     */
    UInt localStart = M_vertexDistribution[M_myPID];
    UInt localEnd   = M_vertexDistribution[M_myPID + 1];

    // this vector contains the weights for the edges of the graph,
    // it is set to null if it is not used.
    std::vector<Int> graphEdgeWeights;

    M_adjacencyGraphKeys.resize(0);
    M_adjacencyGraphKeys.push_back(0);

    UInt sum = 0;

    for (UInt ie = localStart; ie < localEnd; ++ie)
    {
        for (UInt ifacet = 0; ifacet < M_elementFacets; ++ifacet)
        {
            // global ID of the ifacet-th facet in element ie
            UInt facet = M_mesh->localFacetId(ie, ifacet);
            // first adjacent element to face "facet"
            UInt elem = M_mesh->facet(facet).firstAdjacentElementIdentity();
            if (elem == ie)
            {
                elem = M_mesh->facet(facet).secondAdjacentElementIdentity();
            }
            if (elem != NotAnId)
            {
                // this is the list of adjacency
                // for each graph vertex, push back the ID of its neighbors
                M_adjacencyGraphValues.push_back(elem);
                ++sum;
            }
        }
        // this is the list of "keys" to access M_adjacencyGraphValues
        // graph element i has neighbors M_adjacencyGraphValues[ k ],
        // with M_adjacencyGraphKeys[i] <= k < M_adjacencyGraphKeys[i+1]
        M_adjacencyGraphKeys.push_back(sum);
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
    std::vector<Int>  options(3,0);
    options[0] = 1; // means that additional options are actually passed
    options[1] = 3; // level of information to be returned during execution
    options[2] = 1; // random number seed for the ParMETIS routine

    // fraction of vertex weight to be distributed to each subdomain.
    // here we want the subdomains to be of the same size
    std::vector<float> tpwgts(ncon * M_numParts, 1. / M_numParts);
    // imbalance tolerance for each vertex weight
    std::vector<float> ubvec(ncon, 1.05);

    boost::shared_ptr<Epetra_MpiComm> mpiComm
    		= boost::dynamic_pointer_cast <Epetra_MpiComm> (M_comm);
    MPI_Comm MPIcomm = mpiComm->Comm();

    Int nprocs;
    MPI_Comm_size(MPIcomm, &nprocs);

    /*
      (from ParMETIS v 3.1 manual)
      This routine is used to compute a k-way partitioning of a graph
      on p processors using the multilevel k-way multi-constraint
      partitioning algorithm.
    */

    Int numberParts = (Int) M_numParts;

    Int* adjwgtPtr(0);
    if (graphEdgeWeights.size() > 0)
    {
        adjwgtPtr = static_cast<Int*>(&graphEdgeWeights[0]);
    }
    ParMETIS_V3_PartKway(static_cast<Int*>(&M_vertexDistribution[0]),
                         static_cast<Int*>(&M_adjacencyGraphKeys[0]),
                         static_cast<Int*>(&M_adjacencyGraphValues[0]),
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
    M_elementDomains.resize(M_numParts);

    // cycling on locally stored vertices
    for (UInt ii = 0; ii < M_graphVertexLocations.size(); ++ii)
    {
        // here we are associating the vertex global ID to the subdomain ID
        M_elementDomains[ M_graphVertexLocations[ ii ] ].push_back( ii );
    }
}

} // Namespace LifeV

#endif // GRAPH_PARTITION_TOOL_PARMETIS_H
