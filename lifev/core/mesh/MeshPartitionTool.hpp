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
  @brief Class that does flexible mesh partitioning

  @date 16-11-2011
  @author Radu Popescu <radu.popescu@epfl.ch>

  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#ifndef MESH_PARTITION_TOOL_H
#define MESH_PARTITION_TOOL_H 1

#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <Epetra_MpiComm.h>
#include <Teuchos_ParameterList.hpp>

#include <lifev/core/LifeV.hpp>
#include <lifev/core/fem/DOF.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/GraphCutterParMETIS.hpp>
#include <lifev/core/mesh/GraphCutterZoltan.hpp>
#include <lifev/core/mesh/GraphUtil.hpp>
#include <lifev/core/mesh/MeshPartBuilder.hpp>

namespace LifeV
{

using namespace GraphUtil;

/*!
  @brief Class that does flexible mesh partitioning
  @author Radu Popescu radu.popescu@epfl.ch

  This class implements the partitioning of a global mesh and allows choosing
  the graph partitioning tool and the method to build the mesh parts.

  The class is configured through the Teuchos::Parameter list passed to the
  constructor. The following parameters are used:

   graph-lib - std::string - "parmetis" or "zoltan" selects the graph partition
                             library to be used (default "parmetis")
   overlap - UInt - level of overlap for the mesh partition process (default 0)
   offline-mode - bool - mode of operation; offline mode can
                         only be used in serial (default false == online)
   num-parts - Int - (for offline-mode only) sets the number of parts for the
                     mesh partition process (no default value)
   hierarchical - bool - enable hierarchical partitioning mode (default false)
   topology - std::string - value which represents the number of mesh parts per
                            compute node
                            (N == num-parts; topology="m"; N % m == 0)
                            (default "1")

   Notes:

   * the value of the "topology" parameter is given as a string due to a
     requirement of the Zoltan interface
   * When using Zoltan as a graph partition library, additional advanced
     parameters are available. See GraphCutterZoltan.hpp for more information.


*/
template < typename MeshType>
class MeshPartitionTool
{
public:
    //! @name Public Types
    //@{
    typedef MeshType                             mesh_Type;
    typedef GraphCutterBase<mesh_Type>           graphCutter_Type;
    typedef MeshPartBuilder<mesh_Type>           meshPartBuilder_Type;
    typedef boost::shared_ptr<mesh_Type>         meshPtr_Type;
    typedef std::vector<meshPtr_Type>            partMesh_Type;
    typedef boost::shared_ptr<partMesh_Type>     partMeshPtr_Type;
    typedef std::vector<idTablePtr_Type>         vertexPartitionTable_Type;
    typedef boost::shared_ptr<vertexPartitionTable_Type>
    vertexPartitionTablePtr_Type;
    //@}

    //! \name Constructors & Destructors
    //@{
    //! Constructor
    /*!
     * The constructor takes as parameters the global uncut mesh, the Epetra
     * comm of the group or processes involved in the mesh partition process
     * and a Teuchos parameter list.
     *
     * After initializing the object, the constructor will call the private
     * run() method and perform the mesh partition
     *
     * \param mesh - shared pointer to the global uncut mesh
     * \param comm - shared pointer to the Epetra comm object containing the
     *               processes involved in the mesh partition process
     * \param parameters - Teuchos parameter list
    */
    MeshPartitionTool (const meshPtr_Type& mesh,
                       const boost::shared_ptr<Epetra_Comm>& comm,
                       const Teuchos::ParameterList parameters
                       = Teuchos::ParameterList() );

    //! Empty destructor
    ~MeshPartitionTool() {}
    //@}

    //! \name Public Methods
    //! Prints information about the state (data) of the object
    void showMe (std::ostream& output = std::cout) const;
    //@}

    //! \name Get Methods
    //@{
    //! Return the succeess state of the partitioning (true || false)
    bool success() const
    {
        return M_success;
    }
    //! Return a shared pointer to the mesh part (online mode, const ver.)
    const meshPtr_Type& meshPart() const
    {
        return M_meshPart;
    }
    //! Return a shared pointer to the mesh parts
    /*!
     * Return a shared pointer to the mesh parts on offline partition (const)
     */
    const partMeshPtr_Type& allMeshParts() const
    {
        return M_allMeshParts;
    }

    //! Return a shared pointer to the second stage graph parts (for ShyLU-MT)
    /*!
     * Return a shared pointer to the second stage graph parts (for ShyLU-MT)
     * Offline mode
     */
    const vertexPartitionTablePtr_Type& secondStageParts() const
    {
        return M_secondStageParts;
    }

    //! Return a shared pointer to the second stage graph parts (for ShyLU-MT)
    /*!
     * Return a shared pointer to the second stage graph parts (for ShyLU-MT)
     * Online mode
     */
    const idTablePtr_Type& mySecondStageParts() const
    {
        return M_secondStageParts->at (0);
    }
    //@}

private:
    //! \name Private methods
    //@{
    //! This method performs all the steps for the mesh and graph partitioning
    void run();

    //! Initialize M_entityPID
    void fillEntityPID (idTablePtr_Type graph);

    //! Global to local element ID conversion for second stage
    void globalToLocal (const Int curPart);
    //@}

    // Private copy constructor and assignment operator are disabled
    MeshPartitionTool (const MeshPartitionTool&);
    MeshPartitionTool& operator= (const MeshPartitionTool&);

    //! Private Data Members
    //@{
    boost::shared_ptr<Epetra_Comm>             M_comm;
    Int                                        M_myPID;
    Teuchos::ParameterList                     M_parameters;
    meshPtr_Type                               M_originalMesh;
    meshPtr_Type                               M_meshPart;
    partMeshPtr_Type                           M_allMeshParts;
    std::string                                M_graphLib;
    boost::shared_ptr<graphCutter_Type>        M_graphCutter;
    boost::shared_ptr<meshPartBuilder_Type>    M_meshPartBuilder;
    bool                                       M_success;
    bool                                       M_secondStage;
    Int                                        M_secondStageNumParts;
    vertexPartitionTablePtr_Type               M_secondStageParts;

    //! Store ownership for each entity, subdivided by entity type
    typename meshPartBuilder_Type::entityPID_Type M_entityPID;

    //@}
}; // class MeshPartitionToolOnline

//
// IMPLEMENTATION
//

// =================================
// Constructors and destructor
// =================================

template < typename MeshType>
MeshPartitionTool < MeshType >::MeshPartitionTool (
    const meshPtr_Type& mesh,
    const boost::shared_ptr<Epetra_Comm>& comm,
    const Teuchos::ParameterList parameters) :
    M_comm (comm),
    M_myPID (M_comm->MyPID() ),
    M_parameters (parameters),
    M_originalMesh (mesh),
    M_meshPart(),
    M_allMeshParts(),
    M_graphLib (M_parameters.get<std::string> ("graph-lib", "parmetis") ),
    M_meshPartBuilder (new meshPartBuilder_Type (M_originalMesh,
                                                 M_parameters.get<UInt> ("overlap", 0), M_comm) ),
    M_success (false),
    M_secondStage (M_parameters.get<bool> ("second-stage", false) ),
    M_secondStageNumParts (M_parameters.get<Int> ("second-stage-num-parts", 1) ),
    M_secondStageParts (new vertexPartitionTable_Type)
{
    if (! M_graphLib.compare ("parmetis") )
    {
        M_graphCutter.reset (new GraphCutterParMETIS<mesh_Type> (M_originalMesh, M_comm, M_parameters) );
    }
    else if (! M_graphLib.compare ("zoltan") )
    {
        M_graphCutter.reset (new GraphCutterZoltan<mesh_Type> (M_originalMesh, M_comm, M_parameters) );
    }
    else
    {
        std::cout << "Graph partitioner type not defined.\n";
    }

    run();
}

// =================================
// Public methods
// =================================

template < typename MeshType>
void MeshPartitionTool < MeshType >::run()
{
    if (!M_myPID)
    {
        std::cout << "Partitioning mesh graph ..." << std::endl;
    }
    // If the graph partitioning failed, just abort the mesh partitioning.
    // MeshPartitionTool::success() will return false.
    if (M_graphCutter->run() )
    {
        return;
    }

    // Extract the graph from the graphCutter
    idTablePtr_Type graph = M_graphCutter->getGraph();

    // Dispose of the graph partitioner object
    M_graphCutter.reset();

    // Get the current operation mode and number of parts
    bool offlineMode = M_parameters.get<bool> ("offline-mode", false);

    // Do a second stage graph partitioning for ShyLU-MT
    if (M_secondStage)
    {
        if (!M_myPID)
        {
            std::cout << "Performing second stage partitioning ..."
                      << std::endl;
        }
        // MPI comm object for second stage is always MPI_COMM_SELF
#ifdef EPETRA_MPI
        boost::shared_ptr<Epetra_Comm>
        secondStageComm (new Epetra_MpiComm (MPI_COMM_SELF) );
#else
        boost::shared_ptr<Epetra_Comm>
        secondStageComm (new Epetra_SerialComm);
#endif

        Teuchos::ParameterList secondStageParams;
        secondStageParams.set<Int> ("num-parts", M_secondStageNumParts);
        secondStageParams.set<Int> ("my-pid", M_myPID);
        secondStageParams.set<bool> ("verbose", false);
        if (! offlineMode)
        {
            M_secondStageParts->resize (1);
            // For each set of elements in graph perform a second stage partitioning
            const idListPtr_Type currentIds = graph->at (M_myPID);
            partitionGraphParMETIS (currentIds, *M_originalMesh,
                                    secondStageParams,
                                    M_secondStageParts->at (0),
                                    secondStageComm);
        }
        else
        {
            M_secondStageParts->resize (graph->size() );
            // For each set of elements in graph perform a second stage partitioning
            for (Int i = 0; i < graph->size(); ++i)
            {
                const idListPtr_Type& currentIds = graph->at (i);
                partitionGraphParMETIS (currentIds, *M_originalMesh,
                                        secondStageParams,
                                        M_secondStageParts->at (i),
                                        secondStageComm);
            }
        }
    }

    // Fill entity PID
    if (!M_myPID)
    {
        std::cout << "Filling entity PID lists ..." << std::endl;
    }
    fillEntityPID (graph);

    if (!M_myPID)
    {
        std::cout << "Building mesh parts ..." << std::endl;
    }

    if (! offlineMode)
    {
        // Online partitioning
        M_meshPart.reset (new mesh_Type (M_comm) );
        M_meshPart->setIsPartitioned (true);
        M_meshPartBuilder->run (M_meshPart, graph, M_entityPID, M_myPID);

        // Make the global to local element ID conversion for the second stage
        if (M_secondStage)
        {
            globalToLocal (0);
        }

        // Reset the mesh part builder
        M_meshPartBuilder->reset();

        // Mark the partition as successful
        M_success = true;
    }
    else
    {
        // Offline partitioning
        if (M_comm->NumProc() != 1)
        {
            if (!M_myPID)
            {
                std::cout << "Offline partition must be done in serial."
                          << std::endl;
            }
        }
        else
        {
            /*
             * In offline partitioning mode, with overlap, we must make sure
             * that each time the M_meshPartBuilder is run, which modifies the
             * partition graph, it modifies the original graph, not the
             * augmented graph from the previous part.
             */

            Int numParts = M_parameters.get<Int> ("num-parts");
            M_allMeshParts.reset (new partMesh_Type (numParts) );
            for (Int curPart = 0; curPart < numParts; ++curPart)
            {
                // Backup the elements of the current graph part
                idList_Type backup ( * (graph->at (curPart) ) );
                M_allMeshParts->at (curPart).reset (new mesh_Type);
                M_allMeshParts->at (curPart)->setIsPartitioned (true);
                M_meshPartBuilder->run (M_allMeshParts->at (curPart),
                                        graph, M_entityPID,
                                        curPart);

                // At this point (*graph)[curPart] has been modified. Restore
                // to the original state
                * (graph->at (curPart) ) = backup;

                // Make the global to local element ID conversion for the second stage
                if (M_secondStage)
                {
                    globalToLocal (curPart);
                }

                // Reset the mesh part builder
                M_meshPartBuilder->reset();

                // Mark the partition as successful
                M_success = true;
            }
        }
    }

    if (!M_myPID)
    {
        std::cout << "Mesh partition complete." << std::endl;
    }
    // Destroy the mesh part builder to clear memory
    M_meshPartBuilder.reset();
    // Release the pointer to the original uncut mesh
    M_originalMesh.reset();
}

template<typename MeshType>
void
MeshPartitionTool<MeshType>::fillEntityPID (idTablePtr_Type graph)
{
    Int numParts = graph->size();

    // initialize entity PIDs to 0
    M_entityPID.points.resize   ( M_originalMesh->numPoints(),   0 );
    M_entityPID.elements.resize ( M_originalMesh->numElements(), 0 );
    M_entityPID.facets.resize   ( M_originalMesh->numFacets(),   0 );
    M_entityPID.ridges.resize   ( M_originalMesh->numRidges(),   0 );

    // check: parallel algorithm seems to be slower for this
    // p = 0 can be skipped since M_entityPID is already initialized at that value
    for ( Int p = 1; p < numParts; p++ )
    {
        for ( UInt e = 0; e < graph->at (p)->size(); e++ )
        {
            // point block
            for ( UInt k = 0; k < mesh_Type::element_Type::S_numPoints; k++ )
            {
                const ID& pointID = M_originalMesh->element ( graph->at (p)->at (e) ).point ( k ).id();
                // pointPID should be the maximum between the procs that own it
                M_entityPID.points[ pointID ] = std::max ( M_entityPID.points[ pointID ], p );
            }

            // elem block
            const ID& elemID = M_originalMesh->element ( graph->at (p)->at (e) ).id();
            // at his stage each element belongs to a single partition, overlap is not yet done.
            M_entityPID.elements[ elemID ] = p;

            // facet block
            for ( UInt k = 0; k < mesh_Type::element_Type::S_numFacets; k++ )
            {
                const ID& facetID = M_originalMesh->facet ( M_originalMesh->localFacetId ( elemID, k ) ).id();
                // facetPID should be the maximum between the proc that own it
                M_entityPID.facets[ facetID ] = std::max ( M_entityPID.facets[ facetID ], p );
            }

            // ridge block
            for ( UInt k = 0; k < mesh_Type::element_Type::S_numRidges; k++ )
            {
                const ID& ridgeID = M_originalMesh->ridge ( M_originalMesh->localRidgeId ( elemID, k ) ).id();
                // ridgePID should be the maximum between the proc that own it
                M_entityPID.ridges[ ridgeID ] = std::max ( M_entityPID.ridges[ ridgeID ], p );
            }
        }
    }
}

template < typename MeshType>
void
MeshPartitionTool < MeshType >::globalToLocal (const Int curPart)
{
    const std::map<Int, Int>& globalToLocalMap =
        M_meshPartBuilder->globalToLocalElement();
    idTable_Type& currentGraph = * (M_secondStageParts->at (curPart) );

    for (Int i = 0; i < currentGraph.size(); ++i)
    {
        int currentSize = currentGraph[i]->size();
        idList_Type& currentElements = * (currentGraph[i]);
        for (Int j = 0; j < currentSize; ++j)
        {
            currentElements[j] = globalToLocalMap.find (currentElements[j])->second;
        }
    }
}

template < typename MeshType>
void
MeshPartitionTool < MeshType >::showMe (std::ostream& output) const
{
    output << "Sorry, this method is not implemented, yet." << std::endl
           << "We appreciate your interest." << std::endl
           << "Check back in a bit!" << std::endl;
}

} // namespace LifeV

#endif // MESH_PARTITION_TOOL_H
