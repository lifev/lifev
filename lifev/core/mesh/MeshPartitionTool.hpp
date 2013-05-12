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

#ifndef MESH_PARTITION_TOOL_ONLINE_H
#define MESH_PARTITION_TOOL_ONLINE_H 1

#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/shared_ptr.hpp>
#include <Epetra_MpiComm.h>
#include <Teuchos_ParameterList.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>
#include <lifev/core/fem/DOF.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
//#include <lifev/core/mesh/GhostEntityData.hpp>

namespace LifeV
{

/*!
  @brief Class that does flexible mesh partitioning
  @author Radu Popescu radu.popescu@epfl.ch

  This class implements the partitioning of a global mesh and allows choosing
  the graph partitioning tool and the method to build the mesh parts.

  The class template has three parameters: the mesh type, the graph partitioner
  and the mesh part builder. The last two parameters are also class templates,
  parameterized on the mesh type.
*/
template < typename MeshType,
         template <typename> class GraphCutterType,
         template <typename> class MeshPartBuilderType >
class MeshPartitionTool
{
public:
    //! @name Public Types
    //@{
    typedef MeshType                             mesh_Type;
    typedef boost::shared_ptr <
    std::vector<std::vector<Int> > >     graph_Type;
    typedef GraphCutterType<MeshType>            graphCutter_Type;
    typedef MeshPartBuilderType<MeshType>        meshPartBuilder_Type;
    typedef boost::shared_ptr<mesh_Type>         meshPtr_Type;
    typedef std::vector<meshPtr_Type>            partMesh_Type;
    typedef boost::shared_ptr<partMesh_Type>     partMeshPtr_Type;
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
    //! Return a reference to M_ghostDataMap
    //const GhostEntityDataMap_Type&  ghostDataMap() const
    //{
    //  return M_ghostDataMap;
    //}
    //@}

private:
    //! \name Private methods
    //@{
    //! This method performs all the steps for the mesh and graph partitioning
    void run();
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
    boost::shared_ptr<graphCutter_Type>        M_graphCutter;
    boost::shared_ptr<meshPartBuilder_Type>    M_meshPartBuilder;
    bool                                       M_success;
    //@}
}; // class MeshPartitionToolOnline

//
// IMPLEMENTATION
//

// =================================
// Constructors and destructor
// =================================

template < typename MeshType,
         template <typename> class GraphCutterType,
         template <typename> class MeshPartBuilderType >
MeshPartitionTool < MeshType,
                  GraphCutterType,
                  MeshPartBuilderType >::MeshPartitionTool (
                      const meshPtr_Type& mesh,
                      const boost::shared_ptr<Epetra_Comm>& comm,
                      const Teuchos::ParameterList parameters) :
                      M_comm (comm),
                      M_myPID (M_comm->MyPID() ),
                      M_parameters (parameters),
                      M_originalMesh (mesh),
                      M_meshPart(),
                      M_allMeshParts(),
                      M_graphCutter (new graphCutter_Type (M_originalMesh, M_comm, M_parameters) ),
                      M_meshPartBuilder (new meshPartBuilder_Type (M_originalMesh,
                                                                   M_parameters.get<UInt> ("overlap", 0), M_comm) ),
                      M_success (false)
{
    run();
}

// =================================
// Public methods
// =================================

template < typename MeshType,
         template <typename> class GraphCutterType,
         template <typename> class MeshPartBuilderType >
void MeshPartitionTool < MeshType,
     GraphCutterType,
     MeshPartBuilderType >::run()
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

    // Extract the graph from the graphCutter and dispose of the cutter
    graph_Type graph = M_graphCutter->getGraph();
    M_graphCutter.reset();

    if (!M_myPID)
    {
        std::cout << "Building mesh parts ..." << std::endl;
    }

    bool offlineMode = M_parameters.get<bool> ("offline_mode", false);
    if (! offlineMode)
    {
        // Online partitioning
        M_meshPart.reset (new mesh_Type);
        M_meshPartBuilder->run (M_meshPart, graph, M_myPID);

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
            Int numParts = M_parameters.get<Int> ("num_parts");
            /*
             * In offline partitioning mode, with overlap, we must make sure
             * that each time the M_meshPartBuilder is run, which modifies the
             * partition graph, it modifies the original graph, not the
             * augmented graph from the previous part.
             */

            M_allMeshParts.reset (new partMesh_Type (numParts) );
            for (Int curPart = 0; curPart < numParts; ++curPart)
            {
                // Backup the elements of the current graph part
                std::vector<Int> backup ( (*graph) [curPart]);
                M_allMeshParts->at (curPart).reset (new mesh_Type);
                M_meshPartBuilder->run (M_allMeshParts->at (curPart),
                                        graph,
                                        curPart);

                // At this point (*graph)[curPart] has been modified. Restore
                // to the original state
                (*graph) [curPart] = backup;
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

template < typename MeshType,
template <typename> class GraphCutterType,
template <typename> class MeshPartBuilderType >
void
MeshPartitionTool < MeshType,
                  GraphCutterType,
                  MeshPartBuilderType >::showMe (std::ostream& output) const
{
    std::cout << "Sorry, this method is not implemented, yet." << std::endl
    << "We appreciate your interest." << std::endl
    << "Check back in a bit!" << std::endl;
}

} // namespace LifeV

#endif // MESH_PARTITION_TOOL_ONLINE_H
