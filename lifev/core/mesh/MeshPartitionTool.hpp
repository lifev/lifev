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
  @brief Class that does mesh partitioning with flexible graph partitioning

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
  @brief Class that does mesh partitioning with flexible graph partitioning
  @author Radu Popescu radu.popescu@epfl.ch

  This class implements the partitioning of a global mesh using a chosen
  graph partitioning tool. In this way, all the graph operations are
  abstracted. The graph partitioning tool is a member of this class and
  receives the global mesh, computing the redistribution of element GID
  across the given number of partitions.
  Based on this information, the mesh partition tool builds the mesh objects
  corresponding to the new partitions.
*/
template<typename MeshType,
		 template <typename> class GraphCutterType,
		 template <typename> class MeshPartBuilderType >
class MeshPartitionTool
{
public:
    //! @name Public Types
    //@{
    typedef MeshType                             mesh_Type;
    typedef GraphCutterType<MeshType>            graphCutter_Type;
    typedef MeshPartBuilderType<MeshType>        meshPartBuilder_Type;
    typedef boost::shared_ptr<mesh_Type>         meshPtr_Type;
    //! Container for the ghost data
    //typedef std::vector <GhostEntityData>        GhostEntityDataContainer_Type;
    //! Map processor -> container for the ghost data
    //typedef std::map <UInt, GhostEntityDataContainer_Type> GhostEntityDataMap_Type;
    //@}

    //! \name Constructors & Destructors
    //@{
    //! Default empty constructor
    MeshPartitionTool();

    //! Constructor
    /*!
     * TODO: Write description
    */
    MeshPartitionTool (const meshPtr_Type& mesh,
                       const boost::shared_ptr<Epetra_Comm>& comm,
                       const Teuchos::ParameterList parameters
                       	   	   	   = Teuchos::ParameterList());

    //! Empty destructor
    ~MeshPartitionTool() {}
    //@}

    //! \name Public Methods
    //@{
    //! Configures the mesh partitioning tool
    /*!
     * TODO: Write description
    */
    void setup(meshPtr_Type& mesh,
               boost::shared_ptr<Epetra_Comm>& comm,
               Teuchos::ParameterList& parameters);

    //! This method performs all the steps for the mesh and graph partitioning
    void run();

    //! Prints information about the state (data) of the object
    void showMe(std::ostream& output = std::cout) const;
    //@}

    //! \name Get Methods
    //@{
    //! Return a shared pointer to the mesh partition
    const meshPtr_Type& meshPart() const {return M_meshPart;}
    meshPtr_Type& meshPart() {return M_meshPart;}
    //! Return a reference to M_ghostDataMap
    //const GhostEntityDataMap_Type&  ghostDataMap() const
    //{
    //	return M_ghostDataMap;
    //}
    //@}

private:
    // Private copy constructor and assignment operator are disabled
    MeshPartitionTool(const MeshPartitionTool&);
    MeshPartitionTool& operator=(const MeshPartitionTool&);

    //! Private Data Members
    //@{
    boost::shared_ptr<Epetra_Comm>             M_comm;
    Int                                        M_myPID;
    boost::shared_ptr<std::vector<Int> >       M_myElements;
    Teuchos::ParameterList                     M_parameters;
    meshPtr_Type                               M_originalMesh;
    meshPtr_Type                               M_meshPart;
    boost::shared_ptr<graphCutter_Type>        M_graphCutter;
    boost::shared_ptr<meshPartBuilder_Type>    M_meshPartBuilder;
    //GhostEntityDataMap_Type                    M_ghostDataMap;
    //@}
}; // class MeshPartitionToolOnline

//
// IMPLEMENTATION
//

// =================================
// Constructors and destructor
// =================================

template<typename MeshType,
		 template <typename> class GraphCutterType,
		 template <typename> class MeshPartBuilderType>
MeshPartitionTool<MeshType,
				  GraphCutterType,
				  MeshPartBuilderType>::MeshPartitionTool() :
    M_comm(),
    M_myPID(),
    M_parameters(),
    M_originalMesh(),
    M_meshPart(),
    M_graphCutter(),
    M_meshPartBuilder()
{}

template<typename MeshType,
		 template <typename> class GraphCutterType,
		 template <typename> class MeshPartBuilderType>
MeshPartitionTool<MeshType,
				  GraphCutterType,
				  MeshPartBuilderType>::MeshPartitionTool(
		const meshPtr_Type& mesh,
        const boost::shared_ptr<Epetra_Comm>& comm,
        const Teuchos::ParameterList parameters) :
    M_comm(comm),
    M_myPID(M_comm->MyPID()),
    M_parameters(parameters),
    M_originalMesh(mesh),
    M_meshPart(new MeshType),
    M_graphCutter(new graphCutter_Type(M_originalMesh, M_comm, M_parameters)),
	M_meshPartBuilder(new meshPartBuilder_Type(M_originalMesh, M_comm))
{
    run();
}

// =================================
// Public methods
// =================================

template<typename MeshType,
		 template <typename> class GraphCutterType,
		 template <typename> class MeshPartBuilderType>
void
MeshPartitionTool<MeshType,
				  GraphCutterType,
				  MeshPartBuilderType>::setup(
		meshPtr_Type& mesh,
        boost::shared_ptr<Epetra_Comm>& comm,
        Teuchos::ParameterList& parameters)
{
    M_comm = comm;
    M_myPID = M_comm->MyPID();
    M_parameters = parameters;
    M_originalMesh = mesh;
    M_meshPart.reset(new mesh_Type);
    M_graphCutter.reset(new graphCutter_Type(M_originalMesh, M_comm,
    										 M_parameters));
    M_meshPartBuilder.reset(new meshPartBuilder_Type(M_originalMesh, M_comm));
}

template<typename MeshType,
		 template <typename> class GraphCutterType,
		 template <typename> class MeshPartBuilderType>
void MeshPartitionTool<MeshType,
					   GraphCutterType,
					   MeshPartBuilderType>::run()
{
	if (!M_myPID)
	{
		std::cout << "Partitioning mesh graph ..." << std::endl;
	}
    M_graphCutter->run();
    M_myElements = M_graphCutter->getPartition(M_myPID);

	if (!M_myPID)
	{
		std::cout << "Building mesh parts ..." << std::endl;
	}
	M_meshPartBuilder->run(M_meshPart, M_myElements);

	if (!M_myPID)
	{
		std::cout << "Mesh partition complete." << std::endl;
	}
    // Destroy the graph partitioner to clear memory
    M_graphCutter.reset();
    // Destroy the mesh part builder to clear memory
    M_meshPartBuilder.reset();
    // Release the pointer to the original uncut mesh
    M_originalMesh.reset();
}

template<typename MeshType,
		 template <typename> class GraphCutterType,
		 template <typename> class MeshPartBuilderType>
void
MeshPartitionTool<MeshType,
				  GraphCutterType,
				  MeshPartBuilderType>::showMe(std::ostream& output) const
{
    std::cout << "Sorry, this method is not implemented, yet." << std::endl
              << "We appreciate your interest." << std::endl
              << "Check back in a bit!" << std::endl;
}

} // namespace LifeV

#endif // MESH_PARTITION_TOOL_ONLINE_H
