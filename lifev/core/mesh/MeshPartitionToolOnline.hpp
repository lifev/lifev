//@HEADER
/*
*******************************************************************************

Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
Copyright (C) 2010, 2011 EPFL, Politecnico di Milano, Emory University

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
#include <parmetis.h>
#include <Epetra_MpiComm.h>
#include <Teuchos_ParameterList.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifecore/LifeV.hpp>
#include <life/lifefem/DOF.hpp>
#include <life/lifemesh/RegionMesh3D.hpp>
#include <life/lifemesh/GhostEntityData.hpp>

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
template<typename MeshType, template <typename> class GraphPartitionToolType>
class MeshPartitionToolOnline
{
public:
    //! @name Public Types
    //@{
    typedef MeshType mesh_Type;
    typedef GraphPartitionToolType<MeshType>               graphPartitionTool_Type;
    typedef boost::shared_ptr<mesh_Type>                   meshPtr_Type;
    //! Container for the ghost data
    typedef std::vector <GhostEntityData>                  GhostEntityDataContainer_Type;
    //! Map processor -> container for the ghost data
    typedef std::map <UInt, GhostEntityDataContainer_Type> GhostEntityDataMap_Type;
    //@}

    //! \name Constructors & Destructors
    //@{
    //! Default empty constructor
    MeshPartitionToolOnline() {}

    //! Constructor
    /*!
     * TODO: Write description
    */
    MeshPartitionToolOnline (meshPtr_Type& mesh,
                             boost::shared_ptr<Epetra_Comm>& comm,
                             Teuchos::ParameterList& parameters);

    //! Empty destructor
    virtual ~MeshPartitionToolOnline() {}
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

    //! Executes the graph partitioning
    /*!
      Executes the graph partitioning
    */
    void partitionGraph();

    //! Builds the mesh object corresponding to the partition
    /*!
      Builds the mesh object corresponding to the partition
    */
    void buildMeshPartition();

    //! Releases the original unpartitioned mesh
    /*!
      Releases the unpartitioned mesh so that it can be deleted, freeing A LOT of memory
      in some cases.
    */
    void releaseUnpartitionedMesh() {M_originalMesh.reset();}

    //! This method performs all the steps for the mesh and graph partitioning
    void run();

    //! Prints information about the state (data) of the object
    void showMe(std::ostream& output = std::cout) const;
    //@}

    //! \name Get Methods
    //@{
    //! Return a shared pointer to the mesh partition
    const meshPtr_Type& meshPartition() const {return M_meshPartition;}
    //! Return a reference to M_ghostDataMap
    const GhostEntityDataMap_Type&  ghostDataMap() const {return M_ghostDataMap;}
    //@}

private:
    // Private copy constructor and assignment operator are disabled
    MeshPartitionToolOnline(const MeshPartitionToolOnline&);
    MeshPartitionToolOnline& operator=(const MeshPartitionToolOnline&);

    //! Private Methods
    //@{
    //! Construct local mesh
    /*!
      Constructs the data structures for the local mesh partition.
      Updates M_localNodes, M_localEdges, M_localFaces, M_localVolumes,
      M_globalToLocalNode.
    */
    void constructLocalMesh();
    //! Construct nodes
    /*!
      Adds nodes to the partitioned mesh object. Updates M_nBoundaryPoints,
      M_meshPartition.
    */
    void constructNodes();
    //! Construct volumes
    /*!
      Adds volumes to the partitioned mesh object. Updates M_globalToLocalVolume,
      M_meshPartition.
    */
    void constructVolumes();
    //! Construct edges
    /*!
      Adds edges to the partitioned mesh object. Updates M_nBoundaryEdges,
      M_meshPartition.
    */
    void constructEdges();
    //! Construct faces
    /*!
      Adds faces to the partitioned mesh object. Updates M_nBoundaryFaces,
      M_meshPartition.
    */
    void constructFaces();
    //! Final setup of local mesh
    /*!
      Updates the partitioned mesh object data members after adding the mesh
      elements (nodes, edges, faces, volumes).
      Updates M_meshPartition.
    */
    void finalSetup();
    //@}

    //! Private Data Members
    //@{
    boost::shared_ptr<Epetra_Comm>             M_comm;
    Int                                        M_myPID;
    UInt                                       M_nBoundaryPoints;
    UInt                                       M_nBoundaryEdges;
    UInt                                       M_nBoundaryFaces;
    UInt                                       M_elementNodes;
    UInt                                       M_elementFaces;
    UInt                                       M_elementEdges;
    UInt                                       M_faceNodes;
    std::vector<Int>                           M_localNodes;
    std::set<Int>                              M_localEdges;
    std::set<Int>                              M_localFaces;
    std::vector<Int>                           M_localVolumes;
    std::map<Int, Int>                         M_globalToLocalNode;
    std::map<Int, Int>                         M_globalToLocalVolume;
    boost::shared_ptr<std::vector<Int> >       M_myElements;
    Teuchos::ParameterList                     M_parameters;
    meshPtr_Type                               M_originalMesh;
    meshPtr_Type                               M_meshPartition;
    boost::shared_ptr<graphPartitionTool_Type> M_graphPartitionTool;
    GhostEntityDataMap_Type                    M_ghostDataMap;
    //@}
}; // class MeshPartitionToolOnline

//
// IMPLEMENTATION
//

// =================================
// Constructors and destructor
// =================================

template<typename MeshType, template <typename> class GraphPartitionToolType>
MeshPartitionToolOnline<MeshType, GraphPartitionToolType>::MeshPartitionToolOnline(meshPtr_Type& mesh,
                                                                                   boost::shared_ptr<Epetra_Comm>& comm,
                                                                                   Teuchos::ParameterList& parameters) :
    M_comm(comm),
    M_myPID(M_comm->MyPID()),
    M_parameters(parameters),
    M_originalMesh(mesh),
    M_meshPartition(new MeshType),
    M_graphPartitionTool(new graphPartitionTool_Type(M_originalMesh, M_comm, M_parameters))
{
    run();
}

// =================================
// Public methods
// =================================

template<typename MeshType, template <typename> class GraphPartitionToolType>
void MeshPartitionToolOnline<MeshType, GraphPartitionToolType>::setup(meshPtr_Type& mesh,
                                                                      boost::shared_ptr<Epetra_Comm>& comm,
                                                                      Teuchos::ParameterList& parameters)
{
    M_comm = comm;
    M_myPID = M_comm->MyPID();
    M_parameters = parameters;
    M_originalMesh = mesh;
    M_meshPartition.reset(new mesh_Type);
    M_graphPartitionTool.reset(new graphPartitionTool_Type(M_originalMesh, M_comm, M_parameters));
}

template<typename MeshType, template <typename> class GraphPartitionToolType>
void MeshPartitionToolOnline<MeshType, GraphPartitionToolType>::partitionGraph()
{
    M_graphPartitionTool->run();
    M_myElements = M_graphPartitionTool->getPartition(M_myPID);
}

template<typename MeshType, template <typename> class GraphPartitionToolType>
void MeshPartitionToolOnline<MeshType, GraphPartitionToolType>::buildMeshPartition()
{
    M_elementNodes = MeshType::ElementShape::S_numVertices;
    M_elementFaces = MeshType::ElementShape::S_numFaces;
    M_elementEdges = MeshType::ElementShape::S_numEdges;
    M_faceNodes    = MeshType::BElementShape::S_numVertices;

    constructLocalMesh();
    constructNodes();
    constructVolumes();
    constructEdges();

    // new faces can be built only after all local volumes are complete in order to get proper ghost faces data
    M_comm->Barrier();

    constructFaces();

    finalSetup();
}

template<typename MeshType, template <typename> class GraphPartitionToolType>
void MeshPartitionToolOnline<MeshType, GraphPartitionToolType>::run()
{
    partitionGraph();
    buildMeshPartition();
    releaseUnpartitionedMesh();

    // Destroy the graph partitioner to clear memory
    M_graphPartitionTool.reset();
}

template<typename MeshType, template <typename> class GraphPartitionToolType>
void MeshPartitionToolOnline<MeshType, GraphPartitionToolType>::showMe(std::ostream& output) const
{
    std::cout << "Sorry, this method is not implemented, yet." << std::endl
              << "We appreciate your interest." << std::endl
              << "Check back in a bit!" << std::endl;
}

// =================================
// Private methods
// =================================

template<typename MeshType, template <typename> class GraphPartitionToolType>
void MeshPartitionToolOnline<MeshType, GraphPartitionToolType>::constructLocalMesh()
{
    if (!M_myPID)
    {
        std::cout << "Building local mesh ..." << std::endl;
    }

    std::map<Int, Int>::iterator  im;
    std::set<Int>::iterator       is;

    Int count = 0;
    UInt ielem;
    UInt inode;

    count = 0;
    // cycle on local element's ID

    const std::vector<Int>& myElements = *M_myElements;
    for (UInt jj = 0; jj < myElements.size(); ++jj)
    {
        ielem = myElements[jj];
        M_localVolumes.push_back(ielem);

        // cycle on element's nodes
        for (UInt ii = 0; ii < M_elementNodes; ++ii)
        {
            inode = M_originalMesh->volume(ielem).point(ii).id();
            im    = M_globalToLocalNode.find(inode);

            // if the node is not yet present in the list of local nodes, then add it
            if (im == M_globalToLocalNode.end())
            {
                M_globalToLocalNode.insert(std::make_pair(inode, count));
                ++count;
                // store here the global numbering of the node
                M_localNodes.push_back(M_originalMesh->volume(ielem).point(ii).id());
            }
        }

        // cycle on element's edges
        for (UInt ii = 0; ii < M_elementEdges; ++ii)
        {
            // store here the global numbering of the edge
            M_localEdges.insert(M_originalMesh->localEdgeId(ielem, ii));
        }

        // cycle on element's faces
        for (UInt ii = 0; ii < M_elementFaces; ++ii)
        {
            // store here the global numbering of the face
            M_localFaces.insert(M_originalMesh->localFaceId(ielem, ii));
        }
    }
}

template<typename MeshType, template <typename> class GraphPartitionToolType>
void MeshPartitionToolOnline<MeshType, GraphPartitionToolType>::constructNodes()
{
    UInt inode;
    std::vector<Int>::iterator it;

    M_nBoundaryPoints = 0;
    M_meshPartition->pointList.reserve(M_localNodes.size());
    // guessing how many boundary points on this processor.
    M_meshPartition->_bPoints.reserve(M_originalMesh->numBPoints() * M_localNodes.size() /
                                      M_originalMesh->numBPoints());
    inode = 0;
    typename MeshType::point_Type *pp = 0;

    // loop in the list of local nodes:
    // in this loop inode is the local numbering of the points
    for (it = M_localNodes.begin(); it != M_localNodes.end(); ++it, ++inode)
    {
        typename MeshType::point_Type point = 0;

        // create a boundary point in the local mesh, if needed
        bool boundary = M_originalMesh->isBoundaryPoint(*it);
        if (boundary)
        {
            ++M_nBoundaryPoints;
        }

        pp = &(M_meshPartition->addPoint(boundary));
        *pp = M_originalMesh->point( *it );

        pp->setLocalId( inode );

        M_meshPartition->localToGlobalNode().insert(std::make_pair( pp->localId(), pp->id() ));
        M_meshPartition->globalToLocalNode().insert(std::make_pair( pp->id(), pp->localId() ));
    }
}

template<typename MeshType, template <typename> class GraphPartitionToolType>
void MeshPartitionToolOnline<MeshType, GraphPartitionToolType>::constructVolumes()
{
    Int count;
    std::map<Int, Int>::iterator im;
    std::vector<Int>::iterator it;
    count = 0;
    UInt inode;

    typename MeshType::VolumeType * pv = 0;

    M_meshPartition->volumeList.reserve(M_localVolumes.size());

    // loop in the list of local elements
    // CAREFUL! in this loop inode is the global numbering of the points
    // We insert the local numbering of the nodes in the local volume list
    for (it = M_localVolumes.begin(); it != M_localVolumes.end(); ++it, ++count)
    {
        pv = &(M_meshPartition->addVolume());
        *pv = M_originalMesh->volume( *it );
        pv->setLocalId( count );

        M_globalToLocalVolume.insert(std::make_pair( pv->id(), pv->localId() ) );

        for (ID id = 0; id < M_elementNodes; ++id)
        {
            inode = M_originalMesh->volume(*it).point(id).id();
            // im is an iterator to a map element
            // im->first is the key (i. e. the global ID "inode")
            // im->second is the value (i. e. the local ID "count")
            im = M_globalToLocalNode.find(inode);
            pv->setPoint(id, M_meshPartition->pointList( (*im).second ));
        }
    }
}

template<typename MeshType, template <typename> class GraphPartitionToolType>
void MeshPartitionToolOnline<MeshType, GraphPartitionToolType>::constructEdges()
{
    Int count;
    std::map<Int, Int>::iterator im;
    std::set<Int>::iterator is;

    typename MeshType::EdgeType * pe;
    UInt inode;
    count = 0;

    M_nBoundaryEdges = 0;
    M_meshPartition->edgeList.reserve(M_localEdges.size());

    // loop in the list of local edges
    for (is = M_localEdges.begin(); is != M_localEdges.end(); ++is, ++count)
    {
        // create a boundary edge in the local mesh, if needed
        bool boundary = (M_originalMesh->isBoundaryEdge(*is));
        if (boundary)
        {
            // create a boundary edge in the local mesh, if needed
            ++M_nBoundaryEdges;
        }

        pe = &(M_meshPartition->addEdge(boundary));
        *pe = M_originalMesh->edge( *is );

        pe->setLocalId(count);

        for (ID id = 0; id < 2; ++id)
        {
            inode = M_originalMesh->edge(*is).point(id).id();
            // im is an iterator to a map element
            // im->first is the key (i. e. the global ID "inode")
            // im->second is the value (i. e. the local ID "count")
            im = M_globalToLocalNode.find(inode);
            pe->setPoint(id, M_meshPartition->pointList((*im).second));
        }
    }
}

template<typename MeshType, template <typename> class GraphPartitionToolType>
void MeshPartitionToolOnline<MeshType, GraphPartitionToolType>::constructFaces()
{
    Int count;
    std::map<Int, Int>::iterator im;
    std::set<Int>::iterator      is;

    typename MeshType::FaceType * pf = 0;

    UInt inode;
    count = 0;

    M_nBoundaryFaces = 0;
    M_meshPartition->faceList.reserve(M_localFaces.size());

    // loop in the list of local faces
    for (is = M_localFaces.begin(); is != M_localFaces.end(); ++is, ++count)
    {
        // create a boundary face in the local mesh, if needed
        bool boundary = (M_originalMesh->isBoundaryFace(*is));
        if (boundary)
        {
            ++M_nBoundaryFaces;
        }

        Int elem1 = M_originalMesh->face(*is).firstAdjacentElementIdentity();
        Int elem2 = M_originalMesh->face(*is).secondAdjacentElementIdentity();

        // find the mesh elements adjacent to the face
        im =  M_globalToLocalVolume.find(elem1);

        ID localElem1;

        if (im == M_globalToLocalVolume.end())
        {
            localElem1 = NotAnId;
        }
        else
        {
            localElem1 = (*im).second;
        }

        im =  M_globalToLocalVolume.find(elem2);

        ID localElem2;
        if (im == M_globalToLocalVolume.end())
        {
            localElem2 = NotAnId;
        }
        else
        {
            localElem2 = (*im).second;
        }

        pf =  &(M_meshPartition->addFace(boundary));
        *pf = M_originalMesh->face( *is );

        pf->setLocalId( count );

        // true if we are on a subdomain border
        if ( !boundary && ( localElem1 == NotAnId || localElem2 == NotAnId ) )
        {
            // set the flag for faces on the subdomain border
            pf->setFlag( EntityFlags::SUBDOMAIN_INTERFACE );

            // build GhostEntityData
            GhostEntityData ghostFace;
            ghostFace.localFaceId = pf->localId();

            // TODO: make this work. Zoltan's DD didn't cut it. Find another way.
            //// // set the ghostElem to be searched on other subdomains
            //// ID ghostElem = ( localElem1 == NotAnId ) ? elem1 : elem2;
            //// // find which process holds the facing element
            //// Int ghostProc ( M_me );
            //// for ( Int proc = 0; proc < M_comm->NumProc(); proc++ )
            //// {
            ////     if ( proc != M_me )
            ////     {
            ////         std::vector<Int>::const_iterator ghostIt =
            ////                         std::find ( (*M_elementDomains)[ proc ].begin(), (*M_elementDomains)[ proc ].end(), ghostElem );
            ////         if ( ghostIt != ( (*M_elementDomains)[ proc ] ).end() )
            ////         {
            ////             // we have found the proc storing the element
            ////             ghostProc = proc;
            ////             // we can get its local id
            ////             ghostFace.ghostElementLocalId = *ghostIt;
            ////             // TODO: the local face id is the same of the original mesh ?!
            ////             ghostFace.ghostElementPosition = M_originalMesh->face(*is).secondAdjacentElementPosition();
            ////             break;
            ////         }
            ////     }
            //// }
            //// // check that the ghost element is found on another proc ( this test is acceptable only for online partitioning )
            //// ASSERT ( ghostProc != M_me || M_serialMode, "ghost face not found" );
            //// M_ghostDataMap[ ghostProc ].push_back( ghostFace );
        }

        ASSERT((localElem1 != NotAnId)||(localElem2 != NotAnId),"A hanging face in mesh partitioner!");

        if (localElem1 == NotAnId)
         {
             pf->firstAdjacentElementIdentity()  = localElem2;
             pf->firstAdjacentElementPosition()  = M_originalMesh->face(*is).secondAdjacentElementPosition();
             pf->secondAdjacentElementIdentity() = NotAnId;
             pf->secondAdjacentElementPosition() = NotAnId;
             pf->reversePoints();
         }
        else
        {
            pf->firstAdjacentElementIdentity()  = localElem1;
            pf->firstAdjacentElementPosition()  = M_originalMesh->face(*is).firstAdjacentElementPosition();
            pf->secondAdjacentElementIdentity() = localElem2;
            pf->secondAdjacentElementPosition() =
                            localElem2 != NotAnId ?
                                                   M_originalMesh->face(*is).secondAdjacentElementPosition():
                                                   NotAnId;
        }

        for (ID id = 0; id < M_originalMesh->face(*is).S_numLocalVertices; ++id)
        {
            inode = pf->point(id).id();
            im = M_globalToLocalNode.find(inode);
            pf->setPoint(id, M_meshPartition->pointList((*im).second));
        }
    }
    M_meshPartition->setLinkSwitch("HAS_ALL_FACES");
    M_meshPartition->setLinkSwitch("FACES_HAVE_ADIACENCY");
}

template<typename MeshType, template <typename> class GraphPartitionToolType>
void MeshPartitionToolOnline<MeshType, GraphPartitionToolType>::finalSetup()
{
    UInt nVolumes = M_localVolumes.size();
    UInt nNodes   = M_localNodes.size();
    UInt nEdges   = M_localEdges.size();
    UInt nFaces   = M_localFaces.size();

    M_meshPartition->setMaxNumPoints (nNodes, true);
    M_meshPartition->setMaxNumEdges  (nEdges, true);
    M_meshPartition->setMaxNumFaces  (nFaces, true);
    M_meshPartition->setMaxNumVolumes( nVolumes, true);

    M_meshPartition->setMaxNumGlobalPoints (M_originalMesh->numPoints());
    M_meshPartition->setNumGlobalVertices  (M_originalMesh->numPoints());
    M_meshPartition->setMaxNumGlobalEdges  (M_originalMesh->numEdges());
    M_meshPartition->setMaxNumGlobalFaces  (M_originalMesh->numFaces());

    M_meshPartition->setMaxNumGlobalVolumes(M_originalMesh->numVolumes());
    M_meshPartition->setNumBFaces    (M_nBoundaryFaces);

    M_meshPartition->setNumBPoints   (M_nBoundaryPoints);
    M_meshPartition->setNumBEdges    (M_nBoundaryEdges);

    M_meshPartition->setNumVertices (nNodes );
    M_meshPartition->setNumBVertices(M_nBoundaryPoints);

    M_meshPartition->updateElementEdges();

    M_meshPartition->updateElementFaces();
}

}
#endif // MESH_PARTITION_TOOL_ONLINE_H
