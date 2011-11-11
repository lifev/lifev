//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

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
    @file GhostHandler.hpp
    @brief class to manage ghost data across procs

    @author Antonio Cervone <ant.cervone@gmail.com>

    @date 27-10-2011
*/

#ifndef _GHOSTHANDLER_HPP_
#define _GHOSTHANDLER_HPP_

#include <boost/shared_ptr.hpp>

#include <life/lifemesh/NeighborMarker.hpp>
#include <life/lifearray/MapEpetra.hpp>
#include <life/lifearray/VectorEpetra.hpp>

namespace LifeV {

template <typename Mesh>
class GhostHandler
{
public:

    //! @name Public Types
    //@{

    typedef Mesh mesh_Type;
    typedef boost::shared_ptr<mesh_Type> mesh_PtrType;
    typedef Epetra_Comm comm_Type;
    typedef boost::shared_ptr<comm_Type> comm_PtrType;
    typedef std::set<ID> neighborList_Type;
    typedef std::map< ID, neighborList_Type > neighborMap_Type;
    typedef MapEpetra map_Type;
    typedef boost::shared_ptr<map_Type> map_PtrType;

    //@}

    //! @name Constructors & Destructors
    //@{

    //! Constructor
    GhostHandler( mesh_PtrType fullMesh,
                  mesh_PtrType localMesh,
                  map_Type & map,
                  comm_PtrType const & comm );

    //! Destructor
    ~GhostHandler(){}

    //@}

    //! @name Get Methods
    //@{

    mesh_Type const & fullMesh() { return *M_fullMesh; }
    mesh_Type const & localMesh() { return *M_localMesh; }
    map_Type const & map() { return M_map; }
    neighborMap_Type const & nodeNodeNeighborsMap() { return M_nodeNodeNeighborsMap; }
    neighborMap_Type const & nodeElementNeighborsMap() { return M_nodeElementNeighborsMap; }

    //@}

    //! @name General Methods
    //@{

    //! setup
    void setUp();

    //! clean
    void clean();

    //! create node node neighbors on node markers
    void createNodeNeighbors();

    //! create node node neighbors map
    void createNodeNodeNeighborsMap();

    //! create node element neighbors map
    void createNodeElementNeighborsMap();

    //! create ghost map
    map_Type & ghostMapOnElementsP0();

    //! create ghost map
    map_Type & ghostMapOnElementsP1( UInt overlap );

    //! create ghost map
    map_Type & ghostMapOnNodes();

    //! create ghost map
    map_Type & ghostMapOnNodes( UInt overlap );

    //! showMe method
    void showMe( bool const verbose = false, std::ostream & out = std::cout );

    //@}

protected:

    //! @name Ghost Maps
    //@{

    map_PtrType M_ghostMapOnElementsP0;
    map_PtrType M_ghostMapOnElementsP1;
    map_PtrType M_ghostMapOnNodes;

    //@}

    //! @name Protected Members
    //@{

    mesh_PtrType M_fullMesh;
    mesh_PtrType M_localMesh;
    map_Type & M_map;
    comm_PtrType const M_comm;
    UInt M_me;

    neighborMap_Type M_nodeNodeNeighborsMap;
    neighborMap_Type M_nodeElementNeighborsMap;

    const bool M_verbose;

    //@}
};

template <typename Mesh>
GhostHandler<Mesh>::GhostHandler( mesh_PtrType fullMesh,
                                  mesh_PtrType localMesh,
                                  map_Type & map,
                                  comm_PtrType const & comm ):
    M_fullMesh ( fullMesh ),
    M_localMesh ( localMesh ),
    M_map ( map ),
    M_comm ( comm ),
    M_me ( comm->MyPID() ),
    M_nodeNodeNeighborsMap(),
    M_nodeElementNeighborsMap(),
    M_verbose( !M_me )
{
}

template <typename Mesh>
void GhostHandler<Mesh>::setUp()
{
    this->createNodeNodeNeighborsMap();
    this->createNodeElementNeighborsMap();
}

template <typename Mesh>
void GhostHandler<Mesh>::clean()
{
    clearVector( M_nodeNodeNeighborsMap );
    clearVector( M_nodeElementNeighborsMap );
}

//! this routine generates node neighbors for the given mesh
/*! the routine assumes that the mesh is not yet partitioned or reordered
 *  (i.e. the local id and the global id are the same).
 *  if this is not true the method should be changed to use a more
 *  expensive STL find on the mesh points to get the correct point that has
 *  the given global id or construct a globalToLocal map beforehand.
 */
template <typename Mesh>
void GhostHandler<Mesh>::createNodeNeighbors()
{
    // TODO: ASSERT_COMPILE_TIME that MeshType::pointMarker == NeighborMarker
    // this guarantees that the nodeNeighbors structure is available.

    // generate node neighbors by watching edges
    // note: this can be based also on faces or volumes
    for ( UInt ie = 0; ie < M_fullMesh->numEdges(); ie++ )
    {
        ID id0 = M_fullMesh->edge( ie ).point( 0 ).id();
        ID id1 = M_fullMesh->edge( ie ).point( 1 ).id();

        ASSERT ( M_fullMesh->point( id0 ).id() == id0 && M_fullMesh->point( id1 ).id() == id1,
                 "the mesh has been reordered, the point must be found" );

        // fill fullMesh points neighboring
        M_fullMesh->point( id0 ).nodeNeighbors().insert( id1 );
        M_fullMesh->point( id1 ).nodeNeighbors().insert( id0 );
    }

    // update localMesh points
    for ( UInt ip = 0; ip < M_localMesh->numPoints(); ip++ )
    {
        M_localMesh->point( ip ).nodeNeighbors() = M_fullMesh->point( M_localMesh->point( ip ).id() ).nodeNeighbors();
    }
}

template <typename Mesh>
void GhostHandler<Mesh>::createNodeNodeNeighborsMap()
{
    // generate node neighbors by watching edges
    // note: this can be based also on faces or volumes
    for ( UInt ie = 0; ie < M_fullMesh->numEdges(); ie++ )
    {
        ID id0 = M_fullMesh->edge( ie ).point( 0 ).id();
        ID id1 = M_fullMesh->edge( ie ).point( 1 ).id();

        ASSERT ( M_fullMesh->point( id0 ).id() == id0 && M_fullMesh->point( id1 ).id() == id1,
                 "the mesh has been reordered, the point must be found" );

        M_nodeNodeNeighborsMap[ id0 ].insert ( id1 );
        M_nodeNodeNeighborsMap[ id1 ].insert ( id0 );
    }
}

template <typename Mesh>
void GhostHandler<Mesh>::createNodeElementNeighborsMap()
{
    // generate element neighbors by cycling on elements
    for ( UInt ie = 0; ie < M_fullMesh->numElements(); ie++ )
    {
        ASSERT ( M_fullMesh->element( ie ).id() == ie,
                 "the mesh has been reordered, the point must be found" );

        for ( UInt k = 0; k < Mesh::VolumeType::S_numPoints; k++ )
        {
            ID id ( M_fullMesh->element( ie ).point( k ).id() );
            M_nodeElementNeighborsMap[ id ].insert ( ie );
        }
    }
}

template <typename Mesh>
typename GhostHandler<Mesh>::map_Type & GhostHandler<Mesh>::ghostMapOnNodes()
{
    // if the map has already been created, return it
    if ( M_ghostMapOnNodes ) return *M_ghostMapOnNodes;

    if ( M_verbose ) std::cout << " GH- ghostMapOnNodes()" << std::endl;

    // check that the nodeNeighbors have been created
    if ( M_localMesh->point( 0 ).nodeNeighbors().empty()  )
    {
        if ( M_verbose ) std::cerr << "the nodeNeighbors are empty, will be generated now" << std::endl;
        this->createNodeNeighbors();
    }

    // create map
    M_ghostMapOnNodes.reset ( new map_Type() );
    map_Type & ghostMap ( *M_ghostMapOnNodes );

    // use the same Unique map and comm of the original map
    ghostMap.setMap( M_map.map( Unique ), Unique );
    ghostMap.setComm( M_comm );

    // use a set to avoid duplicates
    std::set<Int> myGlobalElementsSet;

    // iterate on local mesh points
    // todo: this can start from the repeated map and add only neighbors for SUBDOMAIN_INTERFACE marked nodes
    for ( UInt k = 0; k < M_localMesh->numPoints(); k++ )
    {
        // iterate on each node neighborhood
        for ( typename mesh_Type::PointMarker::neighborConstIterator_Type neighborIt = M_localMesh->point( k ).nodeNeighbors().begin();
              neighborIt != M_localMesh->point( k ).nodeNeighbors().end(); ++neighborIt )
        {
            myGlobalElementsSet.insert( *neighborIt );
        }
    }

    std::vector<Int> myGlobalElements( myGlobalElementsSet.begin(), myGlobalElementsSet.end() );

    // generate map
    map_Type::map_ptrtype repeatedMap ( new Epetra_Map( -1, myGlobalElements.size(), &myGlobalElements[0], 0, *M_comm ) );
    ghostMap.setMap( repeatedMap, Repeated );

    return *M_ghostMapOnNodes;
}

template <typename Mesh>
typename GhostHandler<Mesh>::map_Type & GhostHandler<Mesh>::ghostMapOnNodes( UInt overlap )
{
    // if the map has already been created, return it
    if ( M_ghostMapOnNodes ) return *M_ghostMapOnNodes;

    if ( M_verbose ) std::cout << " GH- ghostMapOnNodes( o )" << std::endl;

    // check that the nodeNodeNeighborsMap has been created
    if ( M_nodeNodeNeighborsMap.empty()  )
    {
        if ( M_verbose ) std::cerr << "the nodeNodeNeighborsMap is empty, will be generated now" << std::endl;
        this->createNodeNodeNeighborsMap();
    }

    // create map
    M_ghostMapOnNodes.reset ( new map_Type() );
    map_Type & ghostMap ( *M_ghostMapOnNodes );

    // use the same Unique map and comm of the original map
    ghostMap.setMap( M_map.map( Unique ), Unique );
    ghostMap.setComm( M_comm );

    Int*          pointer;
    std::set<Int> myGlobalElementsSet, myOriginalElementsSet;;

    // get all elements from the repeated map
    pointer = M_map.map( Repeated )->MyGlobalElements();
    for ( Int ii = 0; ii < M_map.map( Repeated )->NumMyElements(); ++ii, ++pointer )
    {
        myOriginalElementsSet.insert( *pointer );
    }

    // todo: optimize this!!
    // 1: work only on the boundary
    // 2: copy back only if necessary
    // repeat on actual nodes to expand overlap
    for ( UInt i = 0; i < overlap; i++ )
    {
        // iterate on points adding all neighbors
        for ( std::set<Int>::const_iterator nodeIt = myOriginalElementsSet.begin();
                        nodeIt != myOriginalElementsSet.end(); ++nodeIt )
        {
            // iterate on each node neighborhood
            for ( neighborList_Type::const_iterator neighborIt = M_nodeNodeNeighborsMap[ *nodeIt ].begin();
                            neighborIt != M_nodeNodeNeighborsMap[ *nodeIt ].end(); ++neighborIt )
            {
                myGlobalElementsSet.insert( *neighborIt );
            }
        }
        myOriginalElementsSet = myGlobalElementsSet;
    }

    std::vector<Int> myGlobalElements( myGlobalElementsSet.begin(), myGlobalElementsSet.end() );

    // generate map
    map_Type::map_ptrtype repeatedMap ( new Epetra_Map( -1, myGlobalElements.size(), &myGlobalElements[0], 0, *M_comm ) );
    ghostMap.setMap( repeatedMap, Repeated );

    return *M_ghostMapOnNodes;
}

template <typename Mesh>
typename GhostHandler<Mesh>::map_Type & GhostHandler<Mesh>::ghostMapOnElementsP0()
{
    // if the map has already been created, return it
    if ( M_ghostMapOnElementsP0 ) return *M_ghostMapOnElementsP0;

    if ( M_verbose ) std::cout << " GH- ghostMapOnElementsP0()" << std::endl;

    // create the map
    M_ghostMapOnElementsP0.reset ( new map_Type() );
    map_Type & ghostMap ( *M_ghostMapOnElementsP0 );

    // use the same Unique map and comm of the original map
    ghostMap.setMap( M_map.map( Unique ), Unique );
    ghostMap.setComm( M_comm );

    Int*          pointer;
    std::set<Int> map;

    // get all elements from the repeated map
    pointer = M_map.map( Repeated )->MyGlobalElements();
    for ( Int ii = 0; ii < M_map.map( Repeated )->NumMyElements(); ++ii, ++pointer )
    {
        map.insert( *pointer );
    }

    // add all facing elements
    std::vector<ID> facesOnSubdInt = M_localMesh->faceList.extractElementsWithFlag(
                    EntityFlags::SUBDOMAIN_INTERFACE, &Flag::testOneSet );
    for ( ID faceId = 0; faceId < facesOnSubdInt.size(); faceId++ )
    {
        map.insert( M_localMesh->faceList[ facesOnSubdInt [ faceId ] ].secondAdjacentElementIdentity() );
    }

    // convert unique list to vector to assure continuity in memorization
    std::vector<Int> myGlobalElements ( map.begin(), map.end() );

    // generate map
    map_Type::map_ptrtype repeatedMap ( new Epetra_Map( -1, myGlobalElements.size(), &myGlobalElements[0], 0, *M_comm ) );
    ghostMap.setMap( repeatedMap, Repeated );

    return *M_ghostMapOnElementsP0;
}

template <typename Mesh>
typename GhostHandler<Mesh>::map_Type & GhostHandler<Mesh>::ghostMapOnElementsP1( UInt overlap )
{
    // if the map has already been created, return it
    if ( M_ghostMapOnElementsP1 ) return *M_ghostMapOnElementsP1;

    if ( M_verbose ) std::cout << " GH- ghostMapOnElementsP1()" << std::endl;

    // check that the nodeElementNeighborsMap has been created
    if ( M_nodeElementNeighborsMap.empty()  )
    {
        if ( M_verbose ) std::cerr << "the nodeElementNeighborsMap is empty, will be generated now" << std::endl;
        this->createNodeElementNeighborsMap();
    }

    // create the map
    M_ghostMapOnElementsP1.reset ( new map_Type() );
    map_Type & ghostMap ( *M_ghostMapOnElementsP1 );

    // use the same Unique map and comm of the original map
    ghostMap.setMap( M_map.map( Unique ), Unique );
    ghostMap.setComm( M_comm );

    Int*          pointer;
    std::set<Int> map;

    // get all elements from the repeated map
    pointer = M_map.map( Repeated )->MyGlobalElements();
    for ( Int ii = 0; ii < M_map.map( Repeated )->NumMyElements(); ++ii, ++pointer )
    {
        map.insert( *pointer );
    }

    // add all elements with a node on SUBDOMAIN_INTERFACE
    std::vector<ID> pointsOnSubdInt = M_localMesh->pointList.extractElementsWithFlag(
                    EntityFlags::SUBDOMAIN_INTERFACE, &Flag::testOneSet );
    // must work on global IDs since added elements are not on localMesh
    for ( UInt i = 0; i < pointsOnSubdInt.size(); i++)
        pointsOnSubdInt[i] = M_localMesh->pointList( pointsOnSubdInt[i] ).id();

    std::vector<ID> addedPoints;

    for ( UInt n = 0; n < overlap; n++ )
    {
        for ( std::vector<ID>::const_iterator globalId = pointsOnSubdInt.begin();
                        globalId != pointsOnSubdInt.end(); ++globalId )
        {
            // iterate on each node neighborhood
            for ( neighborList_Type::const_iterator neighborIt = M_nodeElementNeighborsMap[ *globalId ].begin();
                            neighborIt != M_nodeElementNeighborsMap[ *globalId ].end(); ++neighborIt )
            {
                std::pair<std::set<Int>::iterator, bool> isInserted = map.insert( *neighborIt );
                if ( isInserted.second )
                {
                    typename mesh_Type::VolumeType const & elem = M_fullMesh->element ( *neighborIt );
                    for ( UInt elemPoint = 0; elemPoint < mesh_Type::VolumeType::S_numPoints; elemPoint++ )
                    {
                        // TODO exclude already included nodes
                        addedPoints.push_back ( elem.point( elemPoint ).id() );
                    }
                }
            }
        }
        // TODO: this must be done only if overlap > 1
        pointsOnSubdInt = addedPoints;
    }

    // convert unique list to vector to assure continuity in memorization
    std::vector<Int> myGlobalElements ( map.begin(), map.end() );

    // generate map
    map_Type::map_ptrtype repeatedMap ( new Epetra_Map( -1, myGlobalElements.size(), &myGlobalElements[0], 0, *M_comm ) );
    ghostMap.setMap( repeatedMap, Repeated );

    return *M_ghostMapOnElementsP1;
}

template <typename Mesh>
void GhostHandler<Mesh>::showMe( bool const /*verbose*/, std::ostream & out )
{
    out << "GhostHandler<Mesh>showMe" << std::endl;
    out << "M_nodeNodeNeighborsMap" << std::endl;
    for ( neighborMap_Type::const_iterator it = M_nodeNodeNeighborsMap.begin();
                    it != M_nodeNodeNeighborsMap.end(); ++it )
    {
        out << it->first << "> ";
        for ( neighborList_Type::const_iterator nIt = it->second.begin(); nIt != it->second.end(); ++nIt )
        {
            out << *nIt << " ";
        }
        out << std::endl;
    }
    out << "M_nodeElementNeighborsMap" << std::endl;
    for ( neighborMap_Type::const_iterator it = M_nodeElementNeighborsMap.begin();
                    it != M_nodeElementNeighborsMap.end(); ++it )
    {
        out << it->first << "> ";
        for ( neighborList_Type::const_iterator nIt = it->second.begin(); nIt != it->second.end(); ++nIt )
        {
            out << *nIt << " ";
        }
        out << std::endl;
    }
}


}

#endif /* _GHOSTHANDLER_HPP_ */
