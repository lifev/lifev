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

#include <lifev/core/mesh/NeighborMarker.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/GraphCutterBase.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-pedantic"

#ifdef HAVE_HDF5
#include <EpetraExt_HDF5.h>
#endif

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"
#pragma GCC diagnostic warning "-pedantic"

namespace LifeV
{

template <typename Mesh>
class GhostHandler
{
public:

    //! @name Public Types
    //@{

    typedef Mesh mesh_Type;
    typedef boost::shared_ptr<mesh_Type> meshPtr_Type;
    typedef Epetra_Comm comm_Type;
    typedef boost::shared_ptr<comm_Type> commPtr_Type;
    typedef std::set<ID> neighbors_Type;
    typedef std::vector<neighbors_Type> neighborList_Type;
    typedef MapEpetra map_Type;
    typedef boost::shared_ptr<map_Type> mapPtr_Type;
    typedef std::map< UInt, mapPtr_Type > mapList_Type;
    typedef std::vector<std::vector<Int> > graph_Type;
    typedef boost::shared_ptr<graph_Type> graphPtr_Type;
    typedef boost::shared_ptr<GraphCutterBase<Mesh> > graphCutterPtr_Type;
    typedef std::vector<int> flag_Type;

    //@}

    //! @name Constructors & Destructors
    //@{
    //! Constructor
    explicit GhostHandler ( commPtr_Type const& comm );

    //! Constructor
    GhostHandler ( meshPtr_Type fullMesh, commPtr_Type const& comm );

    //! Constructor
    GhostHandler ( meshPtr_Type fullMesh,
                   meshPtr_Type localMesh,
                   mapPtr_Type map,
                   commPtr_Type const& comm );

    //! Destructor
    ~GhostHandler() {}

    //@}

    //! @name Get Methods
    //@{

    //! Full mesh getter
    mesh_Type const& fullMesh()
    {
        return *M_fullMesh;
    }

    //! Local mesh getter
    mesh_Type const& localMesh()
    {
        return *M_localMesh;
    }

    //! Standard map getter
    map_Type const& map()
    {
        return *M_map;
    }

    //! Node to node neighbor map
    neighborList_Type const& nodeNodeNeighborsList()
    {
        ASSERT ( !M_nodeNodeNeighborsList.empty(), "M_nodeNodeNeighborsList is empty" );
        return M_nodeNodeNeighborsList;
    }

    //! Node to element neighbor map
    neighborList_Type const& nodeElementNeighborsList()
    {
        ASSERT ( !M_nodeElementNeighborsList.empty(), "M_nodeElementNeighborsList is empty" );
        return M_nodeElementNeighborsList;
    }

    //@}

    //! @name Set Methods
    //@{

    void setVerbose ( const bool& verbose )
    {
        M_verbose = verbose && ( M_me == 0 );
    }

    //@}

    //! @name General Methods
    //@{

    //! setup
    void setUp();

    //! setup
    void setUp (flag_Type const& flags);

    //! release
    void release();

    //! clean
    void clean();

#ifdef HAVE_HDF5
    //! export
    void exportToHDF5 ( std::string const& fileName = "ghostmap", bool const& truncate = true );

    //! import
    void importFromHDF5 ( std::string const& fileName = "ghostmap" );
#endif // HAVE_HDF5

    //! create node node neighbors on node markers
    void createNodeNeighbors();

    //! create node node neighbors map
    void createNodeNodeNeighborsMap();

    //! create node node neighbors map for some specified flag
    void createNodeNodeNeighborsMap (flag_Type const& flags);

    //! create node node neighbors map, where the neighbors are selected over the first circles
    std::set<ID> createCircleNodeNodeNeighborsMap (UInt Ncircles, UInt GlobalID);

    //! create node node neighbors map, where the neighbors are selected within a certain user-defined radius
    std::set<ID> createNodeNodeNeighborsMapWithinRadius (double Radius, UInt GlobalID);

    //! create node edge neighbors map
    void createNodeEdgeNeighborsMap();

    //! create node element neighbors map
    void createNodeElementNeighborsMap();

    //! create ghost map
    map_Type& ghostMapOnNodes();

    //! create ghost map
    map_Type& ghostMapOnNodes ( UInt overlap );

    //! create ghost map
    map_Type& ghostMapOnEdges ( UInt overlap );

    //! create ghost map
    // ghostMapOnElementsCommonFacet
    map_Type& ghostMapOnElementsP0();

    //! create ghost map
    // ghostMapOnElementsCommonNodes
    map_Type& ghostMapOnElementsP1 ( UInt overlap );

    //! fill entityPID
    void fillEntityPID ( graphPtr_Type elemGraph, std::vector<std::vector<UInt> >& entityPID );

    //! create ghost map
    void ghostMapOnElementsP1 ( graphPtr_Type elemGraph, const std::vector<UInt>& entityPID, UInt overlap );

    //! create ghost map
    void ghostMapOnElementsP1 ( graphPtr_Type graph,
                                const std::vector<UInt>& entityPID,
                                UInt overlap,
                                const UInt partIndex);

    //! create ghost map
    map_Type& ghostMapOnNodesMap ( UInt overlap );

    //! create ghost map
    map_Type& ghostMapOnEdgesMap ( UInt overlap );

    //! create ghost map
    map_Type& ghostMapOnElementsP0Map();

    //! create ghost map
    map_Type& ghostMapOnElementsP1Map ( UInt overlap );

    //! showMe method
    void showMe ( bool const verbose = false, std::ostream& out = std::cout );

    //! Check if the point with markerID pointMarker has to be selected
    bool isInside (ID const&   pointMarker, flag_Type const& flags);

    //@}

protected:

    //! @name Ghost Maps
    //@{

    mapPtr_Type M_ghostMapOnNodes;
    mapPtr_Type M_ghostMapOnEdges;
    mapPtr_Type M_ghostMapOnElementsP0;
    mapPtr_Type M_ghostMapOnElementsP1;
    mapList_Type M_ghostMapOnNodesMap;
    mapList_Type M_ghostMapOnEdgesMap;
    mapList_Type M_ghostMapOnElementsP0Map;
    mapList_Type M_ghostMapOnElementsP1Map;

    //@}

    //! @name Protected Members
    //@{

    meshPtr_Type M_fullMesh;
    meshPtr_Type M_localMesh;
    mapPtr_Type const M_map;
    commPtr_Type const M_comm;
    UInt const M_me;

    neighborList_Type M_nodeNodeNeighborsList;
    neighborList_Type M_nodeEdgeNeighborsList;
    neighborList_Type M_nodeElementNeighborsList;

    bool M_verbose;
    std::ofstream M_debugOut;

    //@}
};

template <typename Mesh>
GhostHandler<Mesh>::GhostHandler ( commPtr_Type const& comm ) :
    M_fullMesh(),
    M_localMesh(),
    M_map(),
    M_comm ( comm ),
    M_me ( comm->MyPID() ),
    M_nodeNodeNeighborsList(),
    M_nodeEdgeNeighborsList(),
    M_nodeElementNeighborsList(),
    M_verbose ( 0 ),
#ifdef HAVE_LIFEV_DEBUG
    M_debugOut ( ( "gh." + ( comm->NumProc() > 1 ? boost::lexical_cast<std::string> ( M_me ) : "s" ) + ".out" ).c_str() )
#else
    M_debugOut()
#endif
{
}

template <typename Mesh>
GhostHandler<Mesh>::GhostHandler ( meshPtr_Type fullMesh,
                                   meshPtr_Type localMesh,
                                   mapPtr_Type map,
                                   commPtr_Type const& comm ) :
    M_fullMesh ( fullMesh ),
    M_localMesh ( localMesh ),
    M_map ( map ),
    M_comm ( comm ),
    M_me ( comm->MyPID() ),
    M_nodeNodeNeighborsList(),
    M_nodeEdgeNeighborsList(),
    M_nodeElementNeighborsList(),
    M_verbose ( 0 ),
#ifdef HAVE_LIFEV_DEBUG
    M_debugOut ( ( "gh." + ( comm->NumProc() > 1 ? boost::lexical_cast<std::string> ( M_me ) : "s" ) + ".out" ).c_str() )
#else
    M_debugOut()
#endif
{
}

template <typename Mesh>
GhostHandler<Mesh>::GhostHandler ( meshPtr_Type fullMesh,
                                   commPtr_Type const& comm ) :
    M_fullMesh ( fullMesh ),
    M_comm ( comm ),
    M_me ( comm->MyPID() ),
    M_nodeNodeNeighborsList(),
    M_nodeEdgeNeighborsList(),
    M_nodeElementNeighborsList(),
    M_verbose ( 0 ),
#ifdef HAVE_LIFEV_DEBUG
    M_debugOut ( ( "gh." + ( comm->NumProc() > 1 ? boost::lexical_cast<std::string> ( M_me ) : "s" ) + ".out" ).c_str() )
#else
    M_debugOut()
#endif
{
}

template <typename Mesh>
void GhostHandler<Mesh>::setUp()
{
    this->createNodeNodeNeighborsMap();
    this->createNodeEdgeNeighborsMap();
    this->createNodeElementNeighborsMap();
}

template <typename Mesh>
void GhostHandler<Mesh>::setUp (flag_Type const& flags)
{
    this->createNodeNodeNeighborsMap (flags);
}

template <typename Mesh>
void GhostHandler<Mesh>::release()
{
    M_fullMesh.reset();
    M_localMesh.reset();
}

template <typename Mesh>
void GhostHandler<Mesh>::clean()
{
    clearVector ( M_nodeNodeNeighborsList );
    clearVector ( M_nodeEdgeNeighborsList );
    clearVector ( M_nodeElementNeighborsList );
}

#ifdef HAVE_HDF5

namespace
{

void writeNeighborMap ( EpetraExt::HDF5& file, neighborList_Type& list, std::string const& name )
{
    // copy the map into vectors
    ASSERT ( list.size() > 0, "the map " + name + " is empty!" )
    std::vector<Int> offsets ( list.size() + 1 );
    std::vector<Int> values;
    Int sum ( 0 );
    for ( UInt i ( 0 ); i < list.size(); i++ )
    {
        sum += list[ i ].size();
        offsets[ i + 1 ] = sum;
        for ( neighbors_Type::const_iterator j ( list[ i ].begin() ); j != list[ i ].end(); ++j )
        {
            values.push_back ( *j );
        }
    }

    // Save the vectors into the file
    file.Write ( name, "offsetSize", static_cast<Int> ( offsets.size() ) );
    file.Write ( name, "offsets", H5T_NATIVE_INT, offsets.size(), &offsets[ 0 ] );
    file.Write ( name, "valueSize", static_cast<Int> ( values.size() ) );
    file.Write ( name, "values", H5T_NATIVE_INT, values.size(), &values[ 0 ] );

    /*
    // DEBUG
    std::cerr << name << std::endl;
    for ( UInt i ( 0 ); i < map.size(); i++ )
    {
        std::cerr << i << "> ";
        for ( neighborList_Type::const_iterator j ( map[ i ].begin() ); j != map[ i ].end(); ++j )
        {
            std::cerr << *j << " ";
        }
        std::cerr <<std::endl;
    }
    */
}

void readNeighborMap ( EpetraExt::HDF5& file, neighborList_Type& list, std::string const& name )
{
    // Read the vectors from the file
    Int offsetSize;
    file.Read ( name, "offsetSize", offsetSize );
    std::vector<Int> offsets ( offsetSize );
    file.Read ( name, "offsets", H5T_NATIVE_INT, offsetSize, &offsets[ 0 ] );
    Int valueSize;
    file.Read ( name, "valueSize", valueSize );
    std::vector<Int> values ( valueSize );
    file.Read ( name, "values", H5T_NATIVE_INT, valueSize, &values[ 0 ] );

    // setup the map
    list.resize ( offsetSize - 1 );
    for ( Int i ( 0 ); i < offsetSize - 1; i++ )
    {
        for ( Int j ( offsets[ i ] ); j < offsets[ i + 1 ]; j++ )
        {
            list[ i ].insert ( values[ j ] );
        }
    }

    /*
    // DEBUG
    std::cerr << name << std::endl;
    for ( UInt i ( 0 ); i < map.size(); i++ )
    {
        std::cerr << i << "> ";
        for ( neighborList_Type::const_iterator j ( map[ i ].begin() ); j != map[ i ].end(); ++j )
        {
            std::cerr << *j << " ";
        }
        std::cerr <<std::endl;
    }
    */
}

}

template <typename Mesh>
void GhostHandler<Mesh>::exportToHDF5 ( std::string const& fileName, bool const& truncate )
{
    EpetraExt::HDF5 HDF5 ( *M_comm );

    if ( truncate )
    {
        // Create and open the file / Truncate and open the file
        HDF5.Create ( ( fileName + ".h5" ).data() );
    }
    else
    {
        // Open an existing file without truncating it
        HDF5.Open ( ( fileName + ".h5" ).data() );
    }

    // Check if the file is created
    if ( !HDF5.IsOpen () )
    {
        std::cerr << "Unable to create " + fileName + ".h5";
        abort();
    }

    writeNeighborMap ( HDF5, M_nodeNodeNeighborsList, "nodeNodeNeighborsMap" );
    writeNeighborMap ( HDF5, M_nodeEdgeNeighborsList, "nodeEdgeNeighborsMap" );
    writeNeighborMap ( HDF5, M_nodeElementNeighborsList, "nodeElementNeighborsMap" );

    // Close the file
    HDF5.Close();
}

template <typename Mesh>
void GhostHandler<Mesh>::importFromHDF5 ( std::string const& fileName )
{
    EpetraExt::HDF5 HDF5 ( *M_comm );

    // Open an existing file
    HDF5.Open ( ( fileName + ".h5" ).data() );

    // Check if the file is created
    if ( !HDF5.IsOpen () )
    {
        std::cerr << "Unable to open " + fileName + ".h5";
        abort();
    }

    readNeighborMap ( HDF5, M_nodeNodeNeighborsList, "nodeNodeNeighborsMap" );
    readNeighborMap ( HDF5, M_nodeEdgeNeighborsList, "nodeEdgeNeighborsMap" );
    readNeighborMap ( HDF5, M_nodeElementNeighborsList, "nodeElementNeighborsMap" );

    // Close the file
    HDF5.Close();
}

#endif // HAVE_HDF5

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
    // @TODO: ASSERT_COMPILE_TIME that MeshType::pointMarker == NeighborMarker
    // this guarantees that the nodeNeighbors structure is available.

    // generate node neighbors by watching edges
    // note: this can be based also on faces or volumes
    for ( UInt ie = 0; ie < M_fullMesh->numEdges(); ie++ )
    {
        ID id0 = M_fullMesh->edge ( ie ).point ( 0 ).id();
        ID id1 = M_fullMesh->edge ( ie ).point ( 1 ).id();

        ASSERT ( M_fullMesh->point ( id0 ).id() == id0 && M_fullMesh->point ( id1 ).id() == id1,
                 "the mesh has been reordered, the point must be found" );

        // fill fullMesh points neighboring
        M_fullMesh->point ( id0 ).nodeNeighbors().insert ( id1 );
        M_fullMesh->point ( id1 ).nodeNeighbors().insert ( id0 );
    }

    // update localMesh points
    for ( UInt ip = 0; ip < M_localMesh->numPoints(); ip++ )
    {
        M_localMesh->point ( ip ).nodeNeighbors() = M_fullMesh->point ( M_localMesh->point ( ip ).id() ).nodeNeighbors();
    }
}

template <typename Mesh>
void GhostHandler<Mesh>::createNodeNodeNeighborsMap()
{
    M_nodeNodeNeighborsList.resize ( M_fullMesh->numGlobalPoints() );
    // generate node neighbors by watching edges
    // note: this can be based also on faces or volumes
    for ( UInt ie = 0; ie < M_fullMesh->numEdges(); ie++ )
    {
        ID id0 = M_fullMesh->edge ( ie ).point ( 0 ).id();
        ID id1 = M_fullMesh->edge ( ie ).point ( 1 ).id();

        ASSERT ( M_fullMesh->point ( id0 ).id() == id0 && M_fullMesh->point ( id1 ).id() == id1,
                 "the mesh has been reordered, the point must be found" );

        M_nodeNodeNeighborsList[ id0 ].insert ( id1 );
        M_nodeNodeNeighborsList[ id1 ].insert ( id0 );
    }

#ifdef HAVE_LIFEV_DEBUG
    M_debugOut << "M_nodeNodeNeighborsList on proc " << M_me << std::endl;
    for ( UInt i = 0; i < M_nodeNodeNeighborsList.size(); i++ )
    {
        M_debugOut << i << ": ";
        for ( neighbors_Type::const_iterator it = M_nodeNodeNeighborsList[ i ].begin();
                it != M_nodeNodeNeighborsList[ i ].end(); ++it )
        {
            M_debugOut << *it << " ";
        }
        M_debugOut << std::endl;
    }
#endif
}


template <typename Mesh>
void GhostHandler<Mesh>::createNodeNodeNeighborsMap (flag_Type const& flags)
{
    M_nodeNodeNeighborsList.resize ( M_fullMesh->numGlobalPoints() );
    // generate node neighbors by watching edges
    // note: this can be based also on faces or volumes
    for ( UInt ie = 0; ie < M_fullMesh->numEdges(); ie++ )
    {
        ID id0 = M_fullMesh->edge ( ie ).point ( 0 ).id();
        ID id1 = M_fullMesh->edge ( ie ).point ( 1 ).id();

        ASSERT ( M_fullMesh->point ( id0 ).id() == id0 && M_fullMesh->point ( id1 ).id() == id1,
                 "the mesh has been reordered, the point must be found" );

        if ( this->isInside (M_fullMesh->edge ( ie ).point ( 1 ).markerID(), flags) )
        {
            M_nodeNodeNeighborsList[ id0 ].insert ( id1 );
        }

        if ( this->isInside (M_fullMesh->edge ( ie ).point ( 0 ).markerID(), flags) )
        {
            M_nodeNodeNeighborsList[ id1 ].insert ( id0 );
        }

        //if(M_fullMesh->edge( ie ).point( 1 ).markerID()==20 || M_fullMesh->edge( ie ).point( 1 ).markerID()==1)
        //    M_nodeNodeNeighborsList[ id0 ].insert( id1 );

        //if(M_fullMesh->edge( ie ).point( 0 ).markerID()==20 || M_fullMesh->edge( ie ).point( 0 ).markerID()==1)
        //    M_nodeNodeNeighborsList[ id1 ].insert( id0 );
    }
}

template <typename Mesh>
bool GhostHandler<Mesh>::isInside (ID const& pointMarker, flag_Type const& flags)
{
    int check = 0;
    for (UInt i = 0; i < flags.size(); ++i)
        if (pointMarker == flags[i])
        {
            ++check;
        }
    return (check > 0) ? true : false;
}

template <typename Mesh>
std::set<ID> GhostHandler<Mesh>::createCircleNodeNodeNeighborsMap (UInt Ncircles, UInt GlobalID)
{
    std::set<ID> Neighbors;
    std::set<ID> New_entries;
    std::set<ID> New_neighbors;

    Neighbors = this->nodeNodeNeighborsList() [GlobalID];

    for (UInt i = 0; i < Ncircles - 1; ++i)
    {
        for (std::set<ID>::iterator it = Neighbors.begin(); it != Neighbors.end(); ++it)
        {
            New_entries = this->nodeNodeNeighborsList() [*it];
            for (std::set<ID>::iterator ii = New_entries.begin(); ii != New_entries.end(); ++ii)
            {
                New_neighbors.insert (*ii);
            }
        }
        Neighbors = New_neighbors;
    }

    return Neighbors;
}


template <typename Mesh>
std::set<ID> GhostHandler<Mesh>::createNodeNodeNeighborsMapWithinRadius (double Radius, UInt GlobalID)
{
    std::set<ID> Neighbors;
    std::set<ID> New_entries;
    std::set<ID> New_neighbors;
    bool isInside = true;
    double d = 0;
    UInt mysize = 0;

    Neighbors = this->nodeNodeNeighborsList() [GlobalID];

    while (isInside)
    {
        for (std::set<ID>::iterator it = Neighbors.begin(); it != Neighbors.end(); ++it)
        {
            New_entries = this->nodeNodeNeighborsList() [*it];
            for (std::set<ID>::iterator ii = New_entries.begin(); ii != New_entries.end(); ++ii)
            {
                d = std::sqrt ( pow ( M_fullMesh->point (*ii).x() - M_fullMesh->point (GlobalID).x() , 2) +
                                pow ( M_fullMesh->point (*ii).y() - M_fullMesh->point (GlobalID).y() , 2) +
                                pow ( M_fullMesh->point (*ii).z() - M_fullMesh->point (GlobalID).z() , 2) );
                if (d < Radius)
                {
                    New_neighbors.insert (*ii);
                }
            }
        }

        Neighbors = New_neighbors;

        if (Neighbors.size() == mysize)
        {
            isInside = false;
        }

        mysize = Neighbors.size();

    }

    return Neighbors;

}

template <typename Mesh>
void GhostHandler<Mesh>::createNodeEdgeNeighborsMap()
{
    M_nodeEdgeNeighborsList.resize ( M_fullMesh->numGlobalPoints() );
    // generate node neighbors by watching edges
    // note: this can be based also on faces or volumes
    for ( UInt ie = 0; ie < M_fullMesh->numEdges(); ie++ )
    {
        ID id0 = M_fullMesh->edge ( ie ).point ( 0 ).id();
        ID id1 = M_fullMesh->edge ( ie ).point ( 1 ).id();

        ASSERT ( M_fullMesh->point ( id0 ).id() == id0 && M_fullMesh->point ( id1 ).id() == id1,
                 "the mesh has been reordered, the point must be found" );

        M_nodeEdgeNeighborsList[ id0 ].insert ( ie );
        M_nodeEdgeNeighborsList[ id1 ].insert ( ie );
    }

#ifdef HAVE_LIFEV_DEBUG
    M_debugOut << "M_nodeEdgeNeighborsList on proc " << M_me << std::endl;
    for ( UInt i = 0; i < M_nodeEdgeNeighborsList.size(); i++ )
    {
        M_debugOut << i << ": ";
        for ( neighbors_Type::const_iterator it = M_nodeEdgeNeighborsList[ i ].begin();
                it != M_nodeEdgeNeighborsList[ i ].end(); ++it )
        {
            M_debugOut << *it << " ";
        }
        M_debugOut << std::endl;
    }
#endif

}

template <typename Mesh>
void GhostHandler<Mesh>::createNodeElementNeighborsMap()
{
    M_nodeElementNeighborsList.resize ( M_fullMesh->numGlobalPoints() );
    // generate element neighbors by cycling on elements
    for ( UInt ie = 0; ie < M_fullMesh->numElements(); ie++ )
    {
        ASSERT ( M_fullMesh->element ( ie ).id() == ie,
                 "the mesh has been reordered, the point must be found" );

        for ( UInt k = 0; k < Mesh::element_Type::S_numPoints; k++ )
        {
            ID id ( M_fullMesh->element ( ie ).point ( k ).id() );
            M_nodeElementNeighborsList[ id ].insert ( ie );
        }
    }
}

template <typename Mesh>
typename GhostHandler<Mesh>::map_Type& GhostHandler<Mesh>::ghostMapOnNodes()
{
    // if the map has already been created, return it
    if ( M_ghostMapOnNodes )
    {
        return *M_ghostMapOnNodes;
    }

    if ( M_verbose )
    {
        std::cout << " GH- ghostMapOnNodes()" << std::endl;
    }

    // check that the nodeNeighbors have been created
    if ( M_localMesh->point ( 0 ).nodeNeighbors().empty()  )
    {
        if ( M_verbose )
        {
            std::cerr << "the nodeNeighbors are empty, will be generated now" << std::endl;
        }
        this->createNodeNeighbors();
    }

    // create map
    M_ghostMapOnNodes.reset ( new map_Type() );
    map_Type& ghostMap ( *M_ghostMapOnNodes );

    // use the same Unique map and comm of the original map
    ghostMap.setMap ( M_map->map ( Unique ), Unique );
    ghostMap.setComm ( M_comm );

    // use a set to avoid duplicates
    std::set<Int> myGlobalElementsSet;

    // iterate on local mesh points
    // @todo: this can start from the repeated map and add only neighbors for SUBDOMAIN_INTERFACE marked nodes
    for ( UInt k = 0; k < M_localMesh->numPoints(); k++ )
    {
        // iterate on each node neighborhood
        for ( typename mesh_Type::PointMarker::neighborConstIterator_Type neighborIt = M_localMesh->point ( k ).nodeNeighbors().begin();
                neighborIt != M_localMesh->point ( k ).nodeNeighbors().end(); ++neighborIt )
        {
            myGlobalElementsSet.insert ( *neighborIt );
        }
    }

    std::vector<Int> myGlobalElements ( myGlobalElementsSet.begin(), myGlobalElementsSet.end() );

    // generate map
    map_Type::map_ptrtype repeatedMap ( new Epetra_Map ( -1, myGlobalElements.size(), &myGlobalElements[0], 0, *M_comm ) );
    ghostMap.setMap ( repeatedMap, Repeated );

    // memorize the map in the list
    //    M_ghostMapOnNodesMap[ 1 ] = M_ghostMapOnNodes;

    return *M_ghostMapOnNodes;
}

template <typename Mesh>
typename GhostHandler<Mesh>::map_Type& GhostHandler<Mesh>::ghostMapOnNodes ( UInt overlap )
{
    // if the map has already been created, return it
    if ( M_ghostMapOnNodes )
    {
        return *M_ghostMapOnNodes;
    }

    if ( M_verbose )
    {
        std::cout << " GH- ghostMapOnNodes( UInt )" << std::endl;
    }

    // check that the nodeNodeNeighborsMap has been created
    if ( M_nodeNodeNeighborsList.empty()  )
    {
        if ( M_verbose )
        {
            std::cerr << "the nodeNodeNeighborsMap is empty, will be generated now" << std::endl;
        }
        this->createNodeNodeNeighborsMap();
    }

    // create map
    M_ghostMapOnNodes.reset ( new map_Type() );
    map_Type& ghostMap ( *M_ghostMapOnNodes );

    // use the same Unique map and comm of the original map
    ghostMap.setMap ( M_map->map ( Unique ), Unique );
    ghostMap.setComm ( M_comm );

    Int*          pointer;
    std::set<Int> myGlobalElementsSet, myOriginalElementsSet;;

    // get all elements from the repeated map
    pointer = M_map->map ( Repeated )->MyGlobalElements();
    for ( Int ii = 0; ii < M_map->map ( Repeated )->NumMyElements(); ++ii, ++pointer )
    {
        myOriginalElementsSet.insert ( *pointer );
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
            for ( neighbors_Type::const_iterator neighborIt = M_nodeNodeNeighborsList[ *nodeIt ].begin();
                    neighborIt != M_nodeNodeNeighborsList[ *nodeIt ].end(); ++neighborIt )
            {
                myGlobalElementsSet.insert ( *neighborIt );
            }
        }
        myOriginalElementsSet = myGlobalElementsSet;
    }

    std::vector<Int> myGlobalElements ( myGlobalElementsSet.begin(), myGlobalElementsSet.end() );

    // generate map
    map_Type::map_ptrtype repeatedMap ( new Epetra_Map ( -1, myGlobalElements.size(), &myGlobalElements[0], 0, *M_comm ) );
    ghostMap.setMap ( repeatedMap, Repeated );

    // memorize the map in the list
    //M_ghostMapOnNodesMap[ overlap ] = M_ghostMapOnNodes;

    return *M_ghostMapOnNodes;
}

template <typename Mesh>
typename GhostHandler<Mesh>::map_Type& GhostHandler<Mesh>::ghostMapOnEdges ( UInt overlap )
{
    // if the map has already been created, return it
    if ( M_ghostMapOnEdges )
    {
        return *M_ghostMapOnEdges;
    }

    if ( M_verbose )
    {
        std::cout << " GH- ghostMapOnEdges()" << std::endl;
    }

    // check that the nodeEdgeNeighborsMap has been created
    if ( M_nodeEdgeNeighborsList.empty()  )
    {
        if ( M_verbose )
        {
            std::cerr << "the nodeEdgeNeighborsMap is empty, will be generated now" << std::endl;
        }
        this->createNodeEdgeNeighborsMap();
    }

    // Modified version 25 Jan 2013
    // loop on local mesh edges
    //std::vector<Int> myGlobalElements;
    //myGlobalElements.reserve( M_localMesh->numEdges() );
    // for ( ID ie = 0; ie < M_localMesh->numEdges(); ie++ )
    // {
    //     myGlobalElements.push_back( M_localMesh->edge( ie ).id() );
    // }

    // Original GhostMapEpetra
    // set up Unique (first) and Repeated edges based on the OWNED flag
    std::pair< std::vector<Int>, std::vector<Int> > myGlobalElements;
    myGlobalElements.first.reserve ( M_localMesh->numEdges() );
    myGlobalElements.second.reserve ( M_localMesh->numEdges() );
    // loop on local mesh edges
    for ( ID ie = 0; ie < M_localMesh->numEdges(); ie++ )
    {
        myGlobalElements.second.push_back ( M_localMesh->edge ( ie ).id() );
        if ( Flag::testOneSet ( M_localMesh->edge ( ie ).flag(), EntityFlags::OWNED ) )
        {
            myGlobalElements.first.push_back ( M_localMesh->edge ( ie ).id() );
        }
    }

    // Modified Version
    // create map
    // M_ghostMapOnEdges.reset ( new map_Type( -1, myGlobalElements.size(), &myGlobalElements[0], M_comm ) );

    //Original
    M_ghostMapOnEdges.reset ( new map_Type ( myGlobalElements, M_comm ) );

    if ( overlap > 0 )
    {
        map_Type& ghostMap ( *M_ghostMapOnEdges );

        std::set<Int> myGlobalElementsSet;
        std::set<Int> addedElementsSet;

        // Modified Version
        // for (  UInt k ( 0 ); k < myGlobalElements.size(); k++ )
        // {
        //     typename mesh_Type::EdgeType const & edge = M_fullMesh->edge ( myGlobalElements[ k ] );
        //     for ( UInt edgePoint = 0; edgePoint < mesh_Type::EdgeType::S_numPoints; edgePoint++ )
        //         addedElementsSet.insert( edge.point( edgePoint ).id() );
        // }
        // ( myGlobalElements.begin(), myGlobalElements.end() );

        //Original
        for (  UInt k ( 0 ); k < myGlobalElements.second.size(); k++ )
        {
            typename mesh_Type::edge_Type const& edge = M_fullMesh->edge ( myGlobalElements.second[ k ] );
            for ( UInt edgePoint = 0; edgePoint < mesh_Type::edge_Type::S_numPoints; edgePoint++ )
            {
                addedElementsSet.insert ( edge.point ( edgePoint ).id() );
            }
        }
        ( myGlobalElements.second.begin(), myGlobalElements.second.end() );


        // @todo: optimize this!!
        // 1: work only on the boundary
        // 2: copy back only if necessary
        // repeat on actual nodes to expand overlap
        for ( UInt i = 0; i < overlap; i++ )
        {
            // iterate on points adding all neighbors
            for ( std::set<Int>::const_iterator nodeIt = addedElementsSet.begin();
                    nodeIt != addedElementsSet.end(); ++nodeIt )
            {
                // iterate on each node neighborhood
                for ( neighbors_Type::const_iterator neighborIt = M_nodeEdgeNeighborsList[ *nodeIt ].begin();
                        neighborIt != M_nodeEdgeNeighborsList[ *nodeIt ].end(); ++neighborIt )
                {
                    std::pair<std::set<Int>::iterator, bool> isInserted = myGlobalElementsSet.insert ( *neighborIt );
                    if ( isInserted.second )
                    {
                        typename mesh_Type::EdgeType const& edge = M_fullMesh->edge ( *neighborIt );
                        for ( UInt edgePoint = 0; edgePoint < mesh_Type::EdgeType::S_numPoints; edgePoint++ )
                        {
                            // TODO exclude already included nodes
                            //                            if ( edge.point( edgePoint ).id() != *nodeIt )
                            addedElementsSet.insert ( edge.point ( edgePoint ).id() );
                        }
                    }
                }
            }
        }

        // Modified Version
        // myGlobalElements.assign( myGlobalElementsSet.begin(), myGlobalElementsSet.end() );

        //Original
        myGlobalElements.second.assign ( myGlobalElementsSet.begin(), myGlobalElementsSet.end() );

        //Modified Version
        // generate map
        //map_Type::map_ptrtype repeatedMap ( new Epetra_Map( -1, myGlobalElements.size(), &myGlobalElements[0], 0, *M_comm ) );

        //Original
        map_Type::map_ptrtype repeatedMap ( new Epetra_Map ( -1, myGlobalElements.second.size(), &myGlobalElements.second[0], 0, *M_comm ) );
        ghostMap.setMap ( repeatedMap, Repeated );
    }

    // memorize the map in the list
    //M_ghostMapOnEdgesMap[ overlap ] = M_ghostMapOnEdges;

    return *M_ghostMapOnEdges;
}

template <typename Mesh>
typename GhostHandler<Mesh>::map_Type& GhostHandler<Mesh>::ghostMapOnElementsP0()
{
    // if the map has already been created, return it
    if ( M_ghostMapOnElementsP0 )
    {
        return *M_ghostMapOnElementsP0;
    }

    if ( M_verbose )
    {
        std::cout << " GH- ghostMapOnElementsP0()" << std::endl;
    }

    // create the map
    M_ghostMapOnElementsP0.reset ( new map_Type() );
    map_Type& ghostMap ( *M_ghostMapOnElementsP0 );

    // use the same Unique map and comm of the original map
    ghostMap.setMap ( M_map->map ( Unique ), Unique );
    ghostMap.setComm ( M_comm );

    Int*          pointer;
    std::set<Int> map;

    // get all elements from the repeated map
    pointer = M_map->map ( Repeated )->MyGlobalElements();
    for ( Int ii = 0; ii < M_map->map ( Repeated )->NumMyElements(); ++ii, ++pointer )
    {
        map.insert ( *pointer );
    }

    // add all facing elements
    typedef typename mesh_Type::facet_Type const* facetPtr_Type;
    std::vector<facetPtr_Type>
    facetsOnSubdInt = M_localMesh->facetList().extractElementsWithFlag (
                          EntityFlags::SUBDOMAIN_INTERFACE, &Flag::testOneSet );
    for ( typename std::vector<facetPtr_Type>::const_iterator facetIt = facetsOnSubdInt.begin();
            facetIt != facetsOnSubdInt.end(); ++facetIt )
    {
        map.insert ( (*facetIt)->secondAdjacentElementIdentity() );
    }

    // convert unique list to vector to assure continuity in memorization
    std::vector<Int> myGlobalElements ( map.begin(), map.end() );

    // generate map
    map_Type::map_ptrtype repeatedMap ( new Epetra_Map ( -1, myGlobalElements.size(), &myGlobalElements[0], 0, *M_comm ) );
    ghostMap.setMap ( repeatedMap, Repeated );

    // memorize the map in the list
    //M_ghostMapOnElementsP0Map[ 1 ] = M_ghostMapOnElementsP0;

    return *M_ghostMapOnElementsP0;
}

template <typename Mesh>
typename GhostHandler<Mesh>::map_Type& GhostHandler<Mesh>::ghostMapOnElementsP1 ( UInt overlap )
{
    // if the map has already been created, return it
    if ( M_ghostMapOnElementsP1 )
    {
        return *M_ghostMapOnElementsP1;
    }

    if ( M_verbose )
    {
        std::cout << " GH- ghostMapOnElementsP1()" << std::endl;
    }

    // check that the nodeElementNeighborsMap has been created
    if ( M_nodeElementNeighborsList.empty()  )
    {
        if ( M_verbose )
        {
            std::cerr << "the nodeElementNeighborsMap is empty, will be generated now" << std::endl;
        }
        this->createNodeElementNeighborsMap();
    }

    // create the map
    M_ghostMapOnElementsP1.reset ( new map_Type() );
    map_Type& ghostMap ( *M_ghostMapOnElementsP1 );

    // use the same Unique map and comm of the original map
    ghostMap.setMap ( M_map->map ( Unique ), Unique );
    ghostMap.setComm ( M_comm );

    Int*          pointer;
    std::set<Int> map;

    // get all elements from the repeated map
    pointer = M_map->map ( Repeated )->MyGlobalElements();
    for ( Int ii = 0; ii < M_map->map ( Repeated )->NumMyElements(); ++ii, ++pointer )
    {
        map.insert ( *pointer );
    }

    // add all elements with a node on SUBDOMAIN_INTERFACE
    typedef typename mesh_Type::point_Type const* pointPtr_Type;
    std::vector<pointPtr_Type>
    pointsOnSubdInt = M_localMesh->pointList.extractElementsWithFlag (
                          EntityFlags::SUBDOMAIN_INTERFACE, &Flag::testOneSet );
    // must work on global IDs since added elements are not on localMesh
    std::vector<ID> pointIDOnSubdInt ( pointsOnSubdInt.size() );
    for ( UInt i = 0; i < pointsOnSubdInt.size(); i++)
    {
        pointIDOnSubdInt[i] = pointsOnSubdInt[i]->id();
    }

    std::vector<ID> addedPoints;

    for ( UInt n = 0; n < overlap; n++ )
    {
        for ( std::vector<ID>::const_iterator globalId = pointIDOnSubdInt.begin();
                globalId != pointIDOnSubdInt.end(); ++globalId )
        {
            // iterate on each node neighborhood
            for ( neighbors_Type::const_iterator neighborIt = M_nodeElementNeighborsList[ *globalId ].begin();
                    neighborIt != M_nodeElementNeighborsList[ *globalId ].end(); ++neighborIt )
            {
                std::pair<std::set<Int>::iterator, bool> isInserted = map.insert ( *neighborIt );
                if ( isInserted.second )
                {
                    typename mesh_Type::element_Type const& elem = M_fullMesh->element ( *neighborIt );
                    for ( UInt elemPoint = 0; elemPoint < mesh_Type::element_Type::S_numPoints; elemPoint++ )
                    {
                        // TODO exclude already included nodes
                        addedPoints.push_back ( elem.point ( elemPoint ).id() );
                    }
                }
            }
        }
        // TODO: this must be done only if overlap > 1
        pointIDOnSubdInt = addedPoints;
    }

    // convert unique list to vector to assure continuity in memorization
    std::vector<Int> myGlobalElements ( map.begin(), map.end() );

    // generate map
    map_Type::map_ptrtype repeatedMap ( new Epetra_Map ( -1, myGlobalElements.size(), &myGlobalElements[0], 0, *M_comm ) );
    ghostMap.setMap ( repeatedMap, Repeated );

    // memorize the map in the list
    M_ghostMapOnElementsP1Map[ overlap ] = M_ghostMapOnElementsP1;

    return *M_ghostMapOnElementsP1;
}

template <typename Mesh>
void GhostHandler<Mesh>::fillEntityPID (
    graphPtr_Type graph,
    std::vector<std::vector<UInt> >& entityPID )
{
    if ( M_verbose )
    {
        std::cout << " GH- fillEntityPID()" << std::endl;
    }

    const UInt numParts = graph->size();

    // initialize pointPID to NumProc
    std::vector<UInt>& pointPID = entityPID[ 3 ];
    std::vector<UInt>& elemPID = entityPID[ 0 ];
    std::vector<UInt>& facetPID = entityPID[ 1 ];
    std::vector<UInt>& ridgePID = entityPID[ 2 ];
    pointPID.resize ( M_fullMesh->numPoints(), numParts );
    elemPID.resize ( M_fullMesh->numElements(), numParts );
    facetPID.resize ( M_fullMesh->numFacets(), numParts );
    ridgePID.resize ( M_fullMesh->numRidges(), numParts );

    // @todo: check if parallel building + comm is faster
    for ( UInt p = 0; p < numParts; p++ )
    {
        std::vector<Int>& currentPart = (*graph)[p];
        for ( UInt e = 0; e < currentPart.size(); e++ )
        {
            // point block
            for ( UInt k = 0; k < mesh_Type::element_Type::S_numPoints; k++ )
            {
                const ID& pointID = M_fullMesh->element ( currentPart[ e ] ).point ( k ).id();
                const UInt& pointCurrentPID = pointPID[ pointID ];
                // pointPID should be the minimum between the proc that own it
                if ( p < pointCurrentPID )
                {
                    pointPID[ pointID ] = p;
                }
            }

            // elem block
            const ID& elemID = M_fullMesh->element ( currentPart[ e ] ).id();
            // elemPID is always at its initialization value
            elemPID[ elemID ] = p;

            // facet block
            for ( UInt k = 0; k < mesh_Type::element_Type::S_numFacets; k++ )
            {
                const ID& facetID = M_fullMesh->facet ( M_fullMesh->localFacetId ( elemID, k ) ).id();
                const UInt& facetCurrentPID = facetPID[ facetID ];
                // facetPID should be the minimum between the proc that own it
                if ( p < facetCurrentPID )
                {
                    facetPID[ facetID ] = p;
                }
            }

            // ridge block
            for ( UInt k = 0; k < mesh_Type::element_Type::S_numRidges; k++ )
            {
                const ID& ridgeID = M_fullMesh->ridge ( M_fullMesh->localRidgeId ( elemID, k ) ).id();
                const UInt& ridgeCurrentPID = ridgePID[ ridgeID ];
                // ridgePID should be the minimum between the proc that own it
                if ( p < ridgeCurrentPID )
                {
                    ridgePID[ ridgeID ] = p;
                }
            }
        }
    }
}

template <typename Mesh>
void GhostHandler<Mesh>::ghostMapOnElementsP1 ( graphPtr_Type elemGraph,
                                                const std::vector<UInt>& pointPID,
                                                UInt overlap )
{
    if ( M_verbose )
    {
        std::cout << " GH- ghostMapOnElementsP1( graph )" << std::endl;
    }

    // check that the nodeElementNeighborsMap has been created
    if ( M_nodeElementNeighborsList.empty()  )
    {
        if ( M_verbose )
        {
            std::cerr << "the nodeElementNeighborsMap is empty, will be generated now" << std::endl;
        }
        this->createNodeElementNeighborsMap();
    }
    if ( M_nodeNodeNeighborsList.empty()  )
    {
        if ( M_verbose )
        {
            std::cerr << "the nodeElementNeighborsMap is empty, will be generated now" << std::endl;
        }
        this->createNodeNodeNeighborsMap();
    }

    std::vector<int>& myElems = (*elemGraph) [M_me];
#ifdef HAVE_LIFEV_DEBUG
    // show own elements
    M_debugOut << "own elements on proc " << M_me << std::endl;
    for ( UInt i = 0; i < myElems.size(); i++ )
    {
        M_debugOut << myElems[ i ] << std::endl;
    }
#endif

    // generate graph of points
    graph_Type pointGraph ( M_comm->NumProc() );

    // @todo: check if parallel building + comm is faster
    for ( UInt p = 0; p < static_cast<UInt> ( M_comm->NumProc() ); p++ )
    {
        std::set<int> localPointsSet;
        for ( UInt e = 0; e < (*elemGraph) [ p ].size(); e++ )
        {
            // point block
            for ( UInt k = 0; k < mesh_Type::element_Type::S_numPoints; k++ )
            {
                const ID& pointID = M_fullMesh->element ( (*elemGraph) [ p ][ e ] ).point ( k ).id();
                localPointsSet.insert ( pointID );
            }
        }
        pointGraph[ p ].assign ( localPointsSet.begin(), localPointsSet.end() );
    }

    std::vector<int> const& myPoints = pointGraph[ M_me ];
#ifdef HAVE_LIFEV_DEBUG
    M_debugOut << "own points on proc " << M_me << std::endl;
    for ( UInt i = 0; i < myPoints.size(); i++ )
    {
        M_debugOut << myPoints[ i ] << std::endl;
    }
#endif

    // initialize a bool vector that tells if an element is in the current partition
    std::vector<bool> isInPartition ( M_fullMesh->numElements(), false );
    for ( UInt e = 0; e < myElems.size(); e++ )
    {
        isInPartition[ myElems[ e ] ] = true;
    }

    // find subdomain interface nodes
    std::set<Int> mySubdIntPoints;
    for ( UInt k = 0; k < myPoints.size(); k++ )
    {
        int const& currentPoint = myPoints[ k ];
        // mark as SUBD_INT point only if the point is owned by current process
        if ( pointPID[ currentPoint ] == M_me )
        {

            // check if all element neighbors are on this proc
            for ( neighbors_Type::const_iterator neighborIt = M_nodeElementNeighborsList[ currentPoint ].begin();
                    neighborIt != M_nodeElementNeighborsList[ currentPoint ].end(); ++neighborIt )
            {
                // add the point if a neighbor is missing
                if ( !isInPartition[ *neighborIt ] )
                {
                    mySubdIntPoints.insert ( currentPoint );
                    break;
                }
            }
        }
    }
#ifdef HAVE_LIFEV_DEBUG
    M_debugOut << "own SUBDOMAIN_INTERFACE points on proc " << M_me << std::endl;
    for ( std::set<Int>::const_iterator i = mySubdIntPoints.begin(); i != mySubdIntPoints.end(); ++i )
    {
        M_debugOut << *i << std::endl;
    }
#endif

    std::vector<int> workingPoints ( mySubdIntPoints.begin(), mySubdIntPoints.end() );
    std::set<int> newPoints;
    std::set<int> augmentedElemsSet ( myElems.begin(), myElems.end() );

    for ( UInt o = 0; o < overlap; o++ )
    {
#ifdef HAVE_LIFEV_DEBUG
        M_debugOut << "workingPoints" << std::endl;
        for ( UInt i = 0; i < workingPoints.size(); i++ )
        {
            M_debugOut << workingPoints[ i ] << std::endl;
        }
#endif

        for ( UInt k = 0; k < workingPoints.size(); k++ )
        {
            const int& currentPoint = workingPoints[ k ];
            // iterate on point neighborhood
            for ( neighbors_Type::const_iterator neighborIt = M_nodeElementNeighborsList[ currentPoint ].begin();
                    neighborIt != M_nodeElementNeighborsList[ currentPoint ].end(); ++neighborIt )
            {
                std::pair<std::set<Int>::iterator, bool> isInserted = augmentedElemsSet.insert ( *neighborIt );
                if ( isInserted.second )
                {
                    // if the element is inserted in the list, we add its points to the ones
                    // to be checked for next overlap value
                    for ( UInt j = 0; j < mesh_Type::element_Type::S_numPoints; j++ )
                    {
                        newPoints.insert ( M_fullMesh->element ( *neighborIt ).point ( j ).id() );
                    }
                }
            }
        }

#ifdef HAVE_LIFEV_DEBUG
        M_debugOut << "augmentedElemsSet" << std::endl;
        for ( std::set<Int>::const_iterator i = augmentedElemsSet.begin(); i != augmentedElemsSet.end(); ++i )
        {
            M_debugOut << *i << std::endl;
        }
#endif

        // clean up newPoints from already analized points
        for ( UInt k = 0; k < workingPoints.size(); k++ )
        {
            newPoints.erase ( newPoints.find ( workingPoints[ k ] ) );
        }
#ifdef HAVE_LIFEV_DEBUG
        M_debugOut << "newPoints" << std::endl;
        for ( std::set<int>::const_iterator i = newPoints.begin(); i != newPoints.end(); ++i )
        {
            M_debugOut << *i << std::endl;
        }
#endif

        // set up workingPoints if we are not exiting
        if ( o + 1 < overlap  )
        {
            workingPoints.assign ( newPoints.begin(), newPoints.end() );
        }
    }

    // assign the augmentedElems to the element graph
    myElems.assign ( augmentedElemsSet.begin(), augmentedElemsSet.end() );
}

template <typename Mesh>
void GhostHandler<Mesh>::ghostMapOnElementsP1 ( graphPtr_Type graph,
                                                const std::vector<UInt>& pointPID,
                                                UInt overlap,
                                                const UInt partIndex)
{
    if ( M_verbose )
    {
        std::cout << " GH- ghostMapOnElementsP1( graph )" << std::endl;
    }

    // check that the nodeElementNeighborsMap has been created
    if ( M_nodeElementNeighborsList.empty()  )
    {
        if ( M_verbose )
        {
            std::cerr << "the nodeElementNeighborsMap is empty, will be generated now" << std::endl;
        }
        this->createNodeElementNeighborsMap();
    }
    if ( M_nodeNodeNeighborsList.empty()  )
    {
        if ( M_verbose )
        {
            std::cerr << "the nodeNodeNeighborsMap is empty, will be generated now" << std::endl;
        }
        this->createNodeNodeNeighborsMap();
    }

    const UInt numParts = graph->size();

    std::vector<int>& myElems = (*graph)[partIndex];
#ifdef HAVE_LIFEV_DEBUG
    // show own elements
    M_debugOut << "own elements on proc " << M_me << std::endl;
    for ( UInt i = 0; i < myElems.size(); i++ )
    {
        M_debugOut << myElems[ i ] << std::endl;
    }
#endif

    // generate graph of points
    graph_Type pointGraph ( numParts );

    // @todo: check if parallel building + comm is faster
    for ( UInt p = 0; p < static_cast<UInt> ( numParts ); p++ )
    {
        std::vector<Int>& currentPart = (*graph)[p];
        std::set<int> localPointsSet;
        for ( UInt e = 0; e < currentPart.size(); e++ )
        {
            // point block
            for ( UInt k = 0; k < mesh_Type::element_Type::S_numPoints; k++ )
            {
                const ID& pointID = M_fullMesh->element ( currentPart[ e ] ).point ( k ).id();
                localPointsSet.insert ( pointID );
            }
        }
        pointGraph[ p ].assign ( localPointsSet.begin(), localPointsSet.end() );
    }

    std::vector<int> const& myPoints = pointGraph[ partIndex ];
#ifdef HAVE_LIFEV_DEBUG
    M_debugOut << "own points on proc " << M_me << std::endl;
    for ( UInt i = 0; i < myPoints.size(); i++ )
    {
        M_debugOut << myPoints[ i ] << std::endl;
    }
#endif

    // initialize a bool vector to know if a point is in current partition
    std::vector<bool> isInPartition ( M_fullMesh->numPoints(), false );
    for ( UInt e = 0; e < (*graph)[partIndex].size(); e++ )
    {
        for ( UInt k = 0; k < mesh_Type::element_Type::S_numPoints; k++ )
        {
            const ID& pointID = M_fullMesh->element ( (*graph)[partIndex][ e ] ).point ( k ).id();
            isInPartition[ pointID ] = true;
        }
    }

    // find subdomain interface nodes
    std::set<Int> mySubdIntPoints;
    for ( UInt k = 0; k < myPoints.size(); k++ )
    {
        int const& currentPoint = myPoints[ k ];
        // mark as SUBD_INT point only if the point is owned by current process
        if ( pointPID[ currentPoint ] == partIndex )
        {

            // check if all neighbors are on this proc
            for ( neighbors_Type::const_iterator neighborIt = M_nodeNodeNeighborsList[ currentPoint ].begin();
                    neighborIt != M_nodeNodeNeighborsList[ currentPoint ].end(); ++neighborIt )
            {
                // add the point if a neighbor is missing
                if ( !isInPartition[ *neighborIt ] )
                {
                    mySubdIntPoints.insert ( currentPoint );
                    break;
                }
            }
        }
    }
#ifdef HAVE_LIFEV_DEBUG
    M_debugOut << "own SUBDOMAIN_INTERFACE points on proc " << M_me << std::endl;
    for ( std::set<Int>::const_iterator i = mySubdIntPoints.begin(); i != mySubdIntPoints.end(); ++i )
    {
        M_debugOut << *i << std::endl;
    }
#endif

    std::vector<int> workingPoints ( mySubdIntPoints.begin(), mySubdIntPoints.end() );
    std::set<int> newPoints;
    std::set<int> augmentedElemsSet ( myElems.begin(), myElems.end() );

    for ( UInt o = 0; o < overlap; o++ )
    {
#ifdef HAVE_LIFEV_DEBUG
        M_debugOut << "workingPoints" << std::endl;
        for ( UInt i = 0; i < workingPoints.size(); i++ )
        {
            M_debugOut << workingPoints[ i ] << std::endl;
        }
#endif

        for ( UInt k = 0; k < workingPoints.size(); k++ )
        {
            const int& currentPoint = workingPoints[ k ];
            // iterate on point neighborhood
            for ( neighbors_Type::const_iterator neighborIt = M_nodeElementNeighborsList[ currentPoint ].begin();
                    neighborIt != M_nodeElementNeighborsList[ currentPoint ].end(); ++neighborIt )
            {
                std::pair<std::set<Int>::iterator, bool> isInserted = augmentedElemsSet.insert ( *neighborIt );
                if ( isInserted.second )
                {
                    // if the element is inserted in the list, we add its points to the ones
                    // to be checked for next overlap value
                    for ( UInt j = 0; j < mesh_Type::element_Type::S_numPoints; j++ )
                    {
                        newPoints.insert ( M_fullMesh->element ( *neighborIt ).point ( j ).id() );
                    }
                }
            }
        }

#ifdef HAVE_LIFEV_DEBUG
        M_debugOut << "augmentedElemsSet" << std::endl;
        for ( std::set<Int>::const_iterator i = augmentedElemsSet.begin(); i != augmentedElemsSet.end(); ++i )
        {
            M_debugOut << *i << std::endl;
        }
#endif

        // clean up newPoints from already analized points
        for ( UInt k = 0; k < workingPoints.size(); k++ )
        {
            newPoints.erase ( newPoints.find ( workingPoints[ k ] ) );
        }
#ifdef HAVE_LIFEV_DEBUG
        M_debugOut << "newPoints" << std::endl;
        for ( std::set<int>::const_iterator i = newPoints.begin(); i != newPoints.end(); ++i )
        {
            M_debugOut << *i << std::endl;
        }
#endif

        // set up workingPoints if we are not exiting
        if ( o + 1 < overlap  )
        {
            workingPoints.assign ( newPoints.begin(), newPoints.end() );
        }
    }

    // assign the augmentedElems to the element graph
    myElems.assign ( augmentedElemsSet.begin(), augmentedElemsSet.end() );
}

template <typename Mesh>
void GhostHandler<Mesh>::showMe ( bool const /*verbose*/, std::ostream& out )
{
    out << "GhostHandler::showMe()" << std::endl;
    out << "M_nodeNodeNeighborsMap" << std::endl;
    for ( UInt i = 0; i < M_nodeNodeNeighborsList.size(); i++ )
    {
        out << i << " > ";
        for ( neighbors_Type::const_iterator nIt = M_nodeNodeNeighborsList[ i ].begin();
                nIt != M_nodeNodeNeighborsList[ i ].end(); ++nIt )
        {
            out << *nIt << " ";
        }
        out << std::endl;
    }
    out << "M_nodeElementNeighborsMap" << std::endl;
    for ( UInt i = 0; i < M_nodeNodeNeighborsList.size(); i++ )
    {
        out << i << " > ";
        for ( neighbors_Type::const_iterator nIt = M_nodeNodeNeighborsList[ i ].begin();
                nIt != M_nodeNodeNeighborsList[ i ].end(); ++nIt )
        {
            out << *nIt << " ";
        }
        out << std::endl;
    }
}


}

#endif /* _GHOSTHANDLER_HPP_ */
