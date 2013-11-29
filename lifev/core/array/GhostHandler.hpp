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

#ifdef HAVE_HDF5
#include <EpetraExt_HDF5.h>
#endif

#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/LifeChronoManager.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/NeighborMarker.hpp>

#ifdef HAVE_LIFEV_DEBUG
//#define LIFEV_GHOSTHANDLER_DEBUG 1
#endif

namespace LifeV
{

template <typename MeshType>
class GhostHandler
{
public:

    //! @name Public Types
    //@{

    typedef MeshType mesh_Type;
    typedef boost::shared_ptr<mesh_Type> meshPtr_Type;
    typedef Epetra_Comm comm_Type;
    typedef boost::shared_ptr<comm_Type> commPtr_Type;
    typedef std::set<ID> neighbors_Type;
    typedef std::vector<neighbors_Type> neighborList_Type;
    typedef MapEpetra map_Type;
    typedef boost::shared_ptr<map_Type> mapPtr_Type;
    typedef std::map< UInt, mapPtr_Type > mapList_Type;
    typedef std::vector<Int> idList_Type;
    typedef std::vector<idList_Type> graph_Type;
    typedef boost::shared_ptr<graph_Type> graphPtr_Type;
    typedef std::vector<markerID_Type> markerIDList_Type;

    //@}

    //! @name Constructors & Destructors
    //@{
    //! Constructor
    /*!
     * @param comm. Communicator
     */
    explicit GhostHandler ( commPtr_Type const& comm );

    //! Constructor
    /*!
     * @param fullMesh. Original mesh, before partitioning
     * @param comm. Communicator
     */
    GhostHandler ( meshPtr_Type fullMesh, commPtr_Type const& comm );

    //! Constructor
    /*!
     * @param fullMesh. Original mesh, before partitioning
     * @param localMesh. Local mesh of the current proc
     * @param map. Original map without overlapping
     * @param comm. Communicator
     */
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

    //! List of node neighbors to a node (identified by the global ID)
    neighborList_Type const& nodeNodeNeighborsList()
    {
        ASSERT ( !M_nodeNodeNeighborsList.empty(), "M_nodeNodeNeighborsList is empty" );
        return M_nodeNodeNeighborsList;
    }

    //! List of edge neighbors to a node (identified by the global ID)
    neighborList_Type const& nodeEdgeNeighborsList()
    {
        ASSERT ( !M_nodeEdgeNeighborsList.empty(), "M_nodeEdgeNeighborsList is empty" );
        return M_nodeEdgeNeighborsList;
    }

    //! List of element neighbors to a node (identified by the global ID)
    neighborList_Type const& nodeElementNeighborsList()
    {
        ASSERT ( !M_nodeElementNeighborsList.empty(), "M_nodeElementNeighborsList is empty" );
        return M_nodeElementNeighborsList;
    }

    //@}

    //! @name Set Methods
    //@{

    //! Set verbosity
    /*!
     * @param verbose
     */
    void setVerbose ( const bool& verbose )
    {
        M_verbose = verbose && ( M_me == 0 );
    }

    //@}

    //! @name General Methods
    //@{

    //! Initialize all neighbors list
    void setUp();

    //! Initialize node neighbors list on the subset identified by the given list of MarkerIDs
    void setUp (markerIDList_Type const& flags);

    //! Release pointers to full and local mesh
    void release();

    //! Clean up neighbor lists
    void clean();

#ifdef HAVE_HDF5
    //! Export neighbor lists to an hdf5 file
    /*!
     * @param fileName. Name of the file to write
     * @param truncate. Must be true when the file already exists on disk
     */
    void exportToHDF5 ( std::string const& fileName = "ghostmap", bool const& truncate = true );

    //! Import neighbor lists to an hdf5 file
    /*!
     * @param fileName. Name of the file to write
     */
    void importFromHDF5 ( std::string const& fileName = "ghostmap" );
#endif // HAVE_HDF5

    //! Create node neighbors to nodes and store them in the NeighborMarker
    void createNodeNeighbors();

    //! Create the list of node neighbors to nodes
    void createNodeNodeNeighborsList();

    //! Create the list of node neighbors to nodes that are in the given list of MarkerIDs
    /*!
     * @param flags. The list of MarkerIDs to restrict to.
     */
    void createNodeNodeNeighborsList (markerIDList_Type const& flags);

    //! Create neighbors to a given point, with a specified number of generations
    /*!
     * @param globalID. ID of the point to be examined.
     * @param nCircles. Number of circles (generations) to consider.
     * @return the set of neighbors global IDs
     */
    std::set<ID> circleNeighbors ( UInt globalID, UInt nCircles = 1 );

    //! Create neighbors to a given point within a specified radius
    /*!
     * @param globalID. ID of the point to be examined.
     * @param radius. The value of the circle radius within which neighbors are included
     * @return the set of neighbors global IDs
     */
    std::set<ID> neighborsWithinRadius ( UInt globalID, Real radius );

    //! Create the list of edge neighbors to the nodes
    void createNodeEdgeNeighborsList();

    //! Create the list of element neighbors to the nodes
    void createNodeElementNeighborsList();

    //! Create an overlapped map on nodes
    /*! Create a map based on nodes, expanding it across suddomain interfaces
     * with overlap 1, using NeighborMarker.
     */
    map_Type& ghostMapOnNodes();

    //! Create an overlapped map on nodes
    /*! Create a map based on nodes, expanding it across suddomain interfaces
     * with generic overlap.
     *  @param overlap. Level of overlap between subdomains
     *  @return the overlapped map
     */
    map_Type& ghostMapOnNodes ( UInt overlap );

    //! Create an overlapped map on edges
    /*! Create a map based on edges, expanding it across suddomain interfaces
     * with generic overlap.
     *  @param overlap. Level of overlap between subdomains
     *  @return the overlapped map
     */
    map_Type& ghostMapOnEdges ( UInt overlap );

    //! Create an overlapped map on elements for Finite Volumes
    /*! Create a map based on elements, expanding it across suddomain interfaces.
     *  The elements added are only those that share a facet with the current subdomain.
     *  This type of map is typically used for Finite Volumes.
     *  @return the overlapped map
     */
    // ghostMapOnElementsCommonFacet
    map_Type& ghostMapOnElementsFV();

    //! Create an overlapped map on elements for Finite Elements
    /*! Create a map based on elements, expanding it across suddomain interfaces.
     *  The elements added are all those that share a point with the current subdomain.
     *  This type of map is typically used for Finite Elements.
     *  @param overlap. Level of overlap between subdomains
     */
    // ghostMapOnElementsCommonNodes
    map_Type& ghostMapOnElementsFE ( UInt overlap );

    //! Extend the subdomains graph of the given overlap.
    /*!
     * This method enriches each subdomain with the closest elements such that
     * the partitions have the required overlap.
     * \param elemGraph. The list of subdomain elements
     * \param entityPID. Info about proc ownership of each mesh entity.
     * \param overlap. Level of overlap between partitions.
     */
    void extendGraphFE ( graphPtr_Type elemGraph, idList_Type const& pointPID, UInt overlap );

    //! Extend the subdomains graph of the given overlap.
    /*!
     * This method enriches each subdomain with the closest elements such that
     * the partitions have the required overlap.
     * \param elemGraph. The list of subdomain elements
     * \param entityPID. Info about proc ownership of each mesh entity.
     * \param overlap. Level of overlap between partitions.
     */
    void extendGraphFE ( const boost::shared_ptr<std::vector<std::vector<Int > > >& elemGraph,
                         idList_Type const& pointPID,
                         UInt overlap,
                         UInt partIndex);

    //! showMe method
    void showMe ( bool const verbose = false, std::ostream& out = std::cout );

    //@}

protected:

    //! @name Ghost Maps
    //@{

    mapPtr_Type M_ghostMapOnNodes;
    mapPtr_Type M_ghostMapOnEdges;
    mapPtr_Type M_ghostMapOnElementsFV;
    mapPtr_Type M_ghostMapOnElementsFE;

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
#ifdef LIFEV_GHOSTHANDLER_DEBUG
    std::ofstream M_debugOut;
#endif
    //@}
};

template <typename MeshType>
GhostHandler<MeshType>::GhostHandler ( commPtr_Type const& comm ) :
    M_fullMesh(),
    M_localMesh(),
    M_map(),
    M_comm ( comm ),
    M_me ( comm->MyPID() ),
    M_nodeNodeNeighborsList(),
    M_nodeEdgeNeighborsList(),
    M_nodeElementNeighborsList(),
    M_verbose ( 0 )
#ifdef LIFEV_GHOSTHANDLER_DEBUG
    , M_debugOut ( ( "gh." + ( comm->NumProc() > 1 ? boost::lexical_cast<std::string> ( M_me ) : "s" ) + ".out" ).c_str() )
#endif
{
}

template <typename MeshType>
GhostHandler<MeshType>::GhostHandler ( meshPtr_Type fullMesh,
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
    M_verbose ( 0 )
#ifdef LIFEV_GHOSTHANDLER_DEBUG
    , M_debugOut ( ( "gh." + ( comm->NumProc() > 1 ? boost::lexical_cast<std::string> ( M_me ) : "s" ) + ".out" ).c_str() )
#endif
{
}

template <typename MeshType>
GhostHandler<MeshType>::GhostHandler ( meshPtr_Type fullMesh,
                                       commPtr_Type const& comm ) :
    M_fullMesh ( fullMesh ),
    M_comm ( comm ),
    M_me ( comm->MyPID() ),
    M_nodeNodeNeighborsList(),
    M_nodeEdgeNeighborsList(),
    M_nodeElementNeighborsList(),
    M_verbose ( 0 )
#ifdef LIFEV_GHOSTHANDLER_DEBUG
    , M_debugOut ( ( "gh." + ( comm->NumProc() > 1 ? boost::lexical_cast<std::string> ( M_me ) : "s" ) + ".out" ).c_str() )
#endif
{
}

template <typename MeshType>
void GhostHandler<MeshType>::setUp()
{
    this->createNodeNodeNeighborsList();
    this->createNodeEdgeNeighborsList();
    this->createNodeElementNeighborsList();
}

template <typename MeshType>
void GhostHandler<MeshType>::setUp (markerIDList_Type const& flags)
{
    this->createNodeNodeNeighborsList (flags);
}

template <typename MeshType>
void GhostHandler<MeshType>::release()
{
    M_fullMesh.reset();
    M_localMesh.reset();
}

template <typename MeshType>
void GhostHandler<MeshType>::clean()
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

    //#ifdef LIFEV_GHOSTHANDLER_DEBUG
    //    std::cerr << name << std::endl;
    //    for ( UInt i ( 0 ); i < map.size(); i++ )
    //    {
    //        std::cerr << i << "> ";
    //        for ( neighborList_Type::const_iterator j ( map[ i ].begin() ); j != map[ i ].end(); ++j )
    //        {
    //            std::cerr << *j << " ";
    //        }
    //        std::cerr <<std::endl;
    //    }
    //#endif
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

    //#ifdef LIFEV_GHOSTHANDLER_DEBUG
    //    std::cerr << name << std::endl;
    //    for ( UInt i ( 0 ); i < map.size(); i++ )
    //    {
    //        std::cerr << i << "> ";
    //        for ( neighborList_Type::const_iterator j ( map[ i ].begin() ); j != map[ i ].end(); ++j )
    //        {
    //            std::cerr << *j << " ";
    //        }
    //        std::cerr <<std::endl;
    //    }
    //#endif
}

}

template <typename MeshType>
void GhostHandler<MeshType>::exportToHDF5 ( std::string const& fileName, bool const& truncate )
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

template <typename MeshType>
void GhostHandler<MeshType>::importFromHDF5 ( std::string const& fileName )
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
template <typename MeshType>
void GhostHandler<MeshType>::createNodeNeighbors()
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

template <typename MeshType>
void GhostHandler<MeshType>::createNodeNodeNeighborsList()
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

#ifdef LIFEV_GHOSTHANDLER_DEBUG
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

namespace
{

inline bool isInside ( markerID_Type const& pointMarker, std::vector<markerID_Type> const& markerIDList )
{
    for ( UInt i = 0; i < markerIDList.size(); ++i)
        if ( pointMarker == markerIDList[i] )
        {
            return true;
        }
    return false;
}

}

template <typename MeshType>
void GhostHandler<MeshType>::createNodeNodeNeighborsList (markerIDList_Type const& flags)
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

        if ( isInside (M_fullMesh->edge ( ie ).point ( 1 ).markerID(), flags) )
        {
            M_nodeNodeNeighborsList[ id0 ].insert ( id1 );
        }

        if ( isInside (M_fullMesh->edge ( ie ).point ( 0 ).markerID(), flags) )
        {
            M_nodeNodeNeighborsList[ id1 ].insert ( id0 );
        }

        //if(M_fullMesh->edge( ie ).point( 1 ).markerID()==20 || M_fullMesh->edge( ie ).point( 1 ).markerID()==1)
        //    M_nodeNodeNeighborsList[ id0 ].insert( id1 );

        //if(M_fullMesh->edge( ie ).point( 0 ).markerID()==20 || M_fullMesh->edge( ie ).point( 0 ).markerID()==1)
        //    M_nodeNodeNeighborsList[ id1 ].insert( id0 );
    }
}

template <typename MeshType>
std::set<ID> GhostHandler<MeshType>::circleNeighbors ( UInt globalID, UInt nCircles )
{
    std::set<ID> neighbors;
    std::set<ID> newEntries;
    std::set<ID> newNeighbors;

    neighbors = this->nodeNodeNeighborsList() [ globalID ];

    for (UInt i = 0; i < nCircles - 1; ++i)
    {
        for (std::set<ID>::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
        {
            newEntries = this->nodeNodeNeighborsList() [*it];
            for (std::set<ID>::iterator ii = newEntries.begin(); ii != newEntries.end(); ++ii)
            {
                newNeighbors.insert (*ii);
            }
        }
        neighbors = newNeighbors;
    }

    return neighbors;
}


template <typename MeshType>
std::set<ID> GhostHandler<MeshType>::neighborsWithinRadius ( UInt globalID, Real radius )
{
    std::set<ID> neighbors;
    std::set<ID> newEntries;
    std::set<ID> newNeighbors;
    bool isInside = true;
    Real d = 0;
    UInt mysize = 0;

    neighbors = this->nodeNodeNeighborsList() [globalID];

    typename mesh_Type::point_Type const& p = M_fullMesh->point (globalID);

    while (isInside)
    {
        for (std::set<ID>::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
        {
            newEntries = this->nodeNodeNeighborsList() [*it];
            for (std::set<ID>::iterator ii = newEntries.begin(); ii != newEntries.end(); ++ii)
            {
                typename mesh_Type::point_Type const& n = M_fullMesh->point (*ii);
                d = std::sqrt ( ( n.x() - p.x() ) * ( n.x() - p.x() ) +
                                ( n.y() - p.y() ) * ( n.y() - p.y() ) +
                                ( n.z() - p.z() ) * ( n.z() - p.z() ) );
                if (d < radius)
                {
                    newNeighbors.insert (*ii);
                }
            }
        }

        neighbors = newNeighbors;

        if (neighbors.size() == mysize)
        {
            isInside = false;
        }

        mysize = neighbors.size();

    }

    return neighbors;

}

template <typename MeshType>
void GhostHandler<MeshType>::createNodeEdgeNeighborsList()
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

#ifdef LIFEV_GHOSTHANDLER_DEBUG
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

template <typename MeshType>
void GhostHandler<MeshType>::createNodeElementNeighborsList()
{
    M_nodeElementNeighborsList.resize ( M_fullMesh->numGlobalPoints() );
    // generate element neighbors by cycling on elements
    for ( UInt ie = 0; ie < M_fullMesh->numElements(); ie++ )
    {
        ASSERT ( M_fullMesh->element ( ie ).id() == ie,
                 "the mesh has been reordered, the point must be found" );

        for ( UInt k = 0; k < mesh_Type::element_Type::S_numPoints; k++ )
        {
            ID id ( M_fullMesh->element ( ie ).point ( k ).id() );
            M_nodeElementNeighborsList[ id ].insert ( ie );
        }
    }
}

template <typename MeshType>
typename GhostHandler<MeshType>::map_Type& GhostHandler<MeshType>::ghostMapOnNodes()
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

    return *M_ghostMapOnNodes;
}

template <typename MeshType>
typename GhostHandler<MeshType>::map_Type& GhostHandler<MeshType>::ghostMapOnNodes ( UInt overlap )
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
            std::cerr << "the nodeNodeNeighborsList is empty, will be generated now" << std::endl;
        }
        this->createNodeNodeNeighborsList();
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

    return *M_ghostMapOnNodes;
}

template <typename MeshType>
typename GhostHandler<MeshType>::map_Type& GhostHandler<MeshType>::ghostMapOnEdges ( UInt overlap )
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
            std::cerr << "the nodeEdgeNeighborsList is empty, will be generated now" << std::endl;
        }
        this->createNodeEdgeNeighborsList();
    }

    // set up Unique (first) and Repeated edges based on the GHOST flag
    MapEpetraData mapData;
    mapData.unique.reserve   ( M_localMesh->numEdges() );
    mapData.repeated.reserve ( M_localMesh->numEdges() );
    // loop on local mesh edges
    for ( ID ie = 0; ie < M_localMesh->numEdges(); ie++ )
    {
        mapData.repeated.push_back ( M_localMesh->edge ( ie ).id() );
        if ( M_localMesh->edge ( ie ).isOwned() )
        {
            mapData.unique.push_back ( M_localMesh->edge ( ie ).id() );
        }
    }

    M_ghostMapOnEdges.reset ( new map_Type ( mapData, M_comm ) );

    if ( overlap > 0 )
    {
        map_Type& ghostMap ( *M_ghostMapOnEdges );

        std::set<Int> myGlobalElementsSet;
        std::set<Int> addedElementsSet;

        for (  UInt k ( 0 ); k < mapData.repeated.size(); k++ )
        {
            typename mesh_Type::edge_Type const& edge = M_fullMesh->edge ( mapData.repeated[ k ] );
            for ( UInt edgePoint = 0; edgePoint < mesh_Type::edge_Type::S_numPoints; edgePoint++ )
            {
                addedElementsSet.insert ( edge.point ( edgePoint ).id() );
            }
        }
        ( mapData.repeated.begin(), mapData.repeated.end() );


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

        mapData.repeated.assign ( myGlobalElementsSet.begin(), myGlobalElementsSet.end() );

        map_Type::map_ptrtype repeatedMap ( new Epetra_Map ( -1,
                                                             mapData.repeated.size(),
                                                             &mapData.repeated[0],
                                                             0,
                                                             *M_comm ) );
        ghostMap.setMap ( repeatedMap, Repeated );
    }

    return *M_ghostMapOnEdges;
}

template <typename MeshType>
typename GhostHandler<MeshType>::map_Type& GhostHandler<MeshType>::ghostMapOnElementsFV()
{
    // if the map has already been created, return it
    if ( M_ghostMapOnElementsFV )
    {
        return *M_ghostMapOnElementsFV;
    }

    if ( M_verbose )
    {
        std::cout << " GH- ghostMapOnElementsFV()" << std::endl;
    }

    // create the map
    M_ghostMapOnElementsFV.reset ( new map_Type() );
    map_Type& ghostMap ( *M_ghostMapOnElementsFV );

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

    return *M_ghostMapOnElementsFV;
}

template <typename MeshType>
typename GhostHandler<MeshType>::map_Type& GhostHandler<MeshType>::ghostMapOnElementsFE ( UInt overlap )
{
    // if the map has already been created, return it
    if ( M_ghostMapOnElementsFE )
    {
        return *M_ghostMapOnElementsFE;
    }

    if ( M_verbose )
    {
        std::cout << " GH- ghostMapOnElementsFE()" << std::endl;
    }

    // check that the nodeElementNeighborsMap has been created
    if ( M_nodeElementNeighborsList.empty()  )
    {
        if ( M_verbose )
        {
            std::cerr << "the nodeElementNeighborsList is empty, will be generated now" << std::endl;
        }
        this->createNodeElementNeighborsList();
    }

    // create the map
    M_ghostMapOnElementsFE.reset ( new map_Type() );
    map_Type& ghostMap ( *M_ghostMapOnElementsFE );

    // use the same Unique map and comm of the original map
    ghostMap.setMap ( M_map->map ( Unique ), Unique );
    ghostMap.setComm ( M_comm );

    Int*          pointer;
    std::set<Int> myGlobalElementsSet;

    // get all elements from the repeated map
    pointer = M_map->map ( Repeated )->MyGlobalElements();
    for ( Int ii = 0; ii < M_map->map ( Repeated )->NumMyElements(); ++ii, ++pointer )
    {
        myGlobalElementsSet.insert ( *pointer );
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
                std::pair<std::set<Int>::iterator, bool> isInserted = myGlobalElementsSet.insert ( *neighborIt );
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
        if ( overlap > 1 )
        {
            pointIDOnSubdInt = addedPoints;
        }
    }

    // convert unique list to vector to assure continuity in memorization
    std::vector<Int> myGlobalElements ( myGlobalElementsSet.begin(), myGlobalElementsSet.end() );

    // generate map
    map_Type::map_ptrtype repeatedMap ( new Epetra_Map ( -1, myGlobalElements.size(), &myGlobalElements[0], 0, *M_comm ) );
    ghostMap.setMap ( repeatedMap, Repeated );

    return *M_ghostMapOnElementsFE;
}

template <typename MeshType>
void GhostHandler<MeshType>::extendGraphFE ( graphPtr_Type elemGraph,
                                             idList_Type const& pointPID,
                                             UInt overlap )
{
    if ( M_verbose )
    {
        std::cout << " GH- ghostMapOnElementsP1( graph )" << std::endl;
    }

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    LifeChronoManager<> timeMgr ( M_comm );
    LifeChrono timeNL;
    timeMgr.add ( "node-element ngbr list", &timeNL );
    timeNL.start();
#endif
    // check that the nodeElementNeighborsMap has been created
    if ( M_nodeElementNeighborsList.empty()  )
    {
        if ( M_verbose )
        {
            std::cerr << "the nodeElementNeighborsList is empty, will be generated now" << std::endl;
        }
        this->createNodeElementNeighborsList();
    }
#ifdef LIFEV_GHOSTHANDLER_DEBUG
    timeNL.stop();
#endif

    std::vector<int>& myElems = (*elemGraph) [M_me];

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    // show own elements
    M_debugOut << "own elements on proc " << M_me << std::endl;
    for ( UInt i = 0; i < myElems.size(); i++ )
    {
        M_debugOut << myElems[ i ] << std::endl;
    }
#endif

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    LifeChrono timePG;
    timeMgr.add ( "point graph", &timePG );
    timePG.start();
#endif
    // generate graph of points
    // check: parallel algorithm seems to be faster for this
    graph_Type pointGraph ( M_comm->NumProc() );

    std::set<int> localPointsSet;
    for ( UInt e = 0; e < (*elemGraph) [ M_me ].size(); e++ )
    {
        // point block
        for ( UInt k = 0; k < mesh_Type::element_Type::S_numPoints; k++ )
        {
            const ID& pointID = M_fullMesh->element ( (*elemGraph) [ M_me ][ e ] ).point ( k ).id();
            localPointsSet.insert ( pointID );
        }
    }
    pointGraph[ M_me ].assign ( localPointsSet.begin(), localPointsSet.end() );

    std::vector<Int> pointGraphSize ( M_comm->NumProc(), -1 );
    pointGraphSize[ M_me ] = pointGraph[ M_me ].size();
    for ( UInt p = 0; p < static_cast<UInt> ( M_comm->NumProc() ); p++ )
    {
        M_comm->Broadcast ( &pointGraphSize[ p ], 1, p );
    }

    for ( UInt p = 0; p < static_cast<UInt> ( M_comm->NumProc() ); p++ )
    {
        pointGraph[ p ].resize ( pointGraphSize[ p ] );
    }

    // communicate other proc point graphs
    for ( UInt p = 0; p < static_cast<UInt> ( M_comm->NumProc() ); p++ )
    {
        M_comm->Broadcast ( &pointGraph[p][0], pointGraph[p].size(), p );
    }

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    timePG.stop();
#endif

    std::vector<int> const& myPoints = pointGraph[ M_me ];

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    M_debugOut << "own points on proc " << M_me << std::endl;
    for ( UInt i = 0; i < myPoints.size(); i++ )
    {
        M_debugOut << myPoints[ i ] << std::endl;
    }
#endif

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    LifeChrono timeSI;
    timeMgr.add ( "find SUBD_INT", &timeSI );
    timeSI.start();
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
#ifdef LIFEV_GHOSTHANDLER_DEBUG
    timeSI.stop();
#endif

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    M_debugOut << "own SUBDOMAIN_INTERFACE points on proc " << M_me << std::endl;
    for ( std::set<Int>::const_iterator i = mySubdIntPoints.begin(); i != mySubdIntPoints.end(); ++i )
    {
        M_debugOut << *i << std::endl;
    }
#endif

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    LifeChrono timeOP;
    timeMgr.add ( "overlapping points", &timeOP );
    timeOP.start();
#endif

    std::vector<int> workingPoints ( mySubdIntPoints.begin(), mySubdIntPoints.end() );
    std::set<int> newPoints;
    std::set<int> augmentedElemsSet ( myElems.begin(), myElems.end() );

    for ( UInt o = 0; o < overlap; o++ )
    {

#ifdef LIFEV_GHOSTHANDLER_DEBUG
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

#ifdef LIFEV_GHOSTHANDLER_DEBUG
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

#ifdef LIFEV_GHOSTHANDLER_DEBUG
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
#ifdef LIFEV_GHOSTHANDLER_DEBUG
    timeOP.stop();
#endif

    // assign the augmentedElems to the element graph
    myElems.assign ( augmentedElemsSet.begin(), augmentedElemsSet.end() );

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    timeMgr.print();
#endif
}

template <typename MeshType>
void GhostHandler<MeshType>::extendGraphFE ( const boost::shared_ptr<std::vector<std::vector<Int> > >& elemGraph,
                                             idList_Type const& pointPID,
                                             UInt overlap,
                                             UInt partIndex)
{
    if ( M_verbose )
    {
        std::cout << " GH- ghostMapOnElementsP1( graph )" << std::endl;
    }

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    LifeChronoManager<> timeMgr ( M_comm );
    LifeChrono timeNL;
    timeMgr.add ( "node-element ngbr list", &timeNL );
    timeNL.start();
#endif
    // check that the nodeElementNeighborsMap has been created
    if ( M_nodeElementNeighborsList.empty()  )
    {
        if ( M_verbose )
        {
            std::cerr << "the nodeElementNeighborsList is empty, will be generated now" << std::endl;
        }
        this->createNodeElementNeighborsList();
    }
#ifdef LIFEV_GHOSTHANDLER_DEBUG
    timeNL.stop();
#endif

    std::vector<int>& myElems = (*elemGraph) [partIndex];

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    // show own elements
    M_debugOut << "own elements on proc " << partIndex << std::endl;
    for ( UInt i = 0; i < myElems.size(); i++ )
    {
        M_debugOut << myElems[ i ] << std::endl;
    }
#endif

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    LifeChrono timePG;
    timeMgr.add ( "point graph", &timePG );
    timePG.start();
#endif
    // generate graph of points
    // check: parallel algorithm seems to be faster for this
    Int numParts = elemGraph->size();
    graph_Type pointGraph ( numParts );

    std::set<int> localPointsSet;
    for ( UInt e = 0; e < (*elemGraph) [ partIndex ].size(); e++ )
    {
        // point block
        for ( UInt k = 0; k < mesh_Type::element_Type::S_numPoints; k++ )
        {
            const ID& pointID = M_fullMesh->element ( (*elemGraph) [ partIndex ][ e ] ).point ( k ).id();
            localPointsSet.insert ( pointID );
        }
    }
    pointGraph[ partIndex ].assign ( localPointsSet.begin(), localPointsSet.end() );

    std::vector<Int> pointGraphSize ( numParts, -1 );
    pointGraphSize[ partIndex ] = pointGraph[ partIndex ].size();
    if (M_comm->NumProc() > 1) {
        for ( UInt p = 0; p < static_cast<UInt> ( M_comm->NumProc() ); p++ )
        {
            M_comm->Broadcast ( &pointGraphSize[ p ], 1, p );
        }

        for ( UInt p = 0; p < static_cast<UInt> ( M_comm->NumProc() ); p++ )
        {
            pointGraph[ p ].resize ( pointGraphSize[ p ] );
        }

        // communicate other proc point graphs
        for ( UInt p = 0; p < static_cast<UInt> ( M_comm->NumProc() ); p++ )
        {
            M_comm->Broadcast ( &pointGraph[p][0], pointGraph[p].size(), p );
        }
    }

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    timePG.stop();
#endif

    std::vector<int> const& myPoints = pointGraph[ partIndex ];

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    M_debugOut << "own points on proc " << partIndex << std::endl;
    for ( UInt i = 0; i < myPoints.size(); i++ )
    {
        M_debugOut << myPoints[ i ] << std::endl;
    }
#endif

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    LifeChrono timeSI;
    timeMgr.add ( "find SUBD_INT", &timeSI );
    timeSI.start();
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
        if ( pointPID[ currentPoint ] == partIndex )
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
#ifdef LIFEV_GHOSTHANDLER_DEBUG
    timeSI.stop();
#endif

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    M_debugOut << "own SUBDOMAIN_INTERFACE points on proc " << partIndex << std::endl;
    for ( std::set<Int>::const_iterator i = mySubdIntPoints.begin(); i != mySubdIntPoints.end(); ++i )
    {
        M_debugOut << *i << std::endl;
    }
#endif

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    LifeChrono timeOP;
    timeMgr.add ( "overlapping points", &timeOP );
    timeOP.start();
#endif

    std::vector<int> workingPoints ( mySubdIntPoints.begin(), mySubdIntPoints.end() );
    std::set<int> newPoints;
    std::set<int> augmentedElemsSet ( myElems.begin(), myElems.end() );

    for ( UInt o = 0; o < overlap; o++ )
    {

#ifdef LIFEV_GHOSTHANDLER_DEBUG
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

#ifdef LIFEV_GHOSTHANDLER_DEBUG
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

#ifdef LIFEV_GHOSTHANDLER_DEBUG
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
#ifdef LIFEV_GHOSTHANDLER_DEBUG
    timeOP.stop();
#endif

    // assign the augmentedElems to the element graph
    myElems.assign ( augmentedElemsSet.begin(), augmentedElemsSet.end() );

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    timeMgr.print();
#endif
}

template <typename MeshType>
void GhostHandler<MeshType>::showMe ( bool const /*verbose*/, std::ostream& out )
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
