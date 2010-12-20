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
    @file
    @brief Mesh reader from mesh2D files

    @author Luca Formaggia <luca.formaggia@polimi.it>
    @contributor Nur Aiman Fadel <nur.fadel@mail.polimi.it>
    @maintainer Nur Aiman Fadel <nur.fadel@mail.polimi.it>

    @date 19-08-1999

    Mesh reader that it is able to read 2D meshes from freefem and gmsh files.<br>
    ImporterMesh2D reads only triangle meshes.<br>
    readMppFiles handles only liner&quad triangles.<br>
 */

#ifndef _IMPORTERMESH2D_HH_
#define _IMPORTERMESH2D_HH_ 1

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/lambda/lambda.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifecore/fortran_wrap.hpp>
#include <life/lifecore/util_string.hpp>

#include <life/lifemesh/regionMesh2D.hpp>

#include <life/lifefilters/mesh_util.hpp>

namespace LifeV
{

// ===================================================
// Macros for Fortran interface
// ===================================================

// Subroutine readmesh2d
SUBROUTINE_F77 F77NAME( readmesh2d ) ( I_F77 & ne, I_F77 & np,
                                       I_F77 & nptot, I_F77 & nb, I_F77 & nps, I_F77 & nx,
                                       I_F77 & ndimn, I_F77 & npe, I_F77 & npb,
                                       I_F77 & npc, I_F77 * iel,
                                       I_F77 & nd, R4_F77 * coor, I_F77 & ndc,
                                       R4_F77 & xmin, R4_F77 & xmax, R4_F77 & ymin,
                                       R4_F77 & ymax, I_F77 * ib, I_F77 & nbd,
                                       I_F77 * ic, I_F77 * bc, I_F77 * ie,
                                       I_F77 * cpl, R4_F77 * xmed, I_F77 & isw,
                                       I_F77 & ierr, FortranCharacterString filename );

// Subroutine read_mesh2d_head(filename,ne,np,nptot,npe,nb,nx,npc,ierr)
SUBROUTINE_F77 F77NAME( readmesh2dhead ) ( I_F77 & ne, I_F77 & np,
                                           I_F77 & nptot, I_F77 & npe, I_F77 & nb,
                                           I_F77 & nps, I_F77 & nx,
                                           I_F77 & npc, I_F77 & ierr, FortranCharacterString filename );

//! importerMesh2D - reads a mesh in mesh2D(LF) format.
/*!
  It reads a gmsh mesh (2D) file and store it in a RegionMesh2D.

  @param mesh, the mesh data structure to fill in.
  @param fileName, the name of the mesh file  to read.
  @param regionFlag, the identifier for the region.
  @return true if everything went fine, false otherwise.
*/
template <typename RegionMesh2D>

bool
importerMesh2D( RegionMesh2D      & mesh, //importerMesh2D
            const std::string & fileName,
            entityFlag_Type     regionFlag )
{
    UInt i, i1, i2, i3;
    UInt nVe, nBVe, nFa, nPo, nBPo, nEd, nBEd;
    UInt p1, p2, p3;

    bool p2meshstored, p2meshwanted;

    typedef typename RegionMesh2D::ElementShape ElementShape;

    if ( ElementShape::S_shape != TRIANGLE )
    {
        std::cerr << "Sorry, importerMesh2D reads only triangle meshes" << std::endl;
         std::abort();
    }

    if ( ElementShape::S_numPoints > 6 )
    {
        std::cerr << "Sorry, readMppFiles handles only liner&quad triangles" << std::endl;
        std::abort();
    }

    FortranCharacterString filename( const_cast<char *>( fileName.c_str() ), fileName.length() );

    I_F77 ne, np, nptot, npe, nb, nx, npc, ierr, nps, ndimn, npb;

    F77NAME( readmesh2dhead ) ( ne, np, nptot, npe, nb, nps, nx, npc, ierr, filename );

    if ( ierr != 0 )
    {
        std::cout << " Error in importerMesh2D: file " << fileName << std::endl
                  << " not accessible or incorrectly formatted" << std::endl;
        std::abort();
    }

    std::cout << " Information on 2D mesh as read from mesh2D file" << std::endl
              << " See mesh2D docs for explanation" << std::endl << std::endl
              << " ne = " << ne << std::endl
              << " np = " << np << std::endl
              << " nptot = " << nptot << std::endl
              << " npe = " << npe << std::endl
              << " nb = " << nb << std::endl
              << " nx = " << nx << std::endl
              << " npc = " << npc << std::endl
              << " nps = " << nps << std::endl
              << " ierr = " << ierr << std::endl;

    // Dimensioning all arrays
    I_F77 nd( npe );
    I_F77 nbd( nps );
    I_F77 ndc( 3 );
    I_F77 isw( 0 );

    FortranMatrix<I_F77> iel( npe, ne );
    FortranMatrix<I_F77> ib( nps, nb );
    FortranMatrix<I_F77> ic( nb );
    FortranMatrix<I_F77> bc( nb );
    FortranMatrix<I_F77> ie( nb );
    FortranMatrix<I_F77> cpl( npc );
    FortranMatrix<R4_F77> coor( 3, nptot );
    FortranMatrix<R4_F77> xmed( 3, nb );
    R4_F77 xmin, xmax, ymin, ymax;

    // Reading
    F77NAME( readmesh2d ) ( ne, np, nptot, nb, nps, nx, ndimn, npe, npb,
                            npc, iel, nd, coor, ndc, xmin, xmax, ymin, ymax, ib, nbd,
                            ic, bc, ie, cpl, xmed, isw, ierr, filename );

    switch ( ierr )
    {
    case 1:
        std::cerr << " Error in importerMesh2D: file incorrectly formatted" << std::endl;
        std::abort();

    case 0:
        std::cout << " File succesfully read" << std::endl;
        break;

    default:
        std::cerr << " Error in importerMesh2D: file incomplete" << std::endl;
        std::abort();
    }
    /*
    I use explicit constructors instead of relying on implicit conversion rules
    This to make things more explicit: mesh2D files are (so far) single precision!
    */

    nFa =  UInt( ne );
    nBEd = UInt( nb );
    nVe =  UInt( np );
    nBVe = UInt( nb );

    /*
    I Assume that the mesh is OK, so the number of boundary vertices coincides
    with the number of boundary sides: mesh checkers have still to be implemented
    */

    // Do I want a P2 mesh?
    p2meshwanted = ( ElementShape::S_numPoints == 6 );

    // Do I have a P2 mesh?
    p2meshstored = ( npe == 6 );

    // I still have not yet implemented converters p1->p2 for 2D meshes
    if ( p2meshwanted && ! p2meshstored )
    {
        std::cerr << " Warning in importerMesh2D:" << std::endl;
        std::cout << "file " << fileName << std::endl
                  << "contains a P1 mesh, while we request a P2 mesh" << std::endl
                  << "Construction of 2D P2 mesh from P1 data not yet implemented" << std::endl;
        std::abort();
    }

    if ( !p2meshwanted )
    {
        nPo  = nVe;
        nBPo = nBVe;
    }

    else
    {
        nPo  = UInt( nptot );
        nBPo = 2 * nBVe;
    }

    std::cout << "Using Euler formula to compute number of internal edges, assuming connected domain"
              << std::endl;
    std::cout << "(More precise technique still to be implemented)" << std::endl;

    nEd = 3 * ( nFa - nBEd ) / 2;

    std::cout << "Number of Vertices = " << nVe  << std::endl
              << "Number of BVertices= " << nBVe << std::endl
              << "Number of Faces    = " << nFa  << std::endl
              << "Number of Edges    = " << nEd  << std::endl
              << "Number of BEdges   = " << nBEd << std::endl
              << "Number of Points   = " << nPo  << std::endl
              << "Number of BPoints  = " << nBPo << std::endl;

    // I store all the points
    mesh.setMaxNumPoints( nPo, true );
    mesh.setNumBPoints( nBPo );
    mesh.numVertices()  = nVe;
    mesh.numBVertices() = nBVe;

    // Only Boundary Edges
    mesh.setMaxNumEdges( nBEd );
    mesh.numEdges() = nEd;

    /*
    Here the REAL number of edges (all of them) even if I store only BEdges.
    */

    mesh.setNumBEdges( nBEd );

    // Triangular faces
    mesh.setMaxNumFaces( nFa, true );

    // Add Marker to mesh
    mesh.setMarker( regionFlag );

    // Now put the whole lot into the RegionMesh2D structure
    typename RegionMesh2D::point_Type * pp = 0;
    typename RegionMesh2D::EdgeType  * pe = 0;
    typename RegionMesh2D::FaceType  * pf = 0;


    // first the vertices
    for ( i = 0; i < nVe; i++ )
    {
        pp = &mesh.addPoint( i < nBVe );
        pp->x() = Real( coor( 0, i ) );
        pp->y() = Real( coor( 1, i ) );
    }

    // now the points
    for ( i = nVe; i < nPo; i++ )
    {
        pp = &mesh.addPoint( i - nVe < nBPo );
        pp->x() = Real( coor( 0, i ) );
        pp->y() = Real( coor( 1, i ) );
    }

    std::cout << "Vertices and Points Created " << std::endl;

    // now the boundary edges
    ID ia1;
    Int test;
    entityFlag_Type ibc;

    for ( i = 0; i < nBEd; i++ )
    {
        pe = &mesh.addEdge( true ); // Only boundary edges.
        p1 = ID( ib( 0, i ) ); // Explicit conversion to ID
        p2 = ID( ib( 1, i ) );
        ibc = entityFlag_Type( bc( i ) ); //Explicit conversion to entity flag

        // Boundary condition marker
        pe->setMarker( ibc );
        pe->setPoint( 1, mesh.point( p1 ) ); // set edge conn.
        pe->setPoint( 2, mesh.point( p2 ) ); // set edge conn.

        if ( p2meshwanted )
        {
        pe->setPoint( 2, mesh.point( ID( ib( 2, i ) ) ) );
        }

        // fix bedge adjacency information
        ia1 = ID( ie( i ) );
        pe->firstAdjacentElementIdentity() = ia1--; // Later I use ia1 to have an offset so I postdecrement it by 1

        i1 = iel( 0, ia1 );
        i2 = iel( 1, ia1 );
        i3 = iel( 2, ia1 );

        test = i1 - p1 + i2 - p2 + i3; // Get the other node

        if ( test == i1 )
            {
            pe->firstAdjacentElementPosition() = 2;
            }

        else if ( test == i2 )
            {
            pe->firstAdjacentElementPosition() = 3;
            }

        else if ( test == i3 )
            {
            pe->firstAdjacentElementPosition() = 1;
            }

        else
            std::cerr << "Information on adjacency of boundary edge is wrong!" << std::endl;
    }
    std::cout << "Boundary Edges Created " << std::endl;

    // Finally the triangular faces!
    for ( i = 0; i < nFa; i++ )
    {
        p1 = ID( iel( 0, i ) );
        p2 = ID( iel( 1, i ) );
        p3 = ID( iel( 2, i ) );
        pf = &( mesh.addFace() ); // Only boundary faces

        pf->setMarker( entityFlag_Type( ibc ) );
        pf->setPoint( 1, mesh.point( p1 ) ); // set face conn.
        pf->setPoint( 2, mesh.point( p2 ) ); // set face conn.
        pf->setPoint( 3, mesh.point( p3 ) ); // set face conn.

        if ( p2meshwanted )
        {
            p1 = ID( iel( 3, i ) );
            p2 = ID( iel( 4, i ) );
            p3 = ID( iel( 5, i ) );
            pf->setPoint( 4, mesh.point( p1 ) ); // set face conn.
            pf->setPoint( 5, mesh.point( p2 ) ); // set face conn.
            pf->setPoint( 6, mesh.point( p3 ) ); // set face conn.
        }
    }
    std::cout << "Triangular Faces Created " << std::endl;

    return ierr == 0;
} // Function importerMesh2D


//! readGmshFile - reads a mesh in GMSH 2D format.
/*!
  It reads a gmsh mesh (2D) file and store it in a RegionMesh2D.

  @param mesh, the mesh data structure to fill in.
  @param fileName, the name of the gmsh mesh file  to read.
  @param regionFlag, the identifier for the region.
  @return true if everything went fine, false otherwise.
*/

template <typename GeoShape, typename MC>

bool
readGmshFile( RegionMesh2D<GeoShape, MC> & mesh,
              const std::string          & fileName,
              entityFlag_Type              regionFlag )
{
    std::ifstream __is ( fileName.c_str() );

    char buffer[256];
    __is >> buffer;

#ifdef DEBUG
    std::debug( 3000 ) << "buffer: "<< buffer << "\n";
#endif

    UInt __n;
    __is >> __n;

#ifdef DEBUG
    std::debug( 3000 ) << "number of nodes: " << __n;
#endif

    // Add Marker to list of Markers
    mesh.setMarker( regionFlag );

    std::vector<Real> __x( 3 * __n );
    std::vector<bool> __isonboundary( __n );
    std::vector<UInt> __whichboundary( __n );

#ifdef DEBUG
    std::debug( 3000 ) << "reading "<< __n << " nodes\n";
#endif

    std::map<int,int> itoii;

    for ( UInt __i = 0; __i < __n; ++__i )
    {
        UInt __ni;
        __is >> __ni
        >> __x[ 3 * __i ]
        >> __x[ 3 * __i + 1 ]
        >> __x[ 3 * __i + 2 ];

        itoii[ __ni - 1 ] = __i;
    }
    __is >> buffer;

#ifdef DEBUG
    std::debug( 3000 ) << "buffer: "<< buffer << "\n";
#endif

    __is >> buffer;

#ifdef DEBUG
    std::debug( 3000 ) << "buffer: "<< buffer << "\n";
#endif

    UInt __nele;
    __is >> __nele;

    typename RegionMesh2D<GeoShape, MC>::EdgeType * pe = 0;
    typename RegionMesh2D<GeoShape, MC>::FaceType * pf = 0;

#ifdef DEBUG
    std::debug( 3000 ) << "number of elements: " << __nele << "\n";
#endif

    std::vector<std::vector<int> > __e( __nele );
    std::vector<int>               __et( __nele );
    std::vector<int>               __etype( __nele );
    std::vector<int>               __gt( 16 );

    __gt.assign( 16, 0 );

    for ( UInt __i = 0; __i < __nele; ++__i )
    {
       Int __ne, __t, __tag, __np, __dummy;
        __is >> __ne
        >> __t
        >> __tag
        >> __dummy
        >> __np;

        ++__gt[ __t ];

        __etype[ __i ] = __t;
        __et[ __i ]    = __tag;

        __e[ __i ].resize( __np );

       Int __p = 0;

        while ( __p != __np )
        {
            __is >> __e[ __i ][ __p ];
            __e[ __i ][ __p ]  = itoii[ __e[ __i ][ __p ] - 1 ];
            __e[ __i ][ __p ] += 1;

            ++__p;
        }
    }

    std::for_each( __gt.begin(), __gt.end(),  std::cout << boost::lambda::_1 << " " );
    std::cout << "\n";

    // Euler formulas
    UInt n_elements = __gt[ 2 ];

    // Only Boundary Edges (in a next version I will allow for different choices)
    mesh.setMaxNumEdges( __gt[ 1 ] );

    // Here the REAL number of edges (all of them)
    mesh.numEdges    () = __gt[ 1 ];
    mesh.setNumBEdges( __gt[ 1 ] );

#ifdef DEBUG
    std::debug( 3000 ) << "number of edges= " << __gt[ 1 ] << "\n";
#endif

    // Only Boundary Faces
    mesh.setMaxNumFaces( n_elements);

    // Here the REAL number of edges (all of them)
    mesh.numFaces() = n_elements;

#ifdef DEBUG
    std::debug( 3000 ) << "number of faces= " << n_elements << "\n";
#endif

    __isonboundary.assign ( __n, false );
    __whichboundary.assign( __n, 0 );

    for ( UInt __i = 0; __i < __nele; ++__i )
    {
        switch ( __etype[ __i ] )
        {
        // triangular faces (linear)
        case 2:
        {
            __isonboundary[ __e[ __i ][ 0 ] - 1 ] = true;
            __isonboundary[ __e[ __i ][ 1 ] - 1 ] = true;
            __isonboundary[ __e[ __i ][ 2 ] - 1 ] = true;

            __whichboundary[ __e[ __i ][ 0 ] - 1 ] = __et[ __i ];
            __whichboundary[ __e[ __i ][ 1 ] - 1 ] = __et[ __i ];
            __whichboundary[ __e[ __i ][ 2 ] - 1 ] = __et[ __i ];
        }
        }
    }
    // add the point to the mesh
    typename RegionMesh2D<GeoShape, MC>::point_Type * pp = 0;

    mesh.setMaxNumPoints( __n, true );
    mesh.setNumVertices (__n);
    mesh.numBVertices   () = std::count( __isonboundary.begin(), __isonboundary.end(), true );
    mesh.setNumBPoints  ( mesh.numBVertices() );

#ifdef DEBUG
    std::debug( 3000 ) << "number of points : " << mesh.numPoints() << "\n";
    std::debug( 3000 ) << "number of boundary points : " << mesh.numBPoints() << "\n";
    std::debug( 3000 ) << "number of vertices : " << mesh.numVertices() << "\n";
    std::debug( 3000 ) << "number of boundary vertices : " << mesh.numBVertices() << "\n";
#endif

    for ( UInt __i = 0; __i < __n; ++__i )
    {
        pp = &mesh.addPoint( __isonboundary[ __i ] );
        pp->setMarker( __whichboundary[ __i ] );
        pp->x() = __x[ 2 * __i ];
        pp->y() = __x[ 2 * __i + 1 ];
    }

    // add the element to the mesh
    for ( UInt __i = 0; __i < __nele; ++__i )
    {
        switch ( __etype[ __i ] )
        {
        // segment(linear)
        case 1:
        {
            pe = &( mesh.addEdge( true ) );
            pe->setMarker( entityFlag_Type( __et[ __i ] ) );
            pe->setPoint( 1, mesh.point( __e[ __i ][ 0 ] ) );
            pe->setPoint( 2, mesh.point( __e[ __i ][ 1 ] ) );
        }
        break;

        // triangular faces (linear)
        case 2:
        {
            pf = &( mesh.addFace() );
            pf->setMarker( entityFlag_Type( __et[ __i ] ) );
            pf->setPoint( 1, mesh.point( __e[ __i ][ 0 ] ) );
            pf->setPoint( 2, mesh.point( __e[ __i ][ 1 ] ) );
            pf->setPoint( 3, mesh.point( __e[ __i ][ 2 ] ) );
        }
        break;

        // quadrangular faces(linear)
        case 3:
        {
            pf = &( mesh.addFace() );
            pf->setMarker( entityFlag_Type( __et[ __i ] ) );
            pf->setPoint( 1, mesh.point( __e[ __i ][ 0 ] ) );
            pf->setPoint( 2, mesh.point( __e[ __i ][ 1 ] ) );
            pf->setPoint( 3, mesh.point( __e[ __i ][ 2 ] ) );
            pf->setPoint( 4, mesh.point( __e[ __i ][ 3 ] ) );
        }
        break;
        }
    }
    return true;
} // Function readGmshFile

//! readFreeFemFile - reads a mesh in FreeFem 2D format.
/*!
read a freefem mesh (2D) file and store it in a RegionMesh2D.

@param mesh, the mesh data structure to fill in.
@param fileName, the name of the freefem mesh file to read.
@param regionFlag, the identifier for the region.
@param bool useless, it will be removed.
@return true if everything went fine, false otherwise.
 */

template <typename GeoShape, typename MC>
bool
readFreeFemFile( RegionMesh2D<GeoShape, MC> & mesh,
                 const std::string          & fileName,
                 entityFlag_Type              regionFlag, bool useless )
{
    BareItemsHandler<BareEdge> _be;
    std::pair<BareEdge, bool> _edge;

    typename RegionMesh2D<GeoShape, MC>::point_Type * pp = 0;
    typename RegionMesh2D<GeoShape, MC>::EdgeType * pe = 0;
    typename RegionMesh2D<GeoShape, MC>::FaceType * pf = 0;

    std::ifstream __is ( fileName.c_str() );

    // first row: how many vertices, triangles, edges
    UInt __nv, __nt, __ne, i1, i2, i3;
    __is >> __nv >> __nt >> __ne;

#ifdef DEBUG
    std::debug( 3000 ) << "number of vertices: "<< __nv << "\n";
    std::debug( 3000 ) << "number of triangles: "<< __nt << "\n";
    std::debug( 3000 ) << "number of edges: "<< __ne << "\n";
#endif

    // first section: read the list of vertices
    // on each row find the two coordinates and the label for each node
    std::vector<Real> __x(2*__nv);
    std::vector<bool> __isonboundary(__nv);
    std::vector<UInt> __whichboundary(__nv);

#ifdef DEBUG
    std::debug( 3000 ) << "reading "<< __nv << " nodes\n";
#endif

    // count the number of nodes on the boundary
    UInt __nbv( 0 );

    // reading vertices
    for ( UInt __i = 0; __i < __nv; ++ __i )
    {
        __is >> __x[ 2 * __i ] >> __x[ 2 * __i + 1 ] >> __whichboundary[ __i ];
        __isonboundary[ __i ] = __whichboundary[ __i ];
        __nbv                += __isonboundary[ __i ];
    }

    // second section: read the list of triangles
    // on each row find the three nodes and the label for each triangle
    std::vector<int> __triangle_nodes( 3 * __nt );
    std::vector<int> __triangle_label( __nt );

#ifdef DEBUG
    std::debug( 3000 ) << "reading "<< __nt << " triangles\n";
#endif

    std::map<UInt,UInt> edge_to_firstAdjacentElementIdentity, edge_to_firstAdjacentElementPosition;

    // reading vertices
    for ( UInt __i = 0; __i < __nt; ++__i )
    {
        __is >> __triangle_nodes[ 3 * __i ]
        >> __triangle_nodes[ 3 * __i + 1 ]
        >> __triangle_nodes[ 3 * __i + 2 ]
        >> __triangle_label[ __i ];

        // dump first the existing edges, to maintain the correct numbering
        // if everything is correct the numbering in the bareedge
        // structure will reflect the actual edge numbering

        std::pair<UInt, bool> _check;

        i1 = __triangle_nodes[ 3 * __i ];
        i2 = __triangle_nodes[ 3 * __i + 1 ];
        i3 = __triangle_nodes[ 3 * __i + 2 ];

        _edge                             = makeBareEdge( i1, i2 );
        _check                            = _be.addIfNotThere( _edge.first );
        edge_to_firstAdjacentElementIdentity[ _check.first ]  = __i + 1;
        edge_to_firstAdjacentElementPosition[ _check.first ] = 1;

        _edge                             = makeBareEdge( i2, i3 );
        _check                            = _be.addIfNotThere( _edge.first );
        edge_to_firstAdjacentElementIdentity[ _check.first ]  = __i + 1;
        edge_to_firstAdjacentElementPosition[ _check.first ] = 2;

        _edge                             = makeBareEdge( i3, i1 );
        _check                            = _be.addIfNotThere( _edge.first );
        edge_to_firstAdjacentElementIdentity[ _check.first ]  = __i + 1;
        edge_to_firstAdjacentElementPosition[ _check.first ] = 3;
    }

    //    ( __triangle[ 3 * i + 2] > __triangle[ 3 * i + 1 ] )

    // third section: read the list of edges
    // NOTE: only boundary edges are stored
    // on each row find the two nodes and the label for each edge
    std::vector<int> __edge_nodes( 2 * __ne );
    std::vector<int> __edge_label( __ne );

    // reading edges

    for ( UInt __i = 0; __i < __ne; ++__i )
    {
        __is >> __edge_nodes[ 2 * __i ] >> __edge_nodes[ 2 * __i + 1 ] >> __edge_label[ __i ];
    }

    // Set mesh properties
    // Add Marker to list of Markers
    mesh.setMarker( regionFlag );

    // Till now I only have information about boundary edges - I don't know the MAX num of edges
    // Euler formula: ne = nv + nt - 1
    mesh.setMaxNumEdges      ( __nv + __nt - 1 );
    mesh.setMaxNumGlobalEdges( __nv + __nt - 1 );

    // Here the REAL number of edges (all of them)
    mesh.numEdges() = __nv + __nt - 1;

    mesh.setNumBEdges        ( __ne );
    mesh.setMaxNumFaces      ( __nt );
    mesh.setMaxNumGlobalFaces( __nt );

    // Here the REAL number of edges (all of them)
    mesh.numFaces() = __nt;

    mesh.setMaxNumPoints      ( __nv, true );
    mesh.setMaxNumGlobalPoints( __nv );

    mesh.setNumVertices       ( __nv );
    mesh.setNumGlobalVertices ( __nv );

    mesh.numBVertices()      = __nbv;
    mesh.setNumBPoints( mesh.numBVertices() );

#ifdef DEBUG
    std::debug( 3000 ) << "number of points : " << mesh.numPoints() << "\n";
    std::debug( 3000 ) << "number of boundary points : " << mesh.numBPoints() << "\n";
    std::debug( 3000 ) << "number of vertices : " << mesh.numVertices() << "\n";
    std::debug( 3000 ) << "number of boundary vertices : " << mesh.numBVertices() << "\n";
#endif

    for ( UInt __i = 0; __i < __nv; ++__i )
    {
        pp = &mesh.addPoint( __isonboundary[ __i ] );
        pp->setMarker( __whichboundary[ __i ] );
        pp->x() = __x[ 2 * __i ];
        pp->y() = __x[ 2 * __i + 1 ];
        pp->setId( __i + 1 );
        pp->setLocalId( __i + 1 );

        mesh.localToGlobalNode().insert(std::make_pair( __i + 1 , __i + 1 ) );
        mesh.globalToLocalNode().insert(std::make_pair( __i + 1 , __i + 1 ) );
    }

    // add the edges to the mesh
    for ( UInt __i = 0; __i < __ne; ++__i )
    {
        pe = &( mesh.addEdge( true ) );
        pe->setMarker( entityFlag_Type( __edge_label[ __i ] ) );
        pe->setPoint( 1, mesh.point( __edge_nodes[ 2 * __i ] ) );
        pe->setPoint( 2, mesh.point( __edge_nodes[ 2 * __i + 1 ] ) );
        _edge = makeBareEdge( __edge_nodes[ 2 * __i ], __edge_nodes[ 2 * __i + 1 ] );
        UInt map_it( _be.id( _edge.first ) );
        pe->firstAdjacentElementIdentity() = edge_to_firstAdjacentElementIdentity[ map_it ];
        pe->firstAdjacentElementPosition() = edge_to_firstAdjacentElementPosition[ map_it ];
    }

    // add the triangles to the mesh
    for ( UInt __i = 0; __i < __nt; ++__i )
    {
        pf = &( mesh.addFace() );
        pf->setId     ( __i + 1 );
        pf->setLocalId( __i + 1);
        pf->setMarker( entityFlag_Type( __triangle_label[ __i ] ) );
        pf->setPoint( 1, mesh.point( __triangle_nodes[ 3 * __i ] ) );
        pf->setPoint( 2, mesh.point( __triangle_nodes[ 3 * __i + 1 ] ) );
        pf->setPoint( 3, mesh.point( __triangle_nodes[ 3 * __i + 2 ] ) );
    }

    return true;
} // Function readFreeFemFile

template <typename RegionMesh2D>
bool __attribute__ ((__deprecated__))
readMesh2d( RegionMesh2D      & mesh,
            const std::string & fileName,
            entityFlag_Type     regionFlag )
{
return importerMesh2D( mesh, fileName, regionFlag );
}

} // Namespace LifeV

#endif /* IMPORTERMESH2D_H */
