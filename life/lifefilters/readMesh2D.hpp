/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#ifndef _READMESH2D_HH_
#define _READMESH2D_HH_
#include "regionMesh2D.hpp"
#include "util_string.hpp"
#include "mesh_util.hpp"
#include "fortran_wrap.hpp"

namespace LifeV
{
/*----------------------------------------------------------------------*
|
|
| #Version 0.1 Experimental 19/8/99. Luca Formaggia
|
| Added support for numbering from 1
! Added markers
|
| #Mesh readers
|
*----------------------------------------------------------------------*/ 
/*****************************************************************
     MESH READER form MESH2D files (unformatted)
******************************************************************/ 
// Macros for Fortran interface

SUBROUTINE_F77 F77NAME( readmesh2d ) ( I_F77 & ne, I_F77 & np,
                                       I_F77 & nptot, I_F77 & nb, I_F77 & nps, I_F77 & nx,
                                       I_F77 & ndimn, I_F77 & npe, I_F77 & npb,
                                       I_F77 & npc, I_F77 * iel,
                                       I_F77 & nd, R4_F77 * coor, I_F77 & ndc,
                                       R4_F77 & xmin, R4_F77 & xmax, R4_F77 & ymin,
                                       R4_F77 & ymax, I_F77 * ib, I_F77 & nbd,
                                       I_F77 * ic, I_F77 * bc, I_F77 * ie,
                                       I_F77 * cpl, R4_F77 * xmed, I_F77 & isw,
                                       I_F77 & ierr, CHARACTER filename );

// subroutine read_mesh2d_head(filename,ne,np,nptot,npe,nb,nx,npc,ierr)

SUBROUTINE_F77 F77NAME( readmesh2dhead ) ( I_F77 & ne, I_F77 & np,
        I_F77 & nptot, I_F77 & npe, I_F77 & nb,
        I_F77 & nps, I_F77 & nx,
        I_F77 & npc, I_F77 & ierr, CHARACTER filename );
//! Reeads a mesh in MEsh2d (LF) format.
template <typename RegionMesh2D>
bool
readMesh2d( RegionMesh2D & mesh, const std::string & fname, EntityFlag regionFlag )
{
    UInt p1, p2, p3;
    UInt nVe, nBVe, nFa, nPo, nBPo, nEd, nBEd;
    UInt i, i1, i2, i3;
    bool p2meshstored, p2meshwanted;
    typedef typename RegionMesh2D::ElementShape ElementShape;

    if ( ElementShape::Shape != TRIANGLE )
    {
        std::cerr << "Sorry, Readmesh2d reads only triangle meshes" << std::endl;
        abort();
    }
    if ( ElementShape::numPoints > 6 )
    {
        std::cerr << "Sorry, ReadMppFiles handles only liner&quad triangles" << std::endl;
        abort();
    }

    CHARACTER filename( const_cast<char *>( fname.c_str() ), fname.length() );
    I_F77 ne, np, nptot, npe, nb, nx, npc, ierr, nps, ndimn, npb;
    F77NAME( readmesh2dhead ) ( ne, np, nptot, npe, nb, nps, nx, npc, ierr, filename );
    if ( ierr != 0 )
    {
        std::cout << " Error in readmesh2d: file " << fname << std::endl;
        std::cout << " not accessible or incorrectly formatted" << std::endl;
        abort();
    }

    std::cout << "INFORMATION ON 2D MESH AS READ FROM MESH2D FILE" << std::endl;
    std::cout << " See mesh2d docs for explanation" << std::endl;
    std::cout << "ne,np,nptot,npe,nb,nx,npc,nps,ierr" << std::endl;
    std::cout << ne << " " << np << " " << nptot << " " << npe << " " << nb << " " << nx << " " << npc << " " << nps << " " << ierr << std::endl;

    // Dimensioning all arrays
    I_F77 nd( npe );
    I_F77 nbd( nps );
    I_F77 ndc( 3 );
    I_F77 isw( 0 );
    FMATRIX<I_F77> iel( npe, ne );
    FMATRIX<I_F77> ib( nps, nb );
    FMATRIX<I_F77> ic( nb );
    FMATRIX<I_F77> bc( nb );
    FMATRIX<I_F77> ie( nb );
    FMATRIX<I_F77> cpl( npc );
    FMATRIX<R4_F77> coor( 3, nptot );
    FMATRIX<R4_F77> xmed( 3, nb );
    R4_F77 xmin, xmax, ymin, ymax;

    // Reading
    F77NAME( readmesh2d ) ( ne, np, nptot, nb, nps, nx, ndimn, npe, npb,
                            npc, iel, nd, coor, ndc, xmin, xmax, ymin, ymax, ib, nbd,
                            ic, bc, ie, cpl, xmed, isw, ierr, filename );
    switch ( ierr )
    {
    case 1:
        std::cerr << " Error in readmesh2d: file incorrectly formatted" << std::endl;
        abort();
    case 0:
        std::cout << " File succesfully read" << std::endl;
        break;
    default:
        std::cerr << " Error in readmesh2d: file incomplete" << std::endl;
        abort();
    }

    // I use explicit constructors instead of relying on implicit conversion rules
    // This to make things more explicit: mesh2d files are (so far) single precision!

    nFa = UInt( ne );
    nBEd = UInt( nb );
    nVe = UInt( np );
    nBVe = UInt( nb );

    /* I Assume that the mesh is OK, so the number of boundary vertices coincides
       with the number of boundary sides: mesh checkers have still to be implemented */ 
    // Do I want a P2 mesh?
    p2meshwanted = ( ElementShape::numPoints == 6 );
    // Do I have a P2 mesh?
    p2meshstored = ( npe == 6 );

    // I still have not yet implemented converters p1->p2 for 2d meshes...
    if ( p2meshwanted && ! p2meshstored )
    {
        std::cerr << " Warning in readmesh2d:" << std::endl;
        std::cout << "file " << fname << std::endl;
        std::cout << "contains a P1 mesh, while we request a P2 mesh" << std::endl;
        std::cout << "Construction of 2D P2 mesh from P1 data not yet implemented" << std::endl;
        abort();
    }

    if ( !p2meshwanted )
    {
        nPo = nVe;
        nBPo = nBVe;
    }
    else
    {
        nPo = UInt( nptot );
        nBPo = 2 * nBVe;
    }
    std::cout << "Using Euler formula to compute number of internal edges, assuming connected domain" << std::endl;
    std::cout << "(More precise technique still to be implemented)" << std::endl;
    nEd = 3 * ( nFa - nBEd ) / 2;
    std::cout << "#Vertices = " << nVe;
    std::cout << "#BVertices= " << nBVe << std::endl;
    std::cout << "#Faces    = " << nFa;
    std::cout << "#Edges    = " << nEd;
    std::cout << "#BEdges   = " << nBEd << std::endl;
    std::cout << "#Points   = " << nPo;
    std::cout << "#BPoints  = " << nBPo << std::endl;

    // I store all Points
    mesh.setMaxNumPoints( nPo, true );
    mesh.setNumBPoints( nBPo );
    mesh.numVertices() = nVe;
    mesh.numBVertices() = nBVe;
    // Only Boundary Edges (in a next version I will allow for different choices)
    mesh.setMaxNumEdges( nBEd );
    mesh.numEdges() = nEd; /* Here the REAL number of edges (all of them)
         even if I store only BEdges.*/
    mesh.setNumBEdges( nBEd );
    //Triangular faces
    mesh.setMaxNumFaces( nFa, true );
    // Add Marker to mesh
    mesh.setMarker( regionFlag );

    // Now put the whole lot into the RegionMesh2D structure
    typename RegionMesh2D::PointType * pp = 0;
    typename RegionMesh2D::EdgeType * pe = 0;
    typename RegionMesh2D::FaceType * pf = 0;


    // first the vertices
    for ( i = 0;i < nVe;i++ )
    {
        pp = &mesh.addPoint( i < nBVe );
        pp->x() = Real( coor( 0, i ) );
        pp->y() = Real( coor( 1, i ) );
    }
    // now the points
    for ( i = nVe;i < nPo;i++ )
    {
        pp = &mesh.addPoint( i - nVe < nBPo );
        pp->x() = Real( coor( 0, i ) );
        pp->y() = Real( coor( 1, i ) );
    }
    std::cout << "Vertices and Points Created " << std::endl;
    // now the boundary edges
    ID ia1;
    long int test;
    EntityFlag ibc;
    for ( i = 0;i < nBEd;i++ )
    {
        pe = &mesh.addEdge( true ); // Only boundary edges.
        p1 = ID( ib( 0, i ) ); // Explicit conversion to ID
        p2 = ID( ib( 1, i ) );
        ibc = EntityFlag( bc( i ) ); //Explicit conversion to entity flag
        // Boundary condition marker
        pe->setMarker( ibc );
        pe->setPoint( 1, mesh.point( p1 ) ); // set edge conn.
        pe->setPoint( 2, mesh.point( p2 ) ); // set edge conn.
        if ( p2meshwanted )
            pe->setPoint( 2, mesh.point( ID( ib( 2, i ) ) ) );
        // fix bedge adjacency information
        ia1 = ID( ie( i ) );
        pe->ad_first() = ia1--; /* Later I use ia1 to have an offset
               so I postdecrement it by 1*/
        i1 = iel( 0, ia1 );
        i2 = iel( 1, ia1 );
        i3 = iel( 2, ia1 );
        test = i1 - p1 + i2 - p2 + i3; // Get the other node
        if ( test == i1 )
            pe->pos_first() = 2;
        else if ( test == i2 )
            pe->pos_first() = 3;
        else if ( test == i3 )
            pe->pos_first() = 1;
        else
            std::cerr << "Information on adjacency of boundary edge is wrong!" << std::endl;
    }
    std::cout << "Boundary Edges Created " << std::endl;

    // Finally the triangular faces!
    for ( i = 0;i < nFa;i++ )
    {
        p1 = ID( iel( 0, i ) );
        p2 = ID( iel( 1, i ) );
        p3 = ID( iel( 2, i ) );
        pf = &( mesh.addFace() ); // Only boundary faces

        pf->setMarker( EntityFlag( ibc ) );
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
}
}
#endif
