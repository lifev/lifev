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
    @brief Contains methods which generate 2D structured meshes.

    @author Iori Guido <guido.iori@mail.polimi.it>
    @contributor -

    @date 23-05-2011

*/

#ifndef STRUCTUREDMESH2D_HPP
#define STRUCTUREDMESH2D_HPP 1

#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshChecks.hpp>

namespace LifeV
{

// Labels for the structured 2D mesh
namespace Structured2DLabel
{
//! Label for the internal entities
const markerID_Type INTERNAL = 0;

//! Label for the bottom boundary edge
const markerID_Type BOTTOM = 1;

//! Label for the left boundary edge
const markerID_Type LEFT = 2;

//! Label for the top boundary edge
const markerID_Type TOP = 3;

//! Label for the right boundary edge
const markerID_Type RIGHT = 4;

//! Label for the top and left boundary corner
const markerID_Type TOP_LEFT = 7;

//! Label for the bottom and right boundary corner
const markerID_Type BOTTOM_RIGHT = 5;

//! Label for the bottom and left boundary corner
const markerID_Type BOTTOM_LEFT = 8;

//! Label for the top and right boundary corner
const markerID_Type TOP_RIGHT = 6;

}

/*!
  @brief This method gives the flags for a rectangle

  @param i_x
  @param i_y
  @param n_x Number of elements along the length
  @param n_y Number of elements along the width
*/
markerID_Type regularMeshPointPosition2D ( const UInt& i_x,
                                           const UInt& i_y,
                                           const UInt& n_x,
                                           const UInt& n_y );


/*!
  @brief This method generate a rectangular structured mesh

  For the square \f$ [0,1]^2 \f$ the internal flag is 0.
  <br>
  For the corners the labels are:
  <ul>
  <li> BOTTOM and RIGHT = 5, i.e. \f$ x = 1 \f$ and \f$ y = 0 \f$ </li>
  <li> TOP and RIGHT = 6, i.e. \f$ x = 1 \f$ and \f$ y = 1 \f$ </li>
  <li> TOP and LEFT = 7, i.e. \f$ x = 0 \f$ and \f$ y = 1 \f$ </li>
  <li> BOTTOM and LEFT = 8, i.e. \f$ x = 0 \f$ and \f$ y = 0 \f$ </li>
  </ul>
  For the edges the labels are:
  <ul>
  <li> LEFT    = 1, i.e. \f$ x = 0 \f$ </li>
  <li> BOTTOM  = 2, i.e. \f$ y = 0 \f$ </li>
  <li> RIGHT   = 3, i.e. \f$ x = 1 \f$ </li>
  <li> TOP     = 4, i.e. \f$ y = 1 \f$ </li>
  </ul>

  @param mesh The mesh that we want to generate
  @param regionFlag Flag of the region
  @param m_x Number of elements along the length
  @param m_y Number of elements along the width
  @param l_x length of the mesh
  @param l_y width of the mesh
  @param verbose Verbose mode enabled/disabled
  @param t_x translation of the mesh along the x-axis
  @param t_y translation of the mesh along the y-axis
*/
template <typename MeshType>
void regularMesh2D ( MeshType& mesh,
                     markerID_Type regionFlag,
                     const UInt& m_x,
                     const UInt& m_y,
                     bool verbose = false,
                     const Real& l_x = 1.0,
                     const Real& l_y = 1.0,
                     const Real& t_x = 0.0,
                     const Real& t_y = 0.0 )
{
    typedef MeshType mesh_Type;
    typedef typename mesh_Type::geoShape_Type geoShape_Type;

    ASSERT ( ( geoShape_Type::S_shape == TRIANGLE )
             || ( geoShape_Type::S_shape == QUAD ),
             "Type of 2d structured mesh not available." );

    // discretization
    const Real dx ( l_x / m_x );
    const Real dy ( l_y / m_y );

    // Number of nodes along the side of the rectangle
    const UInt n_x ( m_x + 1 );
    const UInt n_y ( m_y + 1 );

    // Incremental values in order to get the indices of the nodes
    // Due to the structure, if we add N_i to the number of the
    // current node, we end up with the number of the next point
    // in the i axis.
    const UInt N_x ( 1 );
    const UInt N_y ( n_x );

    // Data about the mesh

    //Total-inside points
    const UInt verticesNumber ( n_x * n_y );
    const UInt boundaryVerticesNumber ( verticesNumber - ( n_x - 2 ) * ( n_y - 2 ) );
    const UInt boundaryEdgesNumber (2 * ( m_x + m_y ) );
    const UInt edgesNumber (
        // Edges that draws the rectangles
        m_x * n_y +
        n_x * m_y +
        // Edges that go accross the rectangles
        m_x * m_y
    );

    // Faces
    const UInt elementsNumber ( 2 * ( m_x * m_y ) );

    // Set the data

    // Note: The vertices are the nodes of the mesh while the points
    //       are the nodes and some new points added for example for
    //       the quadratic triangle

    // About points:
    mesh.setNumBoundaryRidges ( boundaryVerticesNumber );
    mesh.setMaxNumRidges      ( verticesNumber, true );
    mesh.setMaxNumGlobalRidges ( verticesNumber );

    // About vertices:
    mesh.setNumVertices      ( verticesNumber );
    mesh.setNumGlobalVertices ( verticesNumber );
    mesh.setNumBVertices     ( boundaryVerticesNumber );

    // About edges:
    mesh.setNumFacets         ( edgesNumber );
    mesh.setNumBoundaryFacets ( boundaryEdgesNumber );
    mesh.setMaxNumFacets      ( edgesNumber );
    mesh.setMaxNumGlobalFacets ( edgesNumber );

    // About faces:
    mesh.setMaxNumElements      ( elementsNumber );
    mesh.setMaxNumGlobalElements ( elementsNumber );

    mesh.setMaxNumFaces      ( elementsNumber );
    mesh.setMaxNumGlobalFaces ( elementsNumber );
    mesh.setNumFaces         ( elementsNumber );


    mesh.setMarkerID ( regionFlag );

    // Declaration of pointers on the different mesh entities
    typename mesh_Type::ridge_Type*   pointPtr   = 0;
    typename mesh_Type::facet_Type*   edgePtr    = 0;
    typename mesh_Type::element_Type* elementPtr = 0;

    // Build the points of the mesh
    Real xPosition ( 0.0 ), yPosition ( 0.0 ), zPosition ( 0.0 );
    markerID_Type nodeFlag ( 0 );
    UInt nodeID ( 0 );
    UInt P0 ( 0 ), P1 ( 0 ), P2 ( 0 ), P3 ( 0 );

    /*
         P3___P2
          | / |
          |/__|
        P0     P1
    */

    for ( UInt j (0); j < n_y; ++j )
    {
        yPosition = dy * j;

        for ( UInt i (0); i < n_x; ++i )
        {
            xPosition = dx * i;
            nodeFlag = regularMeshPointPosition2D (i, j, n_x, n_y );

            // We create the point
            pointPtr = &mesh.addRidge ( nodeFlag > 0 ); // node flag determines if the point is on boundary


            // We set the point properties
            nodeID = j * N_y + i;
            pointPtr->setId ( nodeID );

            pointPtr->setMarkerID ( nodeFlag );
            pointPtr->x() = xPosition + t_x;
            pointPtr->y() = yPosition + t_y;
            pointPtr->z() = zPosition;
        }
    }

    // Build the faces
    UInt faceID (0);

    for ( UInt j (0); j < m_y; ++j )
    {
        for ( UInt i (0); i < m_x; ++i )
        {
            faceID = ( j * m_x + i ) * 2;

            nodeID = j * N_y + i;
            P0 = nodeID;
            P1 = nodeID + N_x;
            P2 = nodeID + N_x + N_y;
            P3 = nodeID + N_y;

            // Triangle 1
            elementPtr = &mesh.addElement();
            elementPtr->setId ( faceID );
            elementPtr->setPoint ( 0, mesh.point (P1) );
            elementPtr->setPoint ( 1, mesh.point (P2) );
            elementPtr->setPoint ( 2, mesh.point (P0) );
            elementPtr->setMarkerID ( regionFlag );

            // Triangle 2
            elementPtr = &mesh.addElement();
            elementPtr->setId ( faceID + 1 );
            elementPtr->setPoint ( 0, mesh.point (P3) );
            elementPtr->setPoint ( 1, mesh.point (P0) );
            elementPtr->setPoint ( 2, mesh.point (P2) );
            elementPtr->setMarkerID ( regionFlag );
        }
    }

    // add the boundary edges to the mesh
    for ( UInt i = 0; i < boundaryEdgesNumber ; ++i )
    {
        UInt edgeLabel = 0;
        UInt adjID = 0;
        UInt pos = 0;

        if (i < m_x)
        {
            nodeID = i;
            P0 = nodeID;
            P1 = nodeID + 1;
            edgeLabel = Structured2DLabel::BOTTOM; //BOTTOMEDGE
            adjID = 2 * i;
            pos = 2;

        }
        else if (i < m_x + m_y)
        {
            nodeID = (i + 1 - m_x) * n_x - 1;
            P0 = nodeID;
            P1 = nodeID + n_x;
            edgeLabel = Structured2DLabel::RIGHT; //RIGHTEDGE
            adjID =  ( (i - m_x) * m_x + m_x - 1) * 2;
            pos = 0;
        }
        else if (i < 2 * m_x + m_y)
        {
            nodeID = n_x * n_y - 1 - (i - m_x - m_y);
            P0 = nodeID;
            P1 = nodeID - 1;
            edgeLabel = Structured2DLabel::TOP; //TOPEDGE
            adjID = 2 * m_x * m_y - (i - m_x - m_y) * 2 - 1;
            pos = 2;
        }
        else
        {
            nodeID =  n_x * n_y - 1 - m_x - (i - 2 * m_x - m_y) * n_x ;
            P0 = nodeID;
            P1 = nodeID - n_x ;
            edgeLabel = Structured2DLabel::LEFT; //LEFTEDGE
            adjID = - (i - 2 * m_x - m_y) * (2 * m_x) + (m_x * (m_y - 1) ) * 2 + 1;
            pos = 0;
        }


        edgePtr = &mesh.addFacet ( true ) ;
        edgePtr->setId ( mesh.facetList().size() - 1 );
        edgePtr->setMarkerID ( edgeLabel );
        edgePtr->setPoint ( 0, mesh.point ( P0 ) );
        edgePtr->setPoint ( 1, mesh.point ( P1 ) );
        edgePtr->firstAdjacentElementIdentity() = adjID;
        edgePtr->firstAdjacentElementPosition() = pos;
    }

    // edges update
    mesh.updateElementFacets ( true, verbose, edgesNumber );

}


} // Namespace LifeV

#endif /* STRUCTUREDMESH2D_HPP */
