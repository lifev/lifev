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
#include <lifev/core/util/LifeAssert.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshChecks.hpp>
#include <lifev/core/mesh/MeshElement.hpp>
#include <lifev/core/mesh/MeshUtility.hpp>
#include <lifev/core/mesh/Marker.hpp>
//#include <fstream>


namespace LifeV
{

//! @name Methods
//@{

//! This method gives the flags for a rectangle 
/*!
  @param i_x 
  @param i_y 
  @param n_x Number of elements along the length
  @param n_y Number of elements along the width
*/

markerID_Type regularMeshPointPosition2D( const UInt& i_x,
                                     	 	const UInt& i_y,
                                     	 	const UInt& n_x,
                                     	 	const UInt& n_y );


//! This method generate a rectangular structured mesh
/*!
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
template <typename GeoShape, typename MC>
void regularMesh2D( RegionMesh<GeoShape,MC>& mesh,
                    markerID_Type regionFlag,
                    const UInt& m_x,
                    const UInt& m_y,
                    bool verbose=false,
                    const Real& l_x=1.0,
                    const Real& l_y=1.0,
                    const Real& t_x=0.0,
                    const Real& t_y=0.0 )
{
    typedef GeoShape geoShape_Type;
    typedef MC       marker_Type;

	std::cout<<geoShape_Type::S_numRidges<<std::endl;

    ASSERT(geoShape_Type::S_numRidges == 3,"Not creating a P1 structured mesh");
	
    // output stream
	std::ostream& oStr = std::cout;
	
	/*
    std::cout<<"l_x: "<<l_x<<std::endl;
    std::cout<<"l_y: "<<l_y<<std::endl;
    std::cout<<"m_x: "<<m_x<<std::endl;
    std::cout<<"m_y: "<<m_y<<std::endl;
    */
	
    // discretization
    Real dx( l_x / m_x ); std::cout<<"dx: "<<dx<<std::endl;
    Real dy( l_y / m_y ); std::cout<<"dy: "<<dy<<std::endl;

    // Number of nodes along the side of the unit cube
    UInt n_x( m_x + 1 );
    UInt n_y( m_y + 1 );

    // Incremental values in order to get the indices of the nodes
    // Due to the structure, if we add N_i to the number of the
    // current node, we end up with the number of the next point
    // in the i axis.
    UInt N_x( 1 );
    UInt N_y( n_x );

    // Data about the mesh
    UInt verticesNumber( n_x * n_y );
    UInt boundaryVerticesNumber( verticesNumber - ( n_x - 2 ) * ( n_y - 2 ) ); //Total-inside points
    UInt boundaryEdgesNumber(2 * ( m_x + m_y ) );
    UInt edgesNumber(
        // Edges that draws the rectangles
        m_x * n_y +                                
        n_x * m_y +                                 
        // Edges that go accross the rectangles
        m_x * m_y                                   
    );
    
    // Faces 
    UInt elementsNumber( 2*( m_x * m_y ) );

    UInt pointsNumber( 0 );
    UInt boundaryPointsNumber( 0 );
    
    oStr << "Linear Tetra Mesh" << std::endl;
    pointsNumber = verticesNumber;
    boundaryPointsNumber = boundaryVerticesNumber;


    // Set the data
    oStr << "initialization of mesh...";

    // Note: The vertices are the nodes of the mesh while the points
    //       are the nodes and some new points added for example for
    //       the quadratic triangle

    // About points:
    mesh.setNumBoundaryRidges( boundaryPointsNumber );
    mesh.setMaxNumRidges( pointsNumber, true );
    mesh.setMaxNumGlobalRidges( pointsNumber );

    // About vertices:
    mesh.setNumVertices( verticesNumber );
    mesh.setNumGlobalVertices( verticesNumber );
    mesh.setNumBVertices( boundaryVerticesNumber );

    // About edges:
    mesh.setNumFacets( edgesNumber );
    mesh.setnumBoundaryFacets( boundaryEdgesNumber );
    mesh.setMaxNumFacets( edgesNumber );
    mesh.setMaxNumGlobalFacets( edgesNumber );

    // About faces:
    mesh.setMaxNumElements( elementsNumber );
    mesh.setMaxNumGlobalElements( elementsNumber );

    mesh.setMaxNumFaces      ( elementsNumber );
    mesh.setMaxNumGlobalFaces( elementsNumber );
    mesh.setNumFaces         ( elementsNumber );


    mesh.setMarkerID( regionFlag );
    oStr << "done" << std::endl;
        
    // Declaration of pointers on the different mesh entities
    typename RegionMesh<GeoShape,MC>::ridge_Type*   pointPtr  = 0;
    typename RegionMesh<GeoShape,MC>::facet_Type*   edgePtr   = 0;
    typename RegionMesh<GeoShape,MC>::element_Type* elementPtr   = 0;
           
    // Build the points of the mesh
    oStr << "building the points of the mesh...";
    Real xPosition( 0.0 ), yPosition( 0.0 ), zPosition( 0.0 );
    markerID_Type nodeFlag( 0 );
    UInt nodeID( 0 );
    UInt P0( 0 ), P1( 0 ), P2( 0 ), P3( 0 );
    
    /*
         P3___P2
          | / |
          |/__|
        P0     P1	  
    */
    
    for ( UInt j(0); j<n_y; ++j )
	{
		yPosition = dy * j;

		for ( UInt i(0); i<n_x; ++i )
		{
			xPosition = dx * i;
			nodeFlag = regularMeshPointPosition2D(i, j, n_x, n_y );

			// We create the point
			if ( nodeFlag>0 )
			{
				pointPtr = &mesh.addRidge( true ); //it is a boundary point
			}
			else
			{
				pointPtr = &mesh.addRidge( false );
			}
			

			// We set the point properties
			nodeID = j * N_y + i;
			pointPtr->setId( nodeID );
			pointPtr->setLocalId( nodeID );
			//std::cout<<"Info point: "<<pointPtr->MeshVertex::showMe(true)<<std::endl;
			
			//mesh.localToGlobalNode().insert( std::make_pair( nodeID, nodeID) );
			//mesh.globalToLocalNode().insert( std::make_pair( nodeID, nodeID) );

			pointPtr->setMarkerID( nodeFlag );
			pointPtr->x() = xPosition + t_x;
			pointPtr->y() = yPosition + t_y;
			pointPtr->z() = zPosition;
			std::cout<<"Info point: "<<pointPtr->MeshVertex::showMe(true)<<std::endl;
			
		}
	}
    
    oStr << "done" << std::endl;

    // Build the faces
    oStr << "building the elements...";
    UInt faceID(0);  
    
	for ( UInt j(0); j<m_y; ++j )
	{
		for ( UInt i(0); i<m_x; ++i )
		{
			faceID = ( j * m_x + i ) * 2;
				
			nodeID = j * N_y + i;
			P0 = nodeID;
			P1 = nodeID + N_x;
			P2 = nodeID + N_x + N_y;
			P3 = nodeID + N_y;
		
		
		// Triangle 1
		elementPtr = &mesh.addElement();
		elementPtr->setId( faceID );
		elementPtr->setLocalId( faceID );
		elementPtr->setPoint( 0, mesh.point(P1) );
		elementPtr->setPoint( 1, mesh.point(P2) );
		elementPtr->setPoint( 2, mesh.point(P0) );
		elementPtr->setMarkerID( regionFlag );
		//std::cout<<elementPtr->showMe(true);
		
		// Triangle 2
		elementPtr = &mesh.addElement();
		elementPtr->setId( faceID + 1 );
		elementPtr->setLocalId( faceID + 1 );
		elementPtr->setPoint( 0, mesh.point(P3) );
		elementPtr->setPoint( 1, mesh.point(P0) );
		elementPtr->setPoint( 2, mesh.point(P2) );
		elementPtr->setMarkerID( regionFlag );
		//std::cout<<elementPtr->showMe(true);
		
		}
	}
	oStr << "done" << std::endl;
	
	// add the boundary edges to the mesh
	oStr << "building the boundary edges...";
    for ( UInt i = 0; i < boundaryEdgesNumber ; ++i )
    {
		UInt edgeLabel = 0;
		UInt adjID = 0;
		UInt pos = 0;

		//edge labels are:
		// BOTTOM = 1
		// LEFT = 2
		// RIGHT = 3
		// TOP = 4

		if(i < m_x)
		{
			nodeID = i;
			P0 = nodeID;
			P1 = nodeID + 1;
			edgeLabel = 1; //BOTTOMEDGE
			adjID = 2*i;
			pos = 2;

		}else if(i < m_x + m_y)
		{
			nodeID = (i + 1 - m_x)*n_x -1;
			P0 = nodeID;
			P1 = nodeID + n_x;
			edgeLabel = 3; //RIGHTEDGE
			adjID =  ((i - m_x)*m_x + m_x - 1)*2;
			pos = 0;
		}
		else if(i < 2*m_x + m_y)
		{
			nodeID = n_x*n_y - 1 - (i - m_x - m_y);
			P0 = nodeID;
			P1 = nodeID - 1;
			edgeLabel = 4; //TOPEDGE
			adjID = 2*m_x*m_y-(i - m_x - m_y)*2 -1;
			pos = 2;
		}
		else{
			nodeID =  n_x*n_y - 1 - m_x - (i - 2*m_x - m_y)*n_x ;
			P0 = nodeID;
			P1 = nodeID - n_x ;
			edgeLabel = 2; //LEFTEDGE
			adjID = -(i - 2*m_x - m_y)*(2*m_x) + (m_x*(m_y-1))*2 + 1;
			pos = 0;
		}


		std::cout<<"Edge# "<<i<<"  PO="<<P0<<"  P1="<<P1<<"  Adj="<<adjID<<std::endl;
        edgePtr = &mesh.addFacet( true ) ;
        edgePtr->setMarkerID( edgeLabel );
        edgePtr->setPoint( 0, mesh.point( P0 ));
        edgePtr->setPoint( 1, mesh.point( P1 ));
        edgePtr->firstAdjacentElementIdentity() = adjID;
        edgePtr->firstAdjacentElementPosition() = pos;
     }

    oStr << "done" << std::endl;
    
    std::cout << "Numero elementi" << mesh.elementList().size();
    std::cout << "Numero bedge" << boundaryEdgesNumber;
    std::cout << "Numero edge" << edgesNumber;
	// edges update
    mesh.updateElementFacets( true, true, edgesNumber );
    
    // mesh control
    // std::cout<<mesh.check( 0, true, true )<<std::endl;
    // mesh.showMe();
}

//@}


} // Namespace LifeV

#endif /* STRUCTUREDMESH2D_HPP */
