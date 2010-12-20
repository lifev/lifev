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
    @brief Class for connecting the dof of a mesh (3D) and an interface (2D)
         that lives on the boundary of the mesh.

    @author Vincent Martin
    @date 00-02-2003

    This file contains the class which may be used to update and hold
    the connections between the dof of a mesh (3D) and an interface (2D)
    that lives on the boundary of the mesh. The interface is referenced
    by a flag.

    @contributor M.A. Fernandez
                 Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */


#ifndef _DOFINTERFACE3DTO2D_HH
#define _DOFINTERFACE3DTO2D_HH

#include <life/lifefem/dofInterfaceBase.hpp>
#include <life/lifefem/localDofPattern.hpp>
#include <life/lifefem/dof.hpp>
#include <life/lifemesh/markers.hpp>

#include <iostream>
#include <fstream>
#include <map>
#include <list>   //necessary to write vertices in order.
#include <vector>   //necessary to write faces and access them arbitrarily.


namespace LifeV
{
/*!
  \class DofInterface3Dto2D

  Base class which holds the connections of the dof between a 3D mesh and its 2D interface

  The connections may be built by calling the update method.
*/
class DofInterface3Dto2D:
        public DofInterfaceBase
{
public:

    //! @name Constructor & Destructor
    //@{

    //! default constructor
    DofInterface3Dto2D():
            M_refFE1( 0 ),
            M_dof1( 0 )
    {}

    //! Constructor for interfacing Dof of the same type (LocalDofPattern)
    /*!
      \param refFe the part of the reference FE that contains the dof patterns (nbDofPerEdge...)
      \param dof1 the Dof object of the mesh in which we want to make the computations
    */
    DofInterface3Dto2D( const LocalDofPattern& refFE, const Dof& dof1 );

    //@}


    //! @name Methods
    //@{

    void setup(const LocalDofPattern& refFE1, const Dof& dof1);

    //! This method builds the Dof connections at the interface
    /*!
      \param mesh1 the mesh in which we want to make the computations
      \param flag1 the marker of the interface in the mesh1
    */
    template <typename MeshType>
    void update( const MeshType& mesh1, const entityFlag_Type& flag1 );

    //! Creates an Inria medit type mesh for the interface (pseudo 3D)
    /*! Write the 2D Inria mesh of the interface
    referenced by the number M_interfaceFlag.

    It uses a inria "me.hpp" format. (for medit)

    You should have filled the lists of vertices and faces before.
    (Call update once previously).
    */
    template <typename MeshType>
    void generate2DMesh( std::string fname, const MeshType& mesh1 ) const;

    //! removes all unuseful list (all except M_faceList). use it properly!
    void clearLists();

    //! output
    std::ostream& showMe( bool verbose = false, std::ostream& out = std::cout ) const ;

    //@}


    //! @name Operators
    //@{

    //! Returns the identity of the i-th elements in the (finalised) face list
    //! (counting from 0 ' a la C')
    ID operator[] ( const UInt& i ) const;

    //! Assignment operator (we have a vector of DofInterface3Dto2D)
    DofInterface3Dto2D & operator=( const DofInterface3Dto2D& dofi );

    //@}


    //! @name Get Methods
    //@{

    //! Returns the reference of the interface
    const entityFlag_Type& interfaceFlag() const
    {
        return M_interfaceFlag;
    }

    //! true if the lists have been updated.
    const bool& finalized() const
    {
        return M_finalized;
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! Transforms the 3d index of a vertex into its 2d (interface) index.
    //! This is a simple algorithm... Find out something better some day...?

    ID vertex3Dto2D( const ID& idpoint3D ) const;

    //! This method builds the connections between faces at the interface (M_faceList container)
    /*!
      \param mesh1 the mesh in which we want to make the computations
      \param flag1 the marker of the interface in the mesh1
    */
    template <typename MeshType>
    void updateFaceConnections( const MeshType& mesh1, const entityFlag_Type& flag1 );

    //! This method builds the list of vertices at the interface (M_vertexList container)
    /*!
      \param mesh1 the mesh in which we want to make the computations
    */
    template <typename MeshType>
    void updateVertices( const MeshType& mesh1 );

    //! This method builds the connections between Dof at the interface (_locDof container)
    /*!
      \param mesh1 the mesh in which we want to make the computations
      \param dof1 the Dof object of the mesh in which we want to make the computations
    */
    template <typename MeshType>
    void updateDofConnections( const MeshType& mesh1 );

    //@}


    //! reference of the interface
    entityFlag_Type M_interfaceFlag;

    //! LocalDofPattern object used in the mesh in which we want to make the computations
    const LocalDofPattern * M_refFE1;

    //! Dof object of the mesh in which we want to make the computations
    const Dof * M_dof1;

    /*! STL list which holds the connections between faces at the interface
      -> first  : global (3D) face number
      -> second : interface (2D) face number

      The name has changed : it was _elc before. (V.M.)
    */
    std::vector< std::pair<ID, ID> > M_faceList;

    /*! Auxiliary STL list which holds the connections between vertices at the interface
      (vertices appear more than once. (Mathematically, it is a family))
      Empty after calling update
    */
    std::list<ID> M_vertexPerFaceList;

    /*! STL list which holds the connections between vertices at the interface
      -> first  : global (3D) vertex number
      -> second : interface (2D) vertex number
    */
    std::list< std::pair<ID, ID> > M_vertexList;

    //! true if the lists have been updated.
    bool M_finalized;

};


// ===================================================
// Helper
// ===================================================


/**
   \brief useful function to sort a list and remove multiple numbers.
   Helper function for DoF interface
*/
void RemoveMultiple( const std::list<ID> & list0, std::list< std::pair<ID, ID> > & listf );


// ===================================================
// Methods
// ===================================================

template <typename MeshType>
void DofInterface3Dto2D::
update( const MeshType& mesh1, const entityFlag_Type& flag1 )
{

    // Updating face connections at the interface
    updateFaceConnections( mesh1, flag1 );

    // Updating vertex connections at the interface
    updateVertices( mesh1 );

    // _updateEdges( mesh1 );

    // Update of the Dof connections without inperpolation
    updateDofConnections( mesh1 );

    M_interfaceFlag = flag1;

    M_finalized = true; //! the lists are updated
}

template <typename MeshType>
void DofInterface3Dto2D::
generate2DMesh( std::string fname, const MeshType& mesh1 ) const
{

    ASSERT_PRE( M_finalized, "The lists of vertices and faces must be finalized before generating the interface mesh." );

    std::ofstream ofile( fname.c_str() );
    ASSERT( ofile, "Error: Output file cannot be open" );

    ID idpoint3D;
    ID idpoint2D;
    ID idface3D;

    typedef typename MeshType::FaceShape FaceShape;
    UInt numVertexPerFace = FaceShape::S_numVertices;

    ofile << "MeshVersionFormatted 1\n";
    ofile << "Dimension 3\n";
    ofile << std::endl;

    ofile << "Vertices\n";
    ofile << M_vertexList.size() << "\n";

    //! Write the coordinates of the vertices
    for ( std::list< std::pair<ID, ID> >::const_iterator i2D = M_vertexList.begin(); i2D != M_vertexList.end(); ++i2D )
    {
        idpoint3D = i2D->first;
        ofile << mesh1.pointList( idpoint3D ).x() << " "
        << mesh1.pointList( idpoint3D ).y() << " "
        << mesh1.pointList( idpoint3D ).z() << " "
        << mesh1.pointList( idpoint3D ).marker() << std::endl;
    }
    ofile << std::endl;

    switch ( FaceShape::S_shape )
    {
    case FaceShape::QUAD:
        ofile << "Quadrilaterals\n";
        break;
    case FaceShape::TRIANGLE:
        ofile << "Triangles\n";
        break;
    default:
        ERROR_MSG( "This shape is not implemented in myinterface!" );
    }

    ofile << M_faceList.size() << "\n";

    //! Write the face table
    for ( std::vector< std::pair<ID, ID> >::const_iterator i2D = M_faceList.begin(); i2D != M_faceList.end(); ++i2D )
    {
        idface3D = i2D->first;

        for ( ID vertex = 1 ; vertex <= numVertexPerFace ; ++ vertex )
        {
            idpoint3D = mesh1.boundaryFace( idface3D ).point( vertex ).id();
            idpoint2D = vertex3Dto2D( idpoint3D ); //Simple algorithm (of Search in the list...)
            ofile << idpoint2D << " ";
        }
        ofile << mesh1.boundaryFace( idface3D ).marker() << std::endl;
    }
}

// ===================================================
// Private Methods
// ===================================================

template <typename MeshType>
void DofInterface3Dto2D::
updateFaceConnections( const MeshType& mesh1, const entityFlag_Type& flag1 )
{

    UInt numBoundaryFace1 = mesh1.numBFaces(); // Number of boundary faces in mesh1

    entityFlag_Type marker1;

    typedef typename MeshType::FaceShape GeoBShape; // Shape of the faces

    ID fcounter = 1;  //! Face on the interface counter

    //! Loop on boundary faces on mesh1
    for ( ID iBoundaryFace1 = 1; iBoundaryFace1 <= numBoundaryFace1; ++iBoundaryFace1 )
    {

        //! The face marker
        marker1 = mesh1.boundaryFace( iBoundaryFace1 ).marker();

        //! Is the face on the interface?
        if ( marker1 == flag1 )
        {

            std::pair<ID, ID> fp( iBoundaryFace1, fcounter );
            M_faceList.push_back( fp );
            fcounter ++;  //! local face number

        }
    }
    ASSERT( M_faceList.size() == --fcounter,
                "Local face counter and list size do not match (in face loop)." );
}


template <typename MeshType>
void DofInterface3Dto2D::
updateVertices( const MeshType& mesh1 )
{

    typedef typename MeshType::FaceShape GeoBShape; // Shape of the faces

    UInt numVertexPerFace = GeoBShape::S_numVertices;

    // Loop on faces at the interface (matching faces)
    for ( std::vector< std::pair<ID, ID> >::iterator i = M_faceList.begin(); i != M_faceList.end(); ++i )
    {
        for ( ID jVertex = 1 ; jVertex <= numVertexPerFace ; ++ jVertex )
        {
            M_vertexPerFaceList.push_back( mesh1.boundaryFace( i->first ).point( jVertex ).id() );
        }
    }

    RemoveMultiple( M_vertexPerFaceList , M_vertexList );

    //save memory
    M_vertexPerFaceList.clear();
}


template <typename MeshType>
void DofInterface3Dto2D::
updateDofConnections( const MeshType& mesh1 )
{
    typedef typename MeshType::VolumeShape GeoShape;
    typedef typename MeshType::FaceShape GeoBShape;

    UInt nbVertexPerFace = GeoBShape::S_numVertices; // Number of face's vertices
    UInt nbEdgePerFace = GeoBShape::S_numEdges;    // Number of face's edges

    UInt nbDofPerVertex1 = M_refFE1->nbDofPerVertex; // number of Dof per vertices on mesh1
    UInt nbDofPerEdge1 = M_refFE1->nbDofPerEdge;   // number of Dof per edges on mesh1
    UInt nbDofPerFace1 = M_refFE1->nbDofPerFace;   // number of Dof per faces on mesh1

    UInt nbVertexPerElement = GeoShape::S_numVertices; // Number of element's vertices
    UInt nbEdgePerElement = GeoShape::S_numEdges;    // Number of element's edges

    UInt nDofElemV1 = nbVertexPerElement * nbDofPerVertex1; // number of vertex's Dof on a Element on mesh1
    UInt nDofElemE1 = nbEdgePerElement * nbDofPerEdge1; // number of edge's Dof on a Element on mesh1

    ID iElAd1, iVeEl1, iFaEl1, iEdEl1, gDof1;

    ID locDofCounter1 = 1;

    std::map<ID, ID> locDofMap;

    // Loop on faces at the interface (matching faces)
    for ( std::vector< std::pair<ID, ID> >::iterator i = M_faceList.begin(); i != M_faceList.end(); ++i )
    {

        iElAd1 = mesh1.boundaryFace( i->first ).firstAdjacentElementIdentity();  // id of the element adjacent to the face (mesh1)

        iFaEl1 = mesh1.boundaryFace( i->first ).firstAdjacentElementPosition(); // local id of the face in its adjacent element (mesh1)

        // Vertex based Dof on mesh1
        if ( nbDofPerVertex1 )
        {

            // loop on face vertices (mesh1)
            for ( ID iVeFa1 = 1; iVeFa1 <= nbVertexPerFace; ++iVeFa1 )
            {

                iVeEl1 = GeoShape::faceToPoint( iFaEl1, iVeFa1 ); // local vertex number (in element)

                // Loop number of Dof per vertex (mesh1)
                for ( ID l = 1; l <= nbDofPerVertex1; ++l )
                {

                    gDof1 = M_dof1->localToGlobal( iElAd1, ( iVeEl1 - 1 ) * nbDofPerVertex1 + l ); // Global Dof on mesh1

                    std::pair<ID, ID> locDof( gDof1, locDofCounter1);   //! May be : invert the 2 ??
                    locDofMap.insert( locDof ); // Updating the list of dof connections

                    locDofCounter1 ++;  //! local Dof (total dof on the interface)
                }
            }
        }

        // Edge based Dof on mesh1
        if ( nbDofPerEdge1 )
        {

            // loop on face edges (mesh1)
            for ( ID iEdFa1 = 1; iEdFa1 <= nbEdgePerFace; ++iEdFa1 )
            {

                iEdEl1 = GeoShape::faceToEdge( iFaEl1, iEdFa1 ).first; // local edge number (in element)

                // Loop number of Dof per edge (mesh1)
                for ( ID l = 1; l <= nbDofPerEdge1; ++l )
                {

                    gDof1 = M_dof1->localToGlobal( iElAd1, nDofElemV1 + ( iEdEl1 - 1 ) * nbDofPerEdge1 + l ); // Global Dof on mesh1

                    std::pair<ID, ID> locDof( gDof1, locDofCounter1 );
                    locDofMap.insert( locDof ); // Updating the list of dof connections

                    locDofCounter1 ++;  //! local Dof (total dof on the interface)
                }
            }
        }

        // Face based Dof on mesh1
        for ( ID l = 1; l <= nbDofPerFace1; ++l )
        {

            gDof1 = M_dof1->localToGlobal( iElAd1, nDofElemE1 + nDofElemV1 + ( iFaEl1 - 1 ) * nbDofPerFace1 + l ); // Global Dof in mesh1

            std::pair<ID, ID> locDof( gDof1, locDofCounter1 );
            locDofMap.insert( locDof ); // Updating the list of dof connections

            locDofCounter1 ++;  //! local Dof (total dof on the interface)
        }
    }

//    RemoveMultiple(locDofMap);
    UInt ii = 1;
    for ( std::map<ID, ID>::const_iterator it = locDofMap.begin(); it != locDofMap.end(); ++it, ++ii )
    {
//        std::cout << it->first << " " << it->second << ". ";
        std::pair<ID, ID> _locDof(it->first, ii);
        M_localDofMap.insert(_locDof);
    }
}


}
#endif
