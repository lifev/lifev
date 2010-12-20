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
    @brief Degrees of freedom, the class that provides the localtoglobal table

    @author M.A. Fernandez & Luca Formaggia
    @date 00-07-2002

    @contributor Vincent Martin
                 Mohamed Belhadj
                 Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    The numbering of the degrees of freedom follow the order
    Vertices - Edges - Faces - Volumes
    If more than one degree of freedon is associated to a geometrical entity, the
    global numbering associated to two "consecutive" degree of freedom  on the given
    entity is consecutive, and the order is according to the entity orientation.

 */

#ifndef _DOF_HH
#define _DOF_HH

#include <life/lifearray/SimpleVect.hpp>

#include <life/lifecore/life.hpp>

#include <life/lifefem/DOFLocalPattern.hpp>

#include <life/lifemesh/basisElSh.hpp>

#include <algorithm>
#include <map>

namespace LifeV
{
/*! Local-to-global table

This class provides the localtoglobal table that relates the local DOF of
a finite element to its global numbering. It needs a LocalDofPattern in
order to obtain all the necessary information about the local pattern. In
fact it stores a copy of it so to make the local pattern available, if
needed.

It is useless until is has not been set up on a specific RegionMesh. This is accomplished either by
passing the mesh to the constructor, or calling the method Dof::update().

\note The methods bulds the table for ALL degrees of freedom, i.e. it does not handle any essential
boundary condition.

Now the class include also a local-to-global table with DOF grouped by (internal) face that was implemented in the old versions into the dofByFace.hpp and dofByFace.cpp files created by D. A. Di Pietro in 2004
*/
class Dof
{
public:

    //! @name Public Types
    //@{

    //@}


    //! @name Constructor & Destructor
    //@{

    /*! The minimal constructor
      \param fePattern is the LocalDofPattern on which the ref FE is built
      \param Offset: the smallest Dof numbering. It might be used if we want the
      degrees of freedom numbering start from a specific value.
    */
    Dof( const LocalDofPattern& fePattern, UInt offset = 1 );

    //! Copy constructor
    Dof( const Dof & dof2 );

    //! Constructor accepting a mesh as parameter
    /*!
      \param mesh a RegionMesh3D
      \param _fe is the LocalDofPattern on which the ref FE is built
      \param Offset: the smalest Dof numbering. It might be used if we want the
      degrees of freedom numbering start from a specific value.
    */
    template <typename MeshType>
    Dof( MeshType& mesh, const LocalDofPattern& fePattern, UInt offset = 1 );

    //@}


    //! @name Methods
    //@{

    //! Build the localToGlobal table
    /*!
      \param mesh A RegionMesh3D
      Updates the LocaltoGlobal array
    */
    template <typename MeshType>
    void update( MeshType & );

    /*!
      Returns the global numbering of a DOF, given an internal face and the
      local numbering
      \param faceId the internal face ID
      \param localDOF the local DOF numbering (starting from 1)
      \return The global numbering of the DOF
    */
    ID localToGlobalByFace(const ID& faceId, const ID& localDof, bool& exist ) const;

    //! Ouput
    void showMe( std::ostream & out = std::cout, bool verbose = false ) const;
    void showMeByFace(std::ostream& out = std::cout, bool verbose = false) const;

    //@}


    //! @name Set Methods
    //@{

    void setTotalDof(const UInt& totalDof)
    {
        M_totalDof = totalDof;
    }

    //@}


    //! @name Get Methods
    //@{

    //! The total number of Dof
    const UInt& numTotalDof() const
    {
        return M_totalDof;
    }

    //! The number of local Dof (nodes) in the finite element
    const UInt& numLocalDof() const
    {
        return M_elementDofPattern.nbLocalDof();
    }

    //! Return the specified entries of the localToGlobal table
    /*!
      Returns the global numbering of a DOF, given an element and the local numbering
      \param ELId the element ID
      \param localNode the local DOF numbering (starting from 1)
      \return The numbering of the DOF
    */
    const ID& localToGlobal( const ID ElId, const ID localNode) const
    {
        return M_localToGlobal( localNode, ElId );
    }

    //! Number of elements in mesh
    const UInt& numElements() const
    {
        return M_numElement;
    }

    //! Number of faces in the mesh
    const UInt& numFaces() const
    {
        return M_nbFace;
    }

    //! Number of local vertices (in a elment)
    const UInt& numLocalVertices() const
    {
        return M_nbLocalVertex;
    }

    //! Number of local edges (in a elment)
    const UInt& numLocalEdges() const
    {
        return M_nbLocalEdge;
    }

    //! Number of local faces (in a elment)
    const UInt& numLocalFaces() const
    {
        return M_nbLocalFace;
    }

    //! Number of Local DofByFace
    const UInt& numLocalDofByFace() const
    {
        ASSERT_PRE( (M_numLocalDofByFace>0) , "This data are not available for this reference element");
        return M_numLocalDofByFace;
    }

    //! Getter for the localDofPattern
    const LocalDofPattern& localDofPattern() const
    {
        return M_elementDofPattern;
    }


    //@}



private:

    typedef SimpleArray<UInt> Container_Type;

    typedef ID ( *faceToPointPtr_Type )( ID const& localFace, ID const& point );

    //! The pattern of the local degrees of freedom.
    const LocalDofPattern& M_elementDofPattern;

    // Offset for the first degree of freedom numerated
    UInt M_offset;

    // Total number of degrees of freedom
    UInt M_totalDof;

    // Number of elements in the considered mesh
    UInt M_numElement;


    UInt M_nbLocalVertex;
    UInt M_nbLocalEdge;
    UInt M_nbLocalFace;

    // The local to global table
    Container_Type M_localToGlobal;

    // number of faces in the mesh
    UInt M_nbFace;

    // The local to global table based on the faces
    std::vector<std::vector<ID> > M_localToGlobalByFace;

    // the face to the global dof
    std::map<ID,ID> M_globalToLocalByFace;

    // local array that maps the local dof of the
    faceToPointPtr_Type M_faceToPoint;

    // face to the local dof the the element
    UInt M_numLocalDofByFace;

    // Just 5 counters
    UInt M_dofPositionByEntity[ 5 ];
};


// ===================================================
// Constructors & Destructor
// ===================================================

//! Constructor that builds the localToglobal table
template <typename MeshType>
Dof::Dof( MeshType& mesh, const LocalDofPattern& fePattern, UInt offset ) :
        M_elementDofPattern       ( fePattern ),
        M_offset  ( offset ),
        M_totalDof( 0 ),
        M_numElement     ( 0 ),
        M_nbLocalVertex      ( 0 ),
        M_nbLocalEdge      ( 0 ),
        M_nbLocalFace      ( 0 ),
        M_localToGlobal     (),
        M_nbFace( 0 ),
        M_localToGlobalByFace(),
        M_globalToLocalByFace()
{
    //Getting the face
    switch ( fePattern.nbLocalDof() )
    {
    case 2:
        // No M_faceToPoint (it is 1D)
        M_numLocalDofByFace = 1;
        break;
    case 4:
        M_faceToPoint = LinearTetra::faceToPoint;
        M_numLocalDofByFace = 3;
        break;
    case 5:
        M_faceToPoint = LinearTetraBubble::faceToPoint;
        M_numLocalDofByFace = 3;
        break;
    case 10:
        M_faceToPoint = QuadraticTetra::faceToPoint;
        M_numLocalDofByFace = 6;
        break;
    case 8:
        M_faceToPoint = LinearHexa::faceToPoint;
        M_numLocalDofByFace = 4;
        break;
    case 27:
        M_faceToPoint = QuadraticHexa::faceToPoint;
        M_numLocalDofByFace = 27;
        break;
    default:
        std::cout << "Warning: This refFE is not available for the dof by face." << std::endl;
        M_numLocalDofByFace = 0;
        break;
    }

    for ( UInt i = 0; i < 5; ++i )
        M_dofPositionByEntity[ i ] = 0;
    update( mesh );
}

// ===================================================
// Methods
// ===================================================

//! Build the localToGlobal table
template <typename MeshType>
void Dof::update( MeshType& mesh )
{

    typedef typename MeshType::ElementShape GeoShapeType;

    // Some useful local variables, to save some typing
    UInt nbLocalDofPerEdge = M_elementDofPattern.nbDofPerEdge();
    UInt nbLocalDofPerVertex = M_elementDofPattern.nbDofPerVertex();
    UInt nbLocalDofPerFace = M_elementDofPattern.nbDofPerFace();
    UInt nbLocalDofPerVolume = M_elementDofPattern.nbDofPerVolume();

    M_nbLocalVertex = GeoShapeType::S_numVertices;
    M_nbLocalEdge = GeoShapeType::S_numEdges;
    M_nbLocalFace = GeoShapeType::S_numFaces;

    M_numElement = mesh.numElements();

    UInt nbGlobalVolume = mesh.numGlobalVolumes();
    UInt nbGlobalEdge = mesh.numGlobalEdges();
    UInt nbGlobalVertex = mesh.numGlobalVertices();
    UInt nbGlobalFace = mesh.numGlobalFaces();

    UInt i, l, ie;

    // Total number of degree of freedom for each element
    UInt nldof = nbLocalDofPerVolume
        + nbLocalDofPerEdge * M_nbLocalEdge
        + nbLocalDofPerVertex * M_nbLocalVertex
        + nbLocalDofPerFace * M_nbLocalFace;

    // Consistency check
    ASSERT_PRE( nldof == UInt( M_elementDofPattern.nbLocalDof() ), "Something wrong in FE specification" ) ;

    // Global total of degrees of freedom
    M_totalDof = nbGlobalVolume * nbLocalDofPerVolume
        + nbGlobalEdge * nbLocalDofPerEdge
        + nbGlobalVertex * nbLocalDofPerVertex
        + nbGlobalFace * nbLocalDofPerFace;

    // Reshape the container to fit the needs
    M_localToGlobal.reshape( nldof, M_numElement );

    // Make sure the mesh has everything needed
    bool update_edges( nbLocalDofPerEdge != 0 && ! mesh.hasLocalEdges() );
    if ( update_edges )
    {
        mesh.updateElementEdges();
    }

    bool update_faces( nbLocalDofPerFace != 0 && ! mesh.hasLocalFaces() );
    if ( update_faces )
    {
        mesh.updateElementFaces();
    }


    UInt gcount( M_offset );
    UInt lcount;
    UInt lc;

    // Vertex Based Dof
    M_dofPositionByEntity[ 0 ] = gcount;
    if ( nbLocalDofPerVertex > 0 )
        for ( ie = 1; ie <= M_numElement; ++ie )//for each element
        {
            lc = 0;
            for ( i = 1; i <= M_nbLocalVertex; ++i )//for each vertex in the element
                for ( l = 0; l < nbLocalDofPerVertex; ++l )//for each degree of freedom per vertex
                {
                    //                       label of the ith point of the mesh element-1
                    M_localToGlobal( ++lc, ie ) = gcount + ( mesh.element( ie ).point( i ).id() - 1 ) * nbLocalDofPerVertex + l;
                    //M_localToGlobal(++lc, ie) is the global label assigned to the ++lc dof of the element.
                }
        }
    // Edge Based Dof
    gcount += nbLocalDofPerVertex * nbGlobalVertex;//dof per vertex * total # vertices
    lcount = nbLocalDofPerVertex * M_nbLocalVertex;
    M_dofPositionByEntity[ 1 ] = gcount;
    if ( nbLocalDofPerEdge > 0 )
        for ( ie = 1; ie <= M_numElement; ++ie )
        {
            lc = lcount;
            for ( i = 1; i <= M_nbLocalEdge; ++i )
                for ( l = 0; l < nbLocalDofPerEdge; ++l )
                {
                    UInt eID = mesh.edgeList(mesh.localEdgeId(ie, i)).id();
                    M_localToGlobal( ++lc, ie ) = gcount + ( eID - 1 ) * nbLocalDofPerEdge + l;

                }
        }
    // Face  Based Dof

    gcount += nbGlobalEdge * nbLocalDofPerEdge;
    lcount += nbLocalDofPerEdge * M_nbLocalEdge;
    M_dofPositionByEntity[ 2 ] = gcount;
    if ( nbLocalDofPerFace > 0 )
        for ( ie = 1; ie <= M_numElement; ++ie )
        {
            lc = lcount;
#ifdef TWODIM
            // when working in 2D we simply iterate over the elements to have faces
            for ( l = 0; l < nbLocalDofPerFace; ++l )
                M_localToGlobal( ++lc, ie ) = gcount + ( ie - 1 ) * nbLocalDofPerFace + l;
#else // THREEDIM
            for ( i = 1; i <= M_nbLocalFace; ++i )
                for ( l = 0; l < nbLocalDofPerFace; ++l )
                {
                    UInt fID = mesh.faceList( mesh.localFaceId( ie, i ) ).id();
                    M_localToGlobal( ++lc, ie ) = gcount + ( fID - 1 ) * nbLocalDofPerFace + l;
                }
#endif
        }

    // Volume  Based Dof
    gcount += nbGlobalFace * nbLocalDofPerFace;
    lcount += nbLocalDofPerFace * M_nbLocalFace;
    M_dofPositionByEntity[ 3 ] = gcount;
    if ( nbLocalDofPerVolume > 0 )
        for ( ie = 1; ie <= M_numElement; ++ie )
        {
            lc = lcount;
            for ( l = 0; l < nbLocalDofPerVolume; ++l )
            {
                M_localToGlobal( ++lc, ie ) = gcount + (mesh.element( ie ).id() - 1) * nbLocalDofPerVolume + l;
            }
        }
    gcount += nbGlobalVolume * nbLocalDofPerVolume;
    M_dofPositionByEntity[ 4 ] = gcount;
    ASSERT_POS( gcount - M_offset == M_totalDof , "Something wrong in Dof Setup " << gcount << " " << M_offset << " " << M_totalDof ) ;

    if ( update_edges )
        mesh.cleanElementEdges();
    if ( update_faces )
        mesh.cleanElementFaces();

    //UPDATE of the boundary face
    M_nbFace = mesh.numFaces();

    if (M_numLocalDofByFace>0)
    {
        //UInt lfID = 0;
        for (UInt k = 1; k <= M_numElement; k++)
        {
            for (UInt j = 1; j <= M_nbLocalFace; j++)
            {
                ID fID  = mesh.faceList( mesh.localFaceId( k, j ) ).id();
                std::vector<ID> v(M_numLocalDofByFace,0);
                for (UInt i = 1; i <= M_numLocalDofByFace; i++)
                {
                    v[i-1] = M_localToGlobal( M_faceToPoint( j, i ), k );
                }
                M_globalToLocalByFace[fID] = M_localToGlobalByFace.size();
                M_localToGlobalByFace.push_back(v);
            }
        }
    }


}

}
#endif
