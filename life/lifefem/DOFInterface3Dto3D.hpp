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
    @brief Class for interfacing dofs between two 3D meshes

    @author M.A. Fernandez
    @date 00-11-2002

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    This file contains the class which may be used to update and hold the connections between the dof
    on two matching meshes.
 */

#ifndef _DOFINTERFACE3DTO3D_HH
#define _DOFINTERFACE3DTO3D_HH

#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <life/lifefem/DOFInterface.hpp>


#include <life/lifefem/refFE.hpp>
#include <life/lifefem/DOF.hpp>
#include <iostream>
#include <map>
#include <life/lifemesh/markers.hpp>
#include <life/lifefem/CurrentBoundaryFE.hpp>
#include <ext/slist>
#include <cmath>

namespace LifeV
{
/*!
  \class DofInterface3Dto3D

  Base class which holds the connections of the dof in two matching meshes

  In order to hold the interface connections the user must give the
  RefFE elements and Dof used in both meshes.  The connections may be
  built by calling the update method. An interpolate method has been
  provided in order to interpolate data at the interface. Finally
  method getInterfaceDof gives the connections.

*/
class DofInterface3Dto3D:
        public DofInterfaceBase
{
public:
    //! @name Public typedefs
    //@{
    typedef boost::numeric::ublas::vector<Real> Vector;
    //@}

    //! @name Constructor & Destructor
    //@{

    //! default constructor
    DofInterface3Dto3D()
            :
            M_refFE1( 0 ),
            M_dof1( 0 ),
            M_refFE2( 0 ),
            M_dof2( 0 )
    {}

    //! Constructor for interfacing Dof of the same type (RefFE)
    /*!
     \param refFe the reference FE used in both meshes
     \param dof1 the Dof object of the mesh in which we want to make the computations
     \param dof2 the Dof object of the mesh which provides the data at the interface
    */
    DofInterface3Dto3D( const RefFE& refFE, const Dof& dof1, const Dof& dof2 );

    //! Constructor for interfacing Dof of different type (RefFE)
    /*!
      \param refFe1 the reference FE used in the mesh in which we want to make the computations
      \param dof1 the Dof object of the mesh in which we want to make the computations
      \param refFe2 the reference FE used in the mesh which provides the data at the interface
      \param dof2 the Dof object of the mesh which provides the data at the interface
     */
    DofInterface3Dto3D( const RefFE& refFE1, const Dof& dof1, const RefFE& refFE2, const Dof& dof2 );

    //@}


    //! @name Methods
    //@{

    //! Setup method for common finite element
    void setup(const RefFE& refFE, const Dof& dof1, const Dof& dof2 );

    //! Setup method for possibly different finite element
    void setup(const RefFE& refFE1, const Dof& dof1, const RefFE& refFE2, const Dof& dof2 );

    //! This method builds the Dof connections at the interface
    /*!
      \param mesh1 the mesh in which we want to make the computations
      \param flag1 the marker of the interface in the mesh1
      \param mesh2 the mesh which provides de data at the interface
      \param flag2 the marker of the interface in the mesh2
      \param tol tolerance for connecting points of both meshes at the interface
      \param flag3 the marker of a region of interface in the mesh1
      \brief{The parameter flag3 is used in test_meshReorder to export the part of interface determined by flag3 on mesh2.}
     */

    template <typename MeshType>
    void update( MeshType& mesh1,
                 const entityFlag_Type& flag1,
                 MeshType& mesh2,
                 const entityFlag_Type& flag2,
                 const Real& tol,
                 Int const* const flag3 = 0 );


    //! This method interpolate data when using different FE
    /*!
      \param mesh2 the mesh which provides de data at the interface
      \param v the data vector on mesh2 holding dofs of type refFE1
      \param vI the interpolated data vector on mesh2 holding dofs of type refFE2

      \note We should use this method ONLY when the accuracy of the data is less that
            the accuracy of the unknowns on mesh1, i.e., when refFE1 is more accurate
     than refFE2
     */
    template <typename MeshType>
    void interpolate( MeshType& mesh2, const UInt nbComp, const Vector& v, Vector& vI );

    //@}


    //! @name Get Methods
    //@{

    //! This method returns the corresponding dof object when interpolation is used
    const UInt& numTotalDof() const {return M_dof->numTotalDof();}

    //@}

private:

    //!  STL iterator type for the lists
    typedef __gnu_cxx::slist< std::pair<ID, ID> >::iterator Iterator;


    //! @name Private Methods
    //@{

    //! This method builds the connections between faces at the interface (M_faceToFaceConnectionList container)
    /*!
      \param mesh1 the mesh in which we want to make the computations
      \param flag1 the marker of the interface in the mesh1
      \param mesh2 the mesh which provides de data at the interface
      \param flag2 the marker of the interface in the mesh2
      \param tol tolerance for connecting points of both meshes at the interface
     */
    template <typename MeshType>
    void updateFaceConnections( const MeshType& mesh1, const entityFlag_Type& flag1,
                                 const MeshType& mesh2, const entityFlag_Type& flag2, const Real& tol );
    //! This method builds the connections between Dof at the interface (M_dofToDofConnectionList container)
    /*!
      \param mesh1 the mesh in which we want to make the computations
      \param dof1 the Dof object of the mesh in which we want to make the computations
      \param mesh2 the mesh which provides de data at the interface
      \param dof2 the Dof object of the mesh which provides de data at the interface
      \param tol tolerance for connecting points of both meshes at the interface
     */
    template <typename MeshType>
    void updateDofConnections( const MeshType& mesh1, const Dof& dof1,
                                const MeshType& mesh2, const Dof& dof2, const Real& tol,
                                Int const* const flag1 = 0 );

    //@}


    //! RefFE object used in the mesh in which we want to make the computations
    const RefFE* M_refFE1;

    //! Dof object of the mesh in which we want to make the computations
    const Dof* M_dof1;

    //! RefFE object used in the mesh which provides de data at the interface
    const RefFE* M_refFE2;

    //! Dof object of the mesh which provides the data at the interface
    const Dof* M_dof2;

    //! Auxiliary Dof object of the mesh which provides the local to global table
    //! when interpolation is used
    boost::shared_ptr<Dof> M_dof;

    //! STL list which holds the connections between faces at the interface
    __gnu_cxx::slist< std::pair<ID, ID> > M_faceToFaceConnectionList;

    //!  Auxiliary STL list which holds the connections between Dof at the interface
    //! Empty after calling update
    __gnu_cxx::slist< std::pair<ID, ID> > M_dofToDofConnectionList;

};

// ===================================================
// Helpers
// ===================================================

//! Returns true if the vectors v1 and v2 are equal with respect to the tolerance tol (in norm 1)
bool coincide( const KN_<Real>& v1, const KN_<Real>& v2, const Real& tol );


//! Returns true if points (x1,y1,z1) and (x2,y2,z2) are equal with respect to the tolerance tol (in norm 1)
bool coincide( const Real& x1, const Real& y1, const Real& z1, const Real& x2, const Real& y2, const Real& z2, const Real& tol );

// ===================================================
// Methods
// ===================================================

template <typename MeshType>
void DofInterface3Dto3D::update( MeshType& mesh1, const entityFlag_Type& flag1,
                                 MeshType& mesh2, const entityFlag_Type& flag2,
                                 const Real& tol, Int const* const flag3 )
{

    // Updating face connections at the interface
    updateFaceConnections( mesh1, flag1, mesh2, flag2, tol );

    if ( M_refFE1->nbDof() > M_refFE2->nbDof() )
    {
        // Update of the Dof connections when we need interpolation
        M_dof->update( mesh2 ); // Building auxiliary dof
        updateDofConnections( mesh1, *M_dof1, mesh2, *M_dof, tol, flag3 ); // Update of the Dof connections
    }
    else
        // Update of the Dof connections without interpolation
        updateDofConnections( mesh1, *M_dof1, mesh2, *M_dof2, tol, flag3 );
}


template <typename MeshType>
void DofInterface3Dto3D::interpolate( MeshType& mesh2, const UInt nbComp, const Vector& v, Vector& vI )
{

    typedef typename MeshType::VolumeShape GeoShape; // Element shape
    typedef typename MeshType::FaceShape GeoBShape;  // Face Shape

    UInt nbVertexPerFace = GeoBShape::S_numVertices; // Number of face's vertices
    UInt nbEdgePerFace = GeoBShape::S_numEdges;    // Number of face's edges

    UInt nbDofPerVertex = M_refFE1->nbDofPerVertex; // number of Dof per vertices
    UInt nbDofPerEdge = M_refFE1->nbDofPerEdge;   // number of Dof per edges
    UInt nbDofPerFace = M_refFE1->nbDofPerFace;   // number of Dof per faces

    UInt nbVertexPerElement = GeoShape::S_numVertices; // Number of element's vertices
    UInt nbEdgePerElement = GeoShape::S_numEdges;    // Number of element's edges

    UInt nbDofPerElement = M_refFE2->nbDof(); // Number of Dof per element in the lowDof mesh

    UInt nbVertexDofPerElement = nbVertexPerElement * nbDofPerVertex; // number of vertex's Dof on a Element
    UInt nbEdgeDofPerElement = nbEdgePerElement * nbDofPerEdge; // number of edge's Dof on a Element

    ID ibF, iElAd, iFaEl, iVeEl, lDof, iEdEl;

    Real x, y, z, sum;

    KN<Real> vLoc( nbDofPerElement * nbComp );


    // Loop on faces at the interface (matching faces)
    for ( Iterator i = M_faceToFaceConnectionList.begin(); i != M_faceToFaceConnectionList.end(); ++i )
    {

        ibF = i->second; // Face number at the interface

        iElAd = mesh2.boundaryFace( ibF ).firstAdjacentElementIdentity();  // id of the element adjacent to the face
        iFaEl = mesh2.boundaryFace( ibF ).firstAdjacentElementPosition(); // local id of the face in its adjacent element

        // Updating the local dof of the data vector in the adjacent element
        for ( UInt icmp = 0; icmp < nbComp; ++icmp )
            for ( ID idof = 0; idof < nbDofPerElement; ++idof )
                vLoc( icmp * nbDofPerElement + idof ) = v ( icmp * M_dof2->numTotalDof() + M_dof2->localToGlobal( iElAd, idof + 1 ) - 1 );

        // Vertex based Dof
        if ( nbDofPerVertex )
        {

            // loop on face vertices
            for ( ID iVeFa = 1; iVeFa <= nbVertexPerFace; ++iVeFa )
            {

                iVeEl = GeoShape::faceToPoint( iFaEl, iVeFa ); // local vertex number (in element)

                // Loop number of Dof per vertex
                for ( ID l = 1; l <= nbDofPerVertex; ++l )
                {
                    lDof = ( iVeEl - 1 ) * nbDofPerVertex + l; // Local dof in the adjacent Element

                    // Nodal coordinates
                    x = M_refFE1->xi( lDof - 1 );
                    y = M_refFE1->eta( lDof - 1 );
                    z = M_refFE1->zeta( lDof - 1 );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {

                        // Interpolating data at the nodal point
                        sum = 0;
                        for ( ID idof = 0; idof < nbDofPerElement; ++idof )  // Loop on local Dof on the adjacent element
                            sum += vLoc( icmp * nbDofPerElement + idof ) * M_refFE2->phi( idof, x, y, z );

                        // Updating interpolating vector
                        vI ( icmp * M_dof->numTotalDof() + M_dof->localToGlobal( iElAd, lDof ) - 1 ) = sum;
                    }
                }
            }
        }

        // Edge based Dof
        if ( nbDofPerEdge )
        {

            // loop on face edges
            for ( ID iEdFa = 1; iEdFa <= nbEdgePerFace; ++iEdFa )
            {

                iEdEl = GeoShape::faceToEdge( iFaEl, iEdFa ).first; // local edge number (in element)

                // Loop number of Dof per edge
                for ( ID l = 1; l <= nbDofPerEdge; ++l )
                {
                    lDof = nbVertexDofPerElement + ( iEdEl - 1 ) * nbDofPerEdge + l; // Local dof in the adjacent Element

                    // Nodal coordinates
                    x = M_refFE1->xi( lDof - 1 );
                    y = M_refFE1->eta( lDof - 1 );
                    z = M_refFE1->zeta( lDof - 1 );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {

                        // Interpolating data at the nodal point
                        sum = 0;
                        for ( ID idof = 0; idof < nbDofPerElement; ++idof )   // Loop on local Dof on the adjacent element
                            sum += vLoc( icmp * nbDofPerElement + idof ) * M_refFE2->phi( idof, x, y, z );

                        // Updating interpolating vector
                        vI ( icmp * M_dof->numTotalDof() + M_dof->localToGlobal( iElAd, lDof ) - 1 ) = sum;
                    }
                }
            }
        }

        // Loop on number of Dof per face
        for ( ID l = 1; l <= nbDofPerFace; ++l )
        {
            lDof = nbEdgeDofPerElement + nbVertexDofPerElement + ( iFaEl - 1 ) * nbDofPerFace + l; // Local dof in the adjacent Element

            // Nodal coordinates
            x = M_refFE1->xi( lDof - 1 );
            y = M_refFE1->eta( lDof - 1 );
            z = M_refFE1->zeta( lDof - 1 );

            // Loop on data vector components
            for ( UInt icmp = 0; icmp < nbComp; ++icmp )
            {

                // Interpolating data at the nodal point
                sum = 0;
                for ( ID idof = 0; idof < nbDofPerElement; ++idof )  // Loop on local Dof on the adjacent element
                    sum += vLoc( icmp * nbDofPerElement + idof ) * M_refFE2->phi( idof, x, y, z );

                // Updating interpolating vector
                vI ( icmp * M_dof->numTotalDof() + M_dof->localToGlobal( iElAd, lDof ) - 1 ) = sum;
            }
        }
    }
}

// ===================================================
// Private Methods
// ===================================================

template <typename MeshType>
void DofInterface3Dto3D::updateFaceConnections( const MeshType& mesh1, const entityFlag_Type& flag1,
                                                 const MeshType& mesh2, const entityFlag_Type& flag2, const Real& tol )
{

    UInt bdnF1 = mesh1.numBFaces(); // Number of boundary faces in mesh1
    UInt bdnF2 = mesh2.numBFaces(); // Number of boundary faces in mesh2

    entityFlag_Type marker1, marker2;

    typedef typename MeshType::FaceShape GeoBShape; // Shape of the faces

    UInt nbVertexPerFace = GeoBShape::S_numVertices; // Number of face's vertices

    KN<Real> v1( nbVertexPerFace * nDimensions ), v2( nDimensions );

    // Loop on boundary faces on mesh1
    for ( ID ibF1 = 1; ibF1 <= bdnF1; ++ibF1 )
    {

        // The face marker
        marker1 = mesh1.boundaryFace( ibF1 ).marker();

        // Is the face on the interface?
        if ( marker1 == flag1 )
        {

            // Loop on face vertices
            for ( ID iVeFa = 1; iVeFa <= nbVertexPerFace; ++iVeFa )
            {

                // Loop on vertex coordinates
                for ( ID j = 1; j <= nDimensions; ++j )
                    v1[ j - 1 + ( iVeFa - 1 ) * nDimensions ] = mesh1.boundaryFace( ibF1 ).point( iVeFa ).coordinate( j );
            }

            // Loop on boundary faces on mesh2
            for ( ID ibF2 = 1; ibF2 <= bdnF2; ++ibF2 )
            {

                // The face marker
                marker2 = mesh2.boundaryFace( ibF2 ).marker();

                // Is the face on the interface?
                if ( marker2 == flag2 )
                {

                    UInt vertexOk = 0; // Number of matched vertices

                    // Loop on face vertices
                    for ( ID iVeFa = 1; iVeFa <= nbVertexPerFace; ++iVeFa )
                    {

                        // Loop on vertex coordinates
                        for ( ID j = 1; j <= nDimensions; ++j )
                            v2[ j - 1 ] = mesh2.boundaryFace( ibF2 ).point( iVeFa ).coordinate( j );

                        // Loop on face vertices on mesh1
                        for ( ID ivefa = 1; ivefa <= nbVertexPerFace; ++ivefa )
                        {
                            // Do the vertices match?
                            if ( coincide( v1( SubArray( nDimensions, ( ivefa - 1 ) * nDimensions ) ), v2, tol ) )
                            {
                                ++vertexOk;
                                break; // Stop loop on face vertices
                            }
                        }
                    }
                    //! Do the faces match?
                    if ( vertexOk == nbVertexPerFace )
                    {
                        std::pair<ID, ID> elc( ibF1, ibF2 );
                        M_faceToFaceConnectionList.push_front( elc );
                        break; // Stop loop on boundary faces on mesh2
                    }
                }
            }
        }
    }
}


//! This method builds the connections between Dof at the interface (M_dofToDofConnectionList container)
/*!
  \param mesh1 the mesh in which we want to make the computations
  \param dof1 the Dof object of the mesh in which we want to make the computations
  \param mesh2 the mesh which provides the data at the interface
  \param dof2 the Dof object of the mesh which provides the data at the interface
  \param tol tolerance for connecting points of both meshes at the interface
  \param flag1 the marker of a region of interface in the mesh1
  \brief{The parameter flag1 is used in test_meshReorder to export the part of interface determined by flag1 on mesh2.}

*/
template <typename Mesh>
void DofInterface3Dto3D::updateDofConnections( const Mesh& mesh1, const Dof& dof1,
                                                const Mesh& mesh2, const Dof& dof2, const Real& tol, Int const* const flag1)
{

    typedef typename Mesh::VolumeShape GeoShape;
    typedef typename Mesh::FaceShape GeoBShape;

    UInt nbVertexPerFace = GeoBShape::S_numVertices; // Number of face's vertices
    UInt nbEdgePerFace = GeoBShape::S_numEdges;    // Number of face's edges

    UInt nbDofPerVertex1 = M_refFE1->nbDofPerVertex(); // number of Dof per vertices on mesh1
    UInt nbDofPerEdge1 = M_refFE1->nbDofPerEdge();   // number of Dof per edges on mesh1
    UInt nbDofPerFace1 = M_refFE1->nbDofPerFace();   // number of Dof per faces on mesh1

    UInt nbDofPerVertex2 = M_refFE2->nbDofPerVertex(); // number of Dof per vertices on mesh2
    UInt nbDofPerEdge2 = M_refFE2->nbDofPerEdge();   // number of Dof per edges on mesh2
    UInt nbDofPerFace2 = M_refFE2->nbDofPerFace();   // number of Dof per faces on mesh2

    UInt nbVertexPerElement = GeoShape::S_numVertices; // Number of element's vertices
    UInt nbEdgePerElement = GeoShape::S_numEdges;    // Number of element's edges

    UInt nbVertexDofPerFace1 = nbDofPerVertex1 * nbVertexPerFace; // number of vertex's Dof on a face on mesh1
    UInt nbEdgeDofPerFace1 = nbDofPerEdge1 * nbEdgePerFace; // number of edge's Dof on a face on mesh1

    UInt nbVertexDofPerFace2 = nbDofPerVertex2 * nbVertexPerFace; // number of vertex's Dof on a face on mesh2
    UInt nbEdgeDofPerFace2 = nbDofPerEdge2 * nbEdgePerFace; // number of edge's Dof on a face on mesh2

    UInt nbVertexDofPerElement1 = nbVertexPerElement * nbDofPerVertex1; // number of vertex's Dof on a Element on mesh1
    UInt nbEdgeDofPerElement1 = nbEdgePerElement * nbDofPerEdge1; // number of edge's Dof on a Element on mesh1

    UInt nbVertexDofPerElement2 = nbVertexPerElement * nbDofPerVertex2; // number of vertex's Dof on a Element on mesh2
    UInt nbEdgeDofPerElement2 = nbEdgePerElement * nbDofPerEdge2; // number of edge's Dof on a Element on mesh2

    ID iElAd1, iVeEl1, iFaEl1, iEdEl1, iElAd2, iVeEl2, iFaEl2, iEdEl2, lDof1, lDof2, gDof1, gDof2;

    Real x1, x2, y1, y2, z1, z2;

    bool test = false;

    CurrentBdFE feBd1( M_refFE1->boundaryFE(), getGeoMap( mesh1 ).boundaryMap() );
    CurrentBdFE feBd2( M_refFE2->boundaryFE(), getGeoMap( mesh2 ).boundaryMap() );

    // Loop on faces at the interface (matching faces)
    for ( Iterator i = M_faceToFaceConnectionList.begin(); i != M_faceToFaceConnectionList.end(); ++i )
    {

        feBd1.update( mesh1.boundaryFace( i->first ) );  // Updating face information on mesh1
        feBd2.update( mesh2.boundaryFace( i->second ) );  // Updating face information on mesh2

        iElAd1 = mesh1.boundaryFace( i->first ).firstAdjacentElementIdentity();  // id of the element adjacent to the face (mesh1)
        iElAd2 = mesh2.boundaryFace( i->second ).firstAdjacentElementIdentity();  // id of the element adjacent to the face (mesh2)

        iFaEl1 = mesh1.boundaryFace( i->first ).firstAdjacentElementPosition(); // local id of the face in its adjacent element (mesh1)
        iFaEl2 = mesh2.boundaryFace( i->second ).firstAdjacentElementPosition(); // local id of the face in its adjacent element (mesh2)

        // Vertex based Dof on mesh1
        if ( nbDofPerVertex1 )
        {

            // loop on face vertices (mesh1)
            for ( ID iVeFa1 = 1; iVeFa1 <= nbVertexPerFace; ++iVeFa1 )
            {

                iVeEl1 = GeoShape::faceToPoint( iFaEl1, iVeFa1 ); // local vertex number (in element)

                if ( flag1 != 0 && Int(mesh1.boundaryFace(i->first).point(iVeFa1).marker()) != *flag1) continue;

                // Loop number of Dof per vertex (mesh1)
                for ( ID l = 1; l <= nbDofPerVertex1; ++l )
                {
                    lDof1 = ( iVeFa1 - 1 ) * nbDofPerVertex1 + l ; // local Dof
                    feBd1.coorMap( x1, y1, z1, feBd1.refFE.xi( lDof1 - 1 ), feBd1.refFE.eta( lDof1 - 1 ) ); // Nodal coordinates on the current face (mesh1)

                    // loop on face vertices (mesh2)
                    for ( ID iVeFa2 = 1; iVeFa2 <= nbVertexPerFace; ++iVeFa2 )
                    {

                        iVeEl2 = GeoShape::faceToPoint( iFaEl2, iVeFa2 ); // local vertex number (in element)

                        // Loop on number of Dof per vertex (mesh2)
                        for ( ID k = 1; k <= nbDofPerVertex2; ++k )
                        {
                            lDof2 = ( iVeFa2 - 1 ) * nbDofPerVertex2 + k ; // local Dof
                            feBd2.coorMap( x2, y2, z2, feBd2.refFE.xi( lDof2 - 1 ), feBd2.refFE.eta( lDof2 - 1 ) ); // Nodal coordinates on the current face (mesh2)

                            // Do the nodal points match?
                            if ( (test = coincide( x1, y1, z1, x2, y2, z2, tol )) )
                            {
                                gDof1 = dof1.localToGlobal( iElAd1, ( iVeEl1 - 1 ) * nbDofPerVertex1 + l ); // Global Dof on mesh1
                                gDof2 = dof2.localToGlobal( iElAd2, ( iVeEl2 - 1 ) * nbDofPerVertex2 + k ); // Global Dof on mesh2
                                std::pair<ID, ID> locDof( gDof1, gDof2 );
                                M_dofToDofConnectionList.push_front( locDof ); // Updating the list of dof connections
                                break;
                            }
                        }
                        // Exit the loop on face vertices on mesh2?
                        if ( test )
                        {
                            test = false;
                            break;
                        }
                    }
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

                    lDof1 = nbVertexDofPerFace1 + ( iEdFa1 - 1 ) * nbDofPerEdge1 + l ; // local Dof
                    feBd1.coorMap( x1, y1, z1, feBd1.refFE.xi( lDof1 - 1 ), feBd1.refFE.eta( lDof1 - 1 ) ); // Nodal coordinates on the current face (mesh1)

                    // loop on face edges (mesh2)
                    for ( ID iEdFa2 = 1; iEdFa2 <= nbVertexPerFace; ++iEdFa2 )
                    {

                        iEdEl2 = GeoShape::faceToEdge( iFaEl2, iEdFa2 ).first; // local edge number (in element)

                        // Loop number of Dof per edge (mesh1)
                        for ( ID k = 1; k <= nbDofPerEdge2; ++k )
                        {

                            lDof2 = nbVertexDofPerFace2 + ( iEdFa2 - 1 ) * nbDofPerEdge2 + k; // local Dof
                            feBd2.coorMap( x2, y2, z2, feBd2.refFE.xi( lDof2 - 1 ), feBd2.refFE.eta( lDof2 - 1 ) ); // Nodal coordinates on the current face (mesh2)

                            // Do the nodal points match?
                            if ( (test = coincide( x1, y1, z1, x2, y2, z2, tol )) )
                            {
                                gDof1 = dof1.localToGlobal( iElAd1, nbVertexDofPerElement1 + ( iEdEl1 - 1 ) * nbDofPerEdge1 + l ); // Global Dof on mesh1
                                gDof2 = dof2.localToGlobal( iElAd2, nbVertexDofPerElement2 + ( iEdEl2 - 1 ) * nbDofPerEdge2 + k ); // Global Dof on mesh2
                                std::pair<ID, ID> locDof( gDof1, gDof2 );
                                M_dofToDofConnectionList.push_front( locDof ); // Updating the list of dof connections
                                break;
                            }
                        }
                        // Exit the loop on face edges on mesh2?
                        if ( test )
                        {
                            test = false;
                            break;
                        }
                    }
                }
            }
        }

        // Face based Dof on mesh1
        for ( ID l = 1; l <= nbDofPerFace1; ++l )
        {
            lDof1 = nbEdgeDofPerFace1 + nbVertexDofPerFace1 + l; // local Dof
            feBd1.coorMap( x1, y1, z1, feBd1.refFE.xi( lDof1 - 1 ), feBd1.refFE.eta( lDof1 - 1 ) ); // Nodal coordinates on the current face (mesh1)

            for ( ID k = 1; k <= nbDofPerFace2; ++k )
            {
                lDof2 = nbEdgeDofPerFace2 + nbVertexDofPerFace2 + k; // local Dof
                feBd2.coorMap( x2, y2, z2, feBd2.refFE.xi( lDof2 - 1 ), feBd2.refFE.eta( lDof2 - 1 ) ); // Nodal coordinates on the current face (mesh2)

                // Do the nodal points match?
                if ( coincide( x1, y1, z1, x2, y2, z2, tol ) )
                {
                    gDof1 = dof1.localToGlobal( iElAd1, nbEdgeDofPerElement1 + nbVertexDofPerElement1 + ( iFaEl1 - 1 ) * nbDofPerFace1 + l ); // Global Dof in mesh1
                    gDof2 = dof2.localToGlobal( iElAd2, nbEdgeDofPerElement2 + nbVertexDofPerElement2 + ( iFaEl2 - 1 ) * nbDofPerFace2 + k ); // Global Dof in mesh2
                    std::pair<ID, ID> locDof( gDof1, gDof2 );
                    M_dofToDofConnectionList.push_front( locDof ); // Updating the list of dof connections
                    break;
                }
            }
        }
    }

    // Updating the map containter with the connections
    for ( Iterator i = M_dofToDofConnectionList.begin(); i != M_dofToDofConnectionList.end(); ++i )
    {
        M_localDofMap[ i->first ] = i->second;
    }

    // Saving memory
    M_dofToDofConnectionList.clear();
}


}
#endif
