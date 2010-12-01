/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

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
/*!
  \file dof.h
  \brief Degrees of freedom, the class that provides the localtoglobal table
  \version 1.0
  \author M.A. Fernandez & Luca Formaggia
  \date 07/2002

  The numbering of the degrees of freedom follow the order
  Vertices - Edges - Faces - Volumes
  If more than one degree of freedon is associated to a geometrical entity, the
  global numbering associated to two "consecutive" degree of freedom  on the given
  entity is consecutive, and the order is according to the entity orientation.


  \author Modified by Vincent MARTIN & Mohamed BELHADJ
  \date 19/07/2002

  We added the HdivFE element here.
*/

#ifndef _DOF_HH
#define _DOF_HH

#include <life/lifecore/life.hpp>
#include <life/lifearray/SimpleVect.hpp>
#include <life/lifefem/localDofPattern.hpp>
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
    //! Type for the localToGlobal table.
    typedef SimpleArray<UInt> Container;

    //! The pattern of the local degrees of freedom.
    /*! It is exposed so that is is possible to interrogate it directly.*/
    const LocalDofPattern& fe; // be careful : use fe.nbLocalDof (fe.nbDof does not exist !)

    /*! The minimal constructor
      \param _fe is the LocalDofPattern on which the ref FE is built
      \param Offset: the smallest Dof numbering. It might be used if we want the
      degrees of freedom numbering start from a specific value.
    */
    Dof( const LocalDofPattern& _fe, UInt offSet = 1 );

    Dof( const Dof & dof2 );

    //! Constructor accepting a mesh as parameter
    /*!
      \param mesh a RegionMesh3D
      \param _fe is the LocalDofPattern on which the ref FE is built
      \param Offset: the smalest Dof numbering. It might be used if we want the
      degrees of freedom numbering start from a specific value.
    */
    template <typename Mesh>
    Dof( Mesh& mesh, const LocalDofPattern& _fe, UInt offSet = 1 );

    //! Build the localToGlobal table
    /*!
      \param mesh A RegionMesh3D
      Updates the LocaltoGlobal array
    */
    template <typename Mesh>
    void update( Mesh & );

    //! The total number of Dof
    inline UInt numTotalDof() const
    {
        return _totalDof;
    }

    inline void setTotalDof(const UInt totalDof)
    {
        _totalDof = totalDof;
    }

    //! The number of local Dof (nodes) in the finite element
    inline UInt numLocalDof() const
    {
        return fe.nbLocalDof();
    }

    //! Return the specified entries of the localToGlobal table
    /*!
      Returns the global numbering of a DOF, given an element and the local numbering
      \param ELId the element ID
      \param localNode the local DOF numbering (starting from 1)
      \return The numbering of the DOF
    */
    inline ID localToGlobal( const ID ElId, const ID localNode) const
    {
        return _ltg( localNode, ElId );
    }

    //! Number of elements in mesh
    UInt numElements() const
    {
        return _nEl;
    }

    //! Number of faces in the mesh
    UInt numFaces() const
    {
        return _numFaces;
    }

    //! Number of local vertices (in a elment)
    UInt numLocalVertices() const
    {
        return nlv;
    }

    //! Number of local edges (in a elment)
    UInt numLocalEdges() const
    {
        return nle;
    }

    //! Number of local faces (in a elment)
    UInt numLocalFaces() const
    {
        return nlf;
    }

    //Internal data

    //! Number of Local DofByFace
    UInt numLocalDofByFace() const
    {
        ASSERT_PRE( (_numLocalDofByFace>0) , "This data are not available for this reference element");
        return _numLocalDofByFace;
    }

    /*!
      Returns the global numbering of a DOF, given an internal face and the
      local numbering
      \param faceId the internal face ID
      \param localDOF the local DOF numbering (starting from 1)
      \return The global numbering of the DOF
    */
    ID localToGlobalByFace(const ID& faceId, const ID& localDOF, bool& exist ) const;

    //! Ouput
    void showMe( std::ostream & out = std::cout, bool verbose = false ) const;
    void showMeByFace(std::ostream& out = std::cout, bool verbose = false) const;

private:
    typedef ID ( *FTOP )( ID const localFace, ID const point );

    UInt _offset;
    UInt _totalDof;
    UInt _nEl;
    UInt nlv;
    UInt nle;
    UInt nlf;
    Container _ltg; // container is:  typedef SimpleArray<UInt> Container;

    UInt _numFaces;             // number of faces in the mesh

    std::vector<std::vector<ID> > _ltgByFace;
    //Container _ltgByFace;       // connection array that maps the local dof of
    // the face to the global dof
    std::map<ID,ID> _gtlByFace; // connection between global dof and local numbering
    FTOP _fToP;                 // local array that maps the local dof of the
    // face to the local dof the the element
    UInt _numLocalDofByFace;    // number of dof on a face
    UInt _ncount[ 5 ];
};



/********************************************************************
                      IMPLEMENTATIONS
********************************************************************/

//! Constructor that builds the localToglobal table
template <typename Mesh>
Dof::Dof( Mesh& mesh, const LocalDofPattern& _fe, UInt off ) :
        fe       ( _fe ),
        _offset  ( off ),
        _totalDof( 0 ),
        _nEl     ( 0 ),
        nlv      ( 0 ),
        nle      ( 0 ),
        nlf      ( 0 ),
        _ltg     (),
        _numFaces( 0 ),
        _ltgByFace(),
        _gtlByFace()
{
    //Getting the face
    switch ( _fe.nbLocalDof() )
    {
    case 2:
        // No _fToP (it is 1D)
        _numLocalDofByFace = 1;
        break;
    case 4:
        _fToP = LinearTetra::fToP;
        _numLocalDofByFace = 3;
        break;
    case 5:
        _fToP = LinearTetraBubble::fToP;
        _numLocalDofByFace = 3;
        break;
    case 10:
        _fToP = QuadraticTetra::fToP;
        _numLocalDofByFace = 6;
        break;
    case 8:
        _fToP = LinearHexa::fToP;
        _numLocalDofByFace = 4;
        break;
    case 27:
        _fToP = QuadraticHexa::fToP;
        _numLocalDofByFace = 27;
        break;
    default:
        std::cout << "Warning: This refFE is not available for the dof by face." << std::endl;
        _numLocalDofByFace = 0;
        break;
    }

    for ( UInt i = 0; i < 5; ++i )
        _ncount[ i ] = 0;
    update( mesh );
}


//! Build the localToGlobal table
template <typename Mesh>
void Dof::update( Mesh& M )
{

    typedef typename Mesh::ElementShape GeoShape;

    // Some useful local variables, to save some typing
    UInt nldpe = fe.nbDofPerEdge();
    UInt nldpv = fe.nbDofPerVertex();
    UInt nldpf = fe.nbDofPerFace();
    UInt nldpV = fe.nbDofPerVolume();

    nlv = GeoShape::numVertices;
    nle = GeoShape::numEdges;
    nlf = GeoShape::numFaces;

    _nEl = M.numElements();

    UInt nV = M.numGlobalVolumes();
    UInt ne = M.numGlobalEdges();
    UInt nv = M.numGlobalVertices();
    UInt nf = M.numGlobalFaces();

//    std::cout << "Num global Edges = " << ne << std::endl;
//    std::cout << M.numGlobalVolumes()  << std::endl;
//    std::cout << M.numGlobalEdges()    << std::endl;
//    std::cout << M.numGlobalVertices() << std::endl;
//    std::cout << M.numGlobalFaces()    << std::endl;

    UInt i, l, ie;

    UInt nldof = nldpV + nldpe * nle + nldpv * nlv + nldpf * nlf;

    ASSERT_PRE( nldof == UInt( fe.nbLocalDof() ), "Something wrong in FE specification" ) ;

    _totalDof = nV * nldpV + ne * nldpe + nv * nldpv + nf * nldpf;

    _ltg.reshape( nldof, _nEl );

    // Make sure the mesh has everything needed
    bool update_edges( nldpe != 0 && ! M.hasLocalEdges() );
    if ( update_edges )
        M.updateElementEdges();

#ifndef TWODIM
    bool update_faces( nldpf != 0 && ! M.hasLocalFaces() );
    if ( update_faces )
        M.updateElementFaces();
#endif

    //  ASSERT_PRE( !(nldpe !=0 && M.hasLocalEdges()) , "Element edges stuff have not been updated") ;
    //  ASSERT_PRE( !(nldpf !=0 && M.hasLocalFaces()) , "Element faces stuff have not been updated") ;
    //ASSERT_PRE( (nldpe == 0 || M.hasLocalEdges()) , "Element edges stuff have not been updated") ;
    //ASSERT_PRE( (nldpf == 0 || M.hasLocalFaces()) , "Element faces stuff have not been updated") ;


    unsigned int gcount( _offset );
    unsigned int lcount;
    unsigned int lc;

    // Vertex Based Dof
    _ncount[ 0 ] = gcount;
    if ( nldpv > 0 )
        for ( ie = 1; ie <= _nEl; ++ie )//for each element
        {
            lc = 0;
            for ( i = 1; i <= nlv; ++i )//for each vertex in the element
                for ( l = 0; l < nldpv; ++l )//for each degree of freedom per vertex
                {
                    //                       label of the ith point of the mesh element-1
                    _ltg( ++lc, ie ) = gcount + ( M.element( ie ).point( i ).id() - 1 ) * nldpv + l;
                    //_ltg(++lc, ie) is the global label assigned to the ++lc dof of the element.
                }
        }
    // Edge Based Dof
    gcount += nldpv * nv;//dof per vertex * total # vertices
    lcount = nldpv * nlv;
    _ncount[ 1 ] = gcount;
    if ( nldpe > 0 )
        for ( ie = 1; ie <= _nEl; ++ie )
        {
            lc = lcount;
            for ( i = 1; i <= nle; ++i )
                for ( l = 0; l < nldpe; ++l )
                {
                    UInt eID = M.edgeList(M.localEdgeId(ie, i)).id();
                    _ltg( ++lc, ie ) = gcount + ( eID - 1 ) * nldpe + l;

//                     std::cout << eID - 1 << " "
//                               << ie << " "
//                               <<  gcount + ( eID - 1 ) * nldpe + l << " "
//                               << M.localEdgeId(ie, i) << std::endl;

                }
        }
    // Face  Based Dof @@ possibly bugged since

    gcount += ne * nldpe;
    lcount += nldpe * nle;
    _ncount[ 2 ] = gcount;
    if ( nldpf > 0 )
        for ( ie = 1; ie <= _nEl; ++ie )
        {
            lc = lcount;
#ifdef TWODIM
            // when working in 2D we simply iterate over the elements to have faces
            for ( l = 0; l < nldpf; ++l )
                _ltg( ++lc, ie ) = gcount + ( ie - 1 ) * nldpf + l;
#else // THREEDIM
            for ( i = 1; i <= nlf; ++i )
                for ( l = 0; l < nldpf; ++l )
                {
                    UInt fID = M.faceList( M.localFaceId( ie, i ) ).id();
                    _ltg( ++lc, ie ) = gcount + ( fID - 1 ) * nldpf + l;
                }
#endif
        }

    // Volume  Based Dof
    gcount += nf * nldpf;
    lcount += nldpf * nlf;
    _ncount[ 3 ] = gcount;
    if ( nldpV > 0 )
        for ( ie = 1; ie <= _nEl; ++ie )
        {
            lc = lcount;
            for ( l = 0; l < nldpV; ++l )
            {
                _ltg( ++lc, ie ) = gcount + (M.element( ie ).id() - 1) * nldpV + l;
            }
        }
    gcount += nV * nldpV;
    _ncount[ 4 ] = gcount;
    ASSERT_POS( gcount - _offset == _totalDof , "Something wrong in Dof Setup " << gcount << " " << _offset << " " << _totalDof ) ;

    if ( update_edges )
        M.cleanElementEdges();
#ifndef TWODIM
    if ( update_faces )
        M.cleanElementFaces();
#endif

    //UPDATE of the boundary face (add the 13th november 2009)
    _numFaces = M.numFaces();

    if (_numLocalDofByFace>0)
    {
        //UInt lfID = 0;
        for (UInt k = 1; k <= _nEl; k++)
        {
            for (UInt j = 1; j <= nlf; j++)
            {
                ID fID  = M.faceList( M.localFaceId( k, j ) ).id();
                std::vector<ID> v(_numLocalDofByFace,0);
                for (UInt i = 1; i <= _numLocalDofByFace; i++)
                {
                    v[i-1] = _ltg( _fToP( j, i ), k );
                }
                _gtlByFace[fID] = _ltgByFace.size();
                _ltgByFace.push_back(v);
            }
        }
    }


}

}
#endif
