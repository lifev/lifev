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
/*!
  \file ElasticStructureHandler.h
  \author M.A. Fernandez
  \date 06/2003
  \version 1.0

  \brief This file contains an abstract class for Elastic Structures.

*/

#ifndef _ELASTICSTRUCTUREHANDLER_H_
#define _ELASTICSTRUCTUREHANDLER_H_

#include <cmath>
#include <sstream>
#include <cstdlib>

#include "dataElasticStructure.hpp"
#include "dataAztec.hpp"
#include "refFE.hpp"
#include "dof.hpp"
#include "lifeV.hpp"
#include "medit_wrtrs.hpp"
#include "bcCond.hpp"

namespace LifeV
{
/*!
  \class ElasticStructureHandler
*/

template <typename Mesh>
class ElasticStructureHandler:
            public DataElasticStructure<Mesh>
{

public:

    typedef Real ( *Function ) ( const Real&, const Real&, const Real&, const Real&, const ID& );

    //! Constructor
    /*!
      \param data_file GetPot data file
      \param refFE reference FE for the displacement
      \param Qr volumic quadrature rule
      \param bdQr surface quadrature
      \param BCh boundary conditions for the displacement
    */
    ElasticStructureHandler( const GetPot& data_file, const RefFE& refFE,
                             const QuadRule& Qr, const QuadRule& bdQr, BCHandler& BCh );

    //! Sets initial condition for the displacment en velocity
    void initialize( const Function& d0, const Function& w0 );

    //! Update the right  hand side  for time advancing
    /*!
      \param source volumic force
      \param time present time
    */
    virtual void timeAdvance( const Function source, const Real& time ) = 0;

    //! Solve the non-linear problem
    virtual void iterate() = 0;

    //! Returns the displacement vector
    PhysVectUnknown<Vector>& d();

    //! Returns the velocity vector
    PhysVectUnknown<Vector>& w();

    //! Returns the displacement Dof
    const Dof& dDof() const;

    //! Returns the quadrature rule
    CurrentFE& currentFE();

    //! Postprocessing
    void postProcess();

    //! Do nothing destructor
    virtual ~ElasticStructureHandler()
    {}

protected:

    //! Reference FE
    const RefFE& _refFE;

    //! The Dof object
    Dof _dof;

    //! The number of total displacement dofs
    UInt _dim;

    //! Quadrature rule for volumic elementary computations
    const QuadRule& _Qr;

    //! Quadrature rule for elementary computations
    const QuadRule& _bdQr;

    //! Current FE
    CurrentFE _fe;

    //! Current boundary FE
    CurrentBdFE _feBd;

    //! The displacement
    PhysVectUnknown<Vector> _d;

    //! The velocity
    PhysVectUnknown<Vector> _w;

    //! The BC handler
    BCHandler& _BCh;

    //! The actual time
    Real _time;

    //! Aux. var. for PostProc
    UInt _count;
};



//
// IMPLEMENTATION
//


// Constructor
template <typename Mesh>
ElasticStructureHandler<Mesh>::
ElasticStructureHandler( const GetPot& data_file, const RefFE& refFE,
                         const QuadRule& Qr, const QuadRule& bdQr, BCHandler& BCh )
        :
        DataElasticStructure<Mesh>( data_file ),
        _refFE( refFE ),
        _dof( this->_mesh, _refFE ),
        _dim( _dof.numTotalDof() ),
        _Qr( Qr ),
        _bdQr( bdQr ),
        _fe( _refFE, getGeoMap( this->_mesh ), _Qr ),
        _feBd( _refFE.boundaryFE(), getGeoMap( this->_mesh ).boundaryMap(), _bdQr ),
        _d( _dim ),
        _w( _dim ),
        _BCh( BCh ),
        _time( 0 ),
        _count( 0 )
{}


// Returns the displacement vector
template <typename Mesh>
PhysVectUnknown<Vector>&
ElasticStructureHandler<Mesh>::d()
{
    return _d;
}


// Returns the velocity vector
template <typename Mesh>
PhysVectUnknown<Vector>&
ElasticStructureHandler<Mesh>::w()
{
    return _w;
}

// Returns the velocity displacement
template <typename Mesh>
const Dof&
ElasticStructureHandler<Mesh>::dDof() const
{
    return _dof;
}

// Returns the velocity displacement
template <typename Mesh>
CurrentFE&
ElasticStructureHandler<Mesh>::currentFE()
{
    return _fe;
}

// Postprocessing
template <typename Mesh>
void
ElasticStructureHandler<Mesh>::postProcess()
{
    std::ostringstream index;
    std::string name, namedef;

    ++_count;

    if ( fmod( float( _count ), float( this->_verbose ) ) == 0.0 )
    {
        std::cout << "  o-  Post-processing \n";
        index << ( _count / this->_verbose );

        switch ( index.str().size() )
        {
        case 1:
            name = "00" + index.str();
            break;
        case 2:
            name = "0" + index.str();
            break;
        case 3:
            name = index.str();
            break;
        }

        namedef = "defor." + name + ".mesh";
        wr_medit_ascii_scalar( "dep_x." + name + ".bb", _d.giveVec(), this->_mesh.numVertices() );
        wr_medit_ascii_scalar( "dep_y." + name + ".bb", _d.giveVec() + _dim, this->_mesh.numVertices() );
        wr_medit_ascii_scalar( "dep_z." + name + ".bb", _d.giveVec() + 2 * _dim, this->_mesh.numVertices() );
        wr_medit_ascii2( namedef, this->_mesh, _d, this->_factor );
        // wr_medit_ascii_vector("veloc."+name+".bb",_u.giveVec(),this->_mesh.numVertices(),_dim_u);
        system( ( "ln -s " + namedef + " dep_x." + name + ".mesh" ).data() );
        system( ( "ln -s " + namedef + " dep_y." + name + ".mesh" ).data() );
        system( ( "ln -s " + namedef + " dep_z." + name + ".mesh" ).data() );
        // system(("ln -s "+this->_mesh_file+" veloc."+name+".mesh").data());
    }
}


// Sets the initial condition
template <typename Mesh>
void
ElasticStructureHandler<Mesh>::initialize( const Function& d0, const Function& w0 )
{


    // Initialize velocity

    typedef typename Mesh::VolumeShape GeoShape; // Element shape

    UInt nDofpV = _refFE.nbDofPerVertex; // number of Dof per vertex
    UInt nDofpE = _refFE.nbDofPerEdge;   // number of Dof per edge
    UInt nDofpF = _refFE.nbDofPerFace;   // number of Dof per face
    UInt nDofpEl = _refFE.nbDofPerVolume; // number of Dof per Volume

    UInt nElemV = GeoShape::numVertices; // Number of element's vertices
    UInt nElemE = GeoShape::numEdges;    // Number of element's edges
    UInt nElemF = GeoShape::numFaces;    // Number of element's faces

    UInt nDofElemV = nElemV * nDofpV; // number of vertex's Dof on a Element
    UInt nDofElemE = nElemE * nDofpE; // number of edge's Dof on a Element
    UInt nDofElemF = nElemF * nDofpF; // number of face's Dof on a Element

    ID nbComp = _d.nbcomp(); // Number of components of the mesh velocity

    Real x, y, z;

    ID lDof;

    // Loop on elements of the mesh
    for ( ID iElem = 1; iElem <= this->_mesh.numVolumes(); ++iElem )
    {

        _fe.updateJac( this->_mesh.volume( iElem ) );

        // Vertex based Dof
        if ( nDofpV )
        {

            // loop on element vertices
            for ( ID iVe = 1; iVe <= nElemV; ++iVe )
            {

                // Loop number of Dof per vertex
                for ( ID l = 1; l <= nDofpV; ++l )
                {
                    lDof = ( iVe - 1 ) * nDofpV + l; // Local dof in this element

                    // Nodal coordinates
                    _fe.coorMap( x, y, z, _fe.refFE.xi( lDof - 1 ), _fe.refFE.eta( lDof - 1 ), _fe.refFE.zeta( lDof - 1 ) );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {
                        _d( icmp * _dim + _dof.localToGlobal( iElem, lDof ) - 1 ) = d0( 0.0, x, y, z, icmp + 1 );
                        _w( icmp * _dim + _dof.localToGlobal( iElem, lDof ) - 1 ) = w0( 0.0, x, y, z, icmp + 1 );
                    }
                }
            }
        }

        // Edge based Dof
        if ( nDofpE )
        {

            // loop on element edges
            for ( ID iEd = 1; iEd <= nElemE; ++iEd )
            {

                // Loop number of Dof per edge
                for ( ID l = 1; l <= nDofpE; ++l )
                {
                    lDof = nDofElemV + ( iEd - 1 ) * nDofpE + l; // Local dof in the adjacent Element

                    // Nodal coordinates
                    _fe.coorMap( x, y, z, _fe.refFE.xi( lDof - 1 ), _fe.refFE.eta( lDof - 1 ), _fe.refFE.zeta( lDof - 1 ) );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {
                        _d( icmp * _dim + _dof.localToGlobal( iElem, lDof ) - 1 ) = d0( 0.0, x, y, z, icmp + 1 );
                        _w( icmp * _dim + _dof.localToGlobal( iElem, lDof ) - 1 ) = w0( 0.0, x, y, z, icmp + 1 );
                    }
                }
            }
        }

        // Face based Dof
        if ( nDofpF )
        {

            // loop on element faces
            for ( ID iFa = 1; iFa <= nElemF; ++iFa )
            {

                // Loop on number of Dof per face
                for ( ID l = 1; l <= nDofpF; ++l )
                {

                    lDof = nDofElemE + nDofElemV + ( iFa - 1 ) * nDofpF + l; // Local dof in the adjacent Element

                    // Nodal coordinates
                    _fe.coorMap( x, y, z, _fe.refFE.xi( lDof - 1 ), _fe.refFE.eta( lDof - 1 ), _fe.refFE.zeta( lDof - 1 ) );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {
                        _d( icmp * _dim + _dof.localToGlobal( iElem, lDof ) - 1 ) = d0( 0.0, x, y, z, icmp + 1 );
                        _w( icmp * _dim + _dof.localToGlobal( iElem, lDof ) - 1 ) = w0( 0.0, x, y, z, icmp + 1 );
                    }
                }
            }
        }
        // Element based Dof
        // Loop on number of Dof per Element
        for ( ID l = 1; l <= nDofpEl; ++l )
        {
            lDof = nDofElemF + nDofElemE + nDofElemV + l; // Local dof in the Element

            // Nodal coordinates
            _fe.coorMap( x, y, z, _fe.refFE.xi( lDof - 1 ), _fe.refFE.eta( lDof - 1 ), _fe.refFE.zeta( lDof - 1 ) );

            // Loop on data vector components
            for ( UInt icmp = 0; icmp < nbComp; ++icmp )
            {
                _d( icmp * _dim + _dof.localToGlobal( iElem, lDof ) - 1 ) = d0( 0.0, x, y, z, icmp + 1 );
                _w( icmp * _dim + _dof.localToGlobal( iElem, lDof ) - 1 ) = w0( 0.0, x, y, z, icmp + 1 );
            }
        }
    }
}
}

#endif
