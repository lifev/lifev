/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2004-11-11

  Copyright (C) 2004 EPFL

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
/**
   \file interpolate.hpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-11-11
 */
#ifndef __interpolate_H
#define __interpolate_H 1

#include <tab.hpp>
#include <boost/function.hpp>

#warning SHOULD NOT BE USED YET

namespace LifeV
{
template <typename Mesh, typename RefFE, typename CurrFE, typename Dof>
Vector
interpolate( boost::function<Real ( node_type const&, id_type const& )> const& u0,
             Mesh& mesh,
             RefFE& refFE, CurrFE& currFE,
             Dof& dof,
             UInt nbComp = 1 )
{
    typedef typename Mesh::VolumeShape GeoShape; // Element shape

    UInt nDofpV = refFE.nbDofPerVertex; // number of Dof per vertex
    UInt nDofpE = refFE.nbDofPerEdge;   // number of Dof per edge
    UInt nDofpF = refFE.nbDofPerFace;   // number of Dof per face
    UInt nDofpEl = refFE.nbDofPerVolume; // number of Dof per Volume

    UInt nElemV = GeoShape::numVertices; // Number of element's vertices
    UInt nElemE = GeoShape::numEdges;    // Number of element's edges
    UInt nElemF = GeoShape::numFaces;    // Number of element's faces

    UInt nDofElemV = nElemV * nDofpV; // number of vertex's Dof on a Element
    UInt nDofElemE = nElemE * nDofpE; // number of edge's Dof on a Element
    UInt nDofElemF = nElemF * nDofpF; // number of face's Dof on a Element

    UInt size_comp = dof.numTotalDof();
    UInt __size = size_comp * nbComp; // Initialization of the dimension of the vector

    // result of the interpolation of u0 on the space V(mesh+refFE)
    Vector __interp( __size );
    __interp = ZeroVector( __size );

    //Real x, y, z;
    node_type x( 3 );
    id_type lDof;

    // Loop on elements of the mesh
    for ( id_type iElem = 1; iElem <= mesh.numVolumes(); ++iElem )
    {

        currFE.updateJac( mesh.volume( iElem ) );

        // Vertex based Dof
        if ( nDofpV )
        {

            // loop on element vertices
            for ( id_type iVe = 1; iVe <= nElemV; ++iVe )
            {

                // Loop number of Dof per vertex
                for ( id_type l = 1; l <= nDofpV; ++l )
                {
                    lDof = ( iVe - 1 ) * nDofpV + l; // Local dof in this element

                    // Nodal coordinates
                    currFE.coorMap( x[0], x[1], x[2], currFE.refFE.xi( lDof - 1 ),
                                    currFE.refFE.eta( lDof - 1 ),
                                    currFE.refFE.zeta( lDof - 1 ) );

                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {
                        size_t __index =  icmp * size_comp + dof.localToGlobal( iElem, lDof ) - 1;
                        __interp[ __index ] = u0( x, icmp + 1 );
                    }
                }
            }
        }

        // Edge based Dof
        if ( nDofpE )
        {

            // loop on element edges
            for ( id_type iEd = 1; iEd <= nElemE; ++iEd )
            {

                // Loop number of Dof per edge
                for ( id_type l = 1; l <= nDofpE; ++l )
                {
                    // Local dof in the adjacent Element
                    lDof = nDofElemV + ( iEd - 1 ) * nDofpE + l;

                    // Nodal coordinates
                    currFE.coorMap( x[0], x[1], x[2],
                                    currFE.refFE.xi( lDof - 1 ),
                                    currFE.refFE.eta( lDof - 1 ),
                                    currFE.refFE.zeta( lDof - 1 ) );

                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {
                        int __index =  icmp * size_comp + dof.localToGlobal( iElem, lDof ) - 1;
                        __interp[__index] = u0( x, icmp + 1 );
                    }

                }
            }
        }

        // Face based Dof
        if ( nDofpF )
        {

            // loop on element faces
            for ( id_type iFa = 1; iFa <= nElemF; ++iFa )
            {

                // Loop on number of Dof per face
                for ( id_type l = 1; l <= nDofpF; ++l )
                {

                    lDof = nDofElemE + nDofElemV + ( iFa - 1 ) * nDofpF + l; // Local dof in the adjacent Element

                    // Nodal coordinates
                    currFE.coorMap( x[0], x[1], x[2],
                                    currFE.refFE.xi( lDof - 1 ),
                                    currFE.refFE.eta( lDof - 1 ),
                                    currFE.refFE.zeta( lDof - 1 ) );

                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {
                        int __index =  icmp * size_comp + dof.localToGlobal( iElem, lDof ) - 1;
                        __interp[__index] = u0( x, icmp + 1 );
                    }

                }
            }
        }

        if ( nDofpEl )
        {
            // Element based Dof
            // Loop on number of Dof per Element
            for ( id_type l = 1; l <= nDofpEl; ++l )
            {
                lDof = nDofElemF + nDofElemE + nDofElemV + l; // Local dof in the Element

                // Nodal coordinates
                currFE.coorMap( x[0], x[1], x[2],
                                currFE.refFE.xi( lDof - 1 ),
                                currFE.refFE.eta( lDof - 1 ),
                                currFE.refFE.zeta( lDof - 1 ) );
                for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                {
                    int __index =  icmp * size_comp + dof.localToGlobal( iElem, lDof ) - 1;
                    __interp[__index] = u0( x, icmp + 1 );
                }
            }
        }
    }
    return __interp;
}
}
#endif /* __interpolate_H */
