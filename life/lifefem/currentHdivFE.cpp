/* -*- mode: c++ -*-
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
#include <life/lifefem/currentHdivFE.hpp>

namespace LifeV
{
CurrentHdivFE::CurrentHdivFE( const RefHdivFE& _refHdivFE, const GeoMap& _geoMap, const QuadRule& _qr ) :
        nbGeoNode( _geoMap.nbDof ), nbNode( _refHdivFE.nbDof ), nbCoor( _refHdivFE.nbCoor ),
        nbQuadPt( _qr.nbQuadPt ), nbDiag( _refHdivFE.nbDiag() ),
        nbUpper( _refHdivFE.nbUpper() ) , nbPattern( _refHdivFE.nbPattern() ),
        point( nbGeoNode, nbCoor ),
        refHdivFE( _refHdivFE ), geoMap( _geoMap ), qr( _qr ),
        phi( nbNode, nbCoor, nbQuadPt ), phiRef( nbNode, nbCoor, nbQuadPt ),
        divPhi( nbNode, nbQuadPt ),
        jacobian( nbCoor, nbCoor, nbQuadPt ),
        phiGeo( nbGeoNode, nbQuadPt ), dPhiGeo( nbGeoNode, nbCoor, nbQuadPt ),
        weightDet( nbQuadPt ), detJac( nbQuadPt )
{
    CONSTRUCTOR( "CurrentHdivFE" );
    for ( int ig = 0;ig < nbQuadPt;ig++ )
    {
        for ( int i = 0;i < nbNode;i++ )
        {
            divPhi( i, ig ) = refHdivFE.divPhi( i, ig, qr );
            for ( int icoor = 0;icoor < nbCoor;icoor++ )
            {
                phiRef( i, icoor, ig ) = refHdivFE.phi( i, icoor, ig, qr );
            }
        }
        for ( int k = 0;k < nbGeoNode;k++ )
        {
            phiGeo( k, ig ) = geoMap.phi( k, ig, qr );
            for ( int icoor = 0;icoor < nbCoor;icoor++ )
            {
                dPhiGeo( k, icoor, ig ) = geoMap.dPhi( k, icoor, ig, qr );
            }
        }
    }
}

CurrentHdivFE::~CurrentHdivFE()
{
    DESTRUCTOR( "CurrentHdivFE" );
}

//----------------------------------------------------------------------
void CurrentHdivFE::coorMap( Real& x, Real& y, Real& z,
                             const Real & xi, const Real & eta, const Real & zeta ) const
{
    x = y = z = 0.;
    for ( int i = 0;i < nbGeoNode;i++ )
    {
        x += point( i, 0 ) * geoMap.phi( i, xi, eta, zeta );
        y += point( i, 1 ) * geoMap.phi( i, xi, eta, zeta );
#if defined(THREEDIM)

        z += point( i, 2 ) * geoMap.phi( i, xi, eta, zeta );
#endif

    }
}
//----------------------------------------------------------------------
Real CurrentHdivFE::measure() const
{
    Real meas = 0.;
    for ( int ig = 0;ig < nbQuadPt;ig++ )
        meas += weightDet( ig );
    return meas;
}
//----------------------------------------------------------------------
void CurrentHdivFE::_comp_jacobian()
{
    Real fctDer;
    // derivatives of geo map:
    for ( int ig = 0;ig < nbQuadPt;ig++ )
    {
        for ( int icoor = 0;icoor < nbCoor;icoor++ )
        {
            for ( int jcoor = 0;jcoor < nbCoor;jcoor++ )
            {
                fctDer = 0.;
                for ( int j = 0;j < nbGeoNode;j++ )
                {
                    fctDer += point( j, icoor ) * dPhiGeo( j, jcoor, ig );
                }
                jacobian( icoor, jcoor, ig ) = fctDer;
            }
        }
    }
    // determinant on integrations points
#if defined(TWODIM)
    // *** 2D code ***
    Real a, b, c, d;
    for ( int ig = 0;ig < nbQuadPt;ig++ )
    {
        a = jacobian( 0, 0, ig );
        b = jacobian( 0, 1, ig );
        c = jacobian( 1, 0, ig );
        d = jacobian( 1, 1, ig );
        detJac( ig ) = a * d - b * c;
        weightDet( ig ) = detJac( ig ) * qr.weight( ig );
    }
#elif defined(THREEDIM)
    // *** 3D code ***
    Real a, b, c, d, e, f, g, h, i, ei, fh, bi, ch, bf, ce;
    for ( int ig = 0;ig < nbQuadPt;ig++ )
    {
        a = jacobian( 0, 0, ig );
        b = jacobian( 0, 1, ig );
        c = jacobian( 0, 2, ig );
        d = jacobian( 1, 0, ig );
        e = jacobian( 1, 1, ig );
        f = jacobian( 1, 2, ig );
        g = jacobian( 2, 0, ig );
        h = jacobian( 2, 1, ig );
        i = jacobian( 2, 2, ig );
        ei = e * i;
        fh = f * h;
        bi = b * i;
        ch = c * h;
        bf = b * f;
        ce = c * e;
        detJac( ig ) = a * ( ei - fh ) + d * ( ch - bi ) + g * ( bf - ce );
        weightDet( ig ) = detJac( ig ) * qr.weight( ig );
    }
#endif
}
}

