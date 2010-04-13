/*-*- mode: c++ -*-
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
#include <life/lifefem/staticBdFE.hpp>

namespace LifeV
{

/*!
 This constructor is typically used for the currentBdFE (derived class)
*/
StaticBdFE::StaticBdFE( const RefFE& _refFE, const GeoMap& _geoMap,
                        const QuadRule& _qr ) :
    nbGeoNode ( _geoMap.nbDof() ),
    nbNode    ( _refFE.nbDof() ),
    nbCoor    ( _refFE.nbCoor() ),
    nbQuadPt  ( _qr.nbQuadPt() ),
        point     ( nbGeoNode, nbCoor + 1 ),
        refFE     ( _refFE ),
        geoMap    ( _geoMap ), qr( _qr ),
        phi       ( ( int ) nbNode,( int ) nbQuadPt ),
        dPhiRef   ( ( int ) nbNode, ( int ) nbCoor,( int ) nbQuadPt ),
        dPhi      ( ( int ) nbNode, ( int ) nbCoor, ( int ) nbQuadPt ),
        phiGeo    ( ( int ) nbGeoNode, ( int ) nbQuadPt ),
        dPhiGeo   ( ( int ) nbGeoNode, ( int ) nbCoor, ( int ) nbQuadPt ),
        weightMeas( ( int ) nbQuadPt ), meas( ( int ) nbQuadPt ),
        normal    ( ( int ) nbCoor + 1, ( int ) nbQuadPt ),
        tangent   ( ( int ) nbCoor, ( int ) nbCoor + 1, ( int ) nbQuadPt ),
        metric    ( ( int ) nbCoor, ( int ) nbCoor, ( int ) nbQuadPt ),
        quadPt    ( ( int ) nbQuadPt, 3 ),
        invArea   ( 1. )
{
    CONSTRUCTOR( "StaticBdFE" );
    for ( int ig = 0;ig < nbQuadPt;ig++ )
    {
        for ( int i = 0;i < nbNode;i++ )
        {
            //phi( i, ig ) = refFE.phi( i, ig, qr );
            phi( i, ig ) = refFE.phi( i,
                                      qr.quadPointCoor(ig,0),
                                      qr.quadPointCoor(ig,1),
                                      qr.quadPointCoor(ig,2));

            for ( int icoor = 0;icoor < nbCoor;icoor++ )
            {
                //dPhiRef( i, icoor, ig ) = refFE.dPhi( i, icoor, ig, qr );
                dPhiRef( i, icoor, ig ) = refFE.dPhi( i,
                                                      icoor,
                                                      qr.quadPointCoor(ig,0),
                                                      qr.quadPointCoor(ig,1),
                                                      qr.quadPointCoor(ig,2));
            }
        }
        for ( int k = 0;k < nbGeoNode;k++ )
        {
            //phiGeo( k, ig ) = geoMap.phi( k, ig, qr );
            phiGeo( k, ig ) = geoMap.phi( k,
                                          qr.quadPointCoor(ig,0),
                                          qr.quadPointCoor(ig,1),
                                          qr.quadPointCoor(ig,2));
            for ( int icoor = 0;icoor < nbCoor;icoor++ )
            {
                //dPhiGeo( k, icoor, ig ) = geoMap.dPhi( k, icoor, ig, qr );
                dPhiGeo( k, icoor, ig ) = geoMap.dPhi( k, icoor, qr.quadPointCoor(ig,0),
                                                       qr.quadPointCoor(ig,1),
                                                       qr.quadPointCoor(ig,2));
            }
        }
    }
#ifdef TEST_PRE
    _hasQR = true;
#endif

}
/*!
 This constructor is typically used for the currentBdFE (derived class)
 when no quadrature rule is needed (to have just the coordinates of the
 nodes on the current element)
*/
StaticBdFE::StaticBdFE( const RefFE& _refFE, const GeoMap& _geoMap ) :
    nbGeoNode( _geoMap.nbDof() ), nbNode( _refFE.nbDof() ), nbCoor( _refFE.nbCoor() ),
    nbQuadPt( quadRuleDummy.nbQuadPt() ), point( nbGeoNode, nbCoor + 1 ),
        refFE( _refFE ), geoMap( _geoMap ), qr( quadRuleDummy ),
        phi( ( int ) nbNode, ( int ) nbQuadPt ), dPhiRef( ( int ) nbNode, ( int ) nbCoor,
                ( int ) nbQuadPt ),
        dPhi( ( int ) nbNode, ( int ) nbCoor, ( int ) nbQuadPt ),
        phiGeo( ( int ) nbGeoNode, ( int ) nbQuadPt ),
        dPhiGeo( ( int ) nbGeoNode, ( int ) nbCoor, ( int ) nbQuadPt ),
        weightMeas( ( int ) nbQuadPt ), meas( ( int ) nbQuadPt ),
        normal( ( int ) nbCoor + 1, ( int ) nbQuadPt ),
        tangent( ( int ) nbCoor, ( int ) nbCoor + 1, ( int ) nbQuadPt ),
        metric( ( int ) nbCoor, ( int ) nbCoor, ( int ) nbQuadPt ), quadPt( ( int ) nbQuadPt, 3 ),
        invArea( 1. )
{
    CONSTRUCTOR( "StaticBdFE (without quadrature rule)" );
    // be carefull : dPhiRef, dPhiGeo are not initialized
#ifdef TEST_PRE

    _hasQR = false;
#endif
}

/*!
 This constructor is typically used for the static boundary of RefHybridFE
*/
StaticBdFE::StaticBdFE( const RefFE& _refFE, const GeoMap& _geoMap,
                        const QuadRule& _qr, const Real* refcoor,
                        UInt currentid, Real _invarea ) :
        _currentId( currentid ),
        nbGeoNode( _geoMap.nbDof() ), nbNode( _refFE.nbDof() ), nbCoor( _refFE.nbCoor() ),
        nbQuadPt( _qr.nbQuadPt() ), point( nbGeoNode, nbCoor + 1 ),
        refFE( _refFE ), geoMap( _geoMap ), qr( _qr ),
        phi( ( int ) nbNode, ( int ) nbQuadPt ), dPhiRef( ( int ) nbNode, ( int ) nbCoor,
                ( int ) nbQuadPt ),
        dPhi( ( int ) nbNode, ( int ) nbCoor, ( int ) nbQuadPt ),
        phiGeo( ( int ) nbGeoNode, ( int ) nbQuadPt ),
        dPhiGeo( ( int ) nbGeoNode, ( int ) nbCoor, ( int ) nbQuadPt ),
        weightMeas( ( int ) nbQuadPt ), meas( ( int ) nbQuadPt ),
        normal( ( int ) nbCoor + 1, ( int ) nbQuadPt ),
        tangent( ( int ) nbCoor, ( int ) nbCoor + 1, ( int ) nbQuadPt ),
        metric( ( int ) nbCoor, ( int ) nbCoor, ( int ) nbQuadPt ), quadPt( ( int ) nbQuadPt, 3 ),
        invArea( _invarea )
{
    CONSTRUCTOR( "StaticBdFE (with refCoor)" );
    for ( int ig = 0;ig < nbQuadPt;ig++ )
    {
        for ( int i = 0;i < nbNode;i++ )
        {
            //phi( i, ig ) = refFE.phi( i, ig, qr );
            phi( i, ig ) = invArea*refFE.phi( i,               // invArea added here
                                      qr.quadPointCoor(ig,0),
                                      qr.quadPointCoor(ig,1),
                                      qr.quadPointCoor(ig,2));

            for ( int icoor = 0;icoor < nbCoor;icoor++ )
            {
                //dPhiRef( i, icoor, ig ) = refFE.dPhi( i, icoor, ig, qr );
                dPhiRef( i, icoor, ig ) = refFE.dPhi( i,
                                                      icoor,
                                                      qr.quadPointCoor(ig,0),
                                                      qr.quadPointCoor(ig,1),
                                                      qr.quadPointCoor(ig,2));
            }
        }
        for ( int k = 0;k < nbGeoNode;k++ )
        {
            //phiGeo( k, ig ) = geoMap.phi( k, ig, qr );
            phiGeo( k, ig ) = geoMap.phi( k,
                                          qr.quadPointCoor(ig,0),
                                          qr.quadPointCoor(ig,1),
                                          qr.quadPointCoor(ig,2));
            for ( int icoor = 0;icoor < nbCoor;icoor++ )
            {
                //dPhiGeo( k, icoor, ig ) = geoMap.dPhi( k, icoor, ig, qr );
                dPhiGeo( k, icoor, ig ) = geoMap.dPhi( k, icoor, qr.quadPointCoor(ig,0),
                                                       qr.quadPointCoor(ig,1),
                                                       qr.quadPointCoor(ig,2));
            }
        }
    }
   
    for ( int i = 0;i < nbGeoNode;i++ )
    {
        point( i, 0 ) = refcoor[ 3 * i ];
        point( i, 1 ) = refcoor[ 3 * i + 1 ];
        point( i, 2 ) = refcoor[ 3 * i + 2 ];
    }
#ifdef TEST_PRE
    _hasQR = true;
#endif

    _comp_meas_normal();
    _comp_quad_point_coor();
#ifdef TEST_PRE

    _hasMeas = true;
    _hasTangent = true;
    _hasNormal = true;
    _hasQuadPtCoor = true;
    _hasFirstDeriv = false;
#endif
}

StaticBdFE::~StaticBdFE()
{
    DESTRUCTOR( "StaticBdFE" )
}

/*
//----------------------------------------------------------------------
void StaticBdFE::coorMap( Real& x, Real& y, Real& z,
                          const Real & xi, const Real & eta ) const
{
    x = y = z = 0.;
    for ( int i = 0;i < ( int ) nbGeoNode;i++ )
    {
        x += point( i, 0 ) * geoMap.phi( i, xi, eta, 0. );
        y += point( i, 1 ) * geoMap.phi( i, xi, eta, 0. );
#if defined(THREEDIM)

        z += point( i, 2 ) * geoMap.phi( i, xi, eta, 0. );
#endif

    }
}
*/


void StaticBdFE::coorMap( Real& x, Real& y, Real& z,
                          const Real & xi, const Real & eta ) const
{
	Vector coor = ZeroVector(3);
	for( UInt icoor=0; icoor< nDimensions; icoor++)
	        for ( UInt i = 0;i < (UInt)nbGeoNode;i++ )
	        	coor[icoor] += point( i, icoor ) * geoMap.phi( i, xi, eta, 0. );
    x = coor[0]; y=coor[1]; z=coor[2];
}


Real StaticBdFE::measure() const
{
    ASSERT_PRE( _hasMeas, "Needs measure. Call an update function" )

    Real meas = 0.;
    for ( int ig = 0;ig < nbQuadPt;ig++ )
        meas += weightMeas( ig );
    return meas;
}

//----------------------------------------------------------------------
//compute the coordinates of the quadrature points on the current element
//----------------------------------------------------------------------
void StaticBdFE::_comp_quad_point_coor()
{
    ASSERT_PRE( _hasQR, "Needs a quadrature rule" )
    for ( int ig = 0;ig < nbQuadPt;ig++ )
    {
        coorMap( quadPt( ig, 0 ), quadPt( ig, 1 ), quadPt( ig, 2 ), qr.quadPointCoor( ig, 0 ), qr.quadPointCoor( ig, 1 ) );
    }
}

//----------------------------------------------------------------------
//         compute the measure and the local basis on the element
//----------------------------------------------------------------------
void StaticBdFE::_comp_meas()
{
    ASSERT_PRE( _hasQR, "Needs a quadrature rule" )
    Real s, fctDer;

    for ( int ig = 0;ig < nbQuadPt;ig++ )
    {
        for ( int icoor = 0;icoor < (int)nDimensions;icoor++ )
        {
            for ( int k = 0;k < (int)nDimensions-1;k++ )
            {
                fctDer = 0.;
                for ( int i = 0;i < nbGeoNode;i++ )
                {
                    fctDer += point( i, icoor ) * dPhiGeo( i, k, ig );
                }
                // local basis (non-unit tangents) at integration points
                tangent( k, icoor, ig ) = fctDer;
            }
        }
        // metric tensor
        s = 0.;
        for ( int icoor = 0;icoor < (int)nDimensions;icoor++ )
            s += tangent( 0, icoor, ig ) * tangent( 0, icoor, ig );
        metric( 0, 0, ig ) = s;
#if defined(TWODIM)
        meas( ig ) = sqrt( metric( 0, 0, ig ) );
        weightMeas( ig ) = meas( ig ) * qr.weight( ig );
}
#elif defined(THREEDIM)
        s = 0.;
        for ( int icoor = 0;icoor < (int)nDimensions;icoor++ )
            s += tangent( 1, icoor, ig ) * tangent( 1, icoor, ig );
        metric( 1, 1, ig ) = s;
        s = 0.;
        for ( int icoor = 0;icoor < (int)nDimensions;icoor++ )
            s += tangent( 0, icoor, ig ) * tangent( 1, icoor, ig );
        metric( 0, 1, ig ) = metric( 1, 0, ig ) = s;
        meas( ig ) = sqrt( metric( 0, 0, ig ) * metric( 1, 1, ig )
                           - metric( 0, 1, ig ) * metric( 1, 0, ig ) );
        weightMeas( ig ) = meas( ig ) * qr.weight( ig );
        for ( int icoor = 0;icoor < (int)nDimensions;icoor++ )
            for ( int k = 0;k < (int)nDimensions-1;k++ )
                    tangent( k, icoor, ig ) = tangent( k, icoor, ig )
                    / sqrt( metric( k, k, ig ) );
    }
#endif
}
void StaticBdFE::_comp_meas_normal()
/*!
  It computes the measure, the local basis and the unit normal on the element
*/
{
    ASSERT_PRE( _hasQR, "Needs a quadrature rule" )

    Real norm, n1, n2;
    _comp_meas();
#if defined(TWODIM)

    for ( int ig = 0;ig < nbQuadPt;ig++ )
    {
        n1 = - tangent( 0, 1, ig );
        n2 = tangent( 0, 0, ig );
        norm = sqrt( n1 * n1 + n2 * n2 );
        normal( 0, ig ) = n1 / norm;
        normal( 1, ig ) = n2 / norm;
    }
#elif defined(THREEDIM)
    Real n3;
    for ( int ig = 0;ig < nbQuadPt;ig++ )
    {
        n1 = tangent( 0, 1, ig ) * tangent( 1, 2, ig )
             - tangent( 0, 2, ig ) * tangent( 1, 1, ig );
        n2 = tangent( 0, 2, ig ) * tangent( 1, 0, ig )
             - tangent( 0, 0, ig ) * tangent( 1, 2, ig );
        n3 = tangent( 0, 0, ig ) * tangent( 1, 1, ig )
             - tangent( 0, 1, ig ) * tangent( 1, 0, ig );
        norm = sqrt( n1 * n1 + n2 * n2 + n3 * n3 );
        normal( 0, ig ) = n1 / norm;
        normal( 1, ig ) = n2 / norm;
        normal( 2, ig ) = n3 / norm;
    }

#endif
}
}
