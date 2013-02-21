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
    @brief  A class for static boundary finite element

    @author Jean-Frederic Gerbeau
            Vincent Martin
    @date 00-09-2002

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#include <lifev/core/fem/CurrentBoundaryFEBase.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================

CurrentBoundaryFEBase::CurrentBoundaryFEBase ( const ReferenceFE& refFE, const GeometricMap& geoMap,
                                               const QuadratureRule& qr ) :
    M_nbGeoNode ( geoMap.nbDof() ),
    M_nbNode    ( refFE.nbDof() ),
    M_nbCoor    ( refFE.nbCoor() ),
    M_nbQuadPt  ( qr.nbQuadPt() ),
    M_point     ( M_nbGeoNode, nDimensions ),
    refFE     ( refFE ),
    geoMap    ( geoMap ),
    qr ( qr ),
    M_phi       ( ( int ) M_nbNode, ( int ) M_nbQuadPt ),
    M_dPhiRef   ( ( int ) M_nbNode, ( int ) M_nbCoor, ( int ) M_nbQuadPt ),
    M_dPhi      ( ( int ) M_nbNode, ( int ) M_nbCoor, ( int ) M_nbQuadPt ),
    M_phiGeo    ( ( int ) M_nbGeoNode, ( int ) M_nbQuadPt ),
    M_dPhiGeo   ( ( int ) M_nbGeoNode, ( int ) M_nbCoor, ( int ) M_nbQuadPt ),
    M_weightMeas ( ( int ) M_nbQuadPt ),
    M_meas    ( ( int ) M_nbQuadPt ),
    M_normal    ( ( int ) M_nbCoor + 1, ( int ) M_nbQuadPt ),
    M_tangent   ( ( int ) M_nbCoor, ( int ) M_nbCoor + 1, ( int ) M_nbQuadPt ),
    M_metric    ( ( int ) M_nbCoor, ( int ) M_nbCoor, ( int ) M_nbQuadPt ),
    M_quadPt    ( ( int ) M_nbQuadPt, 3 ),
    invArea   ( 1. )
{
    Real x[3] = {0, 0, 0};
    for ( UInt iQuadPt (0); iQuadPt < M_nbQuadPt; iQuadPt++ )
    {
        for ( UInt icoor = 0; icoor < M_nbCoor + 1; icoor++)
        {
            x[icoor] =  qr.quadPointCoor (iQuadPt, icoor);
        }

        for ( UInt iNode (0); iNode < M_nbNode; iNode++ )
        {

            M_phi ( iNode, iQuadPt ) = refFE.phi ( iNode, x[0], x[1], x[2]);

            for ( UInt iCoor (0); iCoor < M_nbCoor; iCoor++ )
            {
                M_dPhiRef ( int (iNode), int (iCoor), int (iQuadPt) ) = refFE.dPhi ( iNode, iCoor, x[0], x[1], x[2]);
            }
        }

        for ( UInt iGeoNode (0); iGeoNode < M_nbGeoNode; iGeoNode++ )
        {
            M_phiGeo ( iGeoNode, iQuadPt ) = geoMap.phi ( iGeoNode, x[0], x[1], x[2]);

            for ( UInt iCoor (0); iCoor < M_nbCoor; iCoor++ )
            {
                M_dPhiGeo ( int (iGeoNode), int (iCoor), int (iQuadPt) ) = geoMap.dPhi ( iGeoNode, iCoor, x[0], x[1], x[2]);
            }
        }
    }
#ifdef TEST_PRE
    M_hasQR = true;
#endif

}

CurrentBoundaryFEBase::CurrentBoundaryFEBase ( const ReferenceFE& refFE, const GeometricMap& geoMap ) :
    M_nbGeoNode ( geoMap.nbDof() ),
    M_nbNode ( refFE.nbDof() ),
    M_nbCoor ( refFE.nbCoor() ),
    M_nbQuadPt ( quadRuleDummy.nbQuadPt() ),
    M_point ( M_nbGeoNode, nDimensions ),
    refFE ( refFE ),
    geoMap ( geoMap ),
    qr ( quadRuleDummy ),
    M_phi ( ( int ) M_nbNode, ( int ) M_nbQuadPt ),
    M_dPhiRef ( ( int ) M_nbNode, ( int ) M_nbCoor, ( int ) M_nbQuadPt ),
    M_dPhi ( ( int ) M_nbNode, ( int ) M_nbCoor, ( int ) M_nbQuadPt ),
    M_phiGeo ( ( int ) M_nbGeoNode, ( int ) M_nbQuadPt ),
    M_dPhiGeo ( ( int ) M_nbGeoNode, ( int ) M_nbCoor, ( int ) M_nbQuadPt ),
    M_weightMeas ( ( int ) M_nbQuadPt ),
    M_meas ( ( int ) M_nbQuadPt ),
    M_normal ( ( int ) M_nbCoor + 1, ( int ) M_nbQuadPt ),
    M_tangent ( ( int ) M_nbCoor, ( int ) M_nbCoor + 1, ( int ) M_nbQuadPt ),
    M_metric ( ( int ) M_nbCoor, ( int ) M_nbCoor, ( int ) M_nbQuadPt ),
    M_quadPt ( ( int ) M_nbQuadPt, 3 ),
    invArea ( 1. )
{
#ifdef TEST_PRE
    M_hasQR = false;
#endif
}

CurrentBoundaryFEBase::CurrentBoundaryFEBase ( const ReferenceFE& refFE, const GeometricMap& geoMap,
                                               const QuadratureRule& qr, const Real* refcoor,
                                               UInt currentid, Real invarea ) :
    M_nbGeoNode ( geoMap.nbDof() ),
    M_nbNode ( refFE.nbDof() ),
    M_nbCoor ( refFE.nbCoor() ),
    M_nbQuadPt ( qr.nbQuadPt() ),
    M_point ( M_nbGeoNode, nDimensions ),
    refFE ( refFE ),
    geoMap ( geoMap ),
    qr ( qr ),
    M_phi ( ( int ) M_nbNode, ( int ) M_nbQuadPt ),
    M_dPhiRef ( ( int ) M_nbNode, ( int ) M_nbCoor, ( int ) M_nbQuadPt ),
    M_dPhi ( ( int ) M_nbNode, ( int ) M_nbCoor, ( int ) M_nbQuadPt ),
    M_phiGeo ( ( int ) M_nbGeoNode, ( int ) M_nbQuadPt ),
    M_dPhiGeo ( ( int ) M_nbGeoNode, ( int ) M_nbCoor, ( int ) M_nbQuadPt ),
    M_weightMeas ( ( int ) M_nbQuadPt ),
    M_meas ( ( int ) M_nbQuadPt ),
    M_normal ( ( int ) M_nbCoor + 1, ( int ) M_nbQuadPt ),
    M_tangent ( ( int ) M_nbCoor, ( int ) M_nbCoor + 1, ( int ) M_nbQuadPt ),
    M_metric ( ( int ) M_nbCoor, ( int ) M_nbCoor, ( int ) M_nbQuadPt ),
    M_quadPt ( ( int ) M_nbQuadPt, 3 ),
    invArea ( invarea ),
    M_currentID ( currentid )
{
    Real x[3] = {0, 0, 0};
    for ( UInt iQuadPt (0); iQuadPt < M_nbQuadPt; iQuadPt++ )
    {
        for ( UInt icoor = 0; icoor < M_nbCoor + 1; icoor++)
        {
            x[icoor] =  qr.quadPointCoor (iQuadPt, icoor);
        }

        for ( UInt iNode (0); iNode < M_nbNode; iNode++ )
        {
            M_phi ( iNode, iQuadPt ) = invArea * refFE.phi ( iNode, x[0], x[1], x[2]);

            for ( UInt iCoor (0); iCoor < M_nbCoor; iCoor++ )
            {
                M_dPhiRef ( int (iNode), int (iCoor), int (iQuadPt) ) = refFE.dPhi ( iNode, iCoor, x[0], x[1], x[2]);
            }

        }
        for ( UInt iGeoNode (0); iGeoNode < M_nbGeoNode; iGeoNode++ )
        {
            M_phiGeo ( iGeoNode, iQuadPt ) = geoMap.phi ( iGeoNode, x[0], x[1], x[2]);

            for ( UInt iCoor (0); iCoor < M_nbCoor; iCoor++ )
            {
                M_dPhiGeo ( int (iGeoNode), int (iCoor), int (iQuadPt) ) = geoMap.dPhi ( iGeoNode, iCoor, x[0], x[1], x[2]);
            }
        }
    }

    for ( UInt iGeoNode (0); iGeoNode < M_nbGeoNode; iGeoNode++ )
    {
        for ( UInt icoor = 0; icoor < nDimensions; icoor++)
        {
            M_point ( iGeoNode, icoor ) = refcoor[ nDimensions * iGeoNode + icoor];
        }
    }

#ifdef TEST_PRE
    M_hasQR = true;
#endif

    computeMeasureNormal();
    computeQuadPointCoordinate();

#ifdef TEST_PRE
    M_hasMeasure = true;
    M_hasTangent = true;
    M_hasNormal = true;
    M_hasQuadPtCoor = true;
    M_hasFirstDerivative = false;
#endif
}

CurrentBoundaryFEBase::~CurrentBoundaryFEBase()
{}

// ===================================================
// Methods
// ===================================================

void CurrentBoundaryFEBase::coorMap ( Real& x, Real& y, Real& z,
                                      const Real& xi, const Real& eta ) const
{
    Vector coor = ZeroVector (3);
    for ( UInt icoor = 0; icoor < M_nbCoor + 1; icoor++)
        for ( UInt i = 0; i < (UInt) M_nbGeoNode; i++ )
        {
            coor[icoor] += M_point ( i, icoor ) * geoMap.phi ( i, xi, eta, 0. );
        }
    x = coor[0];
    y = coor[1];
    z = coor[2];
}


Real CurrentBoundaryFEBase::measure() const
{
    ASSERT_PRE ( M_hasMeasure, "Needs measure. Call an update function" );

    Real meas = 0.;
    for ( UInt iQuadPt (0); iQuadPt < M_nbQuadPt; iQuadPt++ )
    {
        meas += M_weightMeas ( iQuadPt );
    }
    return meas;
}

void CurrentBoundaryFEBase::computeQuadPointCoordinate()
{
    ASSERT_PRE ( M_hasQR, "Needs a quadrature rule" );

    for ( UInt iQuadPt (0); iQuadPt < M_nbQuadPt; iQuadPt++ )
    {
        coorMap ( M_quadPt ( iQuadPt, UInt (0) ), M_quadPt ( iQuadPt, UInt (1) ), M_quadPt ( iQuadPt, UInt (2) ),
                  qr.quadPointCoor ( iQuadPt, 0 ), qr.quadPointCoor ( iQuadPt, 1 ) );
    }
}

void CurrentBoundaryFEBase::computeMeasure()
{
    ASSERT_PRE ( M_hasQR, "Needs a quadrature rule" )
    Real s, fctDer;

    for ( UInt iQuadPt (0); iQuadPt < M_nbQuadPt; iQuadPt++ )
    {
        for ( UInt iCoor (0); iCoor < M_nbCoor + 1; iCoor++ )
        {
            for ( UInt jCoor (0); jCoor < M_nbCoor; jCoor++ )
            {
                fctDer = 0.;
                for ( UInt iGeoNode (0); iGeoNode < M_nbGeoNode; iGeoNode++ )
                {
                    fctDer += M_point ( iGeoNode, iCoor ) * M_dPhiGeo ( int (iGeoNode), int (jCoor), int (iQuadPt) );
                }

                M_tangent ( int (jCoor), int (iCoor), int (iQuadPt) ) = fctDer;
            }
        }
        // metric tensor
        s = 0.;
        for ( UInt iCoor = 0; iCoor < M_nbCoor + 1; iCoor++ )
        {
            s += M_tangent ( int (0), int (iCoor), int (iQuadPt) ) * M_tangent ( int (0), int (iCoor), int (iQuadPt) );
        }
        M_metric ( int (0), int (0), int (iQuadPt) ) = s;

        if (M_nbCoor == 1)
        {
            M_meas ( iQuadPt ) = std::sqrt ( M_metric ( 0, 0, int (iQuadPt) ) );
            M_weightMeas ( iQuadPt ) = M_meas ( iQuadPt ) * qr.weight ( iQuadPt );
        }
        else
        {
            s = 0.;
            for ( UInt iCoor = 0; iCoor < M_nbCoor + 1; iCoor++ )
            {
                s += M_tangent ( int (1), int (iCoor), int (iQuadPt) ) * M_tangent ( int (1), int (iCoor), int (iQuadPt) );
            }

            M_metric ( int (1), int (1), int (iQuadPt) ) = s;
            s = 0.;
            for ( UInt iCoor = 0; iCoor < M_nbCoor + 1; iCoor++ )
            {
                s += M_tangent ( int (0), int (iCoor), int (iQuadPt) ) * M_tangent ( int (1), int (iCoor), int (iQuadPt) );
            }

            M_metric ( int (0), int (1), int (iQuadPt) ) = M_metric ( int (1), int (0), int (iQuadPt) ) = s;
            M_meas ( iQuadPt ) = std::sqrt ( M_metric ( int (0), int (0), int (iQuadPt) ) * M_metric ( int (1), int (1), int (iQuadPt) )
                                             - M_metric ( int (0), int (1), int (iQuadPt) ) * M_metric ( int (1), int (0), int (iQuadPt) ) );
            M_weightMeas ( iQuadPt ) = M_meas ( iQuadPt ) * qr.weight ( iQuadPt );

            for ( UInt iCoor = 0; iCoor < M_nbCoor + 1; iCoor++ )
            {
                for ( UInt jCoor = 0; jCoor < M_nbCoor; jCoor++ )
                {
                    M_tangent ( int (jCoor), int (iCoor), int (iQuadPt) ) = M_tangent ( int (jCoor), int (iCoor), int (iQuadPt) )
                                                                            / std::sqrt ( M_metric ( int (jCoor), int (jCoor), int (iQuadPt) ) );
                }
            }
        }
    }
}


void CurrentBoundaryFEBase::computeMeasureNormal()
{
    ASSERT_PRE ( M_hasQR, "Needs a quadrature rule" )

    Real norm, n1, n2;
    computeMeasure();
    if ( M_nbCoor == 1 )
    {
        for ( UInt ig = 0; ig < M_nbQuadPt; ig++ )
        {
            n1 = M_tangent ( int (0), int (1), int (ig) );
            n2 = - M_tangent ( int (0), int (0), int (ig) );
            norm = std::sqrt ( n1 * n1 + n2 * n2 );
            M_normal ( UInt (0), ig ) = n1 / norm;
            M_normal ( UInt (1), ig ) = n2 / norm;
        }
    }
    else if ( M_nbCoor == 2 )
    {
        Real n3;
        for ( UInt ig = 0; ig < M_nbQuadPt; ig++ )
        {
            n1 = M_tangent ( int (0), int (1), int (ig) ) * M_tangent ( int (1), int (2), int (ig) )
                 - M_tangent ( int (0), int (2), int (ig) ) * M_tangent ( int (1), int (1), int (ig) );
            n2 = M_tangent ( int (0), int (2), int (ig) ) * M_tangent ( int (1), int (0), int (ig) )
                 - M_tangent ( int (0), int (0), int (ig) ) * M_tangent ( int (1), int (2), int (ig) );
            n3 = M_tangent ( int (0), int (0), int (ig) ) * M_tangent ( int (1), int (1), int (ig) )
                 - M_tangent ( int (0), int (1), int (ig) ) * M_tangent ( int (1), int (0), int (ig) );
            norm = std::sqrt ( n1 * n1 + n2 * n2 + n3 * n3 );
            M_normal ( UInt (0), ig ) = n1 / norm;
            M_normal ( UInt (1), ig ) = n2 / norm;
            M_normal ( UInt (2), ig ) = n3 / norm;
        }
    }
    else
    {
        ASSERT ( 0, "Wrong space dimension." );
    }
}

}
