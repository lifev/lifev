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
    @brief File containing the CurrentFE class implementation

    @author Jean-Frederic Gerbeau
    @date 00-04-2002

    @contributor V. Martin
                 Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#include <lifev/core/fem/CurrentFE.hpp>

namespace LifeV
{
CurrentFE::CurrentFE ( const ReferenceFE& refFE, const GeometricMap& geoMap, const QuadratureRule& qr )
    :
    M_nbNode ( refFE.nbDof() ),
    M_nbLocalCoor ( refFE.nbLocalCoor() ),
    M_nbDiag ( refFE.nbDiag() ),
    M_nbUpper ( refFE.nbUpper() ) ,
    M_nbPattern ( refFE.nbPattern() ),

    M_nbGeoNode ( geoMap.nbDof() ),
    M_nbQuadPt ( qr.nbQuadPt() ),

    M_refFE ( &refFE ),
    M_geoMap ( &geoMap),
    M_quadRule (new QuadratureRule (qr) ),


    M_cellNodes (boost::extents[geoMap.nbDof()][nDimensions]),
    M_quadNodes (boost::extents[M_nbQuadPt][nDimensions]),

    M_dphiGeometricMap (boost::extents[M_nbGeoNode][M_nbLocalCoor][M_nbQuadPt]),
    M_jacobian (boost::extents[M_nbLocalCoor][M_nbLocalCoor][M_nbQuadPt]),
    M_detJacobian (boost::extents[M_nbQuadPt]),
    M_wDetJacobian (boost::extents[M_nbQuadPt]),
    M_tInverseJacobian (boost::extents[M_nbLocalCoor][M_nbLocalCoor][M_nbQuadPt]),

    M_phi (boost::extents[M_nbNode][M_refFE->feDim()][M_nbQuadPt]),
    M_dphi (boost::extents[M_nbNode][M_nbLocalCoor][M_nbQuadPt]),
    M_d2phi (boost::extents[M_nbNode][M_nbLocalCoor][M_nbLocalCoor][M_nbQuadPt]),
    M_phiVect (boost::extents[M_nbNode][M_refFE->feDim()][M_nbQuadPt]),

    M_dphiRef (boost::extents[M_nbNode][M_nbLocalCoor][M_nbQuadPt]),
    M_d2phiRef (boost::extents[M_nbNode][M_nbLocalCoor][M_nbLocalCoor][M_nbQuadPt]),
    M_divPhiRef (boost::extents[M_nbNode][M_nbQuadPt]),

    M_cellNodesUpdated (false),
    M_quadNodesUpdated (false),

    M_dphiGeometricMapUpdated (false),
    M_jacobianUpdated (false),
    M_detJacobianUpdated (false),
    M_wDetJacobianUpdated (false),
    M_tInverseJacobianUpdated (false),

    M_phiUpdated (false),
    M_dphiUpdated (false),
    M_d2phiUpdated (false),
    M_phiVectUpdated (false),

    M_dphiRefUpdated (false),
    M_d2phiRefUpdated (false),
    M_divPhiRefUpdated (false)

{
    for ( UInt iterQuad (0); iterQuad < M_nbQuadPt; ++iterQuad )
    {
        for ( UInt iterNode (0); iterNode < M_nbNode; ++iterNode )
        {
            // --- PHI ---
            if (M_refFE->hasPhi() )
            {
                for (UInt iterFEDim (0); iterFEDim < M_refFE->feDim(); ++iterFEDim)
                {
                    M_phi[iterNode][iterFEDim][iterQuad] = M_refFE->phi ( iterNode, iterFEDim,
                                                                          M_quadRule->quadPointCoor (iterQuad) );
                }
            }

            // --- DIV ---
            if (M_refFE->hasDivPhi() )
            {
                M_divPhiRef[iterNode][iterQuad] = M_refFE->divPhi (iterNode,
                                                                   M_quadRule->quadPointCoor (iterQuad) );
            };


            for ( UInt icoor (0); icoor < M_nbLocalCoor; ++icoor )
            {
                // --- DPHI ---
                if (M_refFE->hasDPhi() )
                {
                    M_dphiRef[iterNode][icoor][iterQuad] = M_refFE->dPhi ( iterNode,
                                                                           icoor,
                                                                           M_quadRule->quadPointCoor (iterQuad) );
                };

                // --- D2PHI
                if (M_refFE->hasD2Phi() )
                {
                    for ( UInt jcoor (0); jcoor < M_nbLocalCoor; ++jcoor )
                    {
                        M_d2phiRef[iterNode][icoor][jcoor][iterQuad] = M_refFE->d2Phi ( iterNode,
                                                                                        icoor, jcoor,
                                                                                        M_quadRule->quadPointCoor (iterQuad) );
                    }
                }
            }
        }
        for ( UInt k (0); k < M_nbGeoNode; ++k )
        {
            for ( UInt icoor (0); icoor < M_nbLocalCoor; ++icoor )
            {
                M_dphiGeometricMap[k][icoor][iterQuad] = M_geoMap->dPhi ( k, icoor,
                                                                          M_quadRule->quadPointCoor (iterQuad) );
            }
        }
    }

    M_phiUpdated = true;
    M_dphiRefUpdated = true;
    M_divPhiRefUpdated = true;
    M_d2phiRefUpdated = true;
    M_dphiGeometricMapUpdated = true;
}

CurrentFE::CurrentFE ( const ReferenceFE& refFE, const GeometricMap& geoMap ) :
    M_nbNode                  (refFE.nbDof() ),
    M_nbLocalCoor             (refFE.nbLocalCoor() ),
    M_nbDiag                  (refFE.nbDiag() ),
    M_nbUpper                 (refFE.nbUpper() ),
    M_nbPattern               (refFE.nbPattern() ),

    M_nbGeoNode               (geoMap.nbDof() ),
    M_nbQuadPt                (0),

    M_refFE                   (&refFE),
    M_geoMap                  (&geoMap),
    M_quadRule                (0),

    M_cellNodes               (boost::extents[geoMap.nbDof()][nDimensions]),


    M_cellNodesUpdated        (false),
    M_quadNodesUpdated        (false),
    M_dphiGeometricMapUpdated (false),
    M_jacobianUpdated         (false),
    M_detJacobianUpdated      (false),
    M_wDetJacobianUpdated     (false),
    M_tInverseJacobianUpdated (false),
    M_phiUpdated              (false),
    M_dphiUpdated             (false),
    M_d2phiUpdated            (false),
    M_phiVectUpdated          (false),
    M_dphiRefUpdated          (false),
    M_d2phiRefUpdated         (false),
    M_divPhiRefUpdated        (false)

{
    // Nothing to be done here
}

CurrentFE::CurrentFE (const CurrentFE& fe) :
    M_nbNode                  (fe.M_nbNode),
    M_nbLocalCoor             (fe.M_nbLocalCoor),
    M_nbDiag                  (fe.M_nbDiag),
    M_nbUpper                 (fe.M_nbUpper) ,
    M_nbPattern               (fe.M_nbPattern),
    M_currentId               (fe.M_currentId),
    M_currentLocalId          (fe.M_currentLocalId),
    M_nbGeoNode               (fe.M_nbGeoNode),
    M_nbQuadPt                (fe.M_nbQuadPt),

    M_refFE                   (fe.M_refFE),
    M_geoMap                  (fe.M_geoMap),
    M_quadRule                (new QuadratureRule (*fe.M_quadRule) ),

    M_cellNodes               (fe.M_cellNodes),
    M_quadNodes               (fe.M_quadNodes),
    M_dphiGeometricMap        (fe.M_dphiGeometricMap),
    M_jacobian                (fe.M_jacobian),
    M_detJacobian             (fe.M_detJacobian),
    M_wDetJacobian            (fe.M_wDetJacobian),
    M_tInverseJacobian        (fe.M_tInverseJacobian),
    M_phi                     (fe.M_phi),
    M_dphi                    (fe.M_dphi),
    M_d2phi                   (fe.M_d2phi),
    M_phiVect                 (fe.M_phiVect),
    M_dphiRef                 (fe.M_dphiRef),
    M_d2phiRef                (fe.M_d2phiRef),
    M_divPhiRef               (fe.M_divPhiRef),

    M_cellNodesUpdated        (fe.M_cellNodesUpdated),
    M_quadNodesUpdated        (fe.M_quadNodesUpdated),
    M_dphiGeometricMapUpdated (fe.M_dphiGeometricMapUpdated),
    M_jacobianUpdated         (fe.M_jacobianUpdated),
    M_detJacobianUpdated      (fe.M_detJacobianUpdated),
    M_wDetJacobianUpdated     (fe.M_wDetJacobianUpdated),
    M_tInverseJacobianUpdated (fe.M_tInverseJacobianUpdated),
    M_phiUpdated              (fe.M_phiUpdated),
    M_dphiUpdated             (fe.M_dphiUpdated),
    M_d2phiUpdated            (fe.M_d2phiUpdated),
    M_phiVectUpdated          (fe.M_phiVectUpdated),
    M_dphiRefUpdated          (fe.M_dphiRefUpdated),
    M_d2phiRefUpdated         (fe.M_d2phiRefUpdated),
    M_divPhiRefUpdated        (fe.M_divPhiRefUpdated)
{
    // Nothing to be done here
}

//----------------------------------------------------------------------

// Added by S. Quinodoz
void CurrentFE::coorBackMap (Real x, Real y, Real z,
                             Real& xi, Real& eta, Real& zeta) const
{
    // If the Mapping is not P1, this simple "back mapping" does not
    // work properly, as it suppose that the mapping is P1. If another
    // mapping (higher order, or for different shape like Q1) is needed
    // a Newton algorithm has to be implemented to solve the problem
    // M(x)-p = 0
    // with M the mapping, p the given point and x its coordinates in
    // the reference frame. However, this is more expensive and inexact.

    // ASSERT(geoMap.name.compare("P1")==0," The 'coorBackMap' function should NOT be used in this case! ");

    // In the P1 mapping case, we want to solve the problem
    // P = \lamda_1 V_1 + \lambda_2 V_2 + \lambda_3 V_3
    // where V_i is the ith vertice of the cell because we have then
    // \hat{P} = \lambda_1 \hat{V}_1 + \lambda_2 \hat{V}_2 + \lambda_3 \hat{V}_3
    // where \hat{P} is our point in the reference frame
    //
    // For "simplicity", we use Cramer's formula (no need for an heavy solver)


    if (M_nbLocalCoor == 1)
    {

        // 1D Case

        const Real a11 = M_cellNodes[1][0] - M_cellNodes[0][0];
        const Real b1 = x - M_cellNodes[0][0];

        Real lambda = b1 / a11;

        xi  = lambda * 1; // 1= refFE.xi(1)-refFE.xi(0);
        eta = 0;
        zeta = 0;

    }
    else if ( M_nbLocalCoor == 2)
    {

        // 2D Case

        const Real a11 = M_cellNodes[1][0] - M_cellNodes[0][0];
        const Real a12 = M_cellNodes[2][0] - M_cellNodes[0][0];
        const Real a21 = M_cellNodes[1][1] - M_cellNodes[0][1];
        const Real a22 = M_cellNodes[2][1] - M_cellNodes[0][1];
        const Real b1 = x - M_cellNodes[0][0];
        const Real b2 = y - M_cellNodes[0][1];

        Real D = a11 * a22 - a21 * a12;
        Real lambda_1 = b1 * a22 - b2 * a12;
        Real lambda_2 = a11 * b2 - a21 * b1;

        lambda_1 /= D;
        lambda_2 /= D;

        xi  = lambda_1 * 1 + lambda_2 * 0; // 1=refFE.xi(1)  -refFE.xi(0)  ; 0=refFE.xi(2)  -refFE.xi(0);
        eta = lambda_1 * 0 + lambda_2 * 1; // 0=refFE.eta(1) -refFE.eta(0) ; 1=refFE.eta(2) -refFE.eta(0);
        zeta = 0;

    }
    else if (M_nbLocalCoor == 3)
    {

        // 3D Case

        const Real a11 = M_cellNodes[1][0] - M_cellNodes[0][0];
        const Real a12 = M_cellNodes[2][0] - M_cellNodes[0][0];
        const Real a13 = M_cellNodes[3][0] - M_cellNodes[0][0];
        const Real a21 = M_cellNodes[1][1] - M_cellNodes[0][1];
        const Real a22 = M_cellNodes[2][1] - M_cellNodes[0][1];
        const Real a23 = M_cellNodes[3][1] - M_cellNodes[0][1];
        const Real a31 = M_cellNodes[1][2] - M_cellNodes[0][2];
        const Real a32 = M_cellNodes[2][2] - M_cellNodes[0][2];
        const Real a33 = M_cellNodes[3][2] - M_cellNodes[0][2];
        const Real b1 = x - M_cellNodes[0][0];
        const Real b2 = y - M_cellNodes[0][1];
        const Real b3 = z - M_cellNodes[0][2];

        const Real D = a11 * a22 * a33 + a31 * a12 * a23 + a21 * a32 * a13 - a11 * a32 * a23 - a31 * a22 * a13 - a21 * a12 * a33;
        Real lambda_1 = b1 * a22 * a33 + b3 * a12 * a23 + b2 * a32 * a13 - b1 * a32 * a23 - b3 * a22 * a13 - b2 * a12 * a33;
        Real lambda_2 = a11 * b2 * a33 + a31 * b1 * a23 + a21 * b3 * a13 - a11 * b3 * a23 - a31 * b2 * a13 - a21 * b1 * a33;
        Real lambda_3 = a11 * a22 * b3 + a31 * a12 * b2 + a21 * a32 * b1 - a11 * a32 * b2 - a31 * a22 * b1 - a21 * a12 * b3;

        lambda_1 /= D;
        lambda_2 /= D;
        lambda_3 /= D;

        xi  = lambda_1 * 1 + lambda_2 * 0 + lambda_3 * 0;
        eta = lambda_1 * 0 + lambda_2 * 1 + lambda_3 * 0;
        zeta = lambda_1 * 0 + lambda_2 * 0 + lambda_3 * 1;

    }
    else
    {
        ERROR_MSG ("Impossible dimension to invert coordinates");
    }

}

// Compute the Jacobian at the given point
// this means that jacobian(P1,P2,P3,i,j)
// computes d x_i / d zeta_j
// where x are the global coordinates
//       zeta the reference coordinates

Real CurrentFE::pointJacobian (Real hat_x, Real hat_y, Real hat_z,
                               Int comp_x, Int comp_zeta) const
{
    Real jac (0);
    GeoVector P (3);
    P[0] = hat_x;
    P[1] = hat_y;
    P[2] = hat_z;

    for ( UInt i = 0; i < M_nbGeoNode; i++ )
    {
        jac += M_cellNodes[i][comp_x] * M_geoMap->dPhi ( i , comp_zeta , P );
    };

    return jac;
}

// Compute the determinant of the Jacobian at the given point
Real CurrentFE::pointDetJacobian (Real hat_x, Real hat_y, Real hat_z) const
{
    Real det (0);

    switch (M_nbLocalCoor)
    {
        case 1:
        {
            Real a ( pointJacobian (hat_x, hat_y, hat_z, 0, 0) );

            det = a;

            break;
        }
        case 2:
        {
            Real a ( pointJacobian (hat_x, hat_y, hat_z, 0, 0) );
            Real b ( pointJacobian (hat_x, hat_y, hat_z, 0, 1) );
            Real c ( pointJacobian (hat_x, hat_y, hat_z, 1, 0) );
            Real d ( pointJacobian (hat_x, hat_y, hat_z, 1, 1) );

            det = a * d - b * c;

            break;
        }
        case 3:
        {
            Real a ( pointJacobian (hat_x, hat_y, hat_z, 0, 0) );
            Real b ( pointJacobian (hat_x, hat_y, hat_z, 0, 1) );
            Real c ( pointJacobian (hat_x, hat_y, hat_z, 0, 2) );
            Real d ( pointJacobian (hat_x, hat_y, hat_z, 1, 0) );
            Real e ( pointJacobian (hat_x, hat_y, hat_z, 1, 1) );
            Real f ( pointJacobian (hat_x, hat_y, hat_z, 1, 2) );
            Real g ( pointJacobian (hat_x, hat_y, hat_z, 2, 0) );
            Real h ( pointJacobian (hat_x, hat_y, hat_z, 2, 1) );
            Real i ( pointJacobian (hat_x, hat_y, hat_z, 2, 2) );

            Real ei (e * i);
            Real fh (f * h);
            Real bi (b * i);
            Real ch (c * h);
            Real bf (b * f);
            Real ce (c  * e);

            det = a * (ei - fh) + d * (ch - bi) + g * ( bf - ce);

            break;
        }
        default:
            ERROR_MSG ( "Dimension (M_nbLocalCoor): only 1, 2 or 3!" );
            break;
    }

    return det;
}

Real CurrentFE::pointInverseJacobian (Real hat_x, Real hat_y, Real hat_z,
                                      Int comp_x, Int comp_zeta) const
{
    if ( M_nbLocalCoor == 1 )
    {

        return 1 / pointJacobian (hat_x, hat_y, hat_z, comp_x, comp_zeta);
    }
    else if (M_nbLocalCoor == 2)
    {

        Real a11 = pointJacobian (hat_x, hat_y, hat_z, 0, 0);
        Real a12 = pointJacobian (hat_x, hat_y, hat_z, 0, 1);
        Real a21 = pointJacobian (hat_x, hat_y, hat_z, 1, 0);
        Real a22 = pointJacobian (hat_x, hat_y, hat_z, 1, 1);

        Int total (comp_x + comp_zeta);
        Int mysign (1);
        if (total % 2 == 1)
        {
            mysign = -1;
        }

        return mysign * pointJacobian (hat_x, hat_y, hat_z, comp_zeta, comp_x) / (a11 * a22 - a12 * a21);

    }
    else if (M_nbLocalCoor == 3)
    {

        Real a11 = pointJacobian (hat_x, hat_y, hat_z, 0, 0);
        Real a12 = pointJacobian (hat_x, hat_y, hat_z, 0, 1);
        Real a13 = pointJacobian (hat_x, hat_y, hat_z, 0, 2);
        Real a21 = pointJacobian (hat_x, hat_y, hat_z, 1, 0);
        Real a22 = pointJacobian (hat_x, hat_y, hat_z, 1, 1);
        Real a23 = pointJacobian (hat_x, hat_y, hat_z, 1, 2);
        Real a31 = pointJacobian (hat_x, hat_y, hat_z, 2, 0);
        Real a32 = pointJacobian (hat_x, hat_y, hat_z, 2, 1);
        Real a33 = pointJacobian (hat_x, hat_y, hat_z, 2, 2);

        Real det = a11 * a22 * a33 + a12 * a23 * a31 + a13 * a21 * a32 - a11 * a23 * a32 - a13 * a22 * a31 - a12 * a21 * a33;

        std::vector<std::vector< Real > > cof (3, std::vector<Real> (3, 0) );

        cof[0][0] =  (a22 * a33 - a23 * a32);
        cof[0][1] = - (a21 * a33 - a31 * a23);
        cof[0][2] =  (a21 * a32 - a31 * a22);
        cof[1][0] = - (a12 * a33 - a32 * a13);
        cof[1][1] =  (a11 * a33 - a13 * a31);
        cof[1][2] = - (a11 * a32 - a31 * a12);
        cof[2][0] =  (a12 * a23 - a13 * a22);
        cof[2][1] = - (a11 * a23 - a13 * a21);
        cof[2][2] =  (a11 * a22 - a12 * a21);

        // inverse need the TRANSPOSE!
        return cof[comp_x][comp_zeta] / det;
    }
    else
    {
        ERROR_MSG ( "Dimension (M_nbLocalCoor): only 1, 2 or 3!" );
    };
    return 0.;
}

//----------------------------------------------------------------------

void CurrentFE::update (const std::vector<std::vector<Real> >& pts, flag_Type upFlag)
{
    ASSERT (M_nbQuadPt != 0 || upFlag == UPDATE_ONLY_CELL_NODES, " No quadrature rule defined, cannot update!");

    M_cellNodesUpdated = false;
    if ( (upFlag & UPDATE_ONLY_CELL_NODES) != 0)
    {
        computeCellNodes (pts);
    };

    M_quadNodesUpdated = false;
    if ( (upFlag & UPDATE_ONLY_QUAD_NODES) != 0)
    {
        computeQuadNodes();
    };

    M_jacobianUpdated = false;
    if ( (upFlag & UPDATE_ONLY_JACOBIAN) != 0)
    {
        computeJacobian();
    };

    M_tInverseJacobianUpdated = false;
    if ( (upFlag & UPDATE_ONLY_T_INVERSE_JACOBIAN) != 0)
    {
        computeTInverseJacobian();
    };

    M_detJacobianUpdated = false;
    if ( (upFlag & UPDATE_ONLY_DET_JACOBIAN) != 0)
    {
        computeDetJacobian();
    };

    M_wDetJacobianUpdated = false;
    if ( (upFlag & UPDATE_ONLY_W_DET_JACOBIAN) != 0)
    {
        computeWDetJacobian();
    };

    M_dphiUpdated = false;
    if ( (upFlag & UPDATE_ONLY_DPHI) != 0)
    {
        computeDphi();
    };

    M_d2phiUpdated = false;
    if ( (upFlag & UPDATE_ONLY_D2PHI) != 0)
    {
        computeD2phi();
    };

    M_phiVectUpdated = false;
    if ( (upFlag & UPDATE_ONLY_PHI_VECT) != 0)
    {
        computePhiVect();
    };
}

void CurrentFE::update (const std::vector<GeoVector>& pts, flag_Type upFlag)
{
    std::vector< std::vector <Real > > newpts (pts.size(), std::vector<Real> (pts[0].size() ) );

    for (UInt iPt (0); iPt < pts.size(); ++iPt)
    {
        for (UInt iCoord (0); iCoord < pts[iPt].size(); ++iCoord)
        {
            newpts[iPt][iCoord] = pts[iPt][iCoord];
        }
    }

    update (newpts, upFlag);
}

Real CurrentFE::measure() const
{
    ASSERT (M_wDetJacobianUpdated, "Weighted jacobian determinant is not updated!");

    Real meas (0.0);
    for ( UInt iterQuadNode (0); iterQuadNode < M_nbQuadPt; ++iterQuadNode )
    {
        meas += M_wDetJacobian[iterQuadNode];
    }
    return meas;
}

void CurrentFE::barycenter ( Real& x, Real& y, Real& z ) const
{
    ASSERT (M_cellNodesUpdated, "Cell nodes are not updated!");

    x = 0.0;
    y = 0.0;
    z = 0.0;

    for (UInt iNode (0); iNode < M_nbGeoNode; ++iNode)
    {
        x += M_cellNodes[iNode][0];
        y += M_cellNodes[iNode][1];
        z += M_cellNodes[iNode][2];
    }

    x /= M_nbGeoNode;
    y /= M_nbGeoNode;
    z /= M_nbGeoNode;
}

Real CurrentFE::diameter() const
{
    ASSERT (M_cellNodesUpdated, "Cell nodes are not updated!");

    Real s (0.0);
    Real h (0.0);
    for ( UInt i (0); i < M_nbGeoNode - 1; ++i)
    {
        for ( UInt j (i + 1); j < M_nbGeoNode; ++j )
        {
            s = 0.0;
            for ( Int icoor (0); icoor < (Int) nDimensions; ++icoor )
            {
                s += std::fabs ( M_cellNodes[i][icoor] - M_cellNodes[j][icoor] );
            }
            if ( s > h )
            {
                h = s;
            }
        }
    }

    return h;
}

Real CurrentFE::diameter2() const
{
    ASSERT (M_cellNodesUpdated, "Cell nodes are not updated!");

    Real s (0.0);
    Real d (0.0);
    Real h (0.0);
    for ( UInt i (0); i < M_nbGeoNode - 1; ++i)
    {
        for ( UInt j (i + 1); j < M_nbGeoNode; ++j)
        {
            s = 0.;
            for ( Int icoor (0); icoor < (Int) nDimensions; ++icoor)
            {
                d = ( M_cellNodes[i][icoor] - M_cellNodes[j][icoor] );
                d = d * d;
                s += d;
            }
            s = std::sqrt (s);
            if ( s > h )
            {
                h = s;
            }
        }
    }
    return h;
}

void CurrentFE::coorMap ( Real& x, Real& y, Real& z,
                          Real xi, Real eta, Real zeta ) const
{
    ASSERT (M_cellNodesUpdated, "Cell nodes are not updated!");

    x = 0.0;
    y = 0.0;
    z = 0.0;

    // Note: we assume (already in the constructor) that the 2nd dimension of
    // M_cellNodes is nDimensions, that is 3. If in the future nDimensions!=3
    // this piece of code needs to be modified.

    Real geoMap;
    for (UInt iNode (0); iNode < M_nbGeoNode; ++iNode)
    {
        geoMap = M_geoMap->phi ( iNode, xi, eta, zeta );
        x += M_cellNodes[iNode][0] * geoMap;
        y += M_cellNodes[iNode][1] * geoMap;
        z += M_cellNodes[iNode][2] * geoMap;
    }
}

GeoVector CurrentFE::coorMap (const GeoVector& P) const
{
    ASSERT (M_cellNodesUpdated, "Cell nodes are not updated!");

    GeoVector Pcurrent ( nDimensions, 0. );

    for ( UInt i (0); i < M_nbGeoNode; ++i )
    {
        for (UInt j (0); j < nDimensions; ++j)
        {
            Pcurrent[j] += M_cellNodes[i][j] * M_geoMap->phi (i, P);
        }
    }

    return Pcurrent;
}


void CurrentFE::quadRuleVTKexport ( const std::string& filename) const
{
    ASSERT (M_quadNodesUpdated, "Quad nodes are not updated! No export possible");

    std::ofstream output (filename.c_str() );
    ASSERT (!output.fail(), " Unable to open the file for the export of the quadrature ");

    // Header
    output << "# vtk DataFile Version 3.0" << std::endl;
    output << "LifeV : Quadrature export" << std::endl;
    output << "ASCII" << std::endl;
    output << "DATASET POLYDATA" << std::endl;
    output << "POINTS " << M_nbQuadPt << " float" << std::endl;

    for (UInt i (0); i < M_nbQuadPt; ++i)
    {
        output << M_quadNodes[i][0] << " " <<  M_quadNodes[i][1] << " " <<  M_quadNodes[i][2] << std::endl;
    };

    output << "VERTICES " << M_nbQuadPt << " " << 2 * M_nbQuadPt << std::endl;

    for (UInt i (0); i < M_nbQuadPt; ++i)
    {
        output << 1 << " " << i << std::endl;
    };

    output.close();
}


void CurrentFE::setQuadRule (const QuadratureRule& newQuadRule)
{
    // Set the quadrature
    if (M_quadRule != 0)
    {
        delete M_quadRule;
    };
    M_quadRule = new QuadratureRule (newQuadRule);

    // Adjust the constants related to the quadrature
    M_nbQuadPt =  UInt ( newQuadRule.nbQuadPt() );

    // Resize all the arrays that need it
    M_quadNodes.resize (boost::extents[M_nbQuadPt][3]);

    M_dphiGeometricMap.resize (boost::extents[M_nbGeoNode][M_nbLocalCoor][M_nbQuadPt]);
    M_jacobian.resize (boost::extents[M_nbLocalCoor][M_nbLocalCoor][M_nbQuadPt]);
    M_detJacobian.resize (boost::extents[M_nbQuadPt]);
    M_wDetJacobian.resize (boost::extents[M_nbQuadPt]);
    M_tInverseJacobian.resize (boost::extents[M_nbLocalCoor][M_nbLocalCoor][M_nbQuadPt]);

    M_phi.resize (boost::extents[M_nbNode][M_refFE->feDim()][M_nbQuadPt]);
    M_dphi.resize (boost::extents[M_nbNode][M_nbLocalCoor][M_nbQuadPt]);
    M_d2phi.resize (boost::extents[M_nbNode][M_nbLocalCoor][M_nbLocalCoor][M_nbQuadPt]);
    M_phiVect.resize (boost::extents[M_nbNode][M_refFE->feDim()][M_nbQuadPt]);

    M_dphiRef.resize (boost::extents[M_nbNode][M_nbLocalCoor][M_nbQuadPt]);
    M_d2phiRef.resize (boost::extents[M_nbNode][M_nbLocalCoor][M_nbLocalCoor][M_nbQuadPt]);
    M_divPhiRef.resize (boost::extents[M_nbNode][M_nbQuadPt]);

    // All the updates have to be done again
    M_cellNodesUpdated = false;
    M_quadNodesUpdated = false;

    M_dphiGeometricMapUpdated = false;
    M_jacobianUpdated = false;
    M_detJacobianUpdated = false;
    M_wDetJacobianUpdated = false;
    M_tInverseJacobianUpdated = false;

    M_phiUpdated = false;
    M_dphiUpdated = false;
    M_d2phiUpdated = false;
    M_phiVectUpdated = false;

    M_dphiRefUpdated = false;
    M_d2phiRefUpdated = false;
    M_divPhiRefUpdated = false;

    // Update what can be updated
    for ( UInt iterQuad (0); iterQuad < M_nbQuadPt; ++iterQuad )
    {
        for ( UInt iterNode (0); iterNode < M_nbNode; ++iterNode )
        {
            // --- PHI ---
            if (M_refFE->hasPhi() )
            {
                for (UInt iterFEDim (0); iterFEDim < M_refFE->feDim(); ++iterFEDim)
                {
                    M_phi[iterNode][iterFEDim][iterQuad] = M_refFE->phi ( iterNode, iterFEDim,
                                                                          M_quadRule->quadPointCoor (iterQuad) );
                }
            }

            // --- DIV ---
            if (M_refFE->hasDivPhi() )
            {
                M_divPhiRef[iterNode][iterQuad] = M_refFE->divPhi (iterNode,
                                                                   M_quadRule->quadPointCoor (iterQuad) );
            }

            for ( UInt icoor (0); icoor < M_nbLocalCoor; ++icoor )
            {
                // --- DPHI ---
                if (M_refFE->hasDPhi() )
                {
                    M_dphiRef[iterNode][icoor][iterQuad] = M_refFE->dPhi ( iterNode,
                                                                           icoor,
                                                                           M_quadRule->quadPointCoor (iterQuad) );
                }

                // --- D2PHI ---
                if (M_refFE->hasD2Phi() )
                {
                    for ( UInt jcoor (0); jcoor < M_nbLocalCoor; ++jcoor )
                    {
                        M_d2phiRef[iterNode][icoor][jcoor][iterQuad] = M_refFE->d2Phi ( iterNode,
                                                                                        icoor, jcoor,
                                                                                        M_quadRule->quadPointCoor (iterQuad) );
                    }
                }
            }
        }
        for ( UInt k (0); k < M_nbGeoNode; ++k )
        {
            for ( UInt icoor (0); icoor < M_nbLocalCoor; ++icoor )
            {
                M_dphiGeometricMap[k][icoor][iterQuad] = M_geoMap->dPhi ( k, icoor,
                                                                          M_quadRule->quadPointCoor (iterQuad) );
            }
        }
    }

    M_phiUpdated = true;
    M_dphiRefUpdated = true;
    M_divPhiRefUpdated = true;
    M_d2phiRefUpdated = true;
    M_dphiGeometricMapUpdated = true;
}

void CurrentFE::computeCellNodes ( const std::vector< std::vector<Real> >& pts)
{
    for (UInt iterNode (0); iterNode < M_nbGeoNode; ++iterNode)
        for (UInt iterCoord (0); iterCoord < nDimensions; ++iterCoord)
        {
            M_cellNodes[iterNode][iterCoord] = pts[iterNode][iterCoord];
        }

    M_cellNodesUpdated = true;
}

void CurrentFE::computeQuadNodes()
{
    ASSERT (M_cellNodesUpdated, "Missing update: cellNodes");

    for ( UInt iterQuadNode (0); iterQuadNode < M_nbQuadPt; ++iterQuadNode)
    {
        GeoVector quadNode (coorMap (M_quadRule->quadPointCoor (iterQuadNode) ) );
        for (UInt i = 0; i < quadNode.size(); i++)
        {
            M_quadNodes[iterQuadNode][i] = quadNode[i];
        }
    }
    M_quadNodesUpdated = true;
}

void CurrentFE::computeDphiGeometricMap()
{
    for (UInt iterQuad (0); iterQuad < M_nbQuadPt ; ++iterQuad)
    {
        for (UInt iterNode (0); iterNode < M_nbGeoNode; ++iterNode)
        {
            for (UInt iterCoor (0); iterCoor < M_nbLocalCoor; ++iterCoor)
            {
                M_dphiGeometricMap[iterNode][iterCoor][iterQuad] = M_geoMap->dPhi (iterNode, iterCoor,
                                                                                   M_quadRule->quadPointCoor (iterQuad) );
            }
        }
    }
    M_dphiGeometricMapUpdated = true;
}

void CurrentFE::computeJacobian()
{
    ASSERT (M_cellNodesUpdated, "Missing update: cellNodes");
    ASSERT (M_dphiGeometricMapUpdated, "Missing update: dphiGeometricMap");

    Real partialSum;

    for (UInt iterQuad (0); iterQuad < M_nbQuadPt ; ++iterQuad)
    {
        for (UInt icoord (0); icoord < M_nbLocalCoor; ++icoord)
        {
            for (UInt jcoord (0); jcoord < M_nbLocalCoor; ++jcoord)
            {
                partialSum = 0;
                for (UInt iterNode (0); iterNode < M_nbGeoNode; ++iterNode)
                {
                    partialSum += M_cellNodes[iterNode][icoord] * M_dphiGeometricMap[iterNode][jcoord][iterQuad];
                }
                M_jacobian[icoord][jcoord][iterQuad] = partialSum;
            }
        }
    }
    M_jacobianUpdated = true;
}

void CurrentFE::computeTInverseJacobian()
{
    ASSERT (M_jacobianUpdated, "Missing update: jacobian");

    for (UInt iterQuad (0); iterQuad < M_nbQuadPt ; ++iterQuad)
    {
        switch (M_nbLocalCoor)
        {
            case 1:
            {
                Real a ( M_jacobian[0][0][iterQuad] );

                M_tInverseJacobian[0][0][iterQuad] = 1.0 / a ;

                break;
            }
            case 2:
            {
                Real a ( M_jacobian[0][0][iterQuad] );
                Real b ( M_jacobian[0][1][iterQuad] );
                Real c ( M_jacobian[1][0][iterQuad] );
                Real d ( M_jacobian[1][1][iterQuad] );

                Real det ( a * d - b * c );

                M_tInverseJacobian[0][0][iterQuad] = d / det ;
                M_tInverseJacobian[0][1][iterQuad] = -c / det ;

                M_tInverseJacobian[1][0][iterQuad] =  -b / det ;
                M_tInverseJacobian[1][1][iterQuad] =  a / det ;

                break;
            }
            case 3:
            {
                Real a ( M_jacobian[0][0][iterQuad] );
                Real b ( M_jacobian[0][1][iterQuad] );
                Real c ( M_jacobian[0][2][iterQuad] );
                Real d ( M_jacobian[1][0][iterQuad] );
                Real e ( M_jacobian[1][1][iterQuad] );
                Real f ( M_jacobian[1][2][iterQuad] );
                Real g ( M_jacobian[2][0][iterQuad] );
                Real h ( M_jacobian[2][1][iterQuad] );
                Real i ( M_jacobian[2][2][iterQuad] );

                Real ei ( e * i );
                Real fh ( f * h );
                Real bi ( b * i );
                Real ch ( c * h );
                Real bf ( b * f );
                Real ce ( c * e );
                Real det ( a * ( ei - fh ) + d * ( ch - bi ) + g * ( bf - ce ) );

                M_tInverseJacobian[0][0][iterQuad] = ( ei - fh ) / det ;
                M_tInverseJacobian[0][1][iterQuad] = ( -d * i + f * g ) / det ;
                M_tInverseJacobian[0][2][iterQuad] = ( d * h - e * g ) / det ;

                M_tInverseJacobian[1][0][iterQuad] = ( -bi + ch ) / det ;
                M_tInverseJacobian[1][1][iterQuad] = ( a * i - c * g ) / det ;
                M_tInverseJacobian[1][2][iterQuad] = ( -a * h + b * g ) / det ;

                M_tInverseJacobian[2][0][iterQuad] = ( bf - ce ) / det ;
                M_tInverseJacobian[2][1][iterQuad] = ( -a * f + c * d ) / det ;
                M_tInverseJacobian[2][2][iterQuad] = ( a * e - b * d ) / det ;

                break;
            }
            default:
                ERROR_MSG ( "Dimension (M_nbLocalCoor): only 1, 2 or 3!" );
                break;
        }

    }

    M_tInverseJacobianUpdated = true;
}

void CurrentFE::computeDetJacobian()
{
    ASSERT (M_jacobianUpdated, "Missing update: jacobian");

    for (UInt iterQuad (0); iterQuad < M_nbQuadPt ; ++iterQuad)
    {
        switch (M_nbLocalCoor)
        {
            case 1:
            {
                Real a ( M_jacobian[0][0][iterQuad] );
                Real det ( a);
                M_detJacobian[iterQuad] = det;

                break;
            }
            case 2:
            {
                Real a ( M_jacobian[0][0][iterQuad] );
                Real b ( M_jacobian[0][1][iterQuad] );
                Real c ( M_jacobian[1][0][iterQuad] );
                Real d ( M_jacobian[1][1][iterQuad] );

                Real det ( a * d - b * c);
                M_detJacobian[iterQuad] = det;
                break;
            }
            case 3:
            {
                Real a ( M_jacobian[0][0][iterQuad] );
                Real b ( M_jacobian[0][1][iterQuad] );
                Real c ( M_jacobian[0][2][iterQuad] );
                Real d ( M_jacobian[1][0][iterQuad] );
                Real e ( M_jacobian[1][1][iterQuad] );
                Real f ( M_jacobian[1][2][iterQuad] );
                Real g ( M_jacobian[2][0][iterQuad] );
                Real h ( M_jacobian[2][1][iterQuad] );
                Real i ( M_jacobian[2][2][iterQuad] );

                Real ei (e * i);
                Real fh (f * h);
                Real bi (b * i);
                Real ch (c * h);
                Real bf (b * f);
                Real ce (c * e);

                Real det ( a * (ei - fh) + d * (ch - bi) + g * ( bf - ce) );
                M_detJacobian[iterQuad] = det;
                break;
            }
            default:
                ERROR_MSG ( "Dimension (M_nbLocalCoor): only 1, 2 or 3!" );
                break;
        }
    }
    M_detJacobianUpdated = true;
}

void CurrentFE::computeWDetJacobian()
{
    ASSERT (M_detJacobianUpdated, "Missing update: determinant of the jacobian");

    for (UInt iterQuad (0); iterQuad < M_nbQuadPt ; ++iterQuad)
    {
        M_wDetJacobian[iterQuad] = M_detJacobian[iterQuad] * M_quadRule->weight (iterQuad);
    }
    M_wDetJacobianUpdated = true;
}


void CurrentFE::computeDphiRef()
{
    for (UInt iterQuad (0); iterQuad < M_nbQuadPt ; ++iterQuad)
    {
        for (UInt iterNode (0); iterNode < M_nbNode; ++iterNode)
        {
            for (UInt iterCoor (0); iterCoor < M_nbLocalCoor; ++iterCoor)
            {
                M_dphiRef[iterNode][iterCoor][iterQuad] = M_refFE->dPhi (iterNode, iterCoor,
                                                                         M_quadRule->quadPointCoor (iterQuad) );
            }
        }
    }
    M_dphiRefUpdated = true;
}

void CurrentFE::computeDphi()
{
    ASSERT (M_jacobianUpdated, "Missing update: tInverseJacobian");
    ASSERT (M_dphiRefUpdated, "Missing update: dphiRef");

    Real partialSum (0);

    for (UInt iterQuad (0); iterQuad < M_nbQuadPt ; ++iterQuad)
    {
        for (UInt iterNode (0); iterNode < M_nbNode; ++iterNode)
        {
            for (UInt iterCoor (0); iterCoor < M_nbLocalCoor; ++iterCoor)
            {
                partialSum = 0;
                for (UInt iter (0); iter < M_nbLocalCoor; ++iter)
                {
                    partialSum += M_tInverseJacobian[iterCoor][iter][iterQuad] * M_dphiRef[iterNode][iter][iterQuad];
                }
                M_dphi[iterNode][iterCoor][iterQuad] = partialSum;
            }
        }
    }
    M_dphiUpdated = true;
}


void CurrentFE::computeD2phiRef()
{
    for (UInt iterQuad (0); iterQuad < M_nbQuadPt ; ++iterQuad)
    {
        for (UInt iterNode (0); iterNode < M_nbNode; ++iterNode)
        {
            for (UInt iterCoor (0); iterCoor < M_nbLocalCoor; ++iterCoor)
            {
                for (UInt iterCoor2 (0); iterCoor2 < M_nbLocalCoor; ++iterCoor2)
                {
                    M_d2phiRef[iterNode][iterCoor][iterCoor2][iterQuad] = M_refFE->d2Phi (iterNode, iterCoor, iterCoor2,
                                                                                          M_quadRule->quadPointCoor (iterQuad) );
                }
            }
        }
    }
    M_d2phiRefUpdated = true;
}


void CurrentFE::computeD2phi()
{
    ASSERT (M_jacobianUpdated, "Missing update: tInverseJacobian");
    ASSERT (M_d2phiRefUpdated, "Missing update: d2phiRef");

    Real partialSum (0);
    for ( UInt iterQuad (0); iterQuad < M_nbQuadPt ; ++iterQuad )
    {
        for ( UInt iterNode (0); iterNode < M_nbNode; ++iterNode )
        {
            for ( UInt icoor (0); icoor < M_nbLocalCoor; ++icoor )
            {
                for ( UInt jcoor (0); jcoor < M_nbLocalCoor; ++jcoor )
                {
                    partialSum = 0.;
                    for ( UInt k1 (0); k1 < M_nbLocalCoor; ++k1 )
                    {
                        for ( UInt k2 (0) ; k2 < M_nbLocalCoor; ++k2 )
                        {
                            partialSum += M_tInverseJacobian[icoor][k1][iterQuad]
                                          * M_d2phiRef[iterNode][k1][k2][iterQuad]
                                          * M_tInverseJacobian[jcoor][k2][iterQuad];
                        }
                    }
                    M_d2phi[iterNode][icoor][jcoor][iterQuad] = partialSum;
                }
            }
        }
    }
    M_d2phiUpdated = true;
}


void CurrentFE::computePhiVect()
{
    ASSERT (M_detJacobianUpdated, "Missing update: detJacobian");

    Real sum (0.);
    for ( UInt ig (0); ig < M_nbQuadPt; ++ig )
    {
        for ( UInt idof (0); idof < M_nbNode; ++idof )
        {
            for ( UInt icoor (0); icoor < M_nbLocalCoor; ++icoor )
            {
                sum = 0.;
                for ( UInt jcoor (0); jcoor < M_nbLocalCoor; ++jcoor )
                {
                    sum += M_jacobian[ icoor ][ jcoor ][ ig ] * M_phi[ idof ][ jcoor ][ ig ];
                }
                M_phiVect[ idof ][ icoor ][ ig ] = sum / M_detJacobian[ ig ];
            }
        }
    }

    M_phiVectUpdated = true;
}

} // Namespace LifeV
