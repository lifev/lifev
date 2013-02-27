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

    @refactoring Luca Bertagna <lbertag@emory.edu>
    @date Sept 2012

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#include <lifev/core/fem/CurrentBoundaryFE.hpp>

namespace LifeV
{

// =================================================== //
//               Constructor(s) & Destructor           //
// =================================================== //

CurrentBoundaryFE::CurrentBoundaryFE (const ReferenceFE& refFE, const GeometricMap& geoMap, const QuadratureRule& qr ) :
    CurrentFE               (refFE, geoMap, qr),
    M_phiGeo                (boost::extents[M_nbGeoNode][M_nbQuadPt]),
    M_tangents              (boost::extents[M_nbCoor][nDimensions][M_nbQuadPt]),
    M_normal                (boost::extents[nDimensions][M_nbQuadPt]),
    M_metric                (boost::extents[M_nbCoor][M_nbCoor][M_nbQuadPt]),
    M_detMetric             (boost::extents[M_nbQuadPt]),
    M_inverseMetric         (boost::extents[M_nbCoor][M_nbCoor][M_nbQuadPt]),
    M_wRootDetMetric        (boost::extents[M_nbQuadPt]),
    M_tangentsUpdated       (false),
    M_normalUpdated         (false),
    M_metricUpdated         (false),
    M_detMetricUpdated      (false),
    M_inverseMetricUpdated  (false),
    M_wRootDetMetricUpdated (false)
{
    // Nothing to be done here
}

CurrentBoundaryFE::CurrentBoundaryFE (const ReferenceFE& refFE, const GeometricMap& geoMap) :
    CurrentFE         (refFE, geoMap),
    M_phiGeo                (boost::extents[M_nbGeoNode][M_nbQuadPt]),
    M_tangents              (boost::extents[M_nbCoor][nDimensions][M_nbQuadPt]),
    M_normal                (boost::extents[nDimensions][M_nbQuadPt]),
    M_metric                (boost::extents[M_nbCoor][M_nbCoor][M_nbQuadPt]),
    M_detMetric             (boost::extents[M_nbQuadPt]),
    M_inverseMetric         (boost::extents[M_nbCoor][M_nbCoor][M_nbQuadPt]),
    M_wRootDetMetric        (boost::extents[M_nbQuadPt]),
    M_tangentsUpdated       (false),
    M_normalUpdated         (false),
    M_metricUpdated         (false),
    M_detMetricUpdated      (false),
    M_inverseMetricUpdated  (false),
    M_wRootDetMetricUpdated (false)
{
    // Nothing to be done here
}

CurrentBoundaryFE::CurrentBoundaryFE (const ReferenceFE& refFE, const GeometricMap& geoMap,
                                      const QuadratureRule& qr, const Real* refCoor,
                                      UInt currentId, Real invArea ) :
    CurrentFE               (refFE, geoMap, qr),
    M_phiGeo                (boost::extents[M_nbGeoNode][M_nbQuadPt]),
    M_tangents              (boost::extents[M_nbCoor][nDimensions][M_nbQuadPt]),
    M_normal                (boost::extents[nDimensions][M_nbQuadPt]),
    M_metric                (boost::extents[M_nbCoor][M_nbCoor][M_nbQuadPt]),
    M_detMetric             (boost::extents[M_nbQuadPt]),
    M_inverseMetric         (boost::extents[M_nbCoor][M_nbCoor][M_nbQuadPt]),
    M_wRootDetMetric        (boost::extents[M_nbQuadPt]),
    M_tangentsUpdated       (false),
    M_normalUpdated         (false),
    M_metricUpdated         (false),
    M_detMetricUpdated      (false),
    M_inverseMetricUpdated  (false),
    M_wRootDetMetricUpdated (false)
{
    M_currentId = currentId;

    for (UInt iterQuad (0); iterQuad < M_nbQuadPt; ++iterQuad)
    {
        for (UInt iDof (0); iDof < M_nbNode; ++iDof)
        {
            M_phi[iDof][0][iterQuad] *= invArea;
        }

        for (UInt iGeoNode (0); iGeoNode < M_nbGeoNode; ++iGeoNode)
        {
            M_phiGeo[iGeoNode][iterQuad] = geoMap.phi (iGeoNode, M_quadRule->quadPointCoor (iterQuad) );
        }
    }

    for (UInt iNode (0); iNode < M_nbGeoNode; ++iNode)
        for (UInt iCoor (0); iCoor < nDimensions; ++iCoor)
        {
            M_cellNodes[iNode][iCoor] = refCoor[nDimensions * iNode + iCoor];
        }

    M_cellNodesUpdated = true;

    // This constructor is used for hybrid FE, where info about normals, tangents and metric
    // are needed only on the reference element. So we update them on the given reference coordinates

    computeQuadNodes();
    computeTangents();
    computeNormal();
    computeMetric();
    computeDetMetric();
    computeInverseMetric();
    computeWRootDetMetric();
}

CurrentBoundaryFE::CurrentBoundaryFE (const CurrentBoundaryFE& bdFE) :
    CurrentFE               (bdFE),
    M_phiGeo                (bdFE.M_phiGeo),
    M_tangents              (bdFE.M_tangents),
    M_normal                (bdFE.M_normal),
    M_metric                (bdFE.M_metric),
    M_detMetric             (bdFE.M_detMetric),
    M_inverseMetric         (bdFE.M_inverseMetric),
    M_wRootDetMetric        (bdFE.M_wRootDetMetric),
    M_tangentsUpdated       (bdFE.M_tangentsUpdated),
    M_normalUpdated         (bdFE.M_normalUpdated),
    M_metricUpdated         (bdFE.M_metricUpdated),
    M_detMetricUpdated      (bdFE.M_detMetricUpdated),
    M_inverseMetricUpdated  (bdFE.M_inverseMetricUpdated),
    M_wRootDetMetricUpdated (bdFE.M_wRootDetMetricUpdated)
{
    M_currentId = bdFE.M_currentId;

    if (bdFE.M_quadRule)
    {
        setQuadRule (*bdFE.M_quadRule);
    }
}

CurrentBoundaryFE::~CurrentBoundaryFE()
{
    // Nothing to be done here
}

// =================================================== //
//                        Methods                      //
// =================================================== //

void CurrentBoundaryFE::coorMap (Real& x, Real& y, Real& z, Real xi, Real eta) const
{
    CurrentFE::coorMap (x, y, z, xi, eta, 0.);
}

Real CurrentBoundaryFE::measure() const
{
    ASSERT (M_wRootDetMetricUpdated, "Weighted root of metric determinant is not updated!");

    Real meas = 0.0;
    for (UInt quadNode (0); quadNode < M_nbQuadPt; ++quadNode)
    {
        meas += M_wRootDetMetric[quadNode];
    }
    return meas;
}

void CurrentBoundaryFE::update (const std::vector<std::vector<Real> >& pts, flag_Type upFlag)
{
    CurrentFE::update (pts, upFlag);

    M_tangentsUpdated = false;
    if ( (upFlag & UPDATE_ONLY_TANGENTS) != 0)
    {
        computeTangents();
    };

    M_normalUpdated = false;
    if ( (upFlag & UPDATE_ONLY_NORMALS) != 0)
    {
        computeNormal();
    };

    M_metricUpdated = false;
    if ( (upFlag & UPDATE_ONLY_METRIC) != 0)
    {
        computeMetric();
    };

    M_detMetricUpdated = false;
    if ( (upFlag & UPDATE_ONLY_DET_METRIC) != 0)
    {
        computeDetMetric();
    };

    M_inverseMetricUpdated = false;
    if ( (upFlag & UPDATE_ONLY_INV_METRIC) != 0)
    {
        computeInverseMetric();
    };

    M_wRootDetMetricUpdated = false;
    if ( (upFlag & UPDATE_ONLY_W_ROOT_DET_METRIC) != 0)
    {
        computeWRootDetMetric();
    };
}

// =============================================== //
//                Protected Methods                //
// =============================================== //

void CurrentBoundaryFE::computeTangents()
{
    ASSERT (M_cellNodesUpdated, "Missing update: cellNodes\n");
    ASSERT (M_dphiGeometricMapUpdated, "Missing update: dphiGeometricMap");

    for (UInt iq (0); iq < M_nbQuadPt; ++iq)
        for (UInt iTangent (0); iTangent < M_nbCoor; ++iTangent)
            for (UInt iCoor (0); iCoor < nDimensions; ++iCoor)
            {
                M_tangents[iTangent][iCoor][iq] = 0.0;
                for (UInt iNode (0); iNode < M_nbGeoNode; ++iNode)
                {
                    M_tangents[iTangent][iCoor][iq] += M_cellNodes[iNode][iCoor] * M_dphiGeometricMap[iNode][iTangent][iq];
                }
            }

    M_tangentsUpdated = true;
}

void CurrentBoundaryFE::computeNormal()
{
    ASSERT (M_tangentsUpdated, "Missing update: tangents\n");

    Real norm, n1, n2;

    // So far I can't find a general method for dimension n, so we need to switch on M_nbCoor
    switch (M_nbCoor)
    {
        case 1:
            for (UInt iq (0); iq < M_nbQuadPt; ++iq)
            {
                n1 = M_tangents[0][1][iq];
                n2 = - M_tangents[0][0][iq];

                norm = sqrt (n1 * n1 + n2 * n2);
                ASSERT (norm > 0, "Error! Something went wrong while computing normals.\n");

                M_normal[0][iq] = n1 / norm;
                M_normal[1][iq] = n2 / norm;
            }
            break;
        case 2:
            Real n3;
            for (UInt iq (0); iq < M_nbQuadPt; ++iq)
            {
                n1 = M_tangents[0][1][iq] * M_tangents[1][2][iq] - M_tangents[0][2][iq] * M_tangents[1][1][iq];
                n2 = M_tangents[0][2][iq] * M_tangents[1][0][iq] - M_tangents[0][0][iq] * M_tangents[1][2][iq];
                n3 = M_tangents[0][0][iq] * M_tangents[1][1][iq] - M_tangents[0][1][iq] * M_tangents[1][0][iq];

                norm = sqrt (n1 * n1 + n2 * n2 + n3 * n3);
                ASSERT (norm > 0, "Error! Something went wrong while computing normals.\n");

                M_normal[0][iq] = n1 / norm;
                M_normal[1][iq] = n2 / norm;
                M_normal[2][iq] = n3 / norm;
            }
            break;
        default:
            ASSERT (0, "No rule to compute normal vector for this dimension.");
    }

    M_normalUpdated = true;
}

void CurrentBoundaryFE::computeMetric()
{
    ASSERT (M_tangentsUpdated, "Missing update: tangents\n");

    for (UInt iq (0); iq < M_nbQuadPt; ++iq)
        for (UInt iTangent (0); iTangent < M_nbCoor; ++iTangent)
        {
            // Diagonal part
            M_metric[iTangent][iTangent][iq] = 0.0;
            for (UInt iCoor (0); iCoor < M_nbCoor + 1; ++iCoor)
            {
                M_metric[iTangent][iTangent][iq] += M_tangents[iTangent][iCoor][iq] * M_tangents[iTangent][iCoor][iq];
            }

            // Extra diagonal part
            for (UInt jTangent (0); jTangent < iTangent; ++jTangent)
            {
                M_metric[iTangent][jTangent][iq] = 0.0;
                for (UInt iCoor (0); iCoor < M_nbCoor + 1; ++iCoor)
                {
                    M_metric[iTangent][jTangent][iq] += M_tangents[iTangent][iCoor][iq] * M_tangents[jTangent][iCoor][iq];
                }
                M_metric[jTangent][iTangent][iq] = M_metric[iTangent][jTangent][iq];
            }
        }

    M_metricUpdated = true;
}

void CurrentBoundaryFE::computeDetMetric()
{
    ASSERT (M_metricUpdated, "Missing update: metric\n");

    for (UInt iq (0); iq < M_nbQuadPt; ++iq)
    {
        switch (M_nbCoor)
        {
            case 1:
                M_detMetric[iq] = M_metric[0][0][iq];
                break;
            case 2:
                M_detMetric[iq] = M_metric[0][0][iq] * M_metric[1][1][iq] - M_metric[0][1][iq] * M_metric[1][0][iq];
                break;
            default:
                ASSERT (0, "No rule to compute determinant of the metric for this dimension.\n");
        }
    }
    M_detMetricUpdated = true;
}

void CurrentBoundaryFE::computeInverseMetric()
{
    ASSERT (M_metricUpdated, "Missing update: metric\n");
    ASSERT (M_detMetricUpdated, "Missing update: detMetric\n");

    for (UInt iq (0); iq < M_nbQuadPt; ++iq)
    {
        switch (M_nbCoor)
        {
            case 1:
                M_inverseMetric[0][0][iq] = 1. / M_metric[0][0][iq];
                break;
            case 2:
                // Diagonal part: switch the two diagonal entries and divide by determinant
                M_inverseMetric[0][0][iq] = M_metric[1][1][iq] / M_detMetric[iq];
                M_inverseMetric[1][1][iq] = M_metric[0][0][iq] / M_detMetric[iq];

                // Extradiagonal part: change sign and divide by determinant (matrix is symmetric, so do it once)
                M_inverseMetric[0][1][iq] = M_inverseMetric[1][0][iq] = - M_metric[0][1][iq] / M_detMetric[iq];
                break;
            default:
                ASSERT (0, "No rule to compute the inverse of the metric for this dimension.\n");
        }
    }

    M_inverseMetricUpdated = true;
}

void CurrentBoundaryFE::computeWRootDetMetric()
{
    ASSERT (M_metricUpdated, "Missing update: metric\n");
    ASSERT (M_detMetricUpdated, "Missing update: detMetric\n");

    for (UInt iq (0); iq < M_nbQuadPt; ++iq)
    {
        M_wRootDetMetric[iq] = std::sqrt (M_detMetric[iq]) * M_quadRule->weight (iq);
    }

    M_wRootDetMetricUpdated = true;
}

} // Namespace LifeV
