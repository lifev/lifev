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
    @brief Adapter for the quadrature on a given level set

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 15 Dec 2011

 */

#ifndef LEVELSETBDQRADAPTER_H
#define LEVELSETBDQRADAPTER_H 1

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/fem/QuadratureBoundary.hpp>
#include <lifev/eta/fem/ETCurrentFE.hpp>

#include <boost/shared_ptr.hpp>

namespace LifeV
{


// Cannot inherit from the volume version because it would
// yield the volume weight, not the surface weight.

template< typename FESpaceType, typename VectorType >
class LevelSetBDQRAdapter
{
public:

    //! @name Public Types
    //@{

    typedef boost::shared_ptr<FESpaceType> FESpaceType_Ptr;

    //@}


    //! @name Constructor & Destructor
    //@{

    LevelSetBDQRAdapter ( FESpaceType_Ptr fespace, const VectorType& vect, const QuadratureBoundary& qrbd)
        : M_lsFESpace (fespace),
          M_lsValue (vect, Repeated),
          M_qrBd (qrbd),
          M_currentFE (fespace->refFE(), getGeometricMap (*fespace->mesh() ), qrbd.qr (0) ),
          M_isAdaptedElement (false),
          M_adaptedQrBd (new QuadratureBoundary (qrbd) )
    {};

    LevelSetBDQRAdapter ( const LevelSetBDQRAdapter& lsbdqra)
        : M_lsFESpace ( lsbdqra.M_lsFESpace),
          M_lsValue ( lsbdqra.M_lsValue),
          M_qrBd ( lsbdqra.M_qrBd),
          M_currentFE ( lsbdqra.M_currentFE),
          M_isAdaptedElement ( lsbdqra.M_isAdaptedElement),
          M_adaptedQrBd (new QuadratureBoundary (lsbdqra.M_qrBd) )
    {};

    ~LevelSetBDQRAdapter() {}

    //@}


    //! @name Methods
    //@{

    void update (UInt elementID, UInt localFaceID);

    //@}


    //! @name Set Methods
    //@{

    //@}


    //! @name Get Methods
    //@{

    QuadratureRule adaptedBdQR ( const UInt i) const
    {
        if (this->M_isAdaptedElement)
        {
            return M_adaptedQrBd->qr (i);
        }
        return M_qrBd.qr (i);
    }

    //@}

private:

    LevelSetBDQRAdapter();


    FESpaceType_Ptr M_lsFESpace;

    VectorType M_lsValue;

    QuadratureBoundary M_qrBd;

    ETCurrentFE<FESpaceType::S_spaceDim, 1> M_currentFE;

    bool M_isAdaptedElement;

    boost::shared_ptr<QuadratureBoundary> M_adaptedQrBd;

};

template<typename FESpaceType, typename VectorType>
LevelSetBDQRAdapter<FESpaceType, VectorType>
adapt (boost::shared_ptr<FESpaceType> fespace, const VectorType& vector, const QuadratureBoundary& qrbd)
{
    return LevelSetBDQRAdapter<FESpaceType, VectorType> (fespace, vector, qrbd);
}


template<typename FESpaceType, typename VectorType>
void
LevelSetBDQRAdapter<FESpaceType, VectorType>::
update (UInt elementID, UInt localFaceID)
{
    // They are ordered so that the "mapping" in
    // QuadratureBoundary is right.

    std::vector<UInt> myLID (3);

    if (localFaceID == 0)
    {
        myLID[0] = 0;
        myLID[1] = 1;
        myLID[2] = 2;
    }
    else if (localFaceID == 1)
    {
        myLID[0] = 0;
        myLID[1] = 1;
        myLID[2] = 3;
    }
    else if (localFaceID == 2)
    {
        myLID[0] = 3;
        myLID[1] = 1;
        myLID[2] = 2;
    }
    else
    {
        myLID[0] = 0;
        myLID[1] = 3;
        myLID[2] = 2;
    }

    // Compute the global IDs
    std::vector<UInt> myGID (3);

    for (UInt i (0); i < 3; ++i)
    {
        myGID[i] = M_lsFESpace->dof().localToGlobalMap (elementID, myLID[i]);
    }

    // Check the level set values
    std::vector<UInt> localPositive;
    std::vector<UInt> localNegative;
    std::vector<Real> localPositiveValue;
    std::vector<Real> localNegativeValue;

    for (UInt i (0); i < 3; ++i)
    {
        if (M_lsValue[myGID[i]] >= 0)
        {
            localPositive.push_back (i);
            localPositiveValue.push_back (M_lsValue[myGID[i]]);
        }
        else
        {
            localNegative.push_back (i);
            localNegativeValue.push_back (M_lsValue[myGID[i]]);
        }
    }

    // Check if the element is actually cut by the interface.
    if (localPositive.empty() || localNegative.empty() )
    {
        M_isAdaptedElement = false;
        *M_adaptedQrBd =  M_qrBd;
    }
    else
    {

        M_isAdaptedElement = true;

        std::vector<std::vector<Real> >refCoor (3, std::vector<Real> (2) );
        refCoor[0][0] = 0.0;
        refCoor[0][1] = 0.0;
        refCoor[1][0] = 1.0;
        refCoor[1][1] = 0.0;
        refCoor[2][0] = 0.0;
        refCoor[2][1] = 1.0;

        QuadratureRule qrRef ("none", TRIANGLE, 3, 0, 0);

        // Disinguish the two cases
        if (localPositive.size() == 1)
        {
            // Do nothing
        }
        else
        {
            // Exchange positive and negative values
            localPositive.swap (localNegative);
            localPositiveValue.swap (localNegativeValue);
        }


        // Two negative case
        Real lambda1 = -localPositiveValue[0] / (localNegativeValue[0] - localPositiveValue[0]);

        std::vector < Real > I1 (2);
        I1[0] = refCoor[localPositive[0]][0] + lambda1 * (refCoor[localNegative[0]][0] - refCoor[localPositive[0]][0]);
        I1[1] = refCoor[localPositive[0]][1] + lambda1 * (refCoor[localNegative[0]][1] - refCoor[localPositive[0]][1]);

        Real lambda2 = -localPositiveValue[0] / (localNegativeValue[1] - localPositiveValue[0]);

        std::vector < Real > I2 (2);
        I2[0] = refCoor[localPositive[0]][0] + lambda2 * (refCoor[localNegative[1]][0] - refCoor[localPositive[0]][0]);
        I2[1] = refCoor[localPositive[0]][1] + lambda2 * (refCoor[localNegative[1]][1] - refCoor[localPositive[0]][1]);

        // First triangle
        {
            std::vector<Real> P0 (2);
            P0[0] = refCoor[localPositive[0]][0];
            P0[1] = refCoor[localPositive[0]][1];

            std::vector<Real> v1 (2);
            v1[0] = I1[0] - refCoor[localPositive[0]][0];
            v1[1] = I1[1] - refCoor[localPositive[0]][1];

            std::vector<Real> v2 (2);
            v2[0] = I2[0] - refCoor[localPositive[0]][0];
            v2[1] = I2[1] - refCoor[localPositive[0]][1];

            Real wRel (std::abs (v1[0]*v2[1] - v1[1]*v2[0]) );

            for (UInt iq (0); iq < M_qrBd.qr (0).nbQuadPt(); ++iq)
            {
                Real x (P0[0] + M_qrBd.qr (0).quadPointCoor (iq, 0) *v1[0] + M_qrBd.qr (0).quadPointCoor (iq, 1) *v2[0]);
                Real y (P0[1] + M_qrBd.qr (0).quadPointCoor (iq, 0) *v1[1] + M_qrBd.qr (0).quadPointCoor (iq, 1) *v2[1]);

                qrRef.addPoint (QuadraturePoint (x, y, M_qrBd.qr (0).weight (iq) *wRel) );
            }
        }

        // Second triangle
        {
            std::vector<Real> P0 (2);
            P0[0] = refCoor[localNegative[0]][0];
            P0[1] = refCoor[localNegative[0]][1];

            std::vector<Real> v1 (2);
            v1[0] = I1[0] - refCoor[localNegative[0]][0];
            v1[1] = I1[1] - refCoor[localNegative[0]][1];

            std::vector<Real> v2 (2);
            v2[0] = I2[0] - refCoor[localNegative[0]][0];
            v2[1] = I2[1] - refCoor[localNegative[0]][1];

            Real wRel (std::abs (v1[0]*v2[1] - v1[1]*v2[0]) );

            for (UInt iq (0); iq < M_qrBd.qr (0).nbQuadPt(); ++iq)
            {
                Real x (P0[0] + M_qrBd.qr (0).quadPointCoor (iq, 0) *v1[0] + M_qrBd.qr (0).quadPointCoor (iq, 1) *v2[0]);
                Real y (P0[1] + M_qrBd.qr (0).quadPointCoor (iq, 0) *v1[1] + M_qrBd.qr (0).quadPointCoor (iq, 1) *v2[1]);

                qrRef.addPoint (QuadraturePoint (x, y, M_qrBd.qr (0).weight (iq) *wRel) );
            }
        }

        // Third triangle
        {
            std::vector<Real> P0 (2);
            P0[0] = refCoor[localNegative[0]][0];
            P0[1] = refCoor[localNegative[0]][1];

            std::vector<Real> v1 (2);
            v1[0] = I2[0] - refCoor[localNegative[0]][0];
            v1[1] = I2[1] - refCoor[localNegative[0]][1];

            std::vector<Real> v2 (2);
            v2[0] = refCoor[localNegative[1]][0] - refCoor[localNegative[0]][0];
            v2[1] = refCoor[localNegative[1]][1] - refCoor[localNegative[0]][1];

            Real wRel (std::abs (v1[0]*v2[1] - v1[1]*v2[0]) );

            for (UInt iq (0); iq < M_qrBd.qr (0).nbQuadPt(); ++iq)
            {
                Real x (P0[0] + M_qrBd.qr (0).quadPointCoor (iq, 0) *v1[0] + M_qrBd.qr (0).quadPointCoor (iq, 1) *v2[0]);
                Real y (P0[1] + M_qrBd.qr (0).quadPointCoor (iq, 0) *v1[1] + M_qrBd.qr (0).quadPointCoor (iq, 1) *v2[1]);

                qrRef.addPoint (QuadraturePoint (x, y, M_qrBd.qr (0).weight (iq) *wRel) );
            }
        }

        // Call the make tretra
        M_adaptedQrBd.reset (new QuadratureBoundary (buildTetraBDQR (qrRef) ) );
    }
}


} // Namespace LifeV

#endif /* LEVELSETBDQRADAPTER_H */
