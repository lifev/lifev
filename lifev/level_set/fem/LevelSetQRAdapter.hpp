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
    @brief Adapter for the quadrature on a given level set.

    This class is intended to be used in conjonction with the ETA framework
    (this module must be available to use this class). Based on the values
    of a level function, it computes (via a piecewise linear reconstruction
    of the interface) a quadrature rule which is adapted to the 0 level set.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 15 Dec 2011

 */


#ifndef LEVELSETQRADAPTER_H
#define LEVELSETQRADAPTER_H 1

#include <lifev/core/LifeV.hpp>

#ifdef LIFEV_HAS_ETA

#include <lifev/eta/fem/ETCurrentFE.hpp>
#include <lifev/core/fem/GeometricMap.hpp>
#include <lifev/core/fem/QuadratureRule.hpp>

#include <lifev/eta/fem/QRAdapterBase.hpp>

#include <boost/shared_ptr.hpp>

namespace LifeV
{


/*!
  The class LevelSetQRAdapter

  P1 reconstruction of the interface
  and quadrature splitting

  3D only for the moment and P1

 */
template< typename FESpaceType, typename VectorType >
class LevelSetQRAdapter : public QRAdapterBase< LevelSetQRAdapter<FESpaceType, VectorType> >
{
public:

    //! @name Public Types
    //@{

    typedef boost::shared_ptr<FESpaceType> FESpaceType_Ptr;

    typedef QRAdapterBase< LevelSetQRAdapter<FESpaceType, VectorType> > base_Type;

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Constructor with both FESpace and level set values
    LevelSetQRAdapter (FESpaceType_Ptr fespace, const VectorType& vect, const QuadratureRule& qr);

    //! Copy constructor
    LevelSetQRAdapter (const LevelSetQRAdapter<FESpaceType, VectorType>& lsqra);

    //! Simple destructor
    ~LevelSetQRAdapter() {}

    //@}


    //! @name Operators
    //@{

    //@}


    //! @name Methods
    //@{

    //! Update this structure with the current element
    /*!
      This checks if the element is crossed by the interface
      and if so, computes the adapted quadrature rule.
     */
    void update (UInt elementID);


    //@}


    //! @name Set Methods
    //@{

    //@}


    //! @name Get Methods
    //@{

    //! Is the current element crossed by the interface
    bool isAdaptedElement() const
    {
        return M_isAdaptedElement;
    }

    //! Getter for the non-adapted quadrature
    const QuadratureRule& standardQR() const
    {
        return M_qr;
    }

    //! Getter for the adapted quadrature
    const QuadratureRule& adaptedQR() const
    {
        return *M_adaptedQR;
    }


    //@}

private:

    //! @name Private Methods
    //@{

    /*! No default constructor as there is no default
      constructor in the ETCurrentFE
     */
    LevelSetQRAdapter();

    //@}

    // Finite element space for the level set values
    FESpaceType_Ptr M_lsFESpace;

    // Value of the level set
    VectorType M_lsValue;

    // Quadrature rule
    QuadratureRule M_qr;

    // CurrentFE for the level set

    ETCurrentFE<FESpaceType::space_dim, 1> M_currentFE;

    // Boolean indicating if the element is crossed by the interface
    bool M_isAdaptedElement;

    // Adapted Quadrature Rule
    boost::shared_ptr<QuadratureRule> M_adaptedQR;

};

template<typename FESpaceType, typename VectorType>
LevelSetQRAdapter<FESpaceType, VectorType>
adapt (boost::shared_ptr<FESpaceType> fespace, const VectorType& vector, const QuadratureRule& qr)
{
    return LevelSetQRAdapter<FESpaceType, VectorType> (fespace, vector, qr);
}


template< typename FESpaceType, typename VectorType >
LevelSetQRAdapter<FESpaceType, VectorType>::
LevelSetQRAdapter (FESpaceType_Ptr fespace, const VectorType& vect, const QuadratureRule& qr)
    :
    base_Type(),
    M_lsFESpace (fespace),
    M_lsValue (vect, Repeated),
    M_qr (qr),
    M_currentFE (fespace->refFE(), getGeometricMap (*fespace->mesh() ), qr),
    M_isAdaptedElement (false),
    M_adaptedQR (new QuadratureRule (qr) )
{

    ASSERT ( FESpaceType::field_dim == 1, "Quadrature adaptation only for scalar fields!");
    ASSERT ( FESpaceType::space_dim == 3, "Quadrature adaptation only 3D cases for the moment!");
}

template< typename FESpaceType, typename VectorType >
LevelSetQRAdapter<FESpaceType, VectorType>::
LevelSetQRAdapter (const LevelSetQRAdapter<FESpaceType, VectorType>& lsqra)
    :
    base_Type(),
    M_lsFESpace ( lsqra.M_lsFESpace ),
    M_lsValue ( lsqra.M_lsValue ),
    M_qr ( lsqra.M_qr ),
    M_currentFE ( lsqra.M_currentFE ),
    M_isAdaptedElement ( lsqra.M_isAdaptedElement ),
    M_adaptedQR ( new QuadratureRule (*lsqra.M_adaptedQR) )
{}


template< typename FESpaceType, typename VectorType >
void
LevelSetQRAdapter<FESpaceType, VectorType>::
update (UInt elementID)
{
    //! Check that the level set values are Repeated
    ASSERT ( M_lsValue.mapType() == Repeated, "Internal error: level values are Unique!");
    ASSERT ( M_lsFESpace != 0, "Internal error: empty pointer to the Space");

    // First, check the values of the level set and detect sign changes

    std::vector<UInt> localPositive;
    std::vector<UInt> localNegative;
    std::vector<Real> localLSValue (4, 0.0);

    for (UInt iDof (0); iDof < 4; ++iDof)
    {
        UInt myGlobalID ( M_lsFESpace->dof().localToGlobalMap (elementID, iDof) );
        localLSValue[iDof] = M_lsValue[myGlobalID];

        if (localLSValue[iDof] >= 0)
        {
            localPositive.push_back (iDof);
        }
        else
        {
            localNegative.push_back (iDof);
        }
    }

    // Check if adaptation is needed
    if ( localPositive.empty() || localNegative.empty() )
    {
        // No adaptation
        M_isAdaptedElement = false;
        *M_adaptedQR = M_qr;
    }
    else
    {
        // Adaptation
        M_isAdaptedElement = true;

        // Define the reference tetrahedra
        std::vector< VectorSmall<3> > refPts (4);
        refPts[1][0] = 1;
        refPts[2][1] = 1;
        refPts[3][2] = 1;

        // Define the CurrentFE that will be used to map the
        // quadrature from its reference shape to the part
        // of the tetrahedra
        ETCurrentFE<3, 1> QRMapper (feTetraP1, getGeometricMap (*M_lsFESpace->mesh() ), M_qr);

        // Reset the adapted quadrature
        M_adaptedQR.reset (new QuadratureRule);
        M_adaptedQR->setDimensionShape (3, TETRA);

        // Constant
        UInt QuadratureSize (M_qr.nbQuadPt() );

        // If there's only one negative value
        // the tetrahedra has to be cut into 4
        // pieces
        if (localPositive.size() == 1)
        {
            // Localize the intersections
            std::vector< VectorSmall<3> > intersections (3);

            for (UInt i (0); i < 3; ++i)
            {
                UInt pos (localPositive[0]);
                UInt neg (localNegative[i]);
                Real mu ( localLSValue[pos] / (localLSValue[pos] - localLSValue[neg]) );
                VectorSmall<3> I ( refPts[pos] + mu * ( refPts[neg] - refPts[pos] ) );
                intersections[i] = I;
            }


            // Build the first tetra
            std::vector< VectorSmall<3> > tetra1 (4);
            tetra1[0] = refPts[localPositive[0]];
            tetra1[1] = intersections[0];
            tetra1[2] = intersections[1];
            tetra1[3] = intersections[2];

            // Map the QR
            QRMapper.update ( tetra1, ET_UPDATE_QUAD_NODE | ET_UPDATE_WDET);

            for (UInt iQuadPt (0); iQuadPt < QuadratureSize; ++iQuadPt)
            {
                M_adaptedQR->addPoint (QuadraturePoint ( QRMapper.quadNode (iQuadPt, 0),
                                                         QRMapper.quadNode (iQuadPt, 1),
                                                         QRMapper.quadNode (iQuadPt, 2),
                                                         std::abs (QRMapper.wDet (iQuadPt) ) ) );
            }

            // Build the second tetra
            std::vector< VectorSmall<3> > tetra2 (4);
            tetra2[0] = intersections[0];
            tetra2[1] = intersections[1];
            tetra2[2] = intersections[2];
            tetra2[3] = refPts[localNegative[0]];

            // Map the QR
            QRMapper.update ( tetra2, ET_UPDATE_QUAD_NODE | ET_UPDATE_WDET);

            for (UInt iQuadPt (0); iQuadPt < QuadratureSize; ++iQuadPt)
            {
                M_adaptedQR->addPoint (QuadraturePoint ( QRMapper.quadNode (iQuadPt, 0),
                                                         QRMapper.quadNode (iQuadPt, 1),
                                                         QRMapper.quadNode (iQuadPt, 2),
                                                         std::abs (QRMapper.wDet (iQuadPt) ) ) );
            }

            // Build the third tetra
            std::vector< VectorSmall<3> > tetra3 (4);
            tetra3[0] = refPts[localNegative[0]];
            tetra3[1] = refPts[localNegative[1]];
            tetra3[2] = refPts[localNegative[2]];
            tetra3[3] = intersections[2];

            // Map the QR
            QRMapper.update ( tetra3, ET_UPDATE_QUAD_NODE | ET_UPDATE_WDET);

            for (UInt iQuadPt (0); iQuadPt < QuadratureSize; ++iQuadPt)
            {
                M_adaptedQR->addPoint (QuadraturePoint ( QRMapper.quadNode (iQuadPt, 0),
                                                         QRMapper.quadNode (iQuadPt, 1),
                                                         QRMapper.quadNode (iQuadPt, 2),
                                                         std::abs (QRMapper.wDet (iQuadPt) ) ) );
            }


            // Build the fourth tetra
            std::vector< VectorSmall<3> > tetra4 (4);
            tetra4[0] = refPts[localNegative[0]];
            tetra4[1] = refPts[localNegative[1]];
            tetra4[2] = intersections[1];
            tetra4[3] = intersections[2];

            // Map the QR
            QRMapper.update ( tetra4, ET_UPDATE_QUAD_NODE | ET_UPDATE_WDET);

            for (UInt iQuadPt (0); iQuadPt < QuadratureSize; ++iQuadPt)
            {
                M_adaptedQR->addPoint (QuadraturePoint ( QRMapper.quadNode (iQuadPt, 0),
                                                         QRMapper.quadNode (iQuadPt, 1),
                                                         QRMapper.quadNode (iQuadPt, 2),
                                                         std::abs (QRMapper.wDet (iQuadPt) ) ) );
            }

        }
        else if ( localPositive.size() == 2)
        {

            // Localize the intersections
            std::vector< VectorSmall<3> > intersections (4);

            for (UInt i (0); i < 2; ++i)
            {
                for (UInt j (0); j < 2; ++j)
                {
                    UInt pos (localPositive[i]);
                    UInt neg (localNegative[j]);
                    Real mu ( localLSValue[pos] / (localLSValue[pos] - localLSValue[neg]) );
                    VectorSmall<3> I ( refPts[pos] + mu * ( refPts[neg] - refPts[pos] ) );
                    intersections[i + 2 * j] = I;
                }
            }


            // Build the first tetra
            std::vector< VectorSmall<3> > tetra1 (4);
            tetra1[0] = refPts[localPositive[0]];
            tetra1[1] = intersections[0];
            tetra1[2] = intersections[1];
            tetra1[3] = intersections[2];

            // Map the QR
            QRMapper.update ( tetra1, ET_UPDATE_QUAD_NODE | ET_UPDATE_WDET);

            for (UInt iQuadPt (0); iQuadPt < QuadratureSize; ++iQuadPt)
            {
                M_adaptedQR->addPoint (QuadraturePoint ( QRMapper.quadNode (iQuadPt, 0),
                                                         QRMapper.quadNode (iQuadPt, 1),
                                                         QRMapper.quadNode (iQuadPt, 2),
                                                         std::abs (QRMapper.wDet (iQuadPt) ) ) );
            }

            // Build the second tetra
            std::vector< VectorSmall<3> > tetra2 (4);
            tetra2[0] = intersections[1];
            tetra2[1] = intersections[2];
            tetra2[2] = intersections[3];
            tetra2[3] = refPts[localPositive[0]];

            // Map the QR
            QRMapper.update ( tetra2, ET_UPDATE_QUAD_NODE | ET_UPDATE_WDET);

            for (UInt iQuadPt (0); iQuadPt < QuadratureSize; ++iQuadPt)
            {
                M_adaptedQR->addPoint (QuadraturePoint ( QRMapper.quadNode (iQuadPt, 0),
                                                         QRMapper.quadNode (iQuadPt, 1),
                                                         QRMapper.quadNode (iQuadPt, 2),
                                                         std::abs (QRMapper.wDet (iQuadPt) ) ) );
            }

            // Build the third tetra
            std::vector< VectorSmall<3> > tetra3 (4);
            tetra3[0] = refPts[localPositive[0]];
            tetra3[1] = refPts[localPositive[1]];
            tetra3[2] = intersections[1];
            tetra3[3] = intersections[3];

            // Map the QR
            QRMapper.update ( tetra3, ET_UPDATE_QUAD_NODE | ET_UPDATE_WDET);

            for (UInt iQuadPt (0); iQuadPt < QuadratureSize; ++iQuadPt)
            {
                M_adaptedQR->addPoint (QuadraturePoint ( QRMapper.quadNode (iQuadPt, 0),
                                                         QRMapper.quadNode (iQuadPt, 1),
                                                         QRMapper.quadNode (iQuadPt, 2),
                                                         std::abs (QRMapper.wDet (iQuadPt) ) ) );
            }

            // Build the fourth tetra
            std::vector< VectorSmall<3> > tetra4 (4);
            tetra4[0] = refPts[localNegative[0]];
            tetra4[1] = intersections[0];
            tetra4[2] = intersections[1];
            tetra4[3] = intersections[2];

            // Map the QR
            QRMapper.update ( tetra4, ET_UPDATE_QUAD_NODE | ET_UPDATE_WDET);

            for (UInt iQuadPt (0); iQuadPt < QuadratureSize; ++iQuadPt)
            {
                M_adaptedQR->addPoint (QuadraturePoint ( QRMapper.quadNode (iQuadPt, 0),
                                                         QRMapper.quadNode (iQuadPt, 1),
                                                         QRMapper.quadNode (iQuadPt, 2),
                                                         std::abs (QRMapper.wDet (iQuadPt) ) ) );
            }

            // Build the fifth tetra
            std::vector< VectorSmall<3> > tetra5 (4);
            tetra5[0] = refPts[localNegative[0]];
            tetra5[1] = intersections[1];
            tetra5[2] = intersections[2];
            tetra5[3] = intersections[3];

            // Map the QR
            QRMapper.update ( tetra5, ET_UPDATE_QUAD_NODE | ET_UPDATE_WDET);

            for (UInt iQuadPt (0); iQuadPt < QuadratureSize; ++iQuadPt)
            {
                M_adaptedQR->addPoint (QuadraturePoint ( QRMapper.quadNode (iQuadPt, 0),
                                                         QRMapper.quadNode (iQuadPt, 1),
                                                         QRMapper.quadNode (iQuadPt, 2),
                                                         std::abs (QRMapper.wDet (iQuadPt) ) ) );
            }

            // Build the sixth tetra
            std::vector< VectorSmall<3> > tetra6 (4);
            tetra6[0] = refPts[localNegative[0]];
            tetra6[1] = refPts[localNegative[1]];
            tetra6[2] = intersections[2];
            tetra6[3] = intersections[3];

            // Map the QR
            QRMapper.update ( tetra6, ET_UPDATE_QUAD_NODE | ET_UPDATE_WDET);

            for (UInt iQuadPt (0); iQuadPt < QuadratureSize; ++iQuadPt)
            {
                M_adaptedQR->addPoint (QuadraturePoint ( QRMapper.quadNode (iQuadPt, 0),
                                                         QRMapper.quadNode (iQuadPt, 1),
                                                         QRMapper.quadNode (iQuadPt, 2),
                                                         std::abs (QRMapper.wDet (iQuadPt) ) ) );
            }



        }
        else
        {
            ASSERT ( localPositive.size() == 3, "Internal inconsistency");

            // Localize the intersections
            std::vector< VectorSmall<3> > intersections (3);

            for (UInt i (0); i < 3; ++i)
            {
                UInt pos (localPositive[i]);
                UInt neg (localNegative[0]);
                Real mu ( localLSValue[pos] / (localLSValue[pos] - localLSValue[neg]) );
                VectorSmall<3> I ( refPts[pos] + mu * ( refPts[neg] - refPts[pos] ) );
                intersections[i] = I;
            }


            // Build the first tetra
            std::vector< VectorSmall<3> > tetra1 (4);
            tetra1[0] = refPts[localNegative[0]];
            tetra1[1] = intersections[0];
            tetra1[2] = intersections[1];
            tetra1[3] = intersections[2];

            // Map the QR
            QRMapper.update ( tetra1, ET_UPDATE_QUAD_NODE | ET_UPDATE_WDET);

            for (UInt iQuadPt (0); iQuadPt < QuadratureSize; ++iQuadPt)
            {
                M_adaptedQR->addPoint (QuadraturePoint ( QRMapper.quadNode (iQuadPt, 0),
                                                         QRMapper.quadNode (iQuadPt, 1),
                                                         QRMapper.quadNode (iQuadPt, 2),
                                                         std::abs (QRMapper.wDet (iQuadPt) ) ) );
            }


            // Build the second tetra
            std::vector< VectorSmall<3> > tetra2 (4);
            tetra2[0] = intersections[0];
            tetra2[1] = intersections[1];
            tetra2[2] = intersections[2];
            tetra2[3] = refPts[localPositive[0]];

            // Map the QR
            QRMapper.update ( tetra2, ET_UPDATE_QUAD_NODE | ET_UPDATE_WDET);

            for (UInt iQuadPt (0); iQuadPt < QuadratureSize; ++iQuadPt)
            {
                M_adaptedQR->addPoint (QuadraturePoint ( QRMapper.quadNode (iQuadPt, 0),
                                                         QRMapper.quadNode (iQuadPt, 1),
                                                         QRMapper.quadNode (iQuadPt, 2),
                                                         std::abs (QRMapper.wDet (iQuadPt) ) ) );
            }

            // Build the third tetra
            std::vector< VectorSmall<3> > tetra3 (4);
            tetra3[0] = refPts[localPositive[0]];
            tetra3[1] = refPts[localPositive[1]];
            tetra3[2] = refPts[localPositive[2]];
            tetra3[3] = intersections[2];

            // Map the QR
            QRMapper.update ( tetra3, ET_UPDATE_QUAD_NODE | ET_UPDATE_WDET);

            for (UInt iQuadPt (0); iQuadPt < QuadratureSize; ++iQuadPt)
            {
                M_adaptedQR->addPoint (QuadraturePoint ( QRMapper.quadNode (iQuadPt, 0),
                                                         QRMapper.quadNode (iQuadPt, 1),
                                                         QRMapper.quadNode (iQuadPt, 2),
                                                         std::abs (QRMapper.wDet (iQuadPt) ) ) );
            }

            // Build the fourth tetra
            std::vector< VectorSmall<3> > tetra4 (4);
            tetra4[0] = refPts[localPositive[0]];
            tetra4[1] = refPts[localPositive[1]];
            tetra4[2] = intersections[1];
            tetra4[3] = intersections[2];

            // Map the QR
            QRMapper.update ( tetra4, ET_UPDATE_QUAD_NODE | ET_UPDATE_WDET);

            for (UInt iQuadPt (0); iQuadPt < QuadratureSize; ++iQuadPt)
            {
                M_adaptedQR->addPoint (QuadraturePoint ( QRMapper.quadNode (iQuadPt, 0),
                                                         QRMapper.quadNode (iQuadPt, 1),
                                                         QRMapper.quadNode (iQuadPt, 2),
                                                         std::abs (QRMapper.wDet (iQuadPt) ) ) );
            }
        }
    }
}


} // Namespace LifeV

#endif /* LEVELSETQRADAPTER_H */

#endif /* LIFEV_HAS_ETA */
