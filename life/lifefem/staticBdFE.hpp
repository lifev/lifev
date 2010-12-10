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

#ifndef _StaticBDFE_H
#define _StaticBDFE_H

#include <vector>


#include <life/lifecore/life.hpp>
#include <life/lifefem/geoMap.hpp>
#include <life/lifefem/refFE.hpp>
#include <life/lifefem/geoMap.hpp>
#include <life/lifearray/tab.hpp>
#include <life/lifefem/quadRule.hpp>

namespace LifeV
{
/*!
  \class StaticBdFE
  \brief A class for static boundary finite element
  \author J.-F. Gerbeau & V. Martin
  \date 09/2002

  This class has two purposes:
  \par
  (1) it is a base class for standard boundary element (see CurrentBdFE.h)
  \par
  (2) it is used by refHybridFE as static boundary for a reference element

*/

class StaticBdFE
{

public:

    //! @name Constructor & Destructor
    //@{

    StaticBdFE( const RefFE& refFE, const GeoMap& geoMap );

    StaticBdFE( const RefFE& refFE, const GeoMap& geoMap, const QuadRule& qr );

    //! new optionnal argument for hybrid hdiv fe invArea
    StaticBdFE( const RefFE& refFE, const GeoMap& geoMap, const QuadRule& qr,
                const Real* refcoor, UInt currentid, Real invarea = 1. );

    virtual ~StaticBdFE();

    //@}


    //! @name Methods
    //@{

    /*!
      return the coordinate (x,y,z)= F(xi,eta), (F: geo mappping)
      where (xi,eta) are the coor in the ref element
      and   (x,y,z) the coor in the current element
      (if the code is compiled in 2D mode then z=0 and eta is disgarded)
    */
    void coorMap( Real& x, Real& y, Real& z,
                  const Real & xi, const Real & eta ) const;
    /*!
      return (x,y,z) = the global coordinates of the quadrature point ig
      in the current element. \warning this function is almost obsolete since if
      you call the function updateMeasQuadPt rather than updateMeas
      (for example), the coordinates of the quadrature points have already
      been computed and may be obtained via quadPt(ig,icoor). This is usually
      much less expensive since it avoids many calls to coorQuadPt
    */
    void coorQuadPt( Real& x, Real& y, Real& z, const int ig ) const
    {
        ASSERT_PRE( M_hasQR, "Needs a quadrature rule" )
        coorMap( x, y, z, qr.quadPointCoor( ig, 0 ), qr.quadPointCoor( ig, 1 ) );
    }
    /*!
      Return the measure of the current element
      \warning either updateMeas(...) or updateMeasNormal(...)
      must have been called before
    */
    Real measure() const;
    /*!
      compute the integral over the current boundary element
      \warning either updateMeas(...) or updateMeasNormal(...)
      must have been called before
    */
    template <typename FunctorType>
    Real integral( const FunctorType & f ) const;

    /*!
      compute the integral of f . n over the current boundary element
      \warning  updateMeasNormal(...) must have been called before
    */
    template <typename FunctorType>
    Real integral_n( const FunctorType & f ) const;

    //@}


    //! @name Get Methods
    //@{

#ifdef TEST_PRE
    /*!
      return true if a quadrature rule has been given
      (can ONLY be used if TEST_PRE is defined at compilation time)
    */
    bool hasQR() const
    {
        return M_hasQR;
    }
    /*!
      return true if the measure has been updated
      (can ONLY be used if TEST_PRE is defined at compilation time)
    */
    bool hasMeas() const
    {
        return M_hasMeasure;
    }
    /*!
      return true if the first derivatives have been updated
      (can ONLY be used if TEST_PRE is defined at compilation time)
    */
    bool hasFirstDeriv() const
    {
        return M_hasFirstDerivative;
    }
    /*!
      return true if the tangents have been updated
      (can ONLY be used if TEST_PRE is defined at compilation time)
    */
    bool hasTangent() const
    {
        return M_hasTangent;
    }
    /*!
      return true if the normal has been updated
      (can ONLY be used if TEST_PRE is defined at compilation time)
    */
    bool hasNormal() const
    {
        return M_hasNormal;
    }
    /*!
      return true if the coordinate of the quadrature points have been updated
      (can ONLY be used if TEST_PRE is defined at compilation time)
    */
    bool hasQuadPtCoor() const
    {
        return M_hasQuadPtCoor;
    }
#endif
    /*!
      return the id of the current element (updated with the update* functions)
    */
    UInt currentId() const
    {
        return M_currentID;
    }

    //@}


    //! Number of geometrical nodes
    const UInt      nbGeoNode;

    //! Number of finite element node
    const UInt       nbNode;

    //! Number of coordinates
    const UInt       nbCoor;

    //! Number of quadrature points
    const UInt       nbQuadPt;

    //! The point that define the geometry
    KNM<Real>       point;

    //! The reference finite element
    const RefFE&    refFE;

    //! The geometical mapping
    const GeoMap&   geoMap;

    //! The quadrature rule
    const QuadRule& qr;

    //! Values of the basis functions on quadrature points
    KNM<Real>  phi;

    //! Values of the derivatives of the basis functions on quadrature points on the reference finite element
    KNMK<Real> dPhiRef;

    //! Values of the derivatives of the basis functions on quadrature points on the current finite element
    KNMK<Real> dPhi;

    //! Values of the geometric basis functions on quadrature points
    KNM<Real>  phiGeo;

    //! Values of the derivatives of the geometric basis functions on quadrature points
    KNMK<Real> dPhiGeo;

    //! Values of the weight times the measure on the quadrature points
    KN<Real>   weightMeas;

    //! Values of the measures on the quadrature points
    KN<Real>   meas;

    //! Values of the normal on the quadrature points
    KNM<Real>  normal;

    //! Values of the tangents on the quadrature points
    KNMK<Real> tangent;

    //! Metric tensor on the quadrature points
    KNMK<Real> metric;

    //! Coordinates of the quadrature points on the current element
    KNM<Real>  quadPt;

    //! Inverse of the area
    const Real invArea;


protected:
    void computeMeasure();
    void computeMeasureNormal();
    void computeQuadPointCoordinate();
    UInt M_currentID;

#ifdef TEST_PRE
    bool M_hasQR;
    bool M_hasMeasure;
    bool M_hasFirstDerivative;
    bool M_hasTangent;
    bool M_hasNormal;
    bool M_hasQuadPtCoor;
#endif

};


template <typename FunctorType>
Real
StaticBdFE::
integral( const FunctorType & f ) const
{
    ASSERT_PRE( M_hasMeasure, "integral needs measure. Call an update function" );
    Real integ( 0.0 );
    Real x, y, z;
    for ( UInt iQuadPt(0); iQuadPt < nbQuadPt; ++iQuadPt )
    {
        coorQuadPt( x, y, z, iQuadPt );
        integ += f( x, y, z ) * weightMeas( iQuadPt );
    }
    return integ;
}

template <typename FunctorType>
Real
StaticBdFE::
integral_n( const FunctorType & f ) const
{
    ASSERT_PRE( M_hasNormal, "integral_n needs measure and normal. Call the appropriate update function" );
    Real integ( 0.0 );
    Real x, y, z;
    std::vector<Real> ret( nbCoor + 1 );
    Real tmp;
    for ( UInt iQuadPt(0); iQuadPt < nbQuadPt; ++iQuadPt )
    {
        coorQuadPt( x, y, z, iQuadPt );
        f( x, y, z, &ret.front() );
        tmp = 0;
        for ( UInt d(0); d <= nbCoor; ++d )
        {
            tmp += ret[ d ] * normal( d, iQuadPt );
        }
        integ += tmp * weightMeas( iQuadPt );
    }
    return integ;
}


}
#endif
