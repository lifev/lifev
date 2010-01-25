/* -*- mode: c++ -*-
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
#ifndef _CURRENTFE_H
#define _CURRENTFE_H

#include <life/lifecore/life.hpp>
#include <life/lifefem/geoMap.hpp>
#include <life/lifefem/refFE.hpp>
//#include <life/lifefem/geoMap.hpp>
/*!
  \file currentFE.h
  \brief Structure for the current finite element
*/

namespace LifeV
{
/*!
  \class CurrentFE
  \brief The class for a finite element
  \author J.-F. Gerbeau
  \date 04/2002

     I also factorized some code, and
     postponed the implementation of template methods
     after the class declaration.
  \author Vincent Martin
  \date 07/2004
*/

class CurrentFE
{
private:
    //! compute only the jacobian matrix
    void _comp_jacobian();
    //! call _comp_jacobian() and compute its determinant
    void _comp_jacobian_and_det();
    //! call _comp_jacobian() and compute its inverse and its determinant
    void _comp_inv_jacobian_and_det();
    void _comp_quad_point_coor();

    //! update the definition of the geo points  (VM 07/2004)
    //! contains a switch over nbCoor (= dimension of pb): slow???
    template <class GEOELE>
    void _update_point( const GEOELE& geoele );

    //! compute phiDer
    void _comp_phiDer();
    //! compute the second derivative phiDer2
    void _comp_phiDer2();
    //! compute phiDer and phiDer2
    void _comp_phiDerDer2();


    UInt _currentId;
    UInt _currentLocalId;



#ifdef TEST_PRE

    bool _hasJac;
    bool _hasFirstDeriv;
    bool _hasSecondDeriv;
    bool _hasQuadPtCoor;
#endif
public:
    CurrentFE( const RefFE& _refFE, const GeoMap& _geoMap, const QuadRule& _qr );
private:
    CurrentFE( );
    CurrentFE( const CurrentFE& );
public:
    const int nbGeoNode;
    const int nbNode;
    const int nbCoor;
    const int nbQuadPt;
    const int nbDiag;
    const int nbUpper;
    const int nbPattern;
    KNM<Real> point;
    const RefFE& refFE;
    const GeoMap& geoMap;
    const QuadRule& qr;
    KNM<Real> phi;
    KNMK<Real> dPhiRef;
    KNMKL<Real> dPhiRef2;
    //
    KNMK<Real> phiDer;
    KNMK<Real> jacobian;
    KNMK<Real> tInvJac;
    KNM<Real> phiGeo;
    KNMK<Real> dPhiGeo;
    KN<Real> weightDet;
    KN<Real> detJac;
    KNM<Real> quadPt;  //!< Coordinates of the quadrature points on the current element
    //second derivatives not yet done (four dimensions KNMKL ?)
    KNMKL<Real> phiDer2;
    //--------------------------------------------------------------------------
    /*!
      return the diameter of the element in the 1-norm
    */
    Real diameter() const;
    /*!
      return the diameter of the element in the 2-norm
    */
    Real diameter2() const;
    /*!
      return the id of the current element (updated with the update* functions)
     */
    inline UInt currentId() const
    {
        return _currentId;
    }

    inline UInt currentLocalId() const
    {
        return _currentLocalId;
    }
#ifdef TEST_PRE
    /*!
      return true if the determinant has been updated
      (can ONLY be used if TEST_PRE is defined at compilation time)
     */
    inline bool hasJac() const
    {
        return _hasJac;
    }
    /*!
      return true if the first derivatives have been updated
      (can ONLY be used if TEST_PRE is defined at compilation time)
     */
    inline bool hasFirstDeriv() const
    {
        return _hasFirstDeriv;
    }
    /*!
      return true if the second derivatives have been updated
      (can ONLY be used if TEST_PRE is defined at compilation time)
     */
    inline bool hasSecondDeriv() const
    {
        return _hasSecondDeriv;
    }
    /*!
      return true if the coordinate of the quadrature points have been updated
      (can ONLY be used if TEST_PRE is defined at compilation time)
     */
    inline bool hasQuadPtCoor() const
    {
        return _hasQuadPtCoor;
    }
#endif
    /*!
      compute the coordinate (x,y,z)= F(xi,eta,zeta), (F: geo mappping)
      where (xi,eta,zeta) are the coor in the ref element
      and   (x,y,z) the coor in the current element
      (if the code is compiled in 2D mode then z=0 and zeta is disgarded)
    */
    void coorMap( Real& x, Real& y, Real& z,
                  const Real & xi, const Real & eta, const Real &
                  zeta ) const;
    /*!
      compute the coordinate (xi,eta,zeta)=inv(F)(x,y,z)
    */
    void coorBackMap(const Real& x, const Real& y, const Real& z,
		     Real & xi, Real & eta, Real& zeta) const;

    /*!
      compute the jacobian at a given point : d x_compx / d zeta_compzeta
    */
    Real pointJacobian(const Real& hat_x, const Real& hat_y, const Real& hat_z,
		  int compx, int compzeta) const;

    /*!
      compute the inverse jacobian
     */

    Real pointInverseJacobian(const Real& hat_x, const Real& hat_y, const Real& hat_z,
		  int compx, int compzeta) const;

    /*!
      return the barycenter of the element
     */
    void barycenter( Real& x, Real& y, Real& z );
    /*!  return (x,y,z) = the global coordinates of the quadrature point ig
      in the current element. \warning this function is almost obsolete since if
      you call the function updateFirstDerivQuadPt rather than updateFirstDeriv
      (for example), the coordinates of the quadrature points have already
      been computed and may be obtained via quadPt(ig,icoor). This is usually
      much less expensive since it avoids many calls to coorQuadPt
    */
    inline void coorQuadPt( Real& x, Real& y, Real& z, const int ig ) const
    {
        coorMap( x, y, z, qr.quadPointCoor( ig, 0 ), qr.quadPointCoor( ig, 1 ),
                 qr.quadPointCoor( ig, 2 ) );
    }
    //!  patternFirst(i): row index in the element matrix of the i-th term of the pattern
    inline int patternFirst( int i ) const
    {
        return refFE.patternFirst( i );
    }
    //! patternSecond(i): column index in the element matrix of the i-th term of the pattern
    inline int patternSecond( int i ) const
    {
        return refFE.patternSecond( i );
    }

    /*!
      Return the measure of the current element
    */
    Real measure() const;

    //---------------------------------------
    //! DECLARATION of the update methods
    //---------------------------------------
    /*!
      minimal update: we just identify the id of the current element
    */
    template <class GEOELE>
    void update( const GEOELE& geoele );
    /*!
      compute the arrays detJac, weightDet, jacobian on
      the current element
    */
    template <class GEOELE>
    void updateJac( const GEOELE& geoele );
    /*!
      compute the arrays detJac, weightDet, jacobian and quadPt
      on the current element
    */
    template <class GEOELE>
    void updateJacQuadPt( const GEOELE& geoele );
    /*!
      compute the arrays detJac, weightDet, jacobian,
      tInvJac, phiDer on the current element
    */
    template <class GEOELE>
    void updateFirstDeriv( const GEOELE& geoele );
    /*!
      compute the arrays detJac, weightDet, jacobian,
      tInvJac, phiDer and quadPt on the current element
    */
    template <class GEOELE>
    void updateFirstDerivQuadPt( const GEOELE& geoele );
    /*!
      compute the arrays detJac, weightDet, jacobian,
      tInvJac, phiDer2 on the current element
    */
    template <class GEOELE>
    void updateSecondDeriv( const GEOELE& geoele );
    /*!
      compute the arrays detJac, weightDet, jacobian,
      tInvJac, phiDer2 on the current element
    */
    template <class GEOELE>
    void updateSecondDerivQuadPt( const GEOELE& geoele );
    /*!
      compute the arrays detJac, weightDet, jacobian,
      tInvJac, phiDer, phiDer2 on the current element
    */
    template <class GEOELE>
    void updateFirstSecondDeriv( const GEOELE& geoele );
    /*!
      compute the arrays detJac, weightDet, jacobian,
      tInvJac, phiDer, phiDer2 on the current element
    */
    template <class GEOELE>
    void updateFirstSecondDerivQuadPt( const GEOELE& geoele );

};


//----------------------------------------------------------------------
//! IMPLEMENTATION
//----------------------------------------------------------------------
//! update the definition of the geo points (VM 07/04)
template <class GEOELE>
void CurrentFE::_update_point( const GEOELE& geoele )
{

    ASSERT( nbCoor < 4, "nbCoor must be smaller than 4");

    //! Nicer way to do the same as in the switch
    Real const * pointCoordinate;
    for ( ID i(0); i < nbGeoNode; ++i )
    {
        pointCoordinate = geoele.point( i + 1  ).coor();
        for( ID icoor(0); icoor < nbCoor; ++icoor)
            point( i, icoor ) =  pointCoordinate[icoor];
    }
    return;
    /*
    switch ( nbCoor )
    {
    case 1:    //! 1D
        for ( int i = 0;i < nbGeoNode;i++ )
        {
            point( i, 0 ) = geoele.point( i + 1 ).x();
        }
        break;
    case 2:    //! 2D
        for ( int i = 0;i < nbGeoNode;i++ )
        {
            point( i, 0 ) = geoele.point( i + 1 ).x();
            point( i, 1 ) = geoele.point( i + 1 ).y();
        }
        break;
    case 3:    //! 3D
        for ( int i = 0;i < nbGeoNode;i++ )
        {
            point( i, 0 ) = geoele.point( i + 1 ).x();
            point( i, 1 ) = geoele.point( i + 1 ).y();
            point( i, 2 ) = geoele.point( i + 1 ).z();
        }
        break;
    default:
        std::ostringstream err_msg;
        err_msg << "Dimension " << nbCoor << " not available.";

        ERROR_MSG( err_msg.str().c_str() );
    }
    */
}

//---------------------------------------
//! IMPLEMENTATION of the CurrentFE::update methods
//---------------------------------------

/*!
    minimal update: we just identify the id of the current element
  */
template <class GEOELE>
void CurrentFE::update( const GEOELE& geoele )
{
#ifdef TEST_PRE
    _hasJac = false;
    _hasFirstDeriv = false;
    _hasSecondDeriv = false;
    _hasQuadPtCoor = false;
#endif

    _currentId      = geoele.id();
    _currentLocalId = geoele.localId();
}

/*!
    compute the arrays detJac, weightDet, jacobian on
    the current element
*/
template <class GEOELE>
void CurrentFE::updateJac( const GEOELE& geoele )
{
#ifdef TEST_PRE
    _hasJac = true;
    _hasFirstDeriv = false;
    _hasSecondDeriv = false;
    _hasQuadPtCoor = false;
#endif

    _currentId      = geoele.id();
    _currentLocalId = geoele.localId();
    //! update the definition of the geo points
    _update_point( geoele );
    //! compute the jacobian and its determinant...
    _comp_jacobian_and_det();
}

/*!
    compute the arrays detJac, weightDet, jacobian and quadPt
    on the current element
*/
template <class GEOELE>
void CurrentFE::updateJacQuadPt( const GEOELE& geoele )
{
#ifdef TEST_PRE
    _hasJac = true;
    _hasFirstDeriv = false;
    _hasSecondDeriv = false;
    _hasQuadPtCoor = true;
#endif

    _currentId      = geoele.id();
    _currentLocalId = geoele.localId();
    //! update the definition of the geo points
    _update_point( geoele );
    //! compute the jacobian and its determinant...
    _comp_jacobian_and_det();
    //! and the coordinates of the quadrature points
    _comp_quad_point_coor();
}

/*!
    compute the arrays detJac, weightDet, jacobian,
    tInvJac, phiDer on the current element
*/
template <class GEOELE>
void CurrentFE::updateFirstDeriv( const GEOELE& geoele )
{
#ifdef TEST_PRE
    _hasJac = true;
    _hasFirstDeriv = true;
    _hasSecondDeriv = false;
    _hasQuadPtCoor = false;
#endif

    _currentId      = geoele.id();
    _currentLocalId = geoele.localId();
    //! update the definition of the geo points
    _update_point( geoele );
    //! compute the inverse jacobian...
    _comp_inv_jacobian_and_det();
    //! product InvJac by dPhiRef to compute phiDer
    _comp_phiDer();
}

/*!
    compute the arrays detJac, weightDet, jacobian,
    tInvJac, phiDer and quadPt on the current element
*/
template <class GEOELE>
void CurrentFE::updateFirstDerivQuadPt( const GEOELE& geoele )
{
#ifdef TEST_PRE
    _hasJac = true;
    _hasFirstDeriv = true;
    _hasSecondDeriv = false;
    _hasQuadPtCoor = true;
#endif

    _currentId      = geoele.id();
    _currentLocalId = geoele.localId();
    //! update the definition of the geo points
    _update_point( geoele );
    //! compute the inverse jacobian...
    _comp_inv_jacobian_and_det();
    //! product InvJac by dPhiRef to compute phiDer
    _comp_phiDer();
    //! and the coordinates of the quadrature points
    _comp_quad_point_coor();
}

// A. Veneziani, October 30, 2002
/*!
  compute the arrays detJac, weightDet, jacobian,
  tInvJac, phiDer2 on the current element
*/
template <class GEOELE>
void CurrentFE::updateSecondDeriv( const GEOELE& geoele )
{
#ifdef TEST_PRE
    _hasJac = true;
    _hasFirstDeriv = false;
    _hasSecondDeriv = true;
    _hasQuadPtCoor = false;
#endif

    _currentId      = geoele.id();
    _currentLocalId = geoele.localId();
    //! update the definition of the geo points
    _update_point( geoele );
    //! compute the inverse jacobian...
    _comp_inv_jacobian_and_det();
    //! compute the second derivative
    _comp_phiDer2();
}

/*!
    compute the arrays detJac, weightDet, jacobian,
    tInvJac, phiDer2 on the current element
  */
template <class GEOELE>
void CurrentFE::updateSecondDerivQuadPt( const GEOELE& geoele )
{
#ifdef TEST_PRE
    _hasJac = true;
    _hasFirstDeriv = false;
    _hasSecondDeriv = true;
    _hasQuadPtCoor = true;
#endif

    _currentId      = geoele.id();
    _currentLocalId = geoele.localId();
    //! update the definition of the geo points
    _update_point( geoele );
    //! compute the inverse jacobian...
    _comp_inv_jacobian_and_det();
    //! compute the second derivative
    _comp_phiDer2();
    //! and the coordinates of the quadrature points
    _comp_quad_point_coor();
}
/*!
    compute the arrays detJac, weightDet, jacobian,
    tInvJac, phiDer, phiDer2 on the current element
  */
template <class GEOELE>
void CurrentFE::updateFirstSecondDeriv( const GEOELE& geoele )
{
#ifdef TEST_PRE
    _hasJac = true;
    _hasFirstDeriv = true;
    _hasSecondDeriv = true;
    _hasQuadPtCoor = false;
#endif

    _currentId      = geoele.id();
    _currentLocalId = geoele.localId();
    //! update the definition of the geo points
    _update_point( geoele );
    //! compute the inverse jacobian...
    _comp_inv_jacobian_and_det();
    //! compute phiDer and phiDer2
    _comp_phiDerDer2();
}
/*!
    compute the arrays detJac, weightDet, jacobian,
    tInvJac, phiDer, phiDer2 on the current element
  */
template <class GEOELE>
void CurrentFE::updateFirstSecondDerivQuadPt( const GEOELE& geoele )
{
#ifdef TEST_PRE
    _hasJac = true;
    _hasFirstDeriv = true;
    _hasSecondDeriv = true;
    _hasQuadPtCoor = true;
#endif

    _currentId      = geoele.id();
    _currentLocalId = geoele.localId();
    //! update the definition of the geo points
    _update_point( geoele );
    //! compute the inverse jacobian...
    _comp_inv_jacobian_and_det();
    //! compute phiDer and phiDer2
    _comp_phiDerDer2();
    //! and the coordinates of the quadrature points
    _comp_quad_point_coor();
}
}

#endif
