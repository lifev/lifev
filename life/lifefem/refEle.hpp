/*-*- mode: c++ -*-
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
#ifndef _REFELE_H_INCLUDE
#define _REFELE_H_INCLUDE

#include "lifeV.hpp"
#include "tab.hpp"
#include "basisElSh.hpp"
#include "quadRule.hpp"

/*!
  \file refEle.h
  \brief Base class for RefFE and GeoMap
*/

namespace LifeV
{
/*!
  \class RefEle
  \brief Base class for RefGeo and RefFE
  \author J.-F. Gerbeau
  \date 04/2002

  It contains the basis functions and their values on quadrature points.
  These functions will be used either by RefFE (finite element) or by
  GeoMap (geometrical mapping).
*/



//! Indicates the LOCAL (reference element) coordinates
typedef const Real & cRRef;
typedef Real (* Fct)(cRRef,cRRef ,cRRef);

class RefEle{
protected:
  const SetOfQuadRule* _sqr; //!< pointer on the set of quadrature rules
  const Fct*  _phi; //!< pointer on the basis functions
  const Fct*  _dPhi;//!< pointer on the derivatives of the basis functions
  const Fct*  _d2Phi;//!< pointer on the second derivatives of the basis functions
  const Real* _refCoor;//!< reference coordinates. Order: xi_1,eta_1,zeta_1,xi_2,eta_2,zeta_2,...
public:
  const string name; //!< name of the reference element
  const ReferenceShapes shape; //!< geometrical shape of the element
  const int nbDof;   //!< Total number of degrees of freedom
  const int nbCoor;  //!< Number of local coordinates
private:
  //!   values of the basis functions on all quadrature points
  KN<Real> _phiQuad;
  //! values of the derivatives of the basis functions on all quadrature points
  KN<Real> _dPhiQuad;
 //! values of the second derivatives of the basis functions on all quadrature points
  KN<Real> _d2PhiQuad;
  //! values of the second derivatives of the basis functions on all quadrature points
  KN<int> _idxQuad;//!< _idxQuad[t] = index of the quadrature rules of id t in _phiQuad
  KN<int> _idxDQuad;//!< _idxDQuad[t] = index of the quadrature rules of id t in _dPhiQuad
  KN<int> _idxD2Quad;//!< _idxD2Quad[t] = index of the quadrature rules of id t in _d2PhiQuad
public:
  //! constructor
  RefEle(string _name,ReferenceShapes _shape,int _nbDof,int _nbCoor,const Fct* phi,const Fct* dPhi,
	 const Fct* d2Phi,const Real* _refCoor,const SetOfQuadRule& sqr);
  ~RefEle();
  //! return the first local coordinate of the i-th node of the reference element
  inline Real xi(int i) const {
    ASSERT_BD(i<nbDof)
      return _refCoor[3*i];
  }
  //! return the second local coordinate of the i-th node of the reference element
  inline Real eta(int i) const {
    ASSERT_BD(i<nbDof)
      return _refCoor[3*i+1];
  }
  //! return the third local coordinate of the i-th node of the reference element
  inline Real zeta(int i) const {
    ASSERT_BD(i<nbDof)
      return _refCoor[3*i+2];
  }
  //! return the icoor-th local coordinate of the i-th node of the reference element
  inline Real refCoor(int i,int icoor) const {
    ASSERT_BD(i<nbDof && icoor < nbCoor)
      return _refCoor[3*i+icoor];
  }
  //! return the value of the i-th basis function on point (x,y,z)
  inline double phi(int i,cRRef x,cRRef y,cRRef z) const{
    ASSERT_BD( i<nbDof)
      return _phi[i](x,y,z);
  }
  //!return the value of the i-th basis function on the ig-th point of integration of quadrature rule qr
  inline Real phi(int i,int ig,const QuadRule& qr)const{
    ASSERT_BD(i < nbDof && ig < qr.nbQuadPt)
      return _phiQuad(_idxQuad(qr.id)+ig*nbDof+i);
  }
  //! return the array of the values of the basis functions phi_1[ig],phi_2[ig],... on the integration point ig of the quadrature rule qr.
  /*!This function is written in order to avoid too many calls to function
    phi(i,ig,qr), please check if it improves really efficiency. If yes,
    analogous function should be written for the derivatives.  */
  inline RN_ phiQuadPt(int ig,const QuadRule& qr)const{
    ASSERT_BD(ig < qr.nbQuadPt)
      return _phiQuad(SubArray(nbDof,_idxQuad(qr.id)+ig*nbDof));
  }
  //! return the array of the values of the i-th basis functions phi_i[ig=1],phi_i[ig=2],...on all the integration points of the quadrature rule qr.
  /*! This function is written in order to avoid too many calls to
    function phi(i,ig,qr), please check if it improves really
    efficiency. If yes, analogous function should be written for the
    derivatives.  */
  inline RN_ phiI(int i,const QuadRule& qr)const{
    ASSERT_BD(i < nbDof)
      return _phiQuad(SubArray(qr.nbQuadPt,_idxQuad(qr.id)+i,nbDof));
  }
  //! return the value of the icoor-th derivative of the i-th basis function on point (x,y,z)
  inline double dPhi(int i,int icoor,cRRef x,cRRef y,cRRef z) const{
    ASSERT_BD(i<nbDof && icoor<nbCoor)
      return  _dPhi[i*nbCoor+icoor](x,y,z);
  }
  //! return the value of icoor-th derivative of the i-th basis function on the ig-th point of integration of quadrature rule qr
  inline Real dPhi(int i,int icoor,int ig,const QuadRule& qr)const{
    ASSERT_BD(i < nbDof && ig < qr.nbQuadPt && icoor < nbCoor)
      return _dPhiQuad( _idxDQuad(qr.id) + (ig*nbDof+ i) * nbCoor + icoor );
  }
  //!  return the value of the (icoor,jcoor)-th second derivative of the i-th basis function on point (x,y,z)
  inline double d2Phi(int i,int icoor,int jcoor,cRRef x,cRRef y,cRRef z) const{
    ASSERT_BD(i<nbDof && icoor <nbCoor && jcoor < nbCoor)
      return  _d2Phi[(i*nbCoor+icoor)*nbCoor + jcoor](x,y,z);
  }
  //! return the value of the (icoor,jcoor)-th second derivative of the i-th basis function on the ig-th point of integration of quadrature rule qr
  inline Real d2Phi(int i,int icoor,int jcoor,int ig,const QuadRule& qr)const{
    ASSERT_BD(i < nbDof && ig < qr.nbQuadPt && icoor < nbCoor)
      return _d2PhiQuad( _idxD2Quad(qr.id) + ((ig*nbDof+ i) * nbCoor + icoor)*nbCoor + jcoor );
  }
  void check() const;//!< A simple check function
  friend std::ostream& operator << (std::ostream& f,const RefEle& fe);
};


//======================================================================
//
//                            P1  (1D)
//
//======================================================================
/*
                           1-----2
*/
Real fct1_P1_1D(cRRef x,cRRef,cRRef );
Real fct2_P1_1D(cRRef x,cRRef,cRRef );

Real derfct1_1_P1_1D(cRRef,cRRef,cRRef );
Real derfct2_1_P1_1D(cRRef,cRRef,cRRef );

Real der2fct1_P1_1D(cRRef,cRRef,cRRef);

static const Real refcoor_P1_1D[6] = {0.  ,0.  ,0.,
				      1.  ,0.  ,0.};

static const Fct fct_P1_1D[2] = {fct1_P1_1D,fct2_P1_1D};
static const Fct derfct_P1_1D[2] = {derfct1_1_P1_1D,derfct2_1_P1_1D};
static const Fct der2fct_P1_1D[2] = {der2fct1_P1_1D,der2fct1_P1_1D};

//======================================================================
//
//                            P2  (1D)
//
//======================================================================
/*
                           1--3--2
*/
Real fct1_P2_1D(cRRef x,cRRef,cRRef );
Real fct2_P2_1D(cRRef x,cRRef,cRRef );
Real fct3_P2_1D(cRRef x,cRRef,cRRef );

Real derfct1_1_P2_1D(cRRef x,cRRef  ,cRRef );
Real derfct2_1_P2_1D(cRRef x,cRRef  ,cRRef );
Real derfct3_1_P2_1D(cRRef x,cRRef  ,cRRef );

Real der2fct1_11_P2_1D(cRRef x,cRRef  ,cRRef );
Real der2fct2_11_P2_1D(cRRef x,cRRef  ,cRRef );
Real der2fct3_11_P2_1D(cRRef x,cRRef  ,cRRef );

static const Real refcoor_P2_1D[9] = {0.  ,0.  ,0.,
				      1.  ,0.  ,0.,
				      0.5 ,0.  ,0.};
static const Fct fct_P2_1D[3] = {fct1_P2_1D,fct2_P2_1D,fct3_P2_1D};
static const Fct derfct_P2_1D[3] = {derfct1_1_P2_1D,derfct2_1_P2_1D,derfct3_1_P2_1D};
static const Fct der2fct_P2_1D[3] = {der2fct1_11_P2_1D,der2fct2_11_P2_1D,der2fct3_11_P2_1D};

//======================================================================
//
//                            P0  (2D)
//
//======================================================================
/*

                           |\
                           | \
                           | 1\
                            ---
*/
Real fct1_P0_2D(cRRef ,cRRef ,cRRef );
// First and Second derivatives are both equal (to 0).
Real derfct1_P0_2D(cRRef,cRRef,cRRef );
Real der2fct1_P0_2D(cRRef,cRRef,cRRef);

static const Real refcoor_P0_2D[3] = {1./3.  ,1./3.  ,0.};  // check this : gravity center??

static const Fct fct_P0_2D[1] = {fct1_P0_2D};

static const Fct derfct_P0_2D[2] = {derfct1_P0_2D, derfct1_P0_2D};

static const Fct der2fct_P0_2D[4] = {der2fct1_P0_2D, der2fct1_P0_2D,
				     der2fct1_P0_2D, der2fct1_P0_2D };

//======================================================================
//
//                            P1  (2D)
//
//======================================================================
/*
                           3
                           |\
                           | \
                           |  \
                           1---2
*/
Real fct1_P1_2D(cRRef x,cRRef y,cRRef );
Real fct2_P1_2D(cRRef x,cRRef  ,cRRef );
Real fct3_P1_2D(cRRef  ,cRRef y,cRRef );

Real derfct1_1_P1_2D(cRRef,cRRef,cRRef );
Real derfct1_2_P1_2D(cRRef,cRRef,cRRef );
Real derfct2_1_P1_2D(cRRef,cRRef,cRRef );
Real derfct2_2_P1_2D(cRRef,cRRef,cRRef );
Real derfct3_1_P1_2D(cRRef,cRRef,cRRef );
Real derfct3_2_P1_2D(cRRef,cRRef,cRRef );

// Second derivatives
Real der2fctx_xx_P1_2D(cRRef,cRRef,cRRef);

static const Real refcoor_P1_2D[9] = {0.  ,0.  ,0.,
				      1.  ,0.  ,0.,
				      0.  ,1.  ,0.};

static const Fct fct_P1_2D[3] = {fct1_P1_2D,fct2_P1_2D,fct3_P1_2D};

static const Fct derfct_P1_2D[6] = {derfct1_1_P1_2D,derfct1_2_P1_2D,
				    derfct2_1_P1_2D,derfct2_2_P1_2D,
				    derfct3_1_P1_2D,derfct3_2_P1_2D};
static const Fct der2fct_P1_2D[12] =
{der2fctx_xx_P1_2D,der2fctx_xx_P1_2D,der2fctx_xx_P1_2D,der2fctx_xx_P1_2D,
 der2fctx_xx_P1_2D,der2fctx_xx_P1_2D,der2fctx_xx_P1_2D,der2fctx_xx_P1_2D,
 der2fctx_xx_P1_2D,der2fctx_xx_P1_2D,der2fctx_xx_P1_2D,der2fctx_xx_P1_2D};
//======================================================================
//
//                            P2  (2D)
//
//======================================================================
/*
                           3
                           |\
                           6 5
                           |  \
                           1-4-2
*/
Real fct1_P2_2D(cRRef x,cRRef y,cRRef );
Real fct2_P2_2D(cRRef x,cRRef y,cRRef );
Real fct3_P2_2D(cRRef x,cRRef y,cRRef );
Real fct4_P2_2D(cRRef x,cRRef y,cRRef );
Real fct5_P2_2D(cRRef x,cRRef y,cRRef );
Real fct6_P2_2D(cRRef x,cRRef y,cRRef );

Real derfct1_1_P2_2D(cRRef x,cRRef y,cRRef );
Real derfct1_2_P2_2D(cRRef x,cRRef y,cRRef );
Real derfct2_1_P2_2D(cRRef x,cRRef  ,cRRef );
Real derfct2_2_P2_2D(cRRef  ,cRRef  ,cRRef );
Real derfct3_1_P2_2D(cRRef  ,cRRef  ,cRRef );
Real derfct3_2_P2_2D(cRRef  ,cRRef y,cRRef );
Real derfct4_1_P2_2D(cRRef x,cRRef y,cRRef );
Real derfct4_2_P2_2D(cRRef x,cRRef  ,cRRef );
Real derfct5_1_P2_2D(cRRef  ,cRRef y,cRRef );
Real derfct5_2_P2_2D(cRRef x,cRRef  ,cRRef );
Real derfct6_1_P2_2D(cRRef  ,cRRef y,cRRef );
Real derfct6_2_P2_2D(cRRef x,cRRef y,cRRef );

Real der2fct1_11_P2_2D(cRRef x,cRRef y,cRRef );
Real der2fct1_12_P2_2D(cRRef x,cRRef y,cRRef );
Real der2fct1_21_P2_2D(cRRef x,cRRef y,cRRef );
Real der2fct1_22_P2_2D(cRRef x,cRRef y,cRRef );

Real der2fct2_11_P2_2D(cRRef x,cRRef  ,cRRef );
Real der2fct2_12_P2_2D(cRRef  ,cRRef  ,cRRef );
Real der2fct2_21_P2_2D(cRRef x,cRRef  ,cRRef );
Real der2fct2_22_P2_2D(cRRef  ,cRRef  ,cRRef );

Real der2fct3_11_P2_2D(cRRef  ,cRRef  ,cRRef );
Real der2fct3_12_P2_2D(cRRef  ,cRRef y,cRRef );
Real der2fct3_21_P2_2D(cRRef  ,cRRef  ,cRRef );
Real der2fct3_22_P2_2D(cRRef  ,cRRef y,cRRef );

Real der2fct4_11_P2_2D(cRRef x,cRRef y,cRRef );
Real der2fct4_12_P2_2D(cRRef x,cRRef  ,cRRef );
Real der2fct4_21_P2_2D(cRRef x,cRRef y,cRRef );
Real der2fct4_22_P2_2D(cRRef x,cRRef  ,cRRef );

Real der2fct5_11_P2_2D(cRRef  ,cRRef y,cRRef );
Real der2fct5_12_P2_2D(cRRef x,cRRef  ,cRRef );
Real der2fct5_21_P2_2D(cRRef  ,cRRef y,cRRef );
Real der2fct5_22_P2_2D(cRRef x,cRRef  ,cRRef );

Real der2fct6_11_P2_2D(cRRef  ,cRRef y,cRRef );
Real der2fct6_12_P2_2D(cRRef x,cRRef y,cRRef );
Real der2fct6_21_P2_2D(cRRef  ,cRRef y,cRRef );
Real der2fct6_22_P2_2D(cRRef x,cRRef y,cRRef );

static const Real refcoor_P2_2D[18] = {0.  ,0.  ,0.,
				       1.  ,0.  ,0.,
				       0.  ,1.  ,0.,
				       0.5 ,0.  ,0.,
				       0.5 ,0.5 ,0.,
				       0.  ,0.5 ,0.};

static const Fct fct_P2_2D[6] = {fct1_P2_2D,fct2_P2_2D,fct3_P2_2D,
				 fct4_P2_2D,fct5_P2_2D,fct6_P2_2D};

static const Fct derfct_P2_2D[12] = {derfct1_1_P2_2D,derfct1_2_P2_2D,
				    derfct2_1_P2_2D,derfct2_2_P2_2D,
				    derfct3_1_P2_2D,derfct3_2_P2_2D,
				    derfct4_1_P2_2D,derfct4_2_P2_2D,
				    derfct5_1_P2_2D,derfct5_2_P2_2D,
				    derfct6_1_P2_2D,derfct6_2_P2_2D};
static const Fct der2fct_P2_2D[24] =
{der2fct1_11_P2_2D,der2fct1_12_P2_2D,der2fct1_21_P2_2D,der2fct1_22_P2_2D,
 der2fct2_11_P2_2D,der2fct2_12_P2_2D,der2fct2_21_P2_2D,der2fct2_22_P2_2D,
 der2fct3_11_P2_2D,der2fct3_12_P2_2D,der2fct3_21_P2_2D,der2fct3_22_P2_2D,
 der2fct4_11_P2_2D,der2fct4_12_P2_2D,der2fct4_21_P2_2D,der2fct4_22_P2_2D,
 der2fct5_11_P2_2D,der2fct5_12_P2_2D,der2fct5_21_P2_2D,der2fct5_22_P2_2D,
 der2fct6_11_P2_2D,der2fct6_12_P2_2D,der2fct6_21_P2_2D,der2fct6_22_P2_2D};


//======================================================================
//
//                            Q0  (2D)
//
//======================================================================
/*
                            -------
                           |       |
                           |   1   |
                           |       |
                            -------

*/
Real fct1_Q0_2D(cRRef , cRRef , cRRef );
Real derfct1_Q0_2D(cRRef  , cRRef , cRRef );
// The second derivative is equal to the first : both = 0.
Real der2fct1_Q0_2D(cRRef , cRRef , cRRef );

static const Real refcoor_Q0_2D[3] = {0.5, 0.5, 0. };

static const Fct fct_Q0_2D[1] = {fct1_Q0_2D};

static const Fct derfct_Q0_2D[2] = {derfct1_Q0_2D,derfct1_Q0_2D};

static const Fct der2fct_Q0_2D[4] = {der2fct1_Q0_2D,der2fct1_Q0_2D,
				     der2fct1_Q0_2D,der2fct1_Q0_2D};

//======================================================================
//
//                            Q1  (2D)
//
//======================================================================
/*
                           4-------3
                           |       |
                           |       |
                           |       |
                           1-------2
*/
Real fct1_Q1_2D(cRRef x,cRRef y,cRRef );
Real fct2_Q1_2D(cRRef x,cRRef y,cRRef );
Real fct3_Q1_2D(cRRef x,cRRef y,cRRef );
Real fct4_Q1_2D(cRRef x,cRRef y,cRRef );

Real derfct1_1_Q1_2D(cRRef x,cRRef y,cRRef );
Real derfct1_2_Q1_2D(cRRef x,cRRef y,cRRef );
Real derfct2_1_Q1_2D(cRRef x,cRRef y,cRRef );
Real derfct2_2_Q1_2D(cRRef x,cRRef y,cRRef );
Real derfct3_1_Q1_2D(cRRef x,cRRef y,cRRef );
Real derfct3_2_Q1_2D(cRRef x,cRRef y,cRRef );
Real derfct4_1_Q1_2D(cRRef x,cRRef y,cRRef );
Real derfct4_2_Q1_2D(cRRef x,cRRef y,cRRef );

// Second derivatives
Real der2fctx_xx_Q1_2D(cRRef,cRRef,cRRef);

static const Real refcoor_Q1_2D[12] = {0.  ,0.  ,0.,
				       1.  ,0.  ,0.,
				       1.  ,1.  ,0.,
				       0.  ,1.  ,0.};

static const Fct fct_Q1_2D[4] = {fct1_Q1_2D,fct2_Q1_2D,fct3_Q1_2D,fct4_Q1_2D};

static const Fct derfct_Q1_2D[8] = {derfct1_1_Q1_2D,derfct1_2_Q1_2D,
				    derfct2_1_Q1_2D,derfct2_2_Q1_2D,
				    derfct3_1_Q1_2D,derfct3_2_Q1_2D,
				    derfct4_1_Q1_2D,derfct4_2_Q1_2D};
static const Fct der2fct_Q1_2D[16] =
{der2fctx_xx_Q1_2D,der2fctx_xx_Q1_2D,der2fctx_xx_Q1_2D,der2fctx_xx_Q1_2D,
 der2fctx_xx_Q1_2D,der2fctx_xx_Q1_2D,der2fctx_xx_Q1_2D,der2fctx_xx_Q1_2D,
 der2fctx_xx_Q1_2D,der2fctx_xx_Q1_2D,der2fctx_xx_Q1_2D,der2fctx_xx_Q1_2D,
 der2fctx_xx_Q1_2D,der2fctx_xx_Q1_2D,der2fctx_xx_Q1_2D,der2fctx_xx_Q1_2D};

//======================================================================
//
//                            Q2  (2D)
//
//======================================================================
/*
                           4---7---3
                           |       |
                           8   9   6
                           |       |
                           1---5---2
*/
Real fct1_Q2_2D(cRRef x,cRRef y,cRRef );
Real fct5_Q2_2D(cRRef x,cRRef y,cRRef );
Real fct2_Q2_2D(cRRef x,cRRef y,cRRef );
Real fct6_Q2_2D(cRRef x,cRRef y,cRRef );
Real fct3_Q2_2D(cRRef x,cRRef y,cRRef );
Real fct7_Q2_2D(cRRef x,cRRef y,cRRef );
Real fct4_Q2_2D(cRRef x,cRRef y,cRRef );
Real fct8_Q2_2D(cRRef x,cRRef y,cRRef );
Real fct9_Q2_2D(cRRef x,cRRef y,cRRef );

Real derfct1_1_Q2_2D(cRRef x,cRRef y,cRRef );
Real derfct1_2_Q2_2D(cRRef x,cRRef y,cRRef );
Real derfct5_1_Q2_2D(cRRef x,cRRef y,cRRef );
Real derfct5_2_Q2_2D(cRRef x,cRRef y,cRRef );
Real derfct2_1_Q2_2D(cRRef x,cRRef y,cRRef );
Real derfct2_2_Q2_2D(cRRef x,cRRef y,cRRef );
Real derfct6_1_Q2_2D(cRRef x,cRRef y,cRRef );
Real derfct6_2_Q2_2D(cRRef x,cRRef y,cRRef );
Real derfct3_1_Q2_2D(cRRef x,cRRef y,cRRef );
Real derfct3_2_Q2_2D(cRRef x,cRRef y,cRRef );
Real derfct7_1_Q2_2D(cRRef x,cRRef y,cRRef );
Real derfct7_2_Q2_2D(cRRef x,cRRef y,cRRef );
Real derfct4_1_Q2_2D(cRRef x,cRRef y,cRRef );
Real derfct4_2_Q2_2D(cRRef x,cRRef y,cRRef );
Real derfct8_1_Q2_2D(cRRef x,cRRef y,cRRef );
Real derfct8_2_Q2_2D(cRRef x,cRRef y,cRRef );
Real derfct9_1_Q2_2D(cRRef x,cRRef y,cRRef );
Real derfct9_2_Q2_2D(cRRef x,cRRef y,cRRef );

Real der2fct1_11_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct1_12_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct1_21_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct1_22_Q2_2D(cRRef x,cRRef y,cRRef );

Real der2fct5_11_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct5_12_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct5_21_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct5_22_Q2_2D(cRRef x,cRRef y,cRRef );

Real der2fct2_11_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct2_12_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct2_21_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct2_22_Q2_2D(cRRef x,cRRef y,cRRef );

Real der2fct6_11_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct6_12_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct6_21_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct6_22_Q2_2D(cRRef x,cRRef y,cRRef );

Real der2fct3_11_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct3_12_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct3_21_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct3_22_Q2_2D(cRRef x,cRRef y,cRRef );

Real der2fct7_11_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct7_12_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct7_21_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct7_22_Q2_2D(cRRef x,cRRef y,cRRef );

Real der2fct4_11_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct4_12_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct4_21_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct4_22_Q2_2D(cRRef x,cRRef y,cRRef );

Real der2fct8_11_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct8_12_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct8_21_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct8_22_Q2_2D(cRRef x,cRRef y,cRRef );

Real der2fct9_11_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct9_12_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct9_21_Q2_2D(cRRef x,cRRef y,cRRef );
Real der2fct9_22_Q2_2D(cRRef x,cRRef y,cRRef );


static const Real refcoor_Q2_2D[27] = {0.  ,0.  ,0.,
				       1.  ,0.  ,0.,
				       1.  ,1.  ,0.,
				       0.  ,1.  ,0.,
				       0.5 ,0.  ,0.,
				       1.  ,0.5 ,0.,
				       0.5 ,1.  ,0.,
				       0.  ,0.5 ,0.,
				       0.5 ,0.5 ,0.};

static const Fct fct_Q2_2D[9] = {fct1_Q2_2D,fct2_Q2_2D,fct3_Q2_2D,fct4_Q2_2D,
				 fct5_Q2_2D,fct6_Q2_2D,fct7_Q2_2D,fct8_Q2_2D,
				 fct9_Q2_2D};


static const Fct derfct_Q2_2D[18] = {derfct1_1_Q2_2D,derfct1_2_Q2_2D,
				     derfct2_1_Q2_2D,derfct2_2_Q2_2D,
				     derfct3_1_Q2_2D,derfct3_2_Q2_2D,
				     derfct4_1_Q2_2D,derfct4_2_Q2_2D,
				     derfct5_1_Q2_2D,derfct5_2_Q2_2D,
				     derfct6_1_Q2_2D,derfct6_2_Q2_2D,
				     derfct7_1_Q2_2D,derfct7_2_Q2_2D,
				     derfct8_1_Q2_2D,derfct8_2_Q2_2D,
				     derfct9_1_Q2_2D,derfct9_2_Q2_2D};

static const Fct der2fct_Q2_2D[36] =
{
  der2fct1_11_Q2_2D,der2fct1_12_Q2_2D,der2fct1_21_Q2_2D,der2fct1_22_Q2_2D,
  der2fct2_11_Q2_2D,der2fct2_12_Q2_2D,der2fct2_21_Q2_2D,der2fct2_22_Q2_2D,
  der2fct3_11_Q2_2D,der2fct3_12_Q2_2D,der2fct3_21_Q2_2D,der2fct3_22_Q2_2D,
  der2fct4_11_Q2_2D,der2fct4_12_Q2_2D,der2fct4_21_Q2_2D,der2fct4_22_Q2_2D,
  der2fct5_11_Q2_2D,der2fct5_12_Q2_2D,der2fct5_21_Q2_2D,der2fct5_22_Q2_2D,
  der2fct6_11_Q2_2D,der2fct6_12_Q2_2D,der2fct6_21_Q2_2D,der2fct6_22_Q2_2D,
  der2fct7_11_Q2_2D,der2fct7_12_Q2_2D,der2fct7_21_Q2_2D,der2fct7_22_Q2_2D,
  der2fct8_11_Q2_2D,der2fct8_12_Q2_2D,der2fct8_21_Q2_2D,der2fct8_22_Q2_2D,
  der2fct9_11_Q2_2D,der2fct9_12_Q2_2D,der2fct9_21_Q2_2D,der2fct9_22_Q2_2D
};

//======================================================================
//
//                            P0  (3D)
//
//======================================================================
/*                 
                
               / .  
              /  \.
             /  . \\
            / . 1  \\
           /.       \!
           ----------
*/
Real fct1_P0_3D(cRRef x,cRRef y,cRRef z);

Real derfct1_P0_3D(cRRef,cRRef,cRRef );

// Second derivatives
Real der2fct1_P0_3D(cRRef,cRRef,cRRef);

static const Real refcoor_P0_3D[3] = {0.25  ,0.25  ,0.25};

static const Fct fct_P0_3D[1] = {fct1_P0_3D};

static const Fct derfct_P0_3D[3] = {derfct1_P0_3D, derfct1_P0_3D, derfct1_P0_3D};
static const Fct der2fct_P0_3D[9] =
{derfct1_P0_3D, derfct1_P0_3D, derfct1_P0_3D,
 derfct1_P0_3D, derfct1_P0_3D, derfct1_P0_3D,
 derfct1_P0_3D, derfct1_P0_3D, derfct1_P0_3D};

//======================================================================
//
//                            P1  (3D)
//
//======================================================================
/*
                4
               / .
              /  \.3
             /  . \\
            / .    \\
           /.       \!
         1 ----------2
*/
Real fct1_P1_3D(cRRef x,cRRef y,cRRef z);
Real fct2_P1_3D(cRRef x,cRRef  ,cRRef  );
Real fct3_P1_3D(cRRef  ,cRRef y,cRRef  );
Real fct4_P1_3D(cRRef  ,cRRef  ,cRRef z);

Real derfct1_1_P1_3D(cRRef,cRRef,cRRef );
Real derfct1_2_P1_3D(cRRef,cRRef,cRRef );
Real derfct1_3_P1_3D(cRRef,cRRef,cRRef );
Real derfct2_1_P1_3D(cRRef,cRRef,cRRef );
Real derfct2_2_P1_3D(cRRef,cRRef,cRRef );
Real derfct2_3_P1_3D(cRRef,cRRef,cRRef );
Real derfct3_1_P1_3D(cRRef,cRRef,cRRef );
Real derfct3_2_P1_3D(cRRef,cRRef,cRRef );
Real derfct3_3_P1_3D(cRRef,cRRef,cRRef );
Real derfct4_1_P1_3D(cRRef,cRRef,cRRef );
Real derfct4_2_P1_3D(cRRef,cRRef,cRRef );
Real derfct4_3_P1_3D(cRRef,cRRef,cRRef );

// Second derivatives
Real der2fctx_xx_P1_3D(cRRef,cRRef,cRRef);

static const Real refcoor_P1_3D[12] = {0.  ,0.  ,0.,
				       1.  ,0.  ,0.,
				       0.  ,1.  ,0.,
				       0.  ,0.  ,1.};

static const Fct fct_P1_3D[4] = {fct1_P1_3D,fct2_P1_3D,fct3_P1_3D,fct4_P1_3D};

static const Fct derfct_P1_3D[12] = {derfct1_1_P1_3D,derfct1_2_P1_3D,derfct1_3_P1_3D,
				     derfct2_1_P1_3D,derfct2_2_P1_3D,derfct2_3_P1_3D,
				     derfct3_1_P1_3D,derfct3_2_P1_3D,derfct3_3_P1_3D,
				     derfct4_1_P1_3D,derfct4_2_P1_3D,derfct4_3_P1_3D};
static const Fct der2fct_P1_3D[36] =
{der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,
 der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,
 der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,
 der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,
 der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,
 der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,
 der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,
 der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,der2fctx_xx_P1_3D,der2fctx_xx_P1_3D};

//======================================================================
//
//                            P1 Bubble (3D)
//
//======================================================================
/*
                4
               / .
              /  \.3
             /  . \\
            / . .5 \\
           /.       \!
         1 ----------2
*/
Real fct1_P1bubble_3D(cRRef x,cRRef y,cRRef z);
Real fct2_P1bubble_3D(cRRef x,cRRef  ,cRRef  );
Real fct3_P1bubble_3D(cRRef  ,cRRef y,cRRef  );
Real fct4_P1bubble_3D(cRRef  ,cRRef  ,cRRef z);
Real fct5_P1bubble_3D(cRRef  ,cRRef  ,cRRef z);

Real derfct1_1_P1bubble_3D(cRRef,cRRef,cRRef );
Real derfct1_2_P1bubble_3D(cRRef,cRRef,cRRef );
Real derfct1_3_P1bubble_3D(cRRef,cRRef,cRRef );
Real derfct2_1_P1bubble_3D(cRRef,cRRef,cRRef );
Real derfct2_2_P1bubble_3D(cRRef,cRRef,cRRef );
Real derfct2_3_P1bubble_3D(cRRef,cRRef,cRRef );
Real derfct3_1_P1bubble_3D(cRRef,cRRef,cRRef );
Real derfct3_2_P1bubble_3D(cRRef,cRRef,cRRef );
Real derfct3_3_P1bubble_3D(cRRef,cRRef,cRRef );
Real derfct4_1_P1bubble_3D(cRRef,cRRef,cRRef );
Real derfct4_2_P1bubble_3D(cRRef,cRRef,cRRef );
Real derfct4_3_P1bubble_3D(cRRef,cRRef,cRRef );
Real derfct5_1_P1bubble_3D(cRRef,cRRef,cRRef );
Real derfct5_2_P1bubble_3D(cRRef,cRRef,cRRef );
Real derfct5_3_P1bubble_3D(cRRef,cRRef,cRRef );

// Second derivatives
Real der2fctx_xx_P1bubble_3D(cRRef,cRRef,cRRef);
Real der2fct5_11_P1bubble_3D(cRRef,cRRef,cRRef);
Real der2fct5_12_P1bubble_3D(cRRef,cRRef,cRRef);
Real der2fct5_13_P1bubble_3D(cRRef,cRRef,cRRef);
Real der2fct5_21_P1bubble_3D(cRRef,cRRef,cRRef);
Real der2fct5_22_P1bubble_3D(cRRef,cRRef,cRRef);
Real der2fct5_23_P1bubble_3D(cRRef,cRRef,cRRef);
Real der2fct5_31_P1bubble_3D(cRRef,cRRef,cRRef);
Real der2fct5_32_P1bubble_3D(cRRef,cRRef,cRRef);
Real der2fct5_33_P1bubble_3D(cRRef,cRRef,cRRef);

static const Real refcoor_P1bubble_3D[15] = {0.  ,0.  ,0.,
				             1.  ,0.  ,0.,
				             0.  ,1.  ,0.,
				             0.  ,0.  ,1.,
                                          0.25,0.25,0.25};

static const Fct fct_P1bubble_3D[5] = {fct1_P1bubble_3D,fct2_P1bubble_3D,fct3_P1bubble_3D,fct4_P1bubble_3D,fct5_P1bubble_3D};

static const Fct derfct_P1bubble_3D[15] = {derfct1_1_P1bubble_3D,derfct1_2_P1bubble_3D,derfct1_3_P1bubble_3D,
			                   derfct2_1_P1bubble_3D,derfct2_2_P1bubble_3D,derfct2_3_P1bubble_3D,
				           derfct3_1_P1bubble_3D,derfct3_2_P1bubble_3D,derfct3_3_P1bubble_3D,
				           derfct4_1_P1bubble_3D,derfct4_2_P1bubble_3D,derfct4_3_P1bubble_3D,
                                           derfct5_1_P1bubble_3D,derfct5_2_P1bubble_3D,derfct5_3_P1bubble_3D};
static const Fct der2fct_P1bubble_3D[45] =
{
  der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D,
  der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D,
  der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D,
  der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D,
  der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D,
  der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D,
  der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D,
  der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D,
  der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D,
  der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D,
  der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D,
  der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D,
  der2fct5_11_P1bubble_3D, der2fct5_12_P1bubble_3D, der2fct5_13_P1bubble_3D,
  der2fct5_21_P1bubble_3D, der2fct5_22_P1bubble_3D, der2fct5_23_P1bubble_3D,
  der2fct5_31_P1bubble_3D, der2fct5_32_P1bubble_3D, der2fct5_33_P1bubble_3D
};

//======================================================================
//
//                            P2  (3D)
//
//======================================================================
/*
                4
               / .10
              /  \.3
             8  . 9\
            / 7    \6
           /.       \!
         1 -----5----2
*/
Real fct1_P2_3D(cRRef x,cRRef y,cRRef z);
Real fct2_P2_3D(cRRef x,cRRef y,cRRef z);
Real fct3_P2_3D(cRRef x,cRRef y,cRRef z);
Real fct4_P2_3D(cRRef x,cRRef y,cRRef z);
Real fct5_P2_3D(cRRef x,cRRef y,cRRef z);
Real fct6_P2_3D(cRRef x,cRRef y,cRRef z);
Real fct7_P2_3D(cRRef x,cRRef y,cRRef z);
Real fct8_P2_3D(cRRef x,cRRef y,cRRef z);
Real fct9_P2_3D(cRRef x,cRRef y,cRRef z);
Real fct10_P2_3D(cRRef x,cRRef y,cRRef z);


Real derfct1_1_P2_3D(cRRef x,cRRef y,cRRef z);
Real derfct1_2_P2_3D(cRRef x,cRRef y,cRRef z);
Real derfct1_3_P2_3D(cRRef x,cRRef y,cRRef z);

Real derfct2_1_P2_3D(cRRef x,cRRef y,cRRef z);
Real derfct2_2_P2_3D(cRRef x,cRRef y,cRRef z);
Real derfct2_3_P2_3D(cRRef x,cRRef y,cRRef z);

Real derfct3_1_P2_3D(cRRef x,cRRef y,cRRef z);
Real derfct3_2_P2_3D(cRRef x,cRRef y,cRRef z);
Real derfct3_3_P2_3D(cRRef x,cRRef y,cRRef z);

Real derfct4_1_P2_3D(cRRef x,cRRef y,cRRef z);
Real derfct4_2_P2_3D(cRRef x,cRRef y,cRRef z);
Real derfct4_3_P2_3D(cRRef x,cRRef y,cRRef z);

Real derfct5_1_P2_3D(cRRef x,cRRef y,cRRef z);
Real derfct5_2_P2_3D(cRRef x,cRRef y,cRRef z);
Real derfct5_3_P2_3D(cRRef x,cRRef y,cRRef z);

Real derfct6_1_P2_3D(cRRef x,cRRef y,cRRef z);
Real derfct6_2_P2_3D(cRRef x,cRRef y,cRRef z);
Real derfct6_3_P2_3D(cRRef x,cRRef y,cRRef z);

Real derfct7_1_P2_3D(cRRef x,cRRef y,cRRef z);
Real derfct7_2_P2_3D(cRRef x,cRRef y,cRRef z);
Real derfct7_3_P2_3D(cRRef x,cRRef y,cRRef z);

Real derfct8_1_P2_3D(cRRef x,cRRef y,cRRef z);
Real derfct8_2_P2_3D(cRRef x,cRRef y,cRRef z);
Real derfct8_3_P2_3D(cRRef x,cRRef y,cRRef z);

Real derfct9_1_P2_3D(cRRef x,cRRef y,cRRef z);
Real derfct9_2_P2_3D(cRRef x,cRRef y,cRRef z);
Real derfct9_3_P2_3D(cRRef x,cRRef y,cRRef z);

Real derfct10_1_P2_3D(cRRef x,cRRef y,cRRef z);
Real derfct10_2_P2_3D(cRRef x,cRRef y,cRRef z);
Real derfct10_3_P2_3D(cRRef x,cRRef y,cRRef z);


Real der2fct1_11_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct1_12_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct1_13_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct1_21_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct1_22_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct1_23_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct1_31_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct1_32_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct1_33_P2_3D(cRRef x,cRRef y,cRRef z);

Real der2fct2_11_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct2_12_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct2_13_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct2_21_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct2_22_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct2_23_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct2_31_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct2_32_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct2_33_P2_3D(cRRef x,cRRef y,cRRef z);

Real der2fct3_11_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct3_12_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct3_13_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct3_21_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct3_22_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct3_23_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct3_31_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct3_32_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct3_33_P2_3D(cRRef x,cRRef y,cRRef z);

Real der2fct4_11_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct4_12_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct4_13_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct4_21_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct4_22_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct4_23_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct4_31_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct4_32_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct4_33_P2_3D(cRRef x,cRRef y,cRRef z);

Real der2fct5_11_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct5_12_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct5_13_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct5_21_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct5_22_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct5_23_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct5_31_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct5_32_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct5_33_P2_3D(cRRef x,cRRef y,cRRef z);

Real der2fct6_11_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct6_12_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct6_13_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct6_21_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct6_22_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct6_23_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct6_31_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct6_32_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct6_33_P2_3D(cRRef x,cRRef y,cRRef z);

Real der2fct7_11_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct7_12_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct7_13_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct7_21_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct7_22_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct7_23_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct7_31_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct7_32_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct7_33_P2_3D(cRRef x,cRRef y,cRRef z);

Real der2fct8_11_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct8_12_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct8_13_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct8_21_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct8_22_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct8_23_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct8_31_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct8_32_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct8_33_P2_3D(cRRef x,cRRef y,cRRef z);

Real der2fct9_11_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct9_12_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct9_13_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct9_21_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct9_22_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct9_23_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct9_31_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct9_32_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct9_33_P2_3D(cRRef x,cRRef y,cRRef z);

Real der2fct10_11_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct10_12_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct10_13_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct10_21_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct10_22_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct10_23_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct10_31_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct10_32_P2_3D(cRRef x,cRRef y,cRRef z);
Real der2fct10_33_P2_3D(cRRef x,cRRef y,cRRef z);


static const Real refcoor_P2_3D[30] = {0.  ,0.  ,0. ,
				       1.  ,0.  ,0. ,
				       0.  ,1.  ,0. ,
				       0.  ,0.  ,1. ,
				       0.5 ,0.  ,0. ,
				       0.5, 0.5 ,0. ,
				       0. , 0.5 ,0. ,
				       0. , 0.  ,0.5,
				       0.5, 0.  ,0.5,
				       0. , 0.5 ,0.5};

static const Fct fct_P2_3D[10] = {fct1_P2_3D,fct2_P2_3D,fct3_P2_3D,fct4_P2_3D,
				  fct5_P2_3D,fct6_P2_3D,fct7_P2_3D,fct8_P2_3D,
				  fct9_P2_3D,fct10_P2_3D};

static const Fct derfct_P2_3D[30] = {derfct1_1_P2_3D,derfct1_2_P2_3D,derfct1_3_P2_3D,
				     derfct2_1_P2_3D,derfct2_2_P2_3D,derfct2_3_P2_3D,
				     derfct3_1_P2_3D,derfct3_2_P2_3D,derfct3_3_P2_3D,
				     derfct4_1_P2_3D,derfct4_2_P2_3D,derfct4_3_P2_3D,
				     derfct5_1_P2_3D,derfct5_2_P2_3D,derfct5_3_P2_3D,
				     derfct6_1_P2_3D,derfct6_2_P2_3D,derfct6_3_P2_3D,
				     derfct7_1_P2_3D,derfct7_2_P2_3D,derfct7_3_P2_3D,
				     derfct8_1_P2_3D,derfct8_2_P2_3D,derfct8_3_P2_3D,
				     derfct9_1_P2_3D,derfct9_2_P2_3D,derfct9_3_P2_3D,
				     derfct10_1_P2_3D,derfct10_2_P2_3D,derfct10_3_P2_3D};
/*the perl-script:
  #!/usr/bin/perl
  for($i=1;$i<=10;$i++){
  for($j=1;$j<=3;$j++){
  printf "der2fct$i\_$j"."1\_P2\_3D, der2fct$i\_$j"."2\_P2\_3D, der2fct$i\_$j"."3_P2\_3D,\n";
  }
}
*/
static const Fct der2fct_P2_3D[90] =
{
  der2fct1_11_P2_3D, der2fct1_12_P2_3D, der2fct1_13_P2_3D,
  der2fct1_21_P2_3D, der2fct1_22_P2_3D, der2fct1_23_P2_3D,
  der2fct1_31_P2_3D, der2fct1_32_P2_3D, der2fct1_33_P2_3D,
  der2fct2_11_P2_3D, der2fct2_12_P2_3D, der2fct2_13_P2_3D,
  der2fct2_21_P2_3D, der2fct2_22_P2_3D, der2fct2_23_P2_3D,
  der2fct2_31_P2_3D, der2fct2_32_P2_3D, der2fct2_33_P2_3D,
  der2fct3_11_P2_3D, der2fct3_12_P2_3D, der2fct3_13_P2_3D,
  der2fct3_21_P2_3D, der2fct3_22_P2_3D, der2fct3_23_P2_3D,
  der2fct3_31_P2_3D, der2fct3_32_P2_3D, der2fct3_33_P2_3D,
  der2fct4_11_P2_3D, der2fct4_12_P2_3D, der2fct4_13_P2_3D,
  der2fct4_21_P2_3D, der2fct4_22_P2_3D, der2fct4_23_P2_3D,
  der2fct4_31_P2_3D, der2fct4_32_P2_3D, der2fct4_33_P2_3D,
  der2fct5_11_P2_3D, der2fct5_12_P2_3D, der2fct5_13_P2_3D,
  der2fct5_21_P2_3D, der2fct5_22_P2_3D, der2fct5_23_P2_3D,
  der2fct5_31_P2_3D, der2fct5_32_P2_3D, der2fct5_33_P2_3D,
  der2fct6_11_P2_3D, der2fct6_12_P2_3D, der2fct6_13_P2_3D,
  der2fct6_21_P2_3D, der2fct6_22_P2_3D, der2fct6_23_P2_3D,
  der2fct6_31_P2_3D, der2fct6_32_P2_3D, der2fct6_33_P2_3D,
  der2fct7_11_P2_3D, der2fct7_12_P2_3D, der2fct7_13_P2_3D,
  der2fct7_21_P2_3D, der2fct7_22_P2_3D, der2fct7_23_P2_3D,
  der2fct7_31_P2_3D, der2fct7_32_P2_3D, der2fct7_33_P2_3D,
  der2fct8_11_P2_3D, der2fct8_12_P2_3D, der2fct8_13_P2_3D,
  der2fct8_21_P2_3D, der2fct8_22_P2_3D, der2fct8_23_P2_3D,
  der2fct8_31_P2_3D, der2fct8_32_P2_3D, der2fct8_33_P2_3D,
  der2fct9_11_P2_3D, der2fct9_12_P2_3D, der2fct9_13_P2_3D,
  der2fct9_21_P2_3D, der2fct9_22_P2_3D, der2fct9_23_P2_3D,
  der2fct9_31_P2_3D, der2fct9_32_P2_3D, der2fct9_33_P2_3D,
  der2fct10_11_P2_3D, der2fct10_12_P2_3D, der2fct10_13_P2_3D,
  der2fct10_21_P2_3D, der2fct10_22_P2_3D, der2fct10_23_P2_3D,
  der2fct10_31_P2_3D, der2fct10_32_P2_3D, der2fct10_33_P2_3D
};
//======================================================================
//
//                            P2tilde  (3D)
// NAVIER-STOKES P2 Basis Oriented to the mass lumping
//
//======================================================================
/*
                4
               / .10
              /  \.3
             8  . 9\
            / 7 .11 6
           /.       \!
         1 -----5----2
*/
Real fct1_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real fct2_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real fct3_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real fct4_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real fct5_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real fct6_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real fct7_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real fct8_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real fct9_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real fct10_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real fct11_P2tilde_3D(cRRef x,cRRef y,cRRef z);


Real derfct1_1_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real derfct1_2_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real derfct1_3_P2tilde_3D(cRRef x,cRRef y,cRRef z);

Real derfct2_1_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real derfct2_2_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real derfct2_3_P2tilde_3D(cRRef x,cRRef y,cRRef z);

Real derfct3_1_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real derfct3_2_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real derfct3_3_P2tilde_3D(cRRef x,cRRef y,cRRef z);

Real derfct4_1_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real derfct4_2_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real derfct4_3_P2tilde_3D(cRRef x,cRRef y,cRRef z);

Real derfct5_1_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real derfct5_2_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real derfct5_3_P2tilde_3D(cRRef x,cRRef y,cRRef z);

Real derfct6_1_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real derfct6_2_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real derfct6_3_P2tilde_3D(cRRef x,cRRef y,cRRef z);

Real derfct7_1_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real derfct7_2_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real derfct7_3_P2tilde_3D(cRRef x,cRRef y,cRRef z);

Real derfct8_1_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real derfct8_2_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real derfct8_3_P2tilde_3D(cRRef x,cRRef y,cRRef z);

Real derfct9_1_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real derfct9_2_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real derfct9_3_P2tilde_3D(cRRef x,cRRef y,cRRef z);

Real derfct10_1_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real derfct10_2_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real derfct10_3_P2tilde_3D(cRRef x,cRRef y,cRRef z);

Real derfct11_1_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real derfct11_2_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real derfct11_3_P2tilde_3D(cRRef x,cRRef y,cRRef z);


Real der2fct1_11_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct1_12_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct1_13_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct1_21_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct1_22_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct1_23_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct1_31_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct1_32_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct1_33_P2tilde_3D(cRRef x,cRRef y,cRRef z);

Real der2fct2_11_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct2_12_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct2_13_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct2_21_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct2_22_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct2_23_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct2_31_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct2_32_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct2_33_P2tilde_3D(cRRef x,cRRef y,cRRef z);

Real der2fct3_11_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct3_12_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct3_13_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct3_21_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct3_22_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct3_23_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct3_31_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct3_32_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct3_33_P2tilde_3D(cRRef x,cRRef y,cRRef z);

Real der2fct4_11_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct4_12_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct4_13_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct4_21_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct4_22_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct4_23_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct4_31_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct4_32_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct4_33_P2tilde_3D(cRRef x,cRRef y,cRRef z);

Real der2fct5_11_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct5_12_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct5_13_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct5_21_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct5_22_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct5_23_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct5_31_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct5_32_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct5_33_P2tilde_3D(cRRef x,cRRef y,cRRef z);

Real der2fct6_11_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct6_12_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct6_13_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct6_21_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct6_22_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct6_23_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct6_31_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct6_32_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct6_33_P2tilde_3D(cRRef x,cRRef y,cRRef z);

Real der2fct7_11_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct7_12_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct7_13_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct7_21_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct7_22_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct7_23_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct7_31_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct7_32_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct7_33_P2tilde_3D(cRRef x,cRRef y,cRRef z);

Real der2fct8_11_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct8_12_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct8_13_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct8_21_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct8_22_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct8_23_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct8_31_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct8_32_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct8_33_P2tilde_3D(cRRef x,cRRef y,cRRef z);

Real der2fct9_11_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct9_12_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct9_13_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct9_21_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct9_22_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct9_23_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct9_31_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct9_32_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct9_33_P2tilde_3D(cRRef x,cRRef y,cRRef z);

Real der2fct10_11_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct10_12_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct10_13_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct10_21_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct10_22_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct10_23_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct10_31_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct10_32_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct10_33_P2tilde_3D(cRRef x,cRRef y,cRRef z);

Real der2fct11_11_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct11_12_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct11_13_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct11_21_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct11_22_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct11_23_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct11_31_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct11_32_P2tilde_3D(cRRef x,cRRef y,cRRef z);
Real der2fct11_33_P2tilde_3D(cRRef x,cRRef y,cRRef z);


static const Real refcoor_P2tilde_3D[33] = {0.  ,0.  ,0. ,
				       1.  ,0.   ,0. ,
				       0.  ,1.   ,0. ,
				       0.  ,0.   ,1. ,
				       0.5 ,0.   ,0. ,
				       0.5 , 0.5 , 0. ,
				       0.  , 0.5 , 0. ,
				       0.  , 0.  , 0.5,
				       0.5 , 0.  , 0.5,
				       0.  , 0.5 , 0.5,
				       0.25, 0.25, 0.25};

static const Fct fct_P2tilde_3D[11] =
{fct1_P2tilde_3D,fct2_P2tilde_3D,fct3_P2tilde_3D,fct4_P2tilde_3D,
 fct5_P2tilde_3D,fct6_P2tilde_3D,fct7_P2tilde_3D,fct8_P2tilde_3D,
 fct9_P2tilde_3D,fct10_P2tilde_3D,fct11_P2tilde_3D};

static const Fct derfct_P2tilde_3D[33] =
{derfct1_1_P2tilde_3D,derfct1_2_P2tilde_3D,derfct1_3_P2tilde_3D,
 derfct2_1_P2tilde_3D,derfct2_2_P2tilde_3D,derfct2_3_P2tilde_3D,
 derfct3_1_P2tilde_3D,derfct3_2_P2tilde_3D,derfct3_3_P2tilde_3D,
 derfct4_1_P2tilde_3D,derfct4_2_P2tilde_3D,derfct4_3_P2tilde_3D,
 derfct5_1_P2tilde_3D,derfct5_2_P2tilde_3D,derfct5_3_P2tilde_3D,
 derfct6_1_P2tilde_3D,derfct6_2_P2tilde_3D,derfct6_3_P2tilde_3D,
 derfct7_1_P2tilde_3D,derfct7_2_P2tilde_3D,derfct7_3_P2tilde_3D,
 derfct8_1_P2tilde_3D,derfct8_2_P2tilde_3D,derfct8_3_P2tilde_3D,
 derfct9_1_P2tilde_3D,derfct9_2_P2tilde_3D,derfct9_3_P2tilde_3D,
 derfct10_1_P2tilde_3D,derfct10_2_P2tilde_3D,derfct10_3_P2tilde_3D,
 derfct11_1_P2tilde_3D,derfct11_2_P2tilde_3D,derfct11_3_P2tilde_3D};
/*the perl-script:
  #!/usr/bin/perl
  for($i=1;$i<=10;$i++){
  for($j=1;$j<=3;$j++){
  printf "der2fct$i\_$j"."1\_P2tilde\_3D, der2fct$i\_$j"."2\_P2tilde\_3D, der2fct$i\_$j"."3_P2tilde\_3D,\n";
  }
}
*/
static const Fct der2fct_P2tilde_3D[99] =
{
  der2fct1_11_P2tilde_3D, der2fct1_12_P2tilde_3D, der2fct1_13_P2tilde_3D,
  der2fct1_21_P2tilde_3D, der2fct1_22_P2tilde_3D, der2fct1_23_P2tilde_3D,
  der2fct1_31_P2tilde_3D, der2fct1_32_P2tilde_3D, der2fct1_33_P2tilde_3D,
  der2fct2_11_P2tilde_3D, der2fct2_12_P2tilde_3D, der2fct2_13_P2tilde_3D,
  der2fct2_21_P2tilde_3D, der2fct2_22_P2tilde_3D, der2fct2_23_P2tilde_3D,
  der2fct2_31_P2tilde_3D, der2fct2_32_P2tilde_3D, der2fct2_33_P2tilde_3D,
  der2fct3_11_P2tilde_3D, der2fct3_12_P2tilde_3D, der2fct3_13_P2tilde_3D,
  der2fct3_21_P2tilde_3D, der2fct3_22_P2tilde_3D, der2fct3_23_P2tilde_3D,
  der2fct3_31_P2tilde_3D, der2fct3_32_P2tilde_3D, der2fct3_33_P2tilde_3D,
  der2fct4_11_P2tilde_3D, der2fct4_12_P2tilde_3D, der2fct4_13_P2tilde_3D,
  der2fct4_21_P2tilde_3D, der2fct4_22_P2tilde_3D, der2fct4_23_P2tilde_3D,
  der2fct4_31_P2tilde_3D, der2fct4_32_P2tilde_3D, der2fct4_33_P2tilde_3D,
  der2fct5_11_P2tilde_3D, der2fct5_12_P2tilde_3D, der2fct5_13_P2tilde_3D,
  der2fct5_21_P2tilde_3D, der2fct5_22_P2tilde_3D, der2fct5_23_P2tilde_3D,
  der2fct5_31_P2tilde_3D, der2fct5_32_P2tilde_3D, der2fct5_33_P2tilde_3D,
  der2fct6_11_P2tilde_3D, der2fct6_12_P2tilde_3D, der2fct6_13_P2tilde_3D,
  der2fct6_21_P2tilde_3D, der2fct6_22_P2tilde_3D, der2fct6_23_P2tilde_3D,
  der2fct6_31_P2tilde_3D, der2fct6_32_P2tilde_3D, der2fct6_33_P2tilde_3D,
  der2fct7_11_P2tilde_3D, der2fct7_12_P2tilde_3D, der2fct7_13_P2tilde_3D,
  der2fct7_21_P2tilde_3D, der2fct7_22_P2tilde_3D, der2fct7_23_P2tilde_3D,
  der2fct7_31_P2tilde_3D, der2fct7_32_P2tilde_3D, der2fct7_33_P2tilde_3D,
  der2fct8_11_P2tilde_3D, der2fct8_12_P2tilde_3D, der2fct8_13_P2tilde_3D,
  der2fct8_21_P2tilde_3D, der2fct8_22_P2tilde_3D, der2fct8_23_P2tilde_3D,
  der2fct8_31_P2tilde_3D, der2fct8_32_P2tilde_3D, der2fct8_33_P2tilde_3D,
  der2fct9_11_P2tilde_3D, der2fct9_12_P2tilde_3D, der2fct9_13_P2tilde_3D,
  der2fct9_21_P2tilde_3D, der2fct9_22_P2tilde_3D, der2fct9_23_P2tilde_3D,
  der2fct9_31_P2tilde_3D, der2fct9_32_P2tilde_3D, der2fct9_33_P2tilde_3D,
  der2fct10_11_P2tilde_3D, der2fct10_12_P2tilde_3D, der2fct10_13_P2tilde_3D,
  der2fct10_21_P2tilde_3D, der2fct10_22_P2tilde_3D, der2fct10_23_P2tilde_3D,
  der2fct10_31_P2tilde_3D, der2fct10_32_P2tilde_3D, der2fct10_33_P2tilde_3D,
  der2fct11_11_P2tilde_3D, der2fct11_12_P2tilde_3D, der2fct11_13_P2tilde_3D,
  der2fct11_21_P2tilde_3D, der2fct11_22_P2tilde_3D, der2fct11_23_P2tilde_3D,
  der2fct11_31_P2tilde_3D, der2fct11_32_P2tilde_3D, der2fct11_33_P2tilde_3D
};

//======================================================================
//
//                            Q0  (3D)
//
//======================================================================
/*
                      ________
                     /.      /|
		    / .     / |
		   /_______/  |
		   |  .  1 |  |
		   |  .....|..|
		   | .     | /
		   |.      |/
		   |_______|

*/
Real fct1_Q0_3D(cRRef ,cRRef ,cRRef );
Real derfct1_Q0_3D(cRRef ,cRRef ,cRRef );
// The second derivative is equal to the first : both = 0.
Real der2fct1_Q0_3D(cRRef ,cRRef ,cRRef );

static const Real refcoor_Q0_3D[3] = {0.5  ,0.5  ,0.5};


static const Fct fct_Q0_3D[1] = {fct1_Q0_3D};

static const Fct derfct_Q0_3D[3] = {derfct1_Q0_3D,derfct1_Q0_3D,derfct1_Q0_3D};

static const Fct der2fct_Q0_3D[9] ={
  der2fct1_Q0_3D,der2fct1_Q0_3D,der2fct1_Q0_3D,
  der2fct1_Q0_3D,der2fct1_Q0_3D,der2fct1_Q0_3D,
  der2fct1_Q0_3D,der2fct1_Q0_3D,der2fct1_Q0_3D};

//======================================================================
//
//                            Q1  (3D)
//
//======================================================================
/*
                      8-------7
                     /.      /|
		    / .     / |
		   5_______6  |
		   |  .    |  |
		   |  4....|..3
		   | .     | /
		   |.      |/
		   1_______2
*/
Real fct1_Q1_3D(cRRef x,cRRef y,cRRef z);
Real fct2_Q1_3D(cRRef x,cRRef y,cRRef z);
Real fct3_Q1_3D(cRRef x,cRRef y,cRRef z);
Real fct4_Q1_3D(cRRef x,cRRef y,cRRef z);
Real fct5_Q1_3D(cRRef x,cRRef y,cRRef z);
Real fct6_Q1_3D(cRRef x,cRRef y,cRRef z);
Real fct7_Q1_3D(cRRef x,cRRef y,cRRef z);
Real fct8_Q1_3D(cRRef x,cRRef y,cRRef z);

Real derfct1_1_Q1_3D(cRRef x,cRRef y,cRRef z);
Real derfct1_2_Q1_3D(cRRef x,cRRef y,cRRef z);
Real derfct1_3_Q1_3D(cRRef x,cRRef y,cRRef z);
Real derfct2_1_Q1_3D(cRRef x,cRRef y,cRRef z);
Real derfct2_2_Q1_3D(cRRef x,cRRef y,cRRef z);
Real derfct2_3_Q1_3D(cRRef x,cRRef y,cRRef z);
Real derfct3_1_Q1_3D(cRRef x,cRRef y,cRRef z);
Real derfct3_2_Q1_3D(cRRef x,cRRef y,cRRef z);
Real derfct3_3_Q1_3D(cRRef x,cRRef y,cRRef z);
Real derfct4_1_Q1_3D(cRRef x,cRRef y,cRRef z);
Real derfct4_2_Q1_3D(cRRef x,cRRef y,cRRef z);
Real derfct4_3_Q1_3D(cRRef x,cRRef y,cRRef z);
Real derfct5_1_Q1_3D(cRRef x,cRRef y,cRRef z);
Real derfct5_2_Q1_3D(cRRef x,cRRef y,cRRef z);
Real derfct5_3_Q1_3D(cRRef x,cRRef y,cRRef z);
Real derfct6_1_Q1_3D(cRRef x,cRRef y,cRRef z);
Real derfct6_2_Q1_3D(cRRef x,cRRef y,cRRef z);
Real derfct6_3_Q1_3D(cRRef x,cRRef y,cRRef z);
Real derfct7_1_Q1_3D(cRRef x,cRRef y,cRRef z);
Real derfct7_2_Q1_3D(cRRef x,cRRef y,cRRef z);
Real derfct7_3_Q1_3D(cRRef x,cRRef y,cRRef z);
Real derfct8_1_Q1_3D(cRRef x,cRRef y,cRRef z);
Real derfct8_2_Q1_3D(cRRef x,cRRef y,cRRef z);
Real derfct8_3_Q1_3D(cRRef x,cRRef y,cRRef z);

Real der2fct1_11_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct1_12_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct1_13_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct1_21_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct1_22_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct1_23_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct1_31_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct1_32_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct1_33_Q1_3D(cRRef x,cRRef y,cRRef z);

Real der2fct2_11_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct2_12_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct2_13_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct2_21_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct2_22_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct2_23_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct2_31_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct2_32_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct2_33_Q1_3D(cRRef x,cRRef y,cRRef z);

Real der2fct3_11_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct3_12_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct3_13_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct3_21_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct3_22_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct3_23_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct3_31_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct3_32_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct3_33_Q1_3D(cRRef x,cRRef y,cRRef z);

Real der2fct4_11_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct4_12_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct4_13_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct4_21_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct4_22_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct4_23_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct4_31_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct4_32_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct4_33_Q1_3D(cRRef x,cRRef y,cRRef z);

Real der2fct5_11_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct5_12_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct5_13_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct5_21_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct5_22_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct5_23_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct5_31_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct5_32_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct5_33_Q1_3D(cRRef x,cRRef y,cRRef z);

Real der2fct6_11_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct6_12_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct6_13_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct6_21_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct6_22_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct6_23_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct6_31_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct6_32_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct6_33_Q1_3D(cRRef x,cRRef y,cRRef z);

Real der2fct7_11_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct7_12_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct7_13_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct7_21_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct7_22_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct7_23_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct7_31_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct7_32_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct7_33_Q1_3D(cRRef x,cRRef y,cRRef z);

Real der2fct8_11_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct8_12_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct8_13_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct8_21_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct8_22_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct8_23_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct8_31_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct8_32_Q1_3D(cRRef x,cRRef y,cRRef z);
Real der2fct8_33_Q1_3D(cRRef x,cRRef y,cRRef z);

static const Real refcoor_Q1_3D[24] = {0.  ,0.  ,0. ,
				       1.  ,0.  ,0. ,
				       1.  ,1.  ,0. ,
				       0.  ,1.  ,0. ,
                                       0.  ,0.  ,1. ,
				       1.  ,0.  ,1. ,
				       1.  ,1.  ,1. ,
				       0.  ,1.  ,1.};


static const Fct fct_Q1_3D[8] = {fct1_Q1_3D,fct2_Q1_3D,fct3_Q1_3D,fct4_Q1_3D,fct5_Q1_3D,
				 fct6_Q1_3D,fct7_Q1_3D,fct8_Q1_3D};


static const Fct derfct_Q1_3D[24] = {derfct1_1_Q1_3D,derfct1_2_Q1_3D,derfct1_3_Q1_3D,
				     derfct2_1_Q1_3D,derfct2_2_Q1_3D,derfct2_3_Q1_3D,
				     derfct3_1_Q1_3D,derfct3_2_Q1_3D,derfct3_3_Q1_3D,
				     derfct4_1_Q1_3D,derfct4_2_Q1_3D,derfct4_3_Q1_3D,
				     derfct5_1_Q1_3D,derfct5_2_Q1_3D,derfct5_3_Q1_3D,
				     derfct6_1_Q1_3D,derfct6_2_Q1_3D,derfct6_3_Q1_3D,
				     derfct7_1_Q1_3D,derfct7_2_Q1_3D,derfct7_3_Q1_3D,
				     derfct8_1_Q1_3D,derfct8_2_Q1_3D,derfct8_3_Q1_3D};
/*the perl-script:
  #!/usr/bin/perl
  for($i=1;$i<=10;$i++){
  for($j=1;$j<=3;$j++){
  printf "der2fct$i\_$j"."1\_P2\_3D, der2fct$i\_$j"."2\_P2\_3D, der2fct$i\_$j"."3_P2\_3D,\n";
  }
}
*/
static const Fct der2fct_Q1_3D[72] =
{
  der2fct1_11_Q1_3D, der2fct1_12_Q1_3D, der2fct1_13_Q1_3D,
  der2fct1_21_Q1_3D, der2fct1_22_Q1_3D, der2fct1_23_Q1_3D,
  der2fct1_31_Q1_3D, der2fct1_32_Q1_3D, der2fct1_33_Q1_3D,
  der2fct2_11_Q1_3D, der2fct2_12_Q1_3D, der2fct2_13_Q1_3D,
  der2fct2_21_Q1_3D, der2fct2_22_Q1_3D, der2fct2_23_Q1_3D,
  der2fct2_31_Q1_3D, der2fct2_32_Q1_3D, der2fct2_33_Q1_3D,
  der2fct3_11_Q1_3D, der2fct3_12_Q1_3D, der2fct3_13_Q1_3D,
  der2fct3_21_Q1_3D, der2fct3_22_Q1_3D, der2fct3_23_Q1_3D,
  der2fct3_31_Q1_3D, der2fct3_32_Q1_3D, der2fct3_33_Q1_3D,
  der2fct4_11_Q1_3D, der2fct4_12_Q1_3D, der2fct4_13_Q1_3D,
  der2fct4_21_Q1_3D, der2fct4_22_Q1_3D, der2fct4_23_Q1_3D,
  der2fct4_31_Q1_3D, der2fct4_32_Q1_3D, der2fct4_33_Q1_3D,
  der2fct5_11_Q1_3D, der2fct5_12_Q1_3D, der2fct5_13_Q1_3D,
  der2fct5_21_Q1_3D, der2fct5_22_Q1_3D, der2fct5_23_Q1_3D,
  der2fct5_31_Q1_3D, der2fct5_32_Q1_3D, der2fct5_33_Q1_3D,
  der2fct6_11_Q1_3D, der2fct6_12_Q1_3D, der2fct6_13_Q1_3D,
  der2fct6_21_Q1_3D, der2fct6_22_Q1_3D, der2fct6_23_Q1_3D,
  der2fct6_31_Q1_3D, der2fct6_32_Q1_3D, der2fct6_33_Q1_3D,
  der2fct7_11_Q1_3D, der2fct7_12_Q1_3D, der2fct7_13_Q1_3D,
  der2fct7_21_Q1_3D, der2fct7_22_Q1_3D, der2fct7_23_Q1_3D,
  der2fct7_31_Q1_3D, der2fct7_32_Q1_3D, der2fct7_33_Q1_3D,
  der2fct8_11_Q1_3D, der2fct8_12_Q1_3D, der2fct8_13_Q1_3D,
  der2fct8_21_Q1_3D, der2fct8_22_Q1_3D, der2fct8_23_Q1_3D,
  der2fct8_31_Q1_3D, der2fct8_32_Q1_3D, der2fct8_33_Q1_3D};
}
#endif
