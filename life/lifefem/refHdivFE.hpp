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
#ifndef _REFHDIVFE_H_INCLUDE
#define _REFHDIVFE_H_INCLUDE

#include <life/lifecore/life.hpp>
#include <life/lifearray/tab.hpp>
#include <life/lifemesh/basisElSh.hpp>
#include <life/lifefem/quadRule.hpp>
#include <life/lifefem/localDofPattern.hpp>

/*!
  \file refHdivFE.h
  \brief Classes RefHdivFE and RefHdivHybridFE
*/

namespace LifeV
{

//! Indicates the LOCAL (reference element) coordinates
typedef const Real & cRRef;
typedef Real ( * FCT ) ( cRRef, cRRef , cRRef );

/*!
  \class RefHdivFE
  \brief Class for H(div,*) vectorial functions.
  \author M. Belhadj & V. Martin
  \date 07/2002

  This is a duplication-modification of both RefEle and RefFE
  in order to take into account Mixed Finite Elements, which are
  based on a (RT0 - Q0) like discretization of H(div, .) - L2(.).

  Here follows commentaries from refEle and refFE :

  This class contains the basis functions and their values on quadrature points.

  \par How to add a new reference element:

  in refHdiv.h : you declare the functions you need (fct1_Pipo_2D,
  derfct1_1_Pipo_2D, etc...), the static arrays containing these functions
  and the coordinates of the nodes on the reference element.

  and in defQuadRuleFE.cc : you define these functions (fct1_Pipo_2D, etc...)

*/

class RefHdivFE:
            public LocalDofPattern
{
private:
    const SetOfQuadRule* _sqr; //!< pointer on the set of quadrature rules
    const FCT* _phi; //!< pointer on the basis functions
    const FCT* _divPhi; //!< pointer on the divergence of the basis functions
    const Real* _refCoor; //!< reference coordinates. Order: xi_1,eta_1,zeta_1,xi_2,eta_2,zeta_2,...

    //! values of the basis functions on all quadrature points
    KN<Real> _phiQuad;
    //! values of the divergence of the basis functions on all quadrature points
    KN<Real> _divPhiQuad;
    KN<int> _idxQuad; //!< _idxQuad[t] = index of the quadrature rules of id t in _phiQuad
    KN<int> _idxDQuad; //!< _idxDQuad[t] = index of the quadrature rules of id t in _divPhiQuad
public:
    const std::string name; //!< name of the reference element
    const int type; //!< Type of finite element (FE_P1_2D, ..., see the #define at the beginning of refFE.h)
    const ReferenceShapes shape; //!< geometrical shape of the element
    const int nbDof;   //!< Total number of degrees of freedom
    const int nbCoor;  //!< Number of local coordinates
public:
    //! Constructor of a reference finite element.
    /*!
      Constructor of a reference finite element. The arguments are:

      _name : the name of the f.e.

      _type : the type of the f.e. (FE_P1_2D,... see the #define at the begining of refFE.h)

      _shape : the geometry belongs to enum ReferenceShapes {NONE, POINT, LINE, TRIANGLE, QUAD, HEXA, PRISM, TETRA}; (see basisElSh.h)

      _nbDofPerVertex : the number of degrees of freedom per vertex

      _nbDofPerEdge : the number of degrees of freedom per edge

      _nbDofPerFace : the number of degrees of freedom per face

      _nbDofPerVolume : the number of degrees of freedom per volume

      _nbDof : the total number of d.o.f ( = _nbDofPerVertex * nb vertex + _nbDofPerEdge * nb edges + etc...)

      _nbCoor : number of local coordinates

      phi : the static array containing the basis functions (defined in RefHdivFE.h)

      divPhi : the static array containing the divergence of the basis functions (defined in RefHdivFE.h)

      refCoor : the static array containing the coordinates of the nodes on the reference element (defined in
      RefHdivFE.h)

      sqr : a set of quadrature rule (defined in quadRule.cc)

      _patternType : in most of cases STANDARD_PATTERN, except for elements like P1isoP2
      (to define a new pattern, add a new #define in refFE.h and code it in refFE.cc following the
      example of P1ISOP2_TRIA_PATTERN)
    */
    RefHdivFE( std::string _name, int _type, ReferenceShapes _shape,
               int _nbDofPerVertex, int _nbDofPerEdge, int _nbDofPerFace, int _nbDofPerVolume,
               int _nbDof, int _nbCoor, const FCT* phi, const FCT* divPhi,
               const Real* _refCoor, const SetOfQuadRule& sqr, PatternType _patternType );
    ~RefHdivFE();
    //! return the first local coordinate of the i-th node of the reference element
    inline Real xi( int i ) const
    {
        ASSERT_BD( i < nbDof )
        return _refCoor[ 3 * i ];
    }
    //! return the second local coordinate of the i-th node of the reference element
    inline Real eta( int i ) const
    {
        ASSERT_BD( i < nbDof )
        return _refCoor[ 3 * i + 1 ];
    }
    //! return the third local coordinate of the i-th node of the reference element
    inline Real zeta( int i ) const
    {
        ASSERT_BD( i < nbDof )
        return _refCoor[ 3 * i + 2 ];
    }
    //! return the icoor-th local coordinate of the i-th node of the reference element
    inline Real refCoor( int i, int icoor ) const
    {
        ASSERT_BD( i < nbDof && icoor < nbCoor )
        return _refCoor[ 3 * i + icoor ];
    }
    //! return the value of the i-th basis function on point (x,y,z)
    inline double phi( int i, int icomp, cRRef x, cRRef y, cRRef z ) const
    {
        ASSERT_BD( i < nbDof && icomp < nbCoor )
        //      std::cout << _phi[3*i+icomp](x,y,z) << std::endl;
        return _phi[ nbCoor * i + icomp ] ( x, y, z );
    }
    //!return the value of the i-th basis function on the ig-th point of integration of quadrature rule qr
    inline Real phi( int i, int icomp, int ig, const QuadRule& qr ) const
    {
        ASSERT_BD( i < nbDof && ig < qr.nbQuadPt && icomp < nbCoor )
        return _phiQuad( _idxQuad( qr.id ) + ig * nbDof * nbCoor + i * nbCoor + icomp );
    }
    //! return the array of the values of the basis functions phi_1[ig],phi_2[ig],... on the integration point ig of the quadrature rule qr.
    /*!This function is written in order to avoid too many calls to function
      phi(i,ig,qr), please check if it improves really efficiency. If yes,
      analogous function should be written for the derivatives.  */
    inline RN_ phiQuadPt( int ig, const QuadRule& qr ) const
    {
        ASSERT_BD( ig < qr.nbQuadPt )
        return _phiQuad( SubArray( nbDof, _idxQuad( qr.id ) + ig * nbDof ) );
    }
    //! return the array of the values of the i-th basis functions phi_i[ig=1],phi_i[ig=2],...on all the integration points of the quadrature rule qr.
    /*! This function is written in order to avoid too many calls to
      function phi(i,ig,qr), please check if it improves really
      efficiency. If yes, analogous function should be written for the
      derivatives.  */
    inline RN_ phiI( int i, const QuadRule& qr ) const
    {
        ASSERT_BD( i < nbDof )
        return _phiQuad( SubArray( qr.nbQuadPt, _idxQuad( qr.id ) + i, nbDof ) );
    }
    //! return the value of the divergence of the i-th basis function on point (x,y,z)
    inline double divPhi( int i, cRRef x, cRRef y, cRRef z ) const
    {
        ASSERT_BD( i < nbDof )
        return _divPhi[ i ] ( x, y, z );
    }
    //! return the value of icoor-th divergence of the i-th basis function on the ig-th point of integration of quadrature rule qr
    inline Real divPhi( int i, int ig, const QuadRule& qr ) const
    {
        ASSERT_BD( i < nbDof && ig < qr.nbQuadPt )
        return _divPhiQuad( _idxDQuad( qr.id ) + ig * nbDof + i );
    }
    //!< A simple check function
    void check() const;
    friend std::ostream& operator << ( std::ostream& f, const RefHdivFE& fe );
};

//!======================================================================
//!
//!                           RT0  (3D)
//!
//!======================================================================
/*!

8-------7
/.      /|
/ .     / |
5_______6  |
|  .    |  |
|  4....|..3
| .     | /
|.      |/
1_______2

face 1: 1,4,3,2
face 2: 1,5,8,4
face 3: 1,2,6,5
face 4: 2,3,7,6
face 5: 3,4,8,7
face 6: 5,6,7,8

*/

Real fct1_RT0_1_3D( cRRef x, cRRef y, cRRef z );
Real fct1_RT0_2_3D( cRRef x, cRRef y, cRRef z );
Real fct1_RT0_3_3D( cRRef x, cRRef y, cRRef z );

Real fct2_RT0_1_3D( cRRef x, cRRef y, cRRef z );
Real fct2_RT0_2_3D( cRRef x, cRRef y, cRRef z );
Real fct2_RT0_3_3D( cRRef x, cRRef y, cRRef z );

Real fct3_RT0_1_3D( cRRef x, cRRef y, cRRef z );
Real fct3_RT0_2_3D( cRRef x, cRRef y, cRRef z );
Real fct3_RT0_3_3D( cRRef x, cRRef y, cRRef z );

Real fct4_RT0_1_3D( cRRef x, cRRef y, cRRef z );
Real fct4_RT0_2_3D( cRRef x, cRRef y, cRRef z );
Real fct4_RT0_3_3D( cRRef x, cRRef y, cRRef z );

Real fct5_RT0_1_3D( cRRef x, cRRef y, cRRef z );
Real fct5_RT0_2_3D( cRRef x, cRRef y, cRRef z );
Real fct5_RT0_3_3D( cRRef x, cRRef y, cRRef z );

Real fct6_RT0_1_3D( cRRef x, cRRef y, cRRef z );
Real fct6_RT0_2_3D( cRRef x, cRRef y, cRRef z );
Real fct6_RT0_3_3D( cRRef x, cRRef y, cRRef z );

Real fct1_DIV_RT0_3D( cRRef x, cRRef y, cRRef z );
Real fct2_DIV_RT0_3D( cRRef x, cRRef y, cRRef z );
Real fct3_DIV_RT0_3D( cRRef x, cRRef y, cRRef z );
Real fct4_DIV_RT0_3D( cRRef x, cRRef y, cRRef z );
Real fct5_DIV_RT0_3D( cRRef x, cRRef y, cRRef z );
Real fct6_DIV_RT0_3D( cRRef x, cRRef y, cRRef z );

static const Real refcoor_RT0_3D[ 18 ] =
    {
        0.5 , 0.5 , 0. ,
        0. , 0.5 , 0.5 ,
        0.5 , 0. , 0.5 ,
        1. , 0.5 , 0.5 ,
        0.5 , 1. , 0.5 ,
        0.5 , 0.5 , 1.
    };

static const FCT fct_RT0_3D[ 18 ] =
    {
        fct1_RT0_1_3D, fct1_RT0_2_3D, fct1_RT0_3_3D,
        fct2_RT0_1_3D, fct2_RT0_2_3D, fct2_RT0_3_3D,
        fct3_RT0_1_3D, fct3_RT0_2_3D, fct3_RT0_3_3D,
        fct4_RT0_1_3D, fct4_RT0_2_3D, fct4_RT0_3_3D,
        fct5_RT0_1_3D, fct5_RT0_2_3D, fct5_RT0_3_3D,
        fct6_RT0_1_3D, fct6_RT0_2_3D, fct6_RT0_3_3D
    };

static const FCT fct_DIV_RT0_3D[ 6 ] =
    {
        fct1_DIV_RT0_3D, fct2_DIV_RT0_3D,
        fct3_DIV_RT0_3D, fct4_DIV_RT0_3D,
        fct5_DIV_RT0_3D, fct6_DIV_RT0_3D
    };


//!======================================================================
//!
//!                           RT0  (3D)
//!
//!======================================================================
/*

4
/ .
/  \.3
/  . \\
/ .    \\
/.       \!
1 ----------2

SEE basisElSh.cc   for the ORIENTATION CONVENTIONS
point 1: 0, 0, 0
point 2: 1, 0, 0
point 3: 0, 1, 0
point 4: 0, 0, 1

face 1: 2, 3, 4
face 2: 1, 4, 3
face 3: 1, 2, 4
face 4: 1, 3, 2
*/

Real fct1_RT0_1_3D_TETRA( cRRef x, cRRef y, cRRef z );
Real fct1_RT0_2_3D_TETRA( cRRef x, cRRef y, cRRef z );
Real fct1_RT0_3_3D_TETRA( cRRef x, cRRef y, cRRef z );

Real fct2_RT0_1_3D_TETRA( cRRef x, cRRef y, cRRef z );
Real fct2_RT0_2_3D_TETRA( cRRef x, cRRef y, cRRef z );
Real fct2_RT0_3_3D_TETRA( cRRef x, cRRef y, cRRef z );

Real fct3_RT0_1_3D_TETRA( cRRef x, cRRef y, cRRef z );
Real fct3_RT0_2_3D_TETRA( cRRef x, cRRef y, cRRef z );
Real fct3_RT0_3_3D_TETRA( cRRef x, cRRef y, cRRef z );

Real fct4_RT0_1_3D_TETRA( cRRef x, cRRef y, cRRef z );
Real fct4_RT0_2_3D_TETRA( cRRef x, cRRef y, cRRef z );
Real fct4_RT0_3_3D_TETRA( cRRef x, cRRef y, cRRef z );


Real fct1_DIV_RT0_3D_TETRA( cRRef x, cRRef y, cRRef z );
Real fct2_DIV_RT0_3D_TETRA( cRRef x, cRRef y, cRRef z );
Real fct3_DIV_RT0_3D_TETRA( cRRef x, cRRef y, cRRef z );
Real fct4_DIV_RT0_3D_TETRA( cRRef x, cRRef y, cRRef z );


static const Real refcoor_RT0_3D_TETRA[ 12 ] =
    {
        1. / 3 , 1. / 3. , 0. ,
        1. / 3. , 0. , 1. / 3. ,
        1. / 3. , 1. / 3. , 1. / 3. ,
        0. , 1. / 3. , 1. / 3.
    };

static const FCT fct_RT0_3D_TETRA[ 12 ] =
    {
        fct1_RT0_1_3D_TETRA, fct1_RT0_2_3D_TETRA, fct1_RT0_3_3D_TETRA,
        fct2_RT0_1_3D_TETRA, fct2_RT0_2_3D_TETRA, fct2_RT0_3_3D_TETRA,
        fct3_RT0_1_3D_TETRA, fct3_RT0_2_3D_TETRA, fct3_RT0_3_3D_TETRA,
        fct4_RT0_1_3D_TETRA, fct4_RT0_2_3D_TETRA, fct4_RT0_3_3D_TETRA
    };

static const FCT fct_DIV_RT0_3D_TETRA[ 4 ] =
    {
        fct1_DIV_RT0_3D_TETRA, fct2_DIV_RT0_3D_TETRA,
        fct3_DIV_RT0_3D_TETRA, fct4_DIV_RT0_3D_TETRA
    };


//!======================================================================
//!     DECLARATION OF FINITE ELEMENTS (defined in defQuadRule.cc)
extern const RefHdivFE feHexaRT0;

extern const RefHdivFE feTetraRT0;

}
#endif
