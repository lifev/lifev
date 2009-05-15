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
#ifndef _REFHYBRIDFE_H_INCLUDE
#define _REFHYBRIDFE_H_INCLUDE

#include <life/lifecore/life.hpp>
#include <life/lifearray/tab.hpp>
#include <life/lifemesh/basisElSh.hpp>
#include <life/lifefem/refFE.hpp>
#include <life/lifefem/quadRule.hpp>
#include <life/lifefem/staticBdFE.hpp>
#include <life/lifefem/localDofPattern.hpp>

/*!
  \file refHybridFE.h
  \brief Class RefHybridFE
*/

namespace LifeV
{
//! Indicates the LOCAL (reference element) coordinates
typedef const Real & cRRef;
typedef Real ( * FCT ) ( cRRef, cRRef , cRRef );

/*!
  \class RefHybridFE
  \brief Class for Hybrid functions i.e. defined on the boundary of an element.
  \author V. Martin
  \date 08/2002

  This is a duplication-modification of both RefEle and RefFE
  in order to implement Mixed Hybrid Finite Elements, which are
  based on a (RT0 - Q0) like discretization of H(div, .) - L2(.).

  This class contains a list of boundary elements.
                  Thanks to the Piola transform, the computations are performed
    on the boundary of the REFERENCE Element. But in general, the
    boundary of a 3D Reference element is not a 2D Reference element.
    Example:
    REFERENCE TETRA -> 3 REFERENCE TRIA + 1 EQUILATERAL TRIANGLE...
    REFERENCE PRISM -> 2 TRIA + 3 QUAD...?
    REFERENCE HEXA  -> 6 REFERENCE QUAD.

  \par How to add a new reference element:

  in refHybridFE.h : you declare the functions you need (fct1_Pipo_2D,
  derfct1_1_Pipo_2D, etc...), the static arrays containing these functions
  and the coordinates of the nodes on the reference element.

  and in defQuadRuleFE.cc : you define these functions (fct1_Pipo_2D, etc...)

*/
class RefHybridFE:
            public LocalDofPattern
{
private:
    //! number of boundary elements to be stored
    const UInt _nBdFE;

    //! list holding the stored boundary elements that live on the boundary faces (3D)
    //! (or edges (2D)) of the RefHybridFE element.
    //! The boundary elements of a reference
    //! element are NOT in general reference elements themselves, that is why
    //! we use here the StaticBdFE rather that RefFE
    const StaticBdFE* _bdfeList;

    const Real* _refCoor; //!< reference coordinates. Order: xi_1,eta_1,zeta_1,xi_2,eta_2,zeta_2,...
    //! useless...

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

       refCoor : the static array containing the coordinates of the nodes on the reference element.

       _patternType : in most of cases STANDARD_PATTERN, except for elements like P1isoP2
      (to define a new pattern, add a new #define in refFE.h and code it in refFE.cc following the
      example of P1ISOP2_TRIA_PATTERN)
      ( What's the use of this here?? )
     */
    RefHybridFE( const UInt& nbdfe, const StaticBdFE* bdfelist,
                 std::string _name, int _type, ReferenceShapes _shape,
                 int _nbDofPerVertex, int _nbDofPerEdge, int _nbDofPerFace, int _nbDofPerVolume,
                 int _nbDof, int _nbCoor, const Real* _refCoor, PatternType _patternType = STANDARD_PATTERN );
    ~RefHybridFE();

    //! extracting a BdFE from the face list.
    const StaticBdFE& operator[] ( const Index_t& i ) const;  //to be checked

    //! return the number of boundary elements of the reference element.
    INLINE UInt nBdFE() const
    {
        return _nBdFE;
    }

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

    void check() const; //!< A simple check function
    friend std::ostream& operator << ( std::ostream& f, const RefHybridFE& fe );
};



//======================================================================
//     DECLARATION OF FINITE ELEMENTS (defined in defQuadRule.cc)

extern const RefHybridFE feHexaRT0Hyb;
extern const RefHybridFE feHexaRT0VdotNHyb;

extern const RefHybridFE feHexaRT1Hyb;
extern const RefHybridFE feHexaRT1VdotNHyb;

extern const RefHybridFE feTetraRT0Hyb;
extern const RefHybridFE feTetraRT0VdotNHyb;

//======================================================================
//
//                           RT0 HEXA HYBRID (3D)
//                Element defined on FACES :  Q0 on each QUAD face.
//
//======================================================================
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

SEE basisElSh.cc   for the ORIENTATION CONVENTIONS
   point 1: 0, 0, 0
   point 2: 1, 0, 0
   point 3: 1, 1, 0
   point 4: 1, 0, 0
   point 5: 0, 0, 1
   point 6: 1, 0, 1
   point 7: 1, 1, 1
   point 8: 1, 0, 1

   face 1: 1,4,3,2
   face 2: 1,5,8,4
   face 3: 1,2,6,5
   face 4: 2,3,7,6
   face 5: 3,4,8,7
   face 6: 5,6,7,8

*/

// not really useful(?). to be  removed? in this case, remove also xi, eta, zeta, etc. in the class.
// this info is included in Staticbdfe.
static const Real refcoor_RT0HYB_HEXA[ 18 ] =
    {
        0.5 , 0.5 , 0. ,
        0. , 0.5 , 0.5 ,
        0.5 , 0. , 0.5 ,
        1. , 0.5 , 0.5 ,
        0.5 , 1. , 0.5 ,
        0.5 , 0.5 , 1.
    };



/*!
  \param refcoor_HYB_HEXA_FACE_I
       Coordinates of the vertices of the 6 faces.
       They are used for the definition of the POINT in the bdfe.

       The same arrays can be used for all the Hybrid elements that have the same shape.
       ex. :  refcoor_HYB_HEXA_FACE_I are used for all the HEXA Hybrid elements.
*/
static const Real refcoor_HYB_HEXA_FACE_1[ 12 ] =
    {
        0. , 0. , 0. ,
        0. , 1. , 0. ,
        1. , 1. , 0. ,
        1. , 0. , 0.
    };

static const Real refcoor_HYB_HEXA_FACE_2[ 12 ] =
    {
        0. , 0. , 0. ,
        0. , 0. , 1. ,
        0. , 1. , 1. ,
        0. , 1. , 0.
    };

static const Real refcoor_HYB_HEXA_FACE_3[ 12 ] =
    {
        0. , 0. , 0. ,
        1. , 0. , 0. ,
        1. , 0. , 1. ,
        0. , 0. , 1.
    };

static const Real refcoor_HYB_HEXA_FACE_4[ 12 ] =
    {
        1. , 0. , 0. ,
        1. , 1. , 0. ,
        1. , 1. , 1. ,
        1. , 0. , 1.
    };

static const Real refcoor_HYB_HEXA_FACE_5[ 12 ] =
    {
        1. , 1. , 0. ,
        0. , 1. , 0. ,
        0. , 1. , 1. ,
        1. , 1. , 1.
    };

static const Real refcoor_HYB_HEXA_FACE_6[ 12 ] =
    {
        0. , 0. , 1. ,
        1. , 0. , 1. ,
        1. , 1. , 1. ,
        0. , 1. , 1.
    };


//======================================================================
//
//                           RT1 HEXA HYBRID (3D)
//                Element defined on FACES :  Q1 on each QUAD face.
//
//======================================================================
// not really useful(?). to be removed?
// this info is included in Staticbdfe.
static const Real refcoor_RT1HYB_HEXA[ 72 ] =
    {
        0. , 0. , 0. ,    // face 1
        0. , 1. , 0. ,
        1. , 1. , 0. ,
        1. , 0. , 0. ,

        0. , 0. , 0. ,    // face 2
        0. , 0. , 1. ,
        0. , 1. , 1. ,
        0. , 1. , 0. ,

        0. , 0. , 0. ,    //face 3
        1. , 0. , 0. ,
        1. , 0. , 1. ,
        0. , 0. , 1. ,

        1. , 0. , 0. ,    //face 4
        1. , 1. , 0. ,
        1. , 1. , 1. ,
        1. , 0. , 1. ,

        1. , 1. , 0. ,    //face 5
        0. , 1. , 0. ,
        0. , 1. , 1. ,
        1. , 1. , 1. ,

        0. , 0. , 1. ,    //face 6
        1. , 0. , 1. ,
        1. , 1. , 1. ,
        0. , 1. , 1.
    };


//======================================================================
//
//                           RT0 TETRA HYBRID (3D)
//                Element defined on FACES :  Q0 on each QUAD face.
//
//======================================================================
/*!

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

   face 1: 1, 3, 2
   face 2: 1, 2, 4
   face 3: 2, 3, 4
   face 4: 1, 4, 3


*/

// not really useful(?). to be  removed? in this case, remove also xi, eta, zeta, etc. in the class.
// this info is included in Staticbdfe.
//! for the TETRA : These values are FALSE!!
static const Real refcoor_RT0HYB_TETRA[ 12 ] =
    {
        1. / 3 , 1. / 3. , 0. ,
        1. / 3. , 0. , 1. / 3. ,
        1. / 3. , 1. / 3. , 1. / 3. ,
        0. , 1. / 3. , 1. / 3.
    };


/*!
  \param refcoor_HYB_TETRA_FACE_I
       Coordinates of the vertices of the 6 faces.
       They are used for the definition of the POINT in the bdfe.

       The same arrays can be used for all the Hybrid elements that have the same shape.
       ex. :  refcoor_HYB_TETRA_FACE_I are used for all the TETRA Hybrid elements.
*/
static const Real refcoor_HYB_TETRA_FACE_1[ 9 ] =
    {
        0. , 0. , 0. ,
        0. , 1. , 0. ,
        1. , 0. , 0.
    };

static const Real refcoor_HYB_TETRA_FACE_2[ 9 ] =
    {
        0. , 0. , 0. ,
        1. , 0. , 0. ,
        0. , 0. , 1.
    };

static const Real refcoor_HYB_TETRA_FACE_3[ 9 ] =
    {
        1. , 0. , 0. ,
        0. , 1. , 0. ,
        0. , 0. , 1.
    };

static const Real refcoor_HYB_TETRA_FACE_4[ 9 ] =
    {
        0. , 0. , 0. ,
        0. , 0. , 1. ,
        0. , 1. , 0.
    };

}
#endif
