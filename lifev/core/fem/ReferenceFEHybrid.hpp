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
    @brief Reference finite element for hybrid FEs.

    @author Alessio Fumagalli
            Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 10-05-2010

    @contributor
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef REFFEHYBRID_H
#define REFFEHYBRID_H 1

#include <lifev/core/fem/CurrentFEManifold.hpp>

namespace LifeV
{

/*!
  @class RefHybridFE
  @brief Class for Hybrid functions, i.e. defined on the boundary of an element.
  @author V. Martin
  @date 08/2002

  This is an enrichment of ReferenceFE in order to implement mixed hybrid finite elements,
  which are based on a (RT0 - Q0) like discretization of \f$ H(div, \Omega) - L^2(\Omega) \f$.

  This class contains a list of boundary elements; thanks to the Piola transform, the computations
  are performed on the boundary of the reference element. But in general, the boundary of a 3D reference
  element is not a 2D reference element.
  <BR>
  Example:
  <BR>
  REFERENCE TETRA -> 3 REFERENCE TRIA + 1 EQUILATERAL TRIANGLE...
  <BR>
  REFERENCE PRISM -> 2 TRIA + 3 QUAD...?
  <BR>
  REFERENCE HEXA  -> 6 REFERENCE QUAD.

  @par How to add a new finite element ?

  @li In refFE.h: you add a new finite element flag in FE_TYPE enum with a command like:
  \code
FE_PIPO = a_new_number
  \endcode

  @li In refHybridFE.h: you declare the functions you need (fct1_Pipo_2D,
      derfct1_1_Pipo_2D, etc...), the static arrays containing these functions
      and the coordinates of the nodes on the reference element.

  @li In defQuadRuleFE.cc: you define these functions (fct1_Pipo_2D, etc...)

  @li In refFE.h you declare your finite element:
  \code
extern const RefHybridFE fePipo;
  \endcode

  @li In defQuadRuleFE.cc: you define a list of CurrentFEManifold with a command like:
  \code
  #define NB_BDFE_PIPO
static const CurrentFEManifold BdFE_PIPO_1( feTriaP0, geoLinearTria, quadRuleTria4pt, refcoor_PIPO_1, 0 );
...
  \endcode

  @li In defQuadRuleFE.cc: you define a static array containing all the CurrentFEManifold
  with a command like
  \code
static const CurrentFEManifold HybPIPOList[ NB_BDFE_PIPO ] =
{
     BdFE_PIPO_1, BdFE_PIPO_2,
     ...
};
  \endcode

  @li In defQuadRuleFE.cc: you define your new element with a command like:
  \code
const ReferenceFEHybrid feTriaPipo("Pipo elements on a tetrahedron", FE_PIPO, TETRA,
                             0, 0, 1, 0, 4, 3, NB_BDFE_PIPO, HybPIPOList, refcoor_PIPO, STANDARD_PATTERN );
  \endcode
  See documentation of ReferenceFEHybrid::ReferenceFEHybrid(...) for a precise description of all arguments.
*/
class ReferenceFEHybrid:
    public ReferenceFE
{
public:

    typedef ReferenceFE::function_Type function_Type;

    //! @name Constructor & Destructor
    //@{

    /*!
      Constructor of a reference hybrid finite element. The arguments are:
      @param name Name of the finite element
      @param type Type of the finite element (FE_P1_2D,... see the #define at the begining of refFE.h)
      @param shape Geometry belongs to enum ReferenceShapes {NONE, POINT, LINE, TRIANGLE, QUAD, HEXA, PRISM, TETRA}; (see ElementShapes.h)
      @param nbDofPerVertex Number of degrees of freedom per vertex
      @param nbDofPerEdge Number of degrees of freedom per edge
      @param nbDofPerFace Number of degrees of freedom per face
      @param nbDofPerVolume Number of degrees of freedom per volume
      @param nbDof Total number of degrees of freedom ( = nbDofPerVertex * nb vertex + nbDofPerEdge * nb edges + etc...)
      @param nbLocalCoor Number of local coordinates
      @param refCoor Static array containing the coordinates of the nodes on the reference element
      @param numBoundaryFE Number of static boundary elements
      @param boundaryFEList List of static boundary elements
      @param refCoor Static array containing the coordinates of the nodes on
      the reference element (defined in refEle.h)
      @param patternType In most of cases is STANDARD_PATTERN, except for elements
      like P1isoP2 (to define a new pattern, add a new #define in refFE.h and
      code it in refFE.cc following the example of P1ISOP2_TRIA_PATTERN)
    */
    ReferenceFEHybrid ( std::string        name,
                        FE_TYPE            type,
                        ReferenceShapes    shape,
                        UInt               nbDofPerVertex,
                        UInt               nbDofPerEdge,
                        UInt               nbDofPerFace,
                        UInt               nbDofPerVolume,
                        UInt               nbDof,
                        UInt                       nbLocalCoor,
                        const UInt&        numberBoundaryFE,
                        const CurrentFEManifold**  boundaryFEList,
                        const Real*        refCoor,
                        DofPatternType     patternType = STANDARD_PATTERN );

    //! Destructor.
    ~ReferenceFEHybrid();

    //@}


    //! @name Operators
    //@{

    //! Extracting a CurrentFEManifold from the faces list.
    const CurrentFEManifold& operator[] ( const ID& i ) const
    {
        ASSERT_BD ( i < static_cast<ID> ( M_numberBoundaryFE ) );
        return * (M_boundaryFEList[ i ]);
    }

    //@}


    //! @name Get Methods
    //@{

    //! Return the number of boundary elements of the reference element.
    const UInt& numberBoundaryFE() const
    {
        return M_numberBoundaryFE;
    }

    //@}


private:

    //! No empty constructor.
    ReferenceFEHybrid();

    //! No copy constructor.
    ReferenceFEHybrid ( const ReferenceFEHybrid& );

    //! Number of boundary elements to be stored.
    const UInt M_numberBoundaryFE;

    /*! List holding the stored boundary elements that live on the boundary faces (3D),
        or edges (2D), of the RefHybridFE element. The boundary elements of a reference
        element are not in general reference elements themselves, that is why
        we use here the CurrentFEManifold rather that ReferenceFE. */
    const CurrentFEManifold** M_boundaryFEList;
};



//======================================================================
//     DECLARATION OF FINITE ELEMENTS (defined in defQuadRule.cc)
extern const ReferenceFEHybrid feTriaRT0Hyb;
extern const ReferenceFEHybrid feTriaRT0VdotNHyb;

extern const ReferenceFEHybrid feHexaRT0Hyb;
extern const ReferenceFEHybrid feHexaRT0VdotNHyb;

extern const ReferenceFEHybrid feTetraRT0Hyb;
extern const ReferenceFEHybrid feTetraRT0VdotNHyb;

//======================================================================
//
//                           RT0 TRIA HYBRID (2D)
//                Element defined on SEGMENTS :  P0 on each TRIA face.
//
//======================================================================
/*


*/

static const Real refcoor_RT0HYB_TRIA [ 9 ] =
{
    1. / 2. , 0.      , 0. ,
    1. / 2. , 1. / 2. , 0. ,
    0.      , 1. / 2. , 0.
};


/* refcoor_HYB_TRIA_SEG_I
   Coordinates of the vertices of the 4 segments. They are used for the definition of the POINT in the bdfe.
   The same arrays can be used for all the Hybrid elements that have the same shape.
   E.g. refcoor_HYB_TRIA_SEG_I are used for all the TRIA Hybrid elements. */

static const Real refcoor_HYB_TRIA_SEG_1[ 6 ] =
{
    0. , 0. , 0. ,
    1. , 0. , 0.
};

static const Real refcoor_HYB_TRIA_SEG_2[ 6 ] =
{
    1. , 0. , 0. ,
    0. , 1. , 0. ,
};

static const Real refcoor_HYB_TRIA_SEG_3[ 6 ] =
{
    0. , 1. , 0. ,
    0. , 0. , 0.
};

//======================================================================
//
//                           RT0 HEXA HYBRID (3D)
//                Element defined on FACES :  Q0 on each QUAD face.
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

SEE ElementShapes.cc   for the ORIENTATION CONVENTIONS
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

/* Not really useful(?). to be  removed? in this case, remove also xi, eta, zeta, etc. in the class.
   this info is included in Staticbdfe. */
static const Real refcoor_RT0HYB_HEXA[ 18 ] =
{
    0.5 , 0.5 , 0. ,
    0.  , 0.5 , 0.5 ,
    0.5 , 0.  , 0.5 ,
    1.  , 0.5 , 0.5 ,
    0.5 , 1.  , 0.5 ,
    0.5 , 0.5 , 1.
};



/* refcoor_HYB_HEXA_FACE_I
   Coordinates of the vertices of the 6 faces. They are used for the definition of the POINT in the bdfe.
   The same arrays can be used for all the Hybrid elements that have the same shape.
   E.g. refcoor_HYB_HEXA_FACE_I are used for all the HEXA Hybrid elements. */

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
//                           RT0 TETRA HYBRID (3D)
//                Element defined on FACES :  Q0 on each QUAD face.
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

SEE ElementShapes.cc   for the ORIENTATION CONVENTIONS
   point 1: 0, 0, 0
   point 2: 1, 0, 0
   point 3: 0, 1, 0
   point 4: 0, 0, 1

   face 1: 1, 3, 2
   face 2: 1, 2, 4
   face 3: 2, 3, 4
   face 4: 1, 4, 3


*/

/* Not really useful(?). to be  removed? in this case, remove also xi, eta, zeta, etc. in the class.
   this info is included in Staticbdfe. */
static const Real refcoor_RT0HYB_TETRA[ 12 ] =
{
    1. / 3 , 1. / 3. , 0. ,
    1. / 3. , 0. , 1. / 3. ,
    1. / 3. , 1. / 3. , 1. / 3. ,
    0. , 1. / 3. , 1. / 3.
};


/* refcoor_HYB_TETRA_FACE_I
   Coordinates of the vertices of the 6 faces. They are used for the definition of the POINT in the bdfe.
   The same arrays can be used for all the Hybrid elements that have the same shape.
   E.g. refcoor_HYB_TETRA_FACE_I are used for all the TETRA Hybrid elements. */

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


} // Namespace LifeV

#endif /* REFFEHYBRID_H */
