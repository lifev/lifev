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
    @brief Reference finite element for scalar lagrangian FEs.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 10-05-2010

    @contributor
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */


#ifndef REFFESCALAR_H
#define REFFESCALAR_H 1


#include <life/lifefem/refFE.hpp>

namespace LifeV
{

//! refFEScalar - Short description of the class
/*!

 */
class RefFEScalar
        : public RefFE
{
public:

    //! @name Public Types
    //@{

    typedef RefFE::function_Type function_Type;
    typedef std::vector<Real> (*ValuesToValuesFct) (const std::vector<Real>&);

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Constructor of a reference Lagrangian finite element.
    /*!
      Constructor of a reference finite element. The arguments are:
      @param name  the name of the f.e.
      @param type  the type of the f.e. (FE_P1_2D,... see the #define at the
      begining of refFE.h)
      @param shape  the geometry belongs to enum ReferenceShapes {NONE, POINT,
      LINE, TRIANGLE, QUAD, HEXA, PRISM, TETRA}; (see basisElSh.h)
      @param nbDofPerVertex  the number of degrees of freedom per vertex
      @param nbDofPerEdge  the number of degrees of freedom per edge
      @param nbDofPerFace  the number of degrees of freedom per face
      @param nbDofPerVolume  the number of degrees of freedom per volume
      @param nbDof  the total number of d.o.f ( = nbDofPerVertex * nb vertex +
      nbDofPerEdge * nb edges + etc...)
      @param nbCoor  number of local coordinates
      @param phi  the static array containing the basis functions (defined in
      refEle.h)
      @param dPhi  the static array containing the derivatives of the basis
      functions (defined in refEle.h)
      @param d2Phi  the static array containing the second derivatives of the
      basis functions (defined in refEle.h)
      @param refCoor  the static array containing the coordinates of the nodes on
      the reference element (defined in refEle.h)
      @param patternType  in most of cases STANDARD_PATTERN, except for elements
      like P1isoP2 (to define a new pattern, add a new #define in refFE.h and
      code it in refFE.cc following the example of P1ISOP2_TRIA_PATTERN)
      @param bdRefFE  a pointer on the associated reference finite element on the boundary
    */
    RefFEScalar( std::string          name,
                 FE_TYPE              type,
                 ReferenceShapes      shape,
                 Int                  nbDofPerVertex,
                 Int                  nbDofPerEdge,
                 Int                  nbDofPerFace,
                 Int                  nbDofPerVolume,
                 Int                  nbDof,
                 Int                  nbCoor,
                 const function_Type*           phi,
                 const function_Type*           dPhi,
                 const function_Type*           d2Phi,
                 const Real*          refCoor,
                 DofPatternType       patternType,
                 const RefFE*         bdRefFE,
                 const ValuesToValuesFct nodalToFE);

    ~RefFEScalar(){};

    //@}


    //! @name Methods
    //@{

    std::vector<Real> nodalToFEValues(const std::vector<Real>& nodalValues) const
    {
        ASSERT( nodalValues.size() == nbDof() ,"Number of nodal values does not match with the number of degrees of freedom");
        return M_nodalToFEValues(nodalValues);
    }

    //@}

private:

    const ValuesToValuesFct M_nodalToFEValues;

};

//--------------------------------------------------

extern const RefFEScalar feSegP1;
extern const RefFEScalar feSegP2;

extern const RefFEScalar feTriaP0;
extern const RefFEScalar feTriaP1;
extern const RefFEScalar feTriaP2;

extern const RefFEScalar feQuadQ0;
extern const RefFEScalar feQuadQ1;
extern const RefFEScalar feQuadQ2;

extern const RefFEScalar feTetraP0;
extern const RefFEScalar feTetraP1;
extern const RefFEScalar feTetraP1bubble;
extern const RefFEScalar feTetraP2;
extern const RefFEScalar feTetraP2tilde;

extern const RefFEScalar feHexaQ0;
extern const RefFEScalar feHexaQ1;

} // Namespace LifeV

#endif /* REFFESCALAR_H */
