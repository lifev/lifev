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
 *  @file
 *  @brief File containing a class for the boundary conditions of the 1D model.
 *
 *  @version 1.0
 *  @date 01-28-2006
 *  @author Lucia Mirabella <lucia@mathcs.emory.edu>
 *  @author Tiziano Passerini <tiziano@mathcs.emory.edu>
 *
 *  @version 2.0
 *  @date 20-04-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @contributors Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef OneDimensionalBC_H
#define OneDimensionalBC_H

// LIFEV - MATHCARD
#include <lifemc/lifefem/OneDimensionalBCFunctionDefault.hpp>
#include <lifemc/lifesolver/OneDimensionalData.hpp>

namespace LifeV
{

//! OneDimensionalBC - Class featuring methods to handle boundary conditions.
/*!
 *  @author Lucia Mirabella, Tiziano Passerini
 */
class OneDimensionalBC
{
public:

    //! @name Type definitions
    //@{

    typedef OneDimensionalBCFunctionDefault              bcFunctionDefault_Type;
    typedef boost::shared_ptr< bcFunctionDefault_Type >  bcFunctionDefaultPtr_Type;

    typedef bcFunctionDefault_Type::bcFunction_Type      bcFunction_Type;
    typedef bcFunctionDefault_Type::bcFunctionPtr_Type   bcFunctionPtr_Type;

    typedef bcFunctionDefault_Type::fluxPtr_Type         fluxPtr_Type;
    typedef bcFunctionDefault_Type::sourcePtr_Type       sourcePtr_Type;
    typedef bcFunctionDefault_Type::solution_Type        solution_Type;
    typedef bcFunctionDefault_Type::solutionPtr_Type     solutionPtr_Type;

    typedef bcFunctionDefault_Type::container2D_Type     container2D_Type;

    typedef bcFunctionDefault_Type::bcLine_Type          bcLine_Type;
    typedef bcFunctionDefault_Type::bcSide_Type          bcSide_Type;
    typedef bcFunctionDefault_Type::bcType_Type          bcType_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit OneDimensionalBC( const bcSide_Type& bcSide );

    //! Copy constructor
    /*!
     * @param BC OneDimensionalBC
     */
    explicit OneDimensionalBC( const OneDimensionalBC& BC );

    //! Destructor
    virtual ~OneDimensionalBC() {}

    //@}


    //! @name Methods
    //@{

    //! Apply boundary conditions
    void applyBC( const Real& time, const Real& timeStep, const solution_Type& solution, const fluxPtr_Type& flux, container2D_Type& BC );

    //@}


    //! @name Set Methods
    //@{

    void setType( const bcLine_Type& bcLine, const bcType_Type& bc ) { M_bcType[bcLine] = bc; }

    void setBCFunction( const bcLine_Type& bcLine, const bcFunction_Type& rhs ) { M_bcFunction[bcLine] = rhs; }

    void setInternalFlag( const bool& bcFlag ) { M_isInternal = bcFlag; }

    //@}


    //! @name Get Methods
    //@{

    const bcType_Type& type( const bcLine_Type& bcLine ) { return M_bcType[bcLine]; }

    bcFunction_Type& bcFunction( const bcLine_Type& bcLine ) { return M_bcFunction[bcLine]; }

    const bool& isInternal() { return M_isInternal; }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    OneDimensionalBC& operator=( const OneDimensionalBC& bc );

    //@}


    //! @name Private Methods
    //@{

    //! Compute the matrix and the RHS for the BC 2x2 linear system
    void computeMatrixAndRHS( const Real& time, const Real& timeStep, const fluxPtr_Type& flux, const bcLine_Type& bcLine,
                              const container2D_Type& leftEigenvector1, const container2D_Type& leftEigenvector2,
                              const UInt& dof, Real& rhs );

    //! Solve a 2x2 linear system by the Cramer method (for the boundary systems)
    /*!
     * Matrix A is given by two pairs corresponding to the 2 M_lines.
     * A = [ M_matrixrow_at_line["first" ] ;
     *       M_matrixrow_at_line["second"] ]
     * @return A^{-1} * rhs2d
     */
    container2D_Type solveLinearSystem( const container2D_Type& line1, const container2D_Type& line2, const container2D_Type& rhs ) const;

    //@}

    std::map<bcLine_Type, bcType_Type>            M_bcType;

    bcSide_Type                                   M_bcSide;

    std::map<bcLine_Type, bcFunction_Type>        M_bcFunction;

    bool                                          M_isInternal;

    std::map<bcLine_Type, container2D_Type>       M_bcMatrix;

    container2D_Type                              M_bcRHS;
};

}

#endif //OneDimensionalBC_H
