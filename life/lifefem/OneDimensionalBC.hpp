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

#include <life/lifefunctions/OneDimensionalBCFunctionSolverDefined.hpp>
#include <life/lifesolver/OneDimensionalData.hpp>

namespace LifeV
{

//! OneDimensionalBC - Class featuring methods to handle boundary conditions.
/*!
 *  @author Lucia Mirabella, Tiziano Passerini, Cristiano Malossi
 *
 *  We need to impose 2 boundary condition on each side of the 1D segment.
 *  These class stores the boundary conditions on one side.
 *
 *  \cond \TODO Improve the description of this class. \endcond
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
    typedef bcFunctionDefault_Type::vectorPtrContainer_Type vectorPtrContainer_Type;

    typedef bcFunctionDefault_Type::vector_Type          vector_Type;
    typedef bcFunctionDefault_Type::matrix_Type          matrix_Type;

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
     * @param bc OneDimensionalBC
     */
    explicit OneDimensionalBC( const OneDimensionalBC& bc );

    //! Destructor
    virtual ~OneDimensionalBC() {}

    //@}


    //! @name Methods
    //@{

    //! Apply boundary conditions to the rhs of the Taylor-Galerkin problem
    /*!
     *  @param time the current time.
     *  @param timeStep the time step.
     *  @param solution the solution container.
     *  @param fluxPtr pointer to the flux class.
     *  @param rhs the rhs of the Taylor-Galerking problem.
     */
    void applyBC( const Real& time, const Real& timeStep, const solution_Type& solution,
                  const fluxPtr_Type& fluxPtr, vectorPtrContainer_Type& rhs );

    //! Apply boundary conditions to the rhs of the viscoelastic problem
    /*!
     *  @param fluxPtr pointer to the flux class.
     *  @param matrix the matrix of the viscoelastic problem.
     *  @param rhs the rhs of the viscoelastic problem.
     */
    void applyViscoelasticBC( const fluxPtr_Type& fluxPtr, matrix_Type& matrix, vector_Type& rhs );

    //@}


    //! @name Set Methods
    //@{

    //! Set the type of boundary condition
    /*!
     *  @param bcLine the line of the boundary condition (first or second).
     *  @param bcType the type of the boundary condition (\f$Q\f$, \f$A\f$, \f$P\f$, \f$S\f$, \f$W_1\f$, \f$W_2\f$).
     */
    void setType( const bcLine_Type& bcLine, const bcType_Type& bcType ) { M_bcType[bcLine] = bcType; }

    //! Set the boundary condition function
    /*!
     *  @param bcLine the line of the boundary condition (first or second).
     *  @param bcFunction the boundary condition function.
     */
    void setBCFunction( const bcLine_Type& bcLine, const bcFunction_Type& rhs ) { M_bcFunction[bcLine] = rhs; }

    //@}


    //! @name Get Methods
    //@{

    //! Get the type of boundary condition
    /*!
     *  @param bcLine the line of the boundary condition (first or second).
     *  @return the type of the boundary condition (\f$Q\f$, \f$A\f$, \f$P\f$, \f$S\f$, \f$W_1\f$, \f$W_2\f$).
     */
    const bcType_Type& type( const bcLine_Type& bcLine ) const { return M_bcType.find( bcLine )->second; }

    //! Get the boundary condition function
    /*!
     *  @param bcLine the line of the boundary condition (first or second).
     *  @return the boundary condition function.
     */
    const bcFunction_Type& bcFunction( const bcLine_Type& bcLine ) const { return M_bcFunction.find( bcLine )->second; }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    OneDimensionalBC& operator=( const OneDimensionalBC& bc );

    //@}


    //! @name Private Methods
    //@{

    //! Compute the matrix and the RHS for the BC 2x2 linear system
    /*!
     *  @param time the current time.
     *  @param timeStep the time step.
     *  @param fluxPtr pointer to the flux class.
     *  @param bcLine the line of the boundary condition (first or second).
     *  @param leftEigenvector1 first line of the left eigenvector matrix.
     *  @param leftEigenvector2 second line of the left eigenvector matrix.
     *  @param dof degree of freedom of the boundary condition.
     *  @param bcMatrix the 2x2 matrix problem for the boundary condition computation.
     *  @param bcRHS the rhs of the 2x2 problem for the boundary condition computation.
     */
    void computeMatrixAndRHS( const Real& time, const Real& timeStep, const fluxPtr_Type& fluxPtr, const bcLine_Type& bcLine,
                              const container2D_Type& leftEigenvector1, const container2D_Type& leftEigenvector2,
                              const UInt& dof, std::map<bcLine_Type, container2D_Type>& bcMatrix, Real& bcRHS );

    //! Solve a 2x2 linear system by the Cramer method (for the boundary conditions)
    /*!
     * Matrix A is given by two pairs corresponding to the 2 lines.
     * @param line1 first line of the 2x2 matrix.
     * @param line2 second line of the 2x2 matrix.
     * @param rhs rhs of the 2x2 system.
     * @return solution
     */
    container2D_Type solveLinearSystem( const container2D_Type& line1,
                                        const container2D_Type& line2, const container2D_Type& rhs ) const;

    //@}

    std::map<bcLine_Type, bcType_Type>            M_bcType;

    bcSide_Type                                   M_bcSide;

    std::map<bcLine_Type, bcFunction_Type>        M_bcFunction;
};

}

#endif //OneDimensionalBC_H
