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
 *  @brief File containing a class for the boundary conditions handling of the 1D model.
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
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef OneDFSIBCHandler_H
#define OneDFSIBCHandler_H

#include <lifev/one_d_fsi/solver/OneDFSIDefinitions.hpp>
#include <lifev/one_d_fsi/fem/OneDFSIBC.hpp>

namespace LifeV
{

//! OneDFSIBCHandler - Class featuring methods to handle boundary conditions.
/*!
 *  @author Lucia Mirabella, Tiziano Passerini, Cristiano Malossi
 *  @see Equations and networks of 1-D models \cite FormaggiaLamponi2003
 *  @see Geometrical multiscale coupling of 1-D models \cite Malossi2011Algorithms \cite Malossi2011Algorithms1D
 *
 *  We need to impose 2 boundary condition on each side of the 1D segment.
 *  These boundary conditions are stored in \c OneDFSIBC objects.
 *
 *  \cond \TODO Improve the description of BC. \endcond
 */
class OneDFSIBCHandler
{
public:

    //! @name Type definitions
    //@{

    typedef OneDFSIBC                           bc_Type;
    typedef boost::shared_ptr< bc_Type >        bcPtr_Type;

    typedef bc_Type::bcFunction_Type            bcFunction_Type;
    typedef bc_Type::bcFunctionPtr_Type         bcFunctionPtr_Type;
    typedef bc_Type::bcFunctionSolverDefined_Type     bcFunctionSolverDefined_Type;
    typedef bc_Type::bcFunctionSolverDefinedPtr_Type  bcFunctionSolverDefinedPtr_Type;

    typedef bc_Type::fluxPtr_Type               fluxPtr_Type;
    typedef bc_Type::sourcePtr_Type             sourcePtr_Type;
    typedef bc_Type::solution_Type              solution_Type;
    typedef bc_Type::solutionPtr_Type           solutionPtr_Type;

    typedef bc_Type::vectorPtrContainer_Type    vectorPtrContainer_Type;

    typedef bc_Type::vector_Type                vector_Type;
    typedef bc_Type::matrix_Type                matrix_Type;

    typedef bc_Type::bcLine_Type                bcLine_Type;
    typedef bc_Type::bcSide_Type                bcSide_Type;
    typedef bc_Type::bcType_Type                bcType_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    explicit OneDFSIBCHandler();

    //! Copy constructor
    /*!
     * @param bcHandler OneDFSIBCHandler
     */
    explicit OneDFSIBCHandler ( const OneDFSIBCHandler& bcHandler );

    //! Destructor
    virtual ~OneDFSIBCHandler() {}

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
    void applyBC ( const Real& time, const Real& timeStep, const solution_Type& solution,
                   const fluxPtr_Type& fluxPtr, vectorPtrContainer_Type& rhs );

    //! Apply boundary conditions to the rhs of the viscoelastic problem
    /*!
     *  @param fluxPtr pointer to the flux class.
     *  @param matrix the matrix of the viscoelastic problem.
     *  @param rhs the rhs of the viscoelastic problem.
     */
    void applyViscoelasticBC (const fluxPtr_Type& fluxPtr, matrix_Type& matrix, vector_Type& rhs );

    //@}


    //! @name Set Methods
    //@{

    //! Set a boundary condition
    /*!
     *  @param bcSide the side of the boundary condition (left or right).
     *  @param bcLine the line of the boundary condition (first or second).
     *  @param bcType the type of the boundary condition (\f$Q\f$, \f$A\f$, \f$P\f$, \f$S\f$, \f$W_1\f$, \f$W_2\f$).
     *  @param bcFunction the boundary condition function.
     */
    void setBC ( const bcSide_Type& bcSide, const bcLine_Type& bcLine, const bcType_Type& bcType, const bcFunction_Type& bcFunction );

    //! Set the default boundary conditions
    /*!
     *  This is done only for the boundary conditions that have not been set yet.
     */
    void setDefaultBC();

    //! Set the flux and the source classes for the problem
    /*!
     *  @param fluxPtr pointer to the flux term of the problem.
     *  @param source pointer to the source term of the problem.
     */
    void setFluxSource ( const fluxPtr_Type& fluxPtr, const sourcePtr_Type& sourcePtr );

    //! Set the solution of the problem
    /*!
     *  @param solutionPtr pointer to the solution of the problem.
     */
    void setSolution ( const solutionPtr_Type& solutionPtr );

    //@}


    //! @name Get Methods
    //@{

    //! Get a specific boundary condition
    /*!
     *  @param bcSide the side of the boundary condition (left or right).
     *  @return the pointer to the boundary conditions on a specific side
     */
    const bcPtr_Type& bc ( const bcSide_Type& bcSide ) const
    {
        return M_boundary.find ( bcSide )->second;
    }

    //! Return true if the boundary condition has been set
    /*!
     *  @param bcSide the side of the boundary condition (left or right).
     *  @param bcLine the line of the boundary condition (first or second).
     *  @return true if the boundary condition has been set, false otherwise.
     */
    const bool& bcReady ( const bcSide_Type& bcSide, const bcLine_Type& bcLine ) const
    {
        return M_boundarySet.find ( bcSide )->second.find ( bcLine )->second;
    }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    OneDFSIBCHandler& operator= ( const OneDFSIBCHandler& bcHandler );

    //@}

    std::map< bcSide_Type, bcPtr_Type >                      M_boundary;
    std::map< bcSide_Type, std::map< bcLine_Type, bool > >   M_boundarySet;

    std::vector < bcFunctionSolverDefinedPtr_Type >                M_defaultFunctions;
};

}

#endif //OneDFSIBCHandler_H
