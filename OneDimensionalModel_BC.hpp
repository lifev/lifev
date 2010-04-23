//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief File containing a class for the boundary conditions of the 1D model.
 *
 *  @version 1.0
 *  @author Lucia Mirabella <lucia@mathcs.emory.edu>
 *  @author Tiziano Passerini <tiziano@mathcs.emory.edu>
 *  @date 01-28-2006
 *
 *  @version 2.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 20-04-2010
 */

#ifndef ONEDIMENSIONALMODEL_BC_H
#define ONEDIMENSIONALMODEL_BC_H

// LIFEV - MATHCARD
#include <lifemc/lifefem/OneDimensionalModel_BCFunction.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Definitions.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Data.hpp>

namespace LifeV {

//! OneDimensionalModel_BC - Class featuring methods to handle boundary conditions.
/*!
 *  @author Lucia Mirabella, Tiziano Passerini
 */
class OneDimensionalModel_BC
{
public:

    //! @name Type definitions
    //@{

    typedef boost::shared_ptr<OneDimensionalModel_BCFunction>    OneDimensionalModel_BCFunction_PtrType;

    typedef OneDimensionalModel_BCFunction::Flux_Type            Flux_Type;
    typedef OneDimensionalModel_BCFunction::Flux_PtrType         Flux_PtrType;

    typedef OneDimensionalModel_BCFunction::Data_Type            Data_Type;
    typedef OneDimensionalModel_BCFunction::Mesh_Type            Mesh_Type;

    typedef OneDimensionalModel_BCFunction::FESpace_Type         FESpace_Type;

    typedef OneDimensionalModel_BCFunction::LinearSolver_Type    LinearSolver_Type;
    typedef OneDimensionalModel_BCFunction::Vector_Type          Vector_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    OneDimensionalModel_BC( const std::vector<Vector_Type>& U_thistime,
                            const Flux_PtrType              fluxFun,
                            const Real&                     dimDof,
                            const std::string&              side );

    //! Destructor
    ~OneDimensionalModel_BC() {}

    //@}


    //! @name Methods
    //@{

    //! Compute [A,Q] at the boundary
    Vec2D Uboundary(const ScalVec& U1, const ScalVec& U2) const;

    //! Apply boundary conditions
    void applyBC( const Real& time_val, Vec2D& BC_dir );

    //@}


    //! @name Get Methods
    //@{

    //! Return the boundary Dof
    UInt boundaryDof() const;

    //! careful I need them to be setters! do not remove the reference
    OneDimensionalModel_BCFunction_PtrType& rhs( const std::string& line );

    std::string& variable( const std::string& line );

    Vec2D& matrixrow( const std::string& line );

    bool& isInternal();

    //@}

protected:

    //! @name Protected Methods
    //@{

    //! Impose the chosen boundary condition
    void compute_resBC( const Real& time_val );

    void compute_resBC_line( std::string line, Vec2D left_eigvec,
                             Vec2D U, Vec2D W, Real& rhs );

    //! Solve a 2x2 linear system by the Cramer method (for the boundary systems)
    /*!
     * Matrix A is given by two pairs corresponding to the 2 M_lines.
     * A = [ M_matrixrow_at_line["first" ] ;
     *       M_matrixrow_at_line["second"] ]
     * @return A^{-1} * rhs2d
     */
    Vec2D _solveLinearSyst2x2( const Vec2D& line1, const Vec2D& line2, const Vec2D& rhs2d ) const;

    //@}

    bool                                          M_isInternal;

    std::map<std::string, std::string>            M_variable_at_line;

    std::map<std::string, Vec2D>                  M_matrixrow_at_line;

    std::map<std::string, OneDimensionalModel_BCFunction_PtrType> M_rhs_at_line;

    std::map<std::string, OneDBCStringValue>      M_OneDimensionalModel_BCMapStringValues;

    //! Result of the 2x2 linear system to be solved at each side
    Vec2D                                         M_resBC;

    //! Reference to the solver current unknowns (U)
    const std::vector<Vector_Type>&               M_U_thistime;

    //! Boundary Dof (right or left)
    UInt                                          M_boundaryDof;

    //! Reference to the solver non linear flux functions
    Flux_PtrType                                  M_fluxFun;
};

}

#endif //ONEDIMENSIONALMODEL_BC_H
