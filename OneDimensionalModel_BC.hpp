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
#include <lifemc/lifefem/OneDimensionalModel_BCFunction_Default.hpp>
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

    typedef OneDimensionalModel_BCFunction_Default::BCFunction_Type      BCFunction_Type;
    typedef OneDimensionalModel_BCFunction_Default::BCFunction_PtrType   BCFunction_PtrType;

    typedef OneDimensionalModel_BCFunction_Default                       BCFunction_Default_Type;
    typedef boost::shared_ptr< BCFunction_Default_Type >                 BCFunction_Default_PtrType;

    typedef OneDimensionalModel_BCFunction_Default::Flux_PtrType         Flux_PtrType;
    typedef OneDimensionalModel_BCFunction_Default::Source_PtrType       Source_PtrType;
    typedef OneDimensionalModel_BCFunction_Default::Solution_Type        Solution_Type;
    typedef OneDimensionalModel_BCFunction_Default::Solution_PtrType     Solution_PtrType;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    OneDimensionalModel_BC( const OneD_BCSide& side );

    //! Copy constructor
    /*!
     * @param BC OneDimensionalModel_BC
     */
    OneDimensionalModel_BC( const OneDimensionalModel_BC& BC );

    //! Destructor
    ~OneDimensionalModel_BC() {}

    //@}


    //! @name Methods
    //@{

    //! Compute [A,Q] at the boundary
    Container2D_Type Uboundary( const ScalVec& U1, const ScalVec& U2 ) const;

    //! Apply boundary conditions
    void applyBC( const Real& time, const Real& timeStep, const Solution_Type& solution,
                  const Flux_PtrType& flux, Container2D_Type& BC );

    //@}


    //! @name Set Methods
    //@{

    void setType( const OneD_BCLine& line, const OneD_BC& bc );

    void setBCFunction( const OneD_BCLine& line, const BCFunction_Type& rhs );

    void setInternalFlag( const bool& flag );

    //void setMatrixRow( const OneD_BCLine& line, const Container2D_Type& matrixrow );

    //@}


    //! @name Get Methods
    //@{

    const OneD_BC& type( const OneD_BCLine& line );

    BCFunction_Type& BCFunction( const OneD_BCLine& line );

    const bool& isInternal();

    //@}

private:

    //! @name Private Methods
    //@{

    //! Compute the matrix and the RHS for the BC 2x2 linear system
    void computeMatrixAndRHS( const Real& time, const Real& timeStep, const Flux_PtrType& flux, const OneD_BCLine& line,
                              const Container2D_Type& leftEigenvector1, const Container2D_Type& leftEigenvector2,
                              const Container2D_Type& U, const Container2D_Type& W, const UInt& dof, Real& rhs );

    //! Solve a 2x2 linear system by the Cramer method (for the boundary systems)
    /*!
     * Matrix A is given by two pairs corresponding to the 2 M_lines.
     * A = [ M_matrixrow_at_line["first" ] ;
     *       M_matrixrow_at_line["second"] ]
     * @return A^{-1} * rhs2d
     */
    Container2D_Type solveLinearSystem( const Container2D_Type& line1,
                                        const Container2D_Type& line2,
                                        const Container2D_Type& rhs ) const;

    //@}

    std::map<OneD_BCLine, OneD_BC>                M_bcType;

    OneD_BCSide                                   M_bcSide;

    std::map<OneD_BCLine, BCFunction_Type>        M_bcFunction;

    bool                                          M_isInternal;

    std::map<OneD_BCLine, Container2D_Type>       M_bcMatrix;

    Container2D_Type                              M_bcRHS;
};

}

#endif //ONEDIMENSIONALMODEL_BC_H
