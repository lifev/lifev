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

#ifndef OneDimensionalBCHandler_H
#define OneDimensionalBCHandler_H

// LIFEV - MATHCARD
#include <lifemc/lifesolver/OneDimensionalModel_Definitions.hpp>
#include <lifemc/lifefem/OneDimensionalModel_BC.hpp>

namespace LifeV
{

//! OneDimensionalBCHandler - Class featuring methods to handle boundary conditions.
/*!
 *  @author Lucia Mirabella, Tiziano Passerini
 */
class OneDimensionalBCHandler
{
public:

    //! @name Type definitions
    //@{

    typedef OneDimensionalBC                    bc_Type;
    typedef boost::shared_ptr< bc_Type >        bcPtr_Type;

    typedef bc_Type::bcFunction_Type            bcFunction_Type;
    typedef bc_Type::bcFunctionPtr_Type         bcFunctionPtr_Type;
    typedef bc_Type::bcFunctionDefault_Type     bcFunctionDefault_Type;
    typedef bc_Type::bcFunctionDefaultPtr_Type  bcFunctionDefaultPtr_Type;

    typedef bc_Type::fluxPtr_Type               fluxPtr_Type;
    typedef bc_Type::sourcePtr_Type             sourcePtr_Type;
    typedef bc_Type::solution_Type              solution_Type;
    typedef bc_Type::solutionPtr_Type           solutionPtr_Type;

    typedef bc_Type::container2D_Type           container2D_Type;

    typedef bc_Type::bcLine_Type                bcLine_Type;
    typedef bc_Type::bcSide_Type                bcSide_Type;
    typedef bc_Type::bcType_Type                bcType_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit OneDimensionalBCHandler();

    //! Copy constructor
    /*!
     * @param BCH OneDimensionalBCHandler
     */
    explicit OneDimensionalBCHandler( const OneDimensionalBCHandler& BCH );

    //! Destructor
    virtual ~OneDimensionalBCHandler() {}

    //@}


    //! @name Methods
    //@{

    //! Apply boundary conditions
    void applyBC( const Real& time, const Real& timeStep, const solution_Type& solution, const fluxPtr_Type& flux, container2D_Type& leftBC, container2D_Type& rightBC );

    //@}


    //! @name Set Methods
    //@{

    void setBC( const bcSide_Type& bcSide, const bcLine_Type& line, const bcType_Type& bcType, const bcFunction_Type& BCfunction );

    void setDefaultBC();

    void setSolution( const solutionPtr_Type& solution );

    void setFluxSource( const fluxPtr_Type& flux, const sourcePtr_Type& source );

    void setInternalNode( const bcSide_Type& bcSide ) { M_boundary[bcSide]->setInternalFlag( true ); }

    //@}


    //! @name Get Methods
    //@{

    const bcPtr_Type& bc( const bcSide_Type& bcSide ) { return M_boundary[bcSide]; }

    const bool& bcReady( const bcSide_Type& bcSide, const bcLine_Type& line ) { return M_boundarySet[bcSide][line]; }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    OneDimensionalBCHandler& operator=( const OneDimensionalBCHandler& bcHandler );

    //@}

    std::map< bcSide_Type, bcPtr_Type >                      M_boundary;
    std::map< bcSide_Type, std::map< bcLine_Type, bool > >   M_boundarySet;

    std::vector < bcFunctionDefaultPtr_Type >                M_defaultFunctions;
};

}

#endif //OneDimensionalBCHandler_H
