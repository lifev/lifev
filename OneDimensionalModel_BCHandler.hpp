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
 *  @author Lucia Mirabella <lucia@mathcs.emory.edu>
 *  @author Tiziano Passerini <tiziano@mathcs.emory.edu>
 *  @date 01-28-2006
 *
 *  @version 2.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 20-04-2010
 */

#ifndef ONEDIMENSIONALMODEL_BCHANDLER_H
#define ONEDIMENSIONALMODEL_BCHANDLER_H

// LIFEV - MATHCARD
#include <lifemc/lifesolver/OneDimensionalModel_Definitions.hpp>
#include <lifemc/lifefem/OneDimensionalModel_BC.hpp>

namespace LifeV
{

//! OneDimensionalModel_BCHandler - Class featuring methods to handle boundary conditions.
/*!
 *  @author Lucia Mirabella, Tiziano Passerini
 */
class OneDimensionalModel_BCHandler
{
public:

    //! @name Type definitions
    //@{

    typedef OneDimensionalModel_BC              BC_Type;
    typedef boost::shared_ptr< BC_Type >        BC_PtrType;

    typedef BC_Type::BCFunction_Type            BCFunction_Type;
    typedef BC_Type::BCFunction_PtrType         BCFunction_PtrType;
    typedef BC_Type::BCFunction_Default_Type    BCFunction_Default_Type;
    typedef BC_Type::BCFunction_Default_PtrType BCFunction_Default_PtrType;

    typedef BC_Type::Flux_PtrType               Flux_PtrType;
    typedef BC_Type::Source_PtrType             Source_PtrType;
    typedef BC_Type::Solution_Type              Solution_Type;
    typedef BC_Type::Solution_PtrType           Solution_PtrType;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    OneDimensionalModel_BCHandler();

    //! Copy constructor
    /*!
     * @param BCH OneDimensionalModel_BCHandler
     */
    OneDimensionalModel_BCHandler( const OneDimensionalModel_BCHandler& BCH );

    //! Destructor
    virtual ~OneDimensionalModel_BCHandler() {}

    //@}


    //! @name Methods
    //@{

    //! Apply boundary conditions
    void applyBC ( const Real&             time,
                   const Real&             timeStep,
                   const Solution_Type&    solution,
                   const Flux_PtrType&     flux,
                   container2D_Type& leftBC,
                   container2D_Type& rightBC );

    //@}


    //! @name Set Methods
    //@{

    void setBC( const OneD_BCSide& side, const OneD_BCLine& line,
                const OneD_BC& bcType,   const BCFunction_Type& BCfunction );

    void setDefaultBC();

    void setSolution( const Solution_PtrType& solution );

    void setFluxSource( const Flux_PtrType& flux, const Source_PtrType& source );

    void setInternalNode( const OneD_BCSide& side );

    //@}


    //! @name Get Methods
    //@{

    const BC_PtrType& BC( const OneD_BCSide& side );

    const bool& BCReady( const OneD_BCSide& side, const OneD_BCLine& line );

    //@}

private:

    std::map< OneD_BCSide, BC_PtrType >                      M_boundary;
    std::map< OneD_BCSide, std::map< OneD_BCLine, bool > >   M_boundarySet;

    std::vector < BCFunction_Default_PtrType >               M_defaultFunctions;
};

}

#endif //ONEDIMENSIONALMODEL_BCHANDLER_H
