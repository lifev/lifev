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
 *  @brief File containing the BCInterfaceFunctionSolverDefined class and specializations
 *
 *  @date 23-04-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterfaceFunctionSolverDefined_H
#define BCInterfaceFunctionSolverDefined_H 1

// BCInterface includes
#include <lifev/bc_interface/core/bc/BCInterfaceData.hpp>

#include <lifev/bc_interface/core/function/BCInterfaceFactory.hpp>

namespace LifeV
{

//! BCInterfaceFunctionSolverDefined - Empty class for solver defined specializations.
/*!
 *  @author Cristiano Malossi
 *
 *  This class provides the base interfaces for the implementation of solver defined boundary functions
 *  through template specializations.
 */
template< typename BcHandlerType, typename PhysicalSolverType >
class BCInterfaceFunctionSolverDefined
{
public:

    //! @name Type definitions
    //@{

    typedef BcHandlerType                                                          bcHandler_Type;
    typedef PhysicalSolverType                                                     physicalSolver_Type;

    typedef boost::shared_ptr< physicalSolver_Type >                               physicalSolverPtr_Type;

    typedef BCInterfaceData                                                        data_Type;
    typedef boost::shared_ptr< data_Type >                                         dataPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterfaceFunctionSolverDefined() {}

    //! Destructor
    virtual ~BCInterfaceFunctionSolverDefined() {}

    //@}


    //! @name Methods
    //@{

    //! Copy the stored parameters in the data container
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    template< typename DataPtrType >
    void exportData ( DataPtrType& /*data*/ ) {}

    //! Assign a boundary function to the boundary condition vector base
    /*!
     * @param physicalSolver FSI physical solver,
     * @param base boundary condition base
     */
    template< typename BCBaseType >
    void assignFunction ( BCBaseType& /*base*/ ) {}

    //! Update the solver variables
    void updatePhysicalSolverVariables() {}

    //@}


    //! @name Set methods
    //@{

    //! Set data
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    void setData ( const dataPtr_Type& /*data*/ ) {}

    //! Set the physical solver
    /*!
     * @param physicalSolver physical solver
     */
    void setPhysicalSolver ( const physicalSolverPtr_Type& /*physicalSolver*/ ) {}

    //@}


    //! @name Get methods
    //@{

    //! Detect the correct base type
    /*!
     * @param bcBaseType the type of the base
     */
    baseContainer_Type baseType() const
    {
        return BASEDefault;
    }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    BCInterfaceFunctionSolverDefined ( const BCInterfaceFunctionSolverDefined& function );

    BCInterfaceFunctionSolverDefined& operator= ( const BCInterfaceFunctionSolverDefined& function );

    //@}

};

// ===================================================
// Factory
// ===================================================
//! Factory create function
template< typename BcHandlerType, typename PhysicalSolverType >
inline BCInterfaceFunctionSolverDefined< BcHandlerType, PhysicalSolverType >* createBCInterfaceFunctionSolverDefined()
{
    return new BCInterfaceFunctionSolverDefined< BcHandlerType, PhysicalSolverType > ();
}

} // Namespace LifeV

#endif /* BCInterfaceFunctionSolverDefined_H */
