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
#include <lifev/bc_interface/0D/bc/BCInterfaceData0D.hpp>
#include <lifev/bc_interface/1D/bc/BCInterfaceData1D.hpp>
#include <lifev/bc_interface/3D/bc/BCInterfaceData3D.hpp>

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
template< class PhysicalSolverType >
class BCInterfaceFunctionSolverDefined
{
public:

    //! @name Type definitions
    //@{

    typedef PhysicalSolverType                                    physicalSolver_Type;
    typedef boost::shared_ptr< physicalSolver_Type >              physicalSolverPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    explicit BCInterfaceFunctionSolverDefined() {}

    virtual ~BCInterfaceFunctionSolverDefined() {}

    //@}


    //! @name Methods
    //@{

    //! Copy the stored parameters in the 0D data container
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    void exportData ( BCInterfaceData0D& /*data*/ ) {}

    //! Copy the stored parameters in the 1D data container
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    void exportData ( BCInterfaceData1D& /*data*/ ) {}

    //! Copy the stored parameters in the 3D data container
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    void exportData ( BCInterfaceData3D& /*data*/ ) {}

    //! Detect the correct base type
    /*!
     * @param bcBaseType the type of the base
     */
    baseContainer_Type baseType() const
    {
        return BASEDefault;
    }

    //! Assign a boundary function to the boundary condition vector base
    /*!
     * @param physicalSolver FSI physical solver,
     * @param base boundary condition base
     */
    template< class BCBaseType >
    void assignFunction ( BCBaseType& /*base*/ ) {}

    //! Update the solver variables
    void updatePhysicalSolverVariables() {}


    //@}


    //! @name Set Methods
    //@{

    //! Set data for 0D boundary conditions
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    void setData ( const BCInterfaceData0D& /*data*/ ) {}

    //! Set data for 1D boundary conditions
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    void setData ( const BCInterfaceData1D& /*data*/ ) {}

    //! Set data for 3D boundary conditions
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    void setData ( const BCInterfaceData3D& /*data*/ ) {}

    //! Set the physical solver
    /*!
     * @param physicalSolver physical solver
     */
    void setPhysicalSolver ( const boost::shared_ptr< PhysicalSolverType >& /*physicalSolver*/ ) {}

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    BCInterfaceFunctionSolverDefined ( const BCInterfaceFunctionSolverDefined& function);

    BCInterfaceFunctionSolverDefined& operator= ( const BCInterfaceFunctionSolverDefined& function );

    //@}

};

// ===================================================
// Factory
// ===================================================
//! Factory create function
template< typename PhysicalSolverType >
inline BCInterfaceFunctionSolverDefined< PhysicalSolverType >* createBCInterfaceFunctionSolverDefined()
{
    return new BCInterfaceFunctionSolverDefined< PhysicalSolverType > ();
}

} // Namespace LifeV

#endif /* BCInterfaceFunctionSolverDefined_H */
