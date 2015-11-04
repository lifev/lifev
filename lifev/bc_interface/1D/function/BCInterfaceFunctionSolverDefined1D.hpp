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

#ifndef BCInterfaceFunctionSolverDefined1D_H
#define BCInterfaceFunctionSolverDefined1D_H 1

// OneDFSI includes
#include <lifev/one_d_fsi/solver/OneDFSISolver.hpp>

// BCInterface includes
#include <lifev/bc_interface/1D/bc/BCInterfaceData1D.hpp>
#include <lifev/bc_interface/core/function/BCInterfaceFunctionSolverDefined.hpp>

namespace LifeV
{

//! BCInterfaceFunctionSolverDefined - Template specialization of \c BCInterfaceFunctionSolverDefined for 1D problems
/*!
 *  @author Cristiano Malossi
 *
 *  The BCInterfaceFunctionSolverDefined class provides a general interface between the
 *  \c BCInterface1D and the solver defined boundary conditions of the \c OneDFSISolver.
 *
 *  <b>DETAILS:</b> <BR>
 *  The list of available conditions is described by the \c solverDefinedFunctions enum type.
 *
 *  They are:
 *  <ol>
 *      <li> Riemann;
 *      <li> Compatibility;
 *      <li> Absorbing;
 *      <li> Resistance.
 *  </ol>
 */
template< >
class BCInterfaceFunctionSolverDefined< OneDFSIBCHandler, OneDFSISolver >
{
public:

    //! @name Type definitions
    //@{

    typedef OneDFSIBCHandler                                                       bcHandler_Type;
    typedef OneDFSISolver                                                          physicalSolver_Type;
    typedef std::shared_ptr< physicalSolver_Type >                               physicalSolverPtr_Type;

    typedef bcHandler_Type::bc_Type                                                bc_Type;
    typedef bc_Type::bcFunctionSolverDefinedPtr_Type                               bcFunctionSolverDefinedPtr_Type;

    typedef bc_Type::vectorPtrContainer_Type                                       vectorPtrContainer_Type;

    typedef BCInterfaceData1D                                                      data_Type;
    typedef std::shared_ptr< data_Type >                                         dataPtr_Type;

    typedef bc_Type::fluxPtr_Type                                                  fluxPtr_Type;
    typedef bc_Type::sourcePtr_Type                                                sourcePtr_Type;
    typedef bc_Type::solutionPtr_Type                                              solutionPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterfaceFunctionSolverDefined();

    //! Destructor
    virtual ~BCInterfaceFunctionSolverDefined() {}

    //@}


    //! @name Methods
    //@{

    //! Assign the function to the base
    /*!
     * @param base base of the bc
     */
    void assignFunction ( OneDFSIFunction& base );

    //! Update the solver variables
    void updatePhysicalSolverVariables() {}

    //@}


    //! @name Set Methods
    //@{

    //! Set data
    /*!
     * @param data BC data loaded from GetPot file
     */
    void setData ( const dataPtr_Type& data );

    //! Set flux and source
    /*!
     * @param flux flux object of the 1D model
     * @param source source object of the 1D model
     */
    void setFluxSource ( const fluxPtr_Type& flux, const sourcePtr_Type& source )
    {
        M_function->setFluxSource ( flux, source );
    }

    //! Set solution
    /*!
     * @param solution solution container of the 1D model
     */
    void setSolution ( const solutionPtr_Type& solution )
    {
        M_function->setSolution ( solution );
    }

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
        return BASEFunction1D;
    }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    BCInterfaceFunctionSolverDefined ( const BCInterfaceFunctionSolverDefined& function );

    BCInterfaceFunctionSolverDefined& operator= ( const BCInterfaceFunctionSolverDefined& function );

    //@}

    enum solverDefinedFunctions
    {
        Riemann,
        Compatibility,
        Absorbing,
        Resistance
    };

    solverDefinedFunctions           M_defaultFunction;
    bcFunctionSolverDefinedPtr_Type  M_function;
};

} // Namespace LifeV

#endif /* BCInterfaceFunctionSolverDefined1D_H */
