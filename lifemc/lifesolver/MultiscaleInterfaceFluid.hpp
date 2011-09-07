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
 *  @brief File containing the Multiscale Interface for Fluid problems
 *
 *  @date 31-03-2011
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleInterfaceFluid_H
#define MultiscaleInterfaceFluid_H 1

#include <lifemc/lifesolver/MultiscaleDefinitions.hpp>

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleInterfaceFluid - The multiscale interface for fluid problems
/*!
 *  @author Cristiano Malossi
 *
 *  The MultiscaleInterfaceFluid class provides a general and abstract interface for the
 *  coupling of fluid problems.
 */
class MultiscaleInterfaceFluid
{
public:

    //! @name Type definitions
    //@{

    typedef boost::function< Real ( const Real&, const Real&, const Real&, const Real&, const ID& ) > function_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! The main constructor.
    explicit MultiscaleInterfaceFluid() {}

    //! Destructor
    virtual ~MultiscaleInterfaceFluid() {}

    //@}


    //! @name Multiscale Interface Virtual Methods
    //@{

    //! Impose the flow rate on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @param function boundary condition function
     */
    virtual void imposeBoundaryFlowRate( const bcFlag_Type& flag, const function_Type& function ) = 0;

    //! Impose the flow rate on a specific boundary face of the model as a valve
    /*!
     * @param flag flag of the boundary face
     * @param function boundary condition function
     */
    virtual void imposeBoundaryFlowRateAsValve( const bcFlag_Type& flag, const function_Type& function, const bool& valveIsOpen ) = 0;

    //! Impose the integral of the normal stress on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @param function boundary condition function
     */
    virtual void imposeBoundaryStress( const bcFlag_Type& flag, const function_Type& function ) = 0;

    //! Get the flow rate on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return flow rate value
     */
    virtual Real boundaryFlowRate( const bcFlag_Type& flag ) const = 0;

    //! Get the integral of the normal stress (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @return stress value
     */
    virtual Real boundaryStress( const bcFlag_Type& flag ) const = 0;

    //! Get the variation of the flow rate (on a specific boundary face) using the linear model
    /*!
     * @param flag flag of the boundary face on which quantity should be computed
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the flow rate
     */
    virtual Real boundaryDeltaFlowRate( const bcFlag_Type& flag, bool& solveLinearSystem ) = 0;

    //! Get the variation of the integral of the normal stress (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the stress
     */
    virtual Real boundaryDeltaStress( const bcFlag_Type& flag, bool& solveLinearSystem ) = 0;

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleInterfaceFluid( const MultiscaleInterfaceFluid& interface );

    MultiscaleInterfaceFluid& operator=( const MultiscaleInterfaceFluid& interface );

    //@}
};

} // Namespace multiscale
} // Namespace LifeV

#endif /* MultiscaleInterfaceFluid_H */
