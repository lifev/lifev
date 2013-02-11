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
 *  @brief File containing the BCInterfaceFunction class
 *
 *  @date 06-04-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterfaceFunction_H
#define BCInterfaceFunction_H 1

#include <lifev/bc_interface/fem/BCInterfaceData0D.hpp>
#include <lifev/bc_interface/fem/BCInterfaceData1D.hpp>
#include <lifev/bc_interface/fem/BCInterfaceData3D.hpp>

namespace LifeV
{

//! BCInterfaceFunction - Base class for \c BCInterface boundary functions
/*!
 *  @author Cristiano Malossi
 *
 *  This class provides the interface between the \c BCInterface and the boundary functions.
 *
 *  \cond \TODO It may be a good idea (or not) to use this interface also for the \c BCInterfaceFunctionSolverDefined functions \endcond
 */
template< typename PhysicalSolverType >
class BCInterfaceFunction
{
public:

    //! @name Type definitions
    //@{

    typedef PhysicalSolverType                                                                      physicalSolver_Type;

    typedef boost::function<Real ( const Real& ) >                                                   boundaryFunctionTime_Type;
    typedef boost::function<Real ( const Real&, const Real& ) >                                      boundaryFunctionTimeTimeStep_Type;
    typedef boost::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID& ) > boundaryFunctionTimeSpaceID_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    explicit BCInterfaceFunction() {}

    //! Destructor
    virtual ~BCInterfaceFunction() {}

    //@}


    //! @name Methods
    //@{

    //! Assign the function to the base of the 1D \c BCHandler
    /*!
     * @param base base of the 1D boundary condition
     */
    virtual void assignFunction ( OneDFSIFunction& base ) = 0;

    //! Assign the function to the base of the 3D \c BCHandler
    /*!
     * @param base base of the 3D boundary condition
     */
    virtual void assignFunction ( BCFunctionBase& base ) = 0;

    //! Function of time
    /*!
     * @param t time
     * @return boundary condition value
     */
    virtual Real functionTime ( const Real& t ) = 0;

    //! Function of time and time step
    /*!
     * @param t time
     * @param timeStep time step
     * @return boundary condition value
     */
    virtual Real functionTimeTimeStep ( const Real& t, const Real& timeStep ) = 0;

    //! Function of time and space
    /*!
     * @param t time
     * @param x x coordinate
     * @param y y coordinate
     * @param z z coordinate
     * @param id id of the boundary condition (not used)
     * @return boundary condition value
     */
    virtual Real functionTimeSpace ( const Real& t, const Real& x, const Real& y, const Real& z, const ID& /*id*/) = 0;

    //! Function of time and space with ID
    /*!
     * @param t time
     * @param x x coordinate
     * @param y y coordinate
     * @param z z coordinate
     * @param id id of the boundary condition
     * @return boundary condition value
     */
    virtual Real functionTimeSpaceID ( const Real& t, const Real& x, const Real& y, const Real& z, const ID& id ) = 0;

    //@}


    //! @name Set Methods
    //@{

    //! Set data for 0D boundary conditions
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    virtual void setData ( const BCInterfaceData0D& data ) = 0;


    //! Set data for 1D boundary conditions
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    virtual void setData ( const BCInterfaceData1D& data ) = 0;

    //! Set data for 3D boundary conditions
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    virtual void setData ( const BCInterfaceData3D& data ) = 0;

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    BCInterfaceFunction ( const BCInterfaceFunction& function );

    BCInterfaceFunction& operator= ( const BCInterfaceFunction& function );

    //@}
};

} // Namespace LifeV

#endif /* BCInterfaceFunction_H */
