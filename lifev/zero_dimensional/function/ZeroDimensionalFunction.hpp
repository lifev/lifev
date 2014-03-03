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
 *  @brief File containing the zero dimensional bc function
 *
 *  @date 13-03-2013
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef ZeroDimensionalFunction_H
#define ZeroDimensionalFunction_H

#include <lifev/zero_dimensional/solver/ZeroDimensionalDefinitions.hpp>

namespace LifeV
{

//! ZeroDimensionalFunction - A boundary conditions function for zero-dimensional models.
/*!
 *  @author Cristiano Malossi
 *
 *  This simple class handles the boundary condition functions for zero-dimensional models.
 */
class ZeroDimensionalFunction
{
public:

    //! @name Type definitions and Enumerators
    //@{

    /*! @typedef function_Type */
    //! Type definition for the 0D boundary function
    typedef boost::function<Real ( const Real& ) > function_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    explicit ZeroDimensionalFunction() : M_function() {}

    //! Constructor by function
    /*!
     *  @param function the user defined function
     */
    explicit ZeroDimensionalFunction ( const function_Type& function ) : M_function ( function ) {}

    //! Copy constructor
    /*!
     *  @param bcFunction ZeroDimensionalFunction
     */
    ZeroDimensionalFunction ( const ZeroDimensionalFunction& bcFunction ) : M_function  ( bcFunction.M_function ) {}

    //! Destructor
    virtual ~ZeroDimensionalFunction() {}

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     *  @param bcFunction ZeroDimensionalFunction
     *  @return reference to a copy of the class
     */
    ZeroDimensionalFunction& operator= ( const ZeroDimensionalFunction& bcFunction )
    {
        if ( this != &bcFunction )
        {
            M_function = bcFunction.M_function;
        }

        return *this;
    }

    //! Operator()
    /*!
     *  Evaluate the function.
     *
     *  @param time the current time.
     *  @param timeStep the time step.
     *  @return the value of the function.
     */
    Real operator() ( const Real& time ) const
    {
        return M_function ( time );
    }

    //@}


    //! @name Set Methods
    //@{

    //! Set the function
    /*!
      @param function the user defined function
    */
    void setFunction ( const function_Type& function )
    {
        M_function = function;
    }

    //@}


    //! @name Get Methods
    //@{

    //! Get the function
    /*!
      @return the user defined function
    */
    const function_Type& function() const
    {
        return M_function;
    }

    //@}

private:

    function_Type M_function;
};

}

#endif // ZeroDimensionalFunction_H
