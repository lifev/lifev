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
 *  @brief File containing the interface class for the boundary function of 1D model.
 *
 *  @version 1.0
 *  @date 01-08-2006
 *  @author Lucia Mirabella  <lucia.mirabella@gmail.com>
 *  @author Tiziano Passerini <tiziano.passerini@gmail.com>
 *
 *  @version 2.0
 *  @date 20-04-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef OneDimensionalBCFunction_H
#define OneDimensionalBCFunction_H

#include <life/lifesolver/OneDimensionalDefinitions.hpp>

namespace LifeV
{

//! OneDimensionalBCFunction - Base class for 1D BC Functions.
/*!
 *  @author Lucia Mirabella, Tiziano Passerini, Cristiano Malossi
 *
 *  The 1D boundary condition function is evaluated as a function of the current time and of the time step.
 */
class OneDimensionalBCFunction
{
public:

    //! @name Type definitions and Enumerators
    //@{

    /*! @typedef function_Type */
    //! Type definition for the 1D boundary function
    typedef boost::function<Real ( const Real&, const Real&  )> function_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    explicit OneDimensionalBCFunction() : M_function() {}

    //! Constructor by function
    /*!
     *  @param function the user defined function
     */
    explicit OneDimensionalBCFunction( const function_Type& function ) : M_function( function ) {}

    //! Copy constructor
    /*!
     *  @param bcFunction OneDimensionalBCFunction
     */
    OneDimensionalBCFunction( const OneDimensionalBCFunction& bcFunction ) : M_function  ( bcFunction.M_function ) {}

    //! Destructor
    virtual ~OneDimensionalBCFunction() {}

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     *  @param bcFunction OneDimensionalBCFunction
     *  @return reference to a copy of the class
     */
    OneDimensionalBCFunction& operator=( const OneDimensionalBCFunction& bcFunction )
    {
        if( this != &bcFunction )
            M_function = bcFunction.M_function;

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
    Real operator()( const Real& time, const Real& timeStep = 0. ) const { return M_function( time, timeStep ); }

    //@}


    //! @name Set Methods
    //@{

    //! Set the function
    /*!
      @param function the user defined function
    */
    void setFunction( const function_Type& function ) { M_function = function; }

    //@}


    //! @name Get Methods
    //@{

    //! Get the function
    /*!
      @return the user defined function
    */
    const function_Type& function() const { return M_function; }

    //@}

private:

    function_Type M_function;
};

/*
//! Factory create function
inline OneDimensionalBCFunction*
Create_OneDimensionalModel_BCFunction( const OneDimensionalBCFunction* bcFunction )
{
    return new OneDimensionalBCFunction( (const OneDimensionalBCFunction&)* bcFunction );
}

namespace
{
    static bool registerOneD_BCFunction = FactoryClone_OneDimensionalModel_BCFunction::instance().registerProduct( typeid(OneDimensionalBCFunction), &Create_OneDimensionalModel_BCFunction );
}
*/

}

#endif // OneDimensionalBCFunction_H
