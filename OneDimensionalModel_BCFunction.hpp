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
 *  @author Lucia Mirabella  <lucia.mirabella@gmail.com>
 *  @date 01-08-2006
 *
 *  @version 2.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 20-04-2010
 *
 */

#ifndef ONEDIMENSIONALMODEL_BCFUNCTION_H
#define ONEDIMENSIONALMODEL_BCFUNCTION_H

// LIFEV - MATHCARD
#include <lifemc/lifesolver/OneDimensionalModel_Definitions.hpp>

namespace LifeV
{

//! OneDimensionalModel_BCFunction - Base class for One Dimensional BC Functions.
/*!
 *  @author Lucia Mirabella
 */
class OneDimensionalModel_BCFunction
{
public:

    //! @name Type definitions and Enumerators
    //@{

    typedef boost::function<Real ( const Real&, const Real&  )> Function_Type;

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    OneDimensionalModel_BCFunction();

    //! Constructor by function
    /*!
      @param function the user defined function
    */
    OneDimensionalModel_BCFunction( const Function_Type& function );

    //! Copy constructor
    OneDimensionalModel_BCFunction( const OneDimensionalModel_BCFunction& BCFunction );

    //! Destructor
    virtual ~OneDimensionalModel_BCFunction() {}

    //@}


    //! @name Operators
    //@{

    OneDimensionalModel_BCFunction& operator= ( const OneDimensionalModel_BCFunction& BCFunction );

    Real operator() ( const Real& time, const Real& timeStep = 0. ) const;

    //@}


    //! @name Set Methods
    //@{

    //! Set the function
    /*!
      @param function the user defined function
    */
    void setFunction( const Function_Type& function );

    //@}


    //! @name Get Methods
    //@{

    //! Get the function
    /*!
      @return the user defined function
    */
    const Function_Type& Function() const;

    //@}

private:

    Function_Type M_function;
};

/*
//! Factory create function
inline OneDimensionalModel_BCFunction*
Create_OneDimensionalModel_BCFunction( const OneDimensionalModel_BCFunction* BCFunction )
{
    return new OneDimensionalModel_BCFunction( (const OneDimensionalModel_BCFunction&)* BCFunction );
}

namespace
{
    static bool registerOneD_BCFunction = FactoryClone_OneDimensionalModel_BCFunction::instance().registerProduct( typeid(OneDimensionalModel_BCFunction), &Create_OneDimensionalModel_BCFunction );
}
*/

}

#endif // ONEDIMENSIONALMODEL_BCFUNCTION_H
