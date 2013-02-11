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
 *  @brief File containing the zero dimensional BCHandler
 *
 *  @date 30-03-2011
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef ZeroDimensionalBC_H
#define ZeroDimensionalBC_H 1

#include <lifev/zero_dimensional/solver/ZeroDimensionalDefinitions.hpp>

namespace LifeV
{

//! ZeroDimensionalBC - A boundary condition for zero-dimensional models
/*!
 *  @author Cristiano Malossi
 *
 *  This simple class is a boundary condition for simple zero-dimensional models.
 */
class ZeroDimensionalBC
{
public:

    //! @name Type definitions
    //@{

    typedef boost::function<Real ( const Real&  ) >                        function_Type;
    typedef ZeroDimensionalBCType                                         bcType_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit ZeroDimensionalBC() : M_function(), M_bcType() {}

    ZeroDimensionalBC ( const ZeroDimensionalBC& bc ) :
        M_function ( bc.M_function ),
        M_bcType  ( bc.M_bcType ) {}

    //! Destructor
    virtual ~ZeroDimensionalBC() {}

    //@}


    //! @name Methods
    //@{

    //! Evaluate the bc
    /*!
     * @param time the current time of the simulation
     * @return the bc value
    */
    Real evaluate ( const Real& time ) const
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
    void setBC ( const bcType_Type& bcType, const function_Type& function )
    {
        M_bcType = bcType;
        M_function = function;
    }

    //! Set the function
    /*!
      @param function the user defined function
    */
    void setFunction ( const function_Type& function )
    {
        M_function = function;
    }

    //! Set the type
    /*!
      @param bcType the bc type
    */
    void setBcType ( const bcType_Type& bcType )
    {
        M_bcType = bcType;
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

    //! Get the type
    /*!
      @return the bc type
    */
    const bcType_Type& bcType() const
    {
        return M_bcType;
    }

    //@}

private:

    function_Type               M_function;
    bcType_Type                 M_bcType;
};

} // Namespace LifeV

#endif /* ZeroDimensionalBCHandler_H */
