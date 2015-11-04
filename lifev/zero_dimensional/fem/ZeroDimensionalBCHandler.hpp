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

#ifndef ZeroDimensionalBCHandler_H
#define ZeroDimensionalBCHandler_H 1

#include <lifev/zero_dimensional/fem/ZeroDimensionalBC.hpp>

namespace LifeV
{

//! ZeroDimensionalBCHandler - A boundary conditions handler for zero-dimensional models.
/*!
 *  @author Cristiano Malossi
 *
 *  This simple class handles the boundary conditions for zero-dimensional models.
 */
class ZeroDimensionalBCHandler
{
public:

    //! @name Type definitions
    //@{

    typedef ZeroDimensionalBC                                             bc_Type;
    typedef bc_Type::bcType_Type                                          bcType_Type;
    typedef bc_Type::bcFunction_Type                                      bcFunction_Type;

    typedef markerID_Type                                                 bcFlag_Type;
    typedef std::map< bcFlag_Type, bc_Type >                              bcContainer_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit ZeroDimensionalBCHandler() : M_bc() {}

    explicit ZeroDimensionalBCHandler ( const ZeroDimensionalBCHandler& handler ) : M_bc ( handler.M_bc ) {}

    //! Destructor
    virtual ~ZeroDimensionalBCHandler() {}

    //@}


    //! @name Set Methods
    //@{

    //! Set the type
    /*!
      @param flag the bc flag
      @param bcType the bc type
      @param bcFunction the bc function
    */
    void setBC ( const bcFlag_Type& flag, const bcType_Type& bcType, const bcFunction_Type& bcFunction )
    {
        M_bc[flag].setBC ( bcType, bcFunction );
    }

    //@}


    //! @name Get Methods
    //@{

    //! Get the boundary condition
    /*!
     *  @param side the boundary condition side
     *  @return boundary condition
     */
    const bc_Type& bc ( const bcFlag_Type& flag ) const
    {
        return M_bc.find ( flag )->second;
    }

    //@}

private:

    bcContainer_Type            M_bc;
};

} // Namespace LifeV

#endif /* ZeroDimensionalBCHandler_H */
