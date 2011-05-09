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
 *  @brief File containing the zero dimensional BCInterface
 *
 *  @date 30-03-2011
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterface0D_H
#define BCInterface0D_H 1

// Mathcard includes
#include <lifemc/lifesolver/BCInterface.hpp>

namespace LifeV
{

//! BCInterface0D - A very simple BCInterface for zero-dimensional models
/*!
 *  @author Cristiano Malossi
 *
 *  This simple class provide the BCInterface for zero-dimensional model BCHandler.
 */
template< class BcHandler, class PhysicalSolverType >
class BCInterface0D : public virtual BCInterface< BcHandler, PhysicalSolverType >
{
public:

    //! @name Type definitions
    //@{

    typedef BCInterface< BcHandler, PhysicalSolverType >          bcInterface_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterface0D() : bcInterface_Type() {}

    //! Destructor
    virtual ~BCInterface0D() {}

    //@}


    //! @name Methods
    //@{

    //! Insert the current boundary condition in the BChandler
    void insertBC()
    {
        bcInterface_Type::insertBC();

        addBcToHandler();
    }

    //@}

private:

    void addBcToHandler()
    {
        if ( !this->M_handler.get() )
            this->createHandler();

        this->M_handler->setBC( this->M_data.side(), this->M_data.quantity(), boost::bind( &BCInterfaceFunction<PhysicalSolverType>::functionTime, this->M_vectorFunction.back(), _1 ) );
    }
};

} // Namespace LifeV

#endif /* BCInterface0D_H */
