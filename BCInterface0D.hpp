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

#include <life/lifesolver/BCInterface.hpp>

// Move this to BCInterfaceDefinitions.hpp
#include <lifemc/lifefem/ZeroDimensionalBCHandler.hpp>

namespace LifeV
{

//! BCInterface0D - LifeV interface to load boundary conditions for 0D problems completely from a \c GetPot file
/*!
 *  @author Cristiano Malossi
 *
 *  This class allows to impose boundary conditions for a 0D problem completely from a file.
 *
 *  <b>EXAMPLE - DATA FILE</b> <BR>
 *  In the GetPot data file, \c BCInterface reads a new section: <CODE> [boundary_conditions] </CODE>.
 *
 *  Inside the new section there is a list of boundary conditions which correspond to other sub-section
 *  with the same name, for example: <CODE> list = 'InFlow OutFlow' </CODE>
 *
 *  Each boundary condition has a similar structure. The list of properties depends from the type of the
 *  boundary condition. For example:
 *
 *  <CODE>
 *  [InFlow]                             <BR>
 *  flag                = 0              <BR>
 *  quantity            = Q              <BR>
 *  function            = 'sin(2*pi*t)'  <BR>
 *
 *  [OutFlow]                            <BR>
 *  flag                = 1              <BR>
 *  quantity            = S              <BR>
 *  function            = 0              <BR>
 *  </CODE>
 *
 *  where \c flag, and \c quantity are the classical parameters for a 0D boundary condition.
 *  The string \c function represents the base module and can be replaced by other derived/alternative modules.
 *  The following functions are available (see the related classes for more information):
 *
 *  <ol>
 *      <li> \c function, which is implemented in \c BCInterfaceFunction;
 *      <li> \c functionFile, which is implemented in \c BCInterfaceFunctionFile;
 *      <li> \c functionSolver, which is implemented in \c BCInterfaceFunctionSolver;
 *      <li> \c functionFileSolver, which is implemented in \c BCInterfaceFunctionFileSolver;
 *  </ol>
 *
 *  All the parameters are case sensitive.
 *
 *  See \c BCInterface base class for more details.
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

        this->M_handler->setBC( this->M_data.flag(), this->M_data.quantity(), boost::bind( &BCInterfaceFunction<PhysicalSolverType>::functionTime, this->M_vectorFunction.back(), _1 ) );
    }
};

} // Namespace LifeV

#endif /* BCInterface0D_H */
