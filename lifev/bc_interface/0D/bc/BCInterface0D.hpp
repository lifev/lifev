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
 *  @brief File containing the BCInterface0D class
 *
 *  @date 30-03-2011
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterface0D_H
#define BCInterface0D_H 1

// BCInterface includes
#include <lifev/bc_interface/core/bc/BCInterface.hpp>
#include <lifev/bc_interface/0D/bc/BCInterfaceData0D.hpp>

// Template specializations
#include <lifev/bc_interface/0D/function/BCInterfaceFunctionParser0D.hpp>
#include <lifev/bc_interface/0D/function/BCInterfaceFunctionParserSolver0D.hpp>
#include <lifev/bc_interface/0D/function/BCInterfaceFunctionUserDefined0D.hpp>

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
 *  type0D              = Current        <BR>
 *  function            = 'sin(2*pi*t)'  <BR>
 *
 *  [OutFlow]                            <BR>
 *  flag                = 1              <BR>
 *  type0D              = Voltage        <BR>
 *  function            = 0              <BR>
 *  </CODE>
 *
 *  where \c flag, and \c type0D are the classical parameters for a 0D boundary condition.
 *  The string \c function represents the base module and can be replaced by other derived/alternative modules.
 *  The following functions are available (see the related classes for more information):
 *
 *  <ol>
 *      <li> \c function, which is implemented in \c BCInterfaceFunctionParser;
 *      <li> \c functionFile, which is implemented in \c BCInterfaceFunctionParserFile;
 *      <li> \c functionSolver, which is implemented in \c BCInterfaceFunctionParserSolver;
 *      <li> \c functionFileSolver, which is implemented in \c BCInterfaceFunctionParserFileSolver;
 *      <li> \c functionUD, which is implemented in \c BCInterfaceFunctionUserDefined.
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

    typedef BCInterface< BcHandler, PhysicalSolverType >                bcInterface_Type;

    typedef typename bcInterface_Type::bcHandler_Type                   bcHandler_Type;
    typedef typename bcInterface_Type::bcHandlerPtr_Type                bcHandlerPtr_Type;

    typedef typename bcInterface_Type::physicalSolver_Type              physicalSolver_Type;
    typedef typename bcInterface_Type::physicalSolverPtr_Type           physicalSolverPtr_Type;

    typedef typename bcInterface_Type::factory_Type                     factory_Type;

    typedef typename bcInterface_Type::bcFunctionPtr_Type               bcFunctionPtr_Type;
    typedef typename bcInterface_Type::vectorFunction_Type              vectorFunction_Type;

    typedef typename bcInterface_Type::bcFunctionSolverDefinedPtr_Type  bcFunctionSolverDefinedPtr_Type;
    typedef typename bcInterface_Type::vectorFunctionSolverDefined_Type vectorFunctionSolverDefined_Type;

    typedef BCInterfaceData0D                                           data_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterface0D() : bcInterface_Type(), M_data() {}

    //! Destructor
    virtual ~BCInterface0D() {}

    //@}


    //! @name Methods
    //@{

    //! Read a specific boundary condition from a file and add it to the data container
    /*!
     * @param fileName Name of the data file
     * @param dataSection section in the data file
     * @param name name of the boundary condition
     */
    void readBC ( const std::string& fileName, const std::string& dataSection, const std::string& name )
    {
        M_data.readBC ( fileName, dataSection, name );
    }

    //! Insert the current boundary condition in the BChandler
    void insertBC();

    //@}


    //! @name Get Methods
    //@{

    //! Get the data container
    /*!
     * @return the data container
     */
    data_Type& dataContainer()
    {
        return M_data;
    }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    BCInterface0D ( const BCInterface0D& interface0D );

    BCInterface0D& operator= ( const BCInterface0D& interface0D );

    //@}


    //! @name Private Methods
    //@{

    template< class BCBaseType >
    void addBcToHandler( BCBaseType& base );

    //@}

    // Data
    data_Type                       M_data;
};

// ===================================================
// Methods
// ===================================================
template< class BcHandler, class PhysicalSolverType >
inline void
BCInterface0D< BcHandler, PhysicalSolverType >::insertBC()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5020 ) << "BCInterface0D::insertBC\n";
#endif

    ZeroDimensionalFunction base;

    switch ( M_data.base().second )
    {
        case BCIFunctionParser:
        case BCIFunctionParserFile:
        case BCIFunctionParserSolver:
        case BCIFunctionParserFileSolver:
        case BCIFunctionUserDefined:
        {
            factory_Type factory;
            this->M_vectorFunction.push_back ( factory.createFunctionParser ( M_data ) );
            this->M_vectorFunction.back()->assignFunction ( base );

            addBcToHandler( base );

            return;
        }

        default:

            std::cout << " !!! Error: " << M_data.base().first << " is not valid in BCInterface0D !!!" << std::endl;
            break;
    }
}



// ===================================================
// Private Methods
// ===================================================
template< class BcHandler, class PhysicalSolverType > template< class BCBaseType >
inline void
BCInterface0D< BcHandler, PhysicalSolverType >::addBcToHandler( BCBaseType& base )
{
    if ( !this->M_handler.get() ) // If BCHandler has not been created yet, we do it now
    {
        this->createHandler();
    }
    //boost::bind ( &BCInterfaceFunction<BcHandler, PhysicalSolverType>::functionTime, this->M_vectorFunction.back(), _1 )
    this->M_handler->setBC ( M_data.flag(), M_data.type(), base );
}

} // Namespace LifeV

#endif /* BCInterface0D_H */
