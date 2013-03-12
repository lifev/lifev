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
 *  @brief File containing the BCInterfaceFunctionParserFileSolver class
 *
 *  @date 26-08-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterfaceFunctionParserFileSolver_H
#define BCInterfaceFunctionParserFileSolver_H 1

#include <lifev/bc_interface/core/function/BCInterfaceFunctionParserFile.hpp>
#include <lifev/bc_interface/core/function/BCInterfaceFunctionParserSolver.hpp>

namespace LifeV
{

//! BCInterfaceFunctionParserFileSolver - LifeV boundary condition function file wrapper for \c BCInterface
/*!
 *  @author Cristiano Malossi
 *
 *  This class is an interface between the \c BCInterface, the \c Parser, and and a general
 *  LifeV physical solver (such as \c OseenSolver or \c FSISolver). It allows to construct LifeV
 *  functions type for boundary conditions, using a \c GetPot file containing a function string and a
 *  table of discrete data (for example a discrete Flux or Pressure depending on time).
 *  The function string can contain physical solver parameters.
 *
 *  See \c BCInterfaceFunctionParser, \c BCInterfaceFunctionParserFile, and \c BCInterfaceFunctionParserSolver classes for more details.
 */
template< class PhysicalSolverType >
class BCInterfaceFunctionParserFileSolver: public virtual BCInterfaceFunctionParserFile< PhysicalSolverType > ,
    public virtual BCInterfaceFunctionParserSolver< PhysicalSolverType >
{
public:

    //! @name Type definitions
    //@{

    typedef PhysicalSolverType                                                    physicalSolver_Type;
    typedef BCInterfaceFunction< physicalSolver_Type >                            function_Type;
    typedef BCInterfaceFunctionParser< physicalSolver_Type >                      functionParser_Type;
    typedef BCInterfaceFunctionParserFile< physicalSolver_Type >                  functionParserFile_Type;
    typedef BCInterfaceFunctionParserSolver< physicalSolver_Type >                functionParserSolver_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterfaceFunctionParserFileSolver();

    //! Destructor
    virtual ~BCInterfaceFunctionParserFileSolver() {}

    //@}


    //! @name Set Methods
    //@{

    //! Set data for 0D boundary conditions
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    virtual void setData ( const BCInterfaceData0D& data );

    //! Set data for 1D boundary conditions
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    virtual void setData ( const BCInterfaceData1D& data );

    //! Set data for 3D boundary conditions
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    virtual void setData ( const BCInterfaceData3D& data );

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    BCInterfaceFunctionParserFileSolver ( const BCInterfaceFunctionParserFileSolver& function );

    BCInterfaceFunctionParserFileSolver& operator= ( const BCInterfaceFunctionParserFileSolver& function );

    //@}

};

// ===================================================
// Factory
// ===================================================
//! Factory create function
template< typename PhysicalSolverType >
inline BCInterfaceFunctionParser< PhysicalSolverType >* createBCInterfaceFunctionParserFileSolver()
{
    return new BCInterfaceFunctionParserFileSolver< PhysicalSolverType > ();
}

// ===================================================
// Constructors
// ===================================================
template< class PhysicalSolverType >
BCInterfaceFunctionParserFileSolver< PhysicalSolverType >::BCInterfaceFunctionParserFileSolver() :
    function_Type                (),
    functionParser_Type          (),
    functionParserFile_Type      (),
    functionParserSolver_Type    ()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5024 ) << "BCInterfaceFunctionFileSolver::BCInterfaceFunctionFileSolver()" << "\n";
#endif

}

// ===================================================
// Set Methods
// ===================================================
template< class PhysicalSolverType >
void
BCInterfaceFunctionParserFileSolver< PhysicalSolverType >::setData ( const BCInterfaceData0D& data )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5024 ) << "BCInterfaceFunctionFileSolver::setData" << "\n";
#endif
    functionParserFile_Type::setData ( data );

    //functionParserSolver_Type::setData( data ); Cannot call directly, because it call again BCInterfaceFunctionParser::setup( data )
    functionParserSolver_Type::M_flag = data.flag();

    functionParserSolver_Type::createAccessList ( data );
}

template< class PhysicalSolverType >
void
BCInterfaceFunctionParserFileSolver< PhysicalSolverType >::setData ( const BCInterfaceData1D& data )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5024 ) << "BCInterfaceFunctionFileSolver::setData" << "\n";
#endif
    functionParserFile_Type::setData ( data );

    //functionParserSolver_Type::setData( data ); Cannot call directly, because it call again BCInterfaceFunctionParser::setup( data )
    functionParserSolver_Type::M_side = data.side();

    functionParserSolver_Type::createAccessList ( data );
}

template< class PhysicalSolverType >
void
BCInterfaceFunctionParserFileSolver< PhysicalSolverType >::setData ( const BCInterfaceData3D& data )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5024 ) << "BCInterfaceFunctionFileSolver::setData" << "\n";
#endif
    functionParserFile_Type::setData ( data );

    //functionParserSolver_Type::setData( data ); Cannot call directly, because it call again BCInterfaceFunctionParser::setup( data )
    functionParserSolver_Type::M_flag = data.flag();

    functionParserSolver_Type::createAccessList ( data );
}

} // Namespace LifeV

#endif /* BCInterfaceFunctionParserFileSolver_H */
