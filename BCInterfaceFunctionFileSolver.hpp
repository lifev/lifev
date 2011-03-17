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
 *  @brief File containing the BCInterfaceFunctionFileSolver class
 *
 *  @date 26-08-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterfaceFunctionFileSolver_H
#define BCInterfaceFunctionFileSolver_H 1

#include <lifemc/lifesolver/BCInterfaceFunctionFile.hpp>
#include <lifemc/lifesolver/BCInterfaceFunctionSolver.hpp>

namespace LifeV
{

//! BCInterfaceFunctionFileSolver - LifeV bcFunction wrapper for BCInterface (with Operators)
/*!
 *  @author Cristiano Malossi
 *
 *  This class is an interface between BCInterface, SpiritParser and and a general
 *  LifeV operator (such as Oseen or FSI). It allows to construct LifeV
 *  functions type for boundary conditions, using a GetPot file containing a function string and a
 *  table of discrete data (for example a discrete Flux or Pressure depending on time).
 *  The function string can contain Operator parameters.
 */
template< class PhysicalSolverType >
class BCInterfaceFunctionFileSolver: public virtual BCInterfaceFunctionFile< PhysicalSolverType > ,
                                       public virtual BCInterfaceFunctionSolver< PhysicalSolverType >
{
public:

    //! @name Type definitions
    //@{

    typedef PhysicalSolverType                                                    physicalSolver_Type;
    typedef BCInterfaceData                                                       data_Type;
    typedef BCInterfaceFunction< physicalSolver_Type >                            function_Type;
    typedef BCInterfaceFunctionFile< physicalSolver_Type >                        functionFile_Type;
    typedef BCInterfaceFunctionSolver< physicalSolver_Type >                      functionSolver_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterfaceFunctionFileSolver();

    //! Constructor
    /*!
     * @param data BC data loaded from GetPot file
     */
    explicit BCInterfaceFunctionFileSolver( const data_Type& data );

    //! Destructor
    virtual ~BCInterfaceFunctionFileSolver() {}

    //@}


    //! @name Set Methods
    //@{

    //! Set data
    /*!
     * @param data BC data loaded from GetPot file
     */
    virtual void setData( const data_Type& data );

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    BCInterfaceFunctionFileSolver( const BCInterfaceFunctionFileSolver& function );

    BCInterfaceFunctionFileSolver& operator=( const BCInterfaceFunctionFileSolver& function );

    //@}

};

// ===================================================
// Factory
// ===================================================
//! Factory create function
template< typename PhysicalSolverType >
inline BCInterfaceFunction< PhysicalSolverType >* createBCInterfaceFunctionFileSolver()
{
    return new BCInterfaceFunctionFileSolver< PhysicalSolverType > ();
}

// ===================================================
// Constructors
// ===================================================
template< class PhysicalSolverType >
BCInterfaceFunctionFileSolver< PhysicalSolverType >::BCInterfaceFunctionFileSolver() :
        function_Type          (),
        functionFile_Type      (),
        functionSolver_Type    ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5024 ) << "BCInterfaceFunctionFileSolver::BCInterfaceFunctionFileSolver()" << "\n";
#endif

}

template< class PhysicalSolverType >
BCInterfaceFunctionFileSolver< PhysicalSolverType >::BCInterfaceFunctionFileSolver( const data_Type& data ) :
        function_Type          (),
        functionFile_Type      (),
        functionSolver_Type    ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5024 ) << "BCInterfaceFunctionFileSolver::BCInterfaceFunctionFileSolver( data )" << "\n";
#endif

    this->setData( data );
}

// ===================================================
// Set Methods
// ===================================================
template< class PhysicalSolverType >
void
BCInterfaceFunctionFileSolver< PhysicalSolverType >::setData( const data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5024 ) << "BCInterfaceFunctionFileSolver::setData" << "\n";
#endif
    functionFile_Type::setData( data );

    //functionSolver_Type::setData( data ); Cannot call directly, because it call again BCInterfaceFunction::setup( data )
    functionSolver_Type::M_side = data.side();
    functionSolver_Type::M_flag = data.flag();

    functionSolver_Type::createAccessList( data );
}

} // Namespace LifeV

#endif /* BCInterfaceFunctionFileSolver_H */
