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
 *  @brief File containing the BCInterface3DFunctionFileSolver class
 *
 *  @date 26-08-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterface3DFunctionFileSolver_H
#define BCInterface3DFunctionFileSolver_H 1

#include <lifemc/lifesolver/BCInterface3DFunctionFile.hpp>
#include <lifemc/lifesolver/BCInterface3DFunctionSolver.hpp>

namespace LifeV
{

//! BCInterface3DFunctionFileSolver - LifeV bcFunction wrapper for BCInterface (with Operators)
/*!
 *  @author Cristiano Malossi
 *
 *  This class is an interface between BCInterface, SpiritParser and and a general
 *  LifeV operator (such as Oseen or FSIOperator). It allows to construct LifeV
 *  functions type for boundary conditions, using a GetPot file containing a function string and a
 *  table of discrete data (for example a discrete Flux or Pressure depending on time).
 *  The function string can contain Operator parameters.
 */
template< class PhysicalSolverType >
class BCInterface3DFunctionFileSolver: public virtual BCInterface3DFunctionFile< PhysicalSolverType > ,
                                       public virtual BCInterface3DFunctionSolver< PhysicalSolverType >
{
public:

    //! @name Type definitions
    //@{

    typedef PhysicalSolverType                                                    physicalSolver_Type;
    typedef BCInterface3DData                                                     data_Type;
    typedef BCInterface3DFunction< physicalSolver_Type >                          function_Type;
    typedef BCInterface3DFunctionFile< physicalSolver_Type >                      functionFile_Type;
    typedef BCInterface3DFunctionSolver< physicalSolver_Type >                    functionSolver_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterface3DFunctionFileSolver();

    //! Constructor
    /*!
     * @param data BC data loaded from GetPot file
     */
    explicit BCInterface3DFunctionFileSolver( const data_Type& data );

    //! Destructor
    virtual ~BCInterface3DFunctionFileSolver() {}

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

    BCInterface3DFunctionFileSolver( const BCInterface3DFunctionFileSolver& function );

    BCInterface3DFunctionFileSolver& operator=( const BCInterface3DFunctionFileSolver& function );

    //@}

};

// ===================================================
// Factory
// ===================================================
//! Factory create function
template< typename PhysicalSolverType >
inline BCInterface3DFunction< PhysicalSolverType >* createBCInterface3DFunctionFileSolver()
{
    return new BCInterface3DFunctionFileSolver< PhysicalSolverType > ();
}

// ===================================================
// Constructors
// ===================================================
template< class PhysicalSolverType >
BCInterface3DFunctionFileSolver< PhysicalSolverType >::BCInterface3DFunctionFileSolver() :
        function_Type          (),
        functionFile_Type      (),
        functionSolver_Type    ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5024 ) << "BCInterface3DFunctionFileSolver::BCInterface3DFunctionFileSolver()" << "\n";
#endif

}

template< class PhysicalSolverType >
BCInterface3DFunctionFileSolver< PhysicalSolverType >::BCInterface3DFunctionFileSolver( const data_Type& data ) :
        function_Type          (),
        functionFile_Type      (),
        functionSolver_Type    ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5024 ) << "BCInterface3DFunctionFileSolver::BCInterface3DFunctionFileSolver( data )" << "\n";
#endif

    this->setData( data );
}

// ===================================================
// Set Methods
// ===================================================
template< class PhysicalSolverType >
void
BCInterface3DFunctionFileSolver< PhysicalSolverType >::setData( const data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5024 ) << "BCInterface3DFunctionFileSolver::setData" << "\n";
#endif
    functionFile_Type::setData( data );

    //functionSolver_Type::setData( data ); Cannot call directly, because it call again BCInterface3DFunction::setup( data )
    functionSolver_Type::M_flag = data.flag();

    functionSolver_Type::createAccessList( data );
}

} // Namespace LifeV

#endif /* BCInterface3DFunctionFileSolver_H */
