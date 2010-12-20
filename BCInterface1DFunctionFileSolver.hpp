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
 *  @brief File containing the BCInterface1DFunctionFileSolver class
 *
 *  @date 10-05-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterface1DFunctionFileSolver_H
#define BCInterface1DFunctionFileSolver_H 1

#include <lifemc/lifesolver/BCInterface1DFunctionFile.hpp>
#include <lifemc/lifesolver/BCInterface1DFunctionSolver.hpp>

namespace LifeV
{

//! BCInterface1DFunctionFileSolver - LifeV bcFunction wrapper for BCInterface1D (with Operators)
/*!
 *  @author Cristiano Malossi
 *
 *  This class is an interface between BCInterface1D, SpiritParser and and a general
 *  LifeV operator (such as Oseen or FSI). It allows to construct LifeV
 *  functions type for boundary conditions, using a GetPot file containing a function string and a
 *  table of discrete data (for example a discrete Flux or Pressure depending on time).
 *  The function string can contain Operator parameters.
 */
template< class PhysicalSolverType >
class BCInterface1DFunctionFileSolver: public virtual BCInterface1DFunctionFile< PhysicalSolverType > ,
        public virtual BCInterface1DFunctionSolver< PhysicalSolverType >
{
public:

    //! @name Type definitions
    //@{

    typedef PhysicalSolverType                                                    physicalSolver_Type;
    typedef BCInterface1DData                                                     data_Type;
    typedef BCInterface1DFunction< physicalSolver_Type >                          function_Type;
    typedef BCInterface1DFunctionFile< physicalSolver_Type >                      functionFile_Type;
    typedef BCInterface1DFunctionSolver< physicalSolver_Type >                    functionSolver_Type;
    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterface1DFunctionFileSolver();

    //! Constructor
    /*!
     * @param data BC data loaded from GetPot file
     */
    explicit BCInterface1DFunctionFileSolver( const data_Type& data );

    //! Destructor
    virtual ~BCInterface1DFunctionFileSolver() {}

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

    BCInterface1DFunctionFileSolver( const BCInterface1DFunctionFileSolver& function );

    BCInterface1DFunctionFileSolver& operator=( const BCInterface1DFunctionFileSolver& function );

    //@}

};

// ===================================================
// Factory
// ===================================================
//! Factory create function
template< typename PhysicalSolverType >
inline BCInterface1DFunction< PhysicalSolverType >* createBCInterface1DFunctionFileSolver()
{
    return new BCInterface1DFunctionFileSolver< PhysicalSolverType > ();
}

// ===================================================
// Constructors
// ===================================================
template< class PhysicalSolverType >
BCInterface1DFunctionFileSolver< PhysicalSolverType >::BCInterface1DFunctionFileSolver() :
        function_Type          (),
        functionFile_Type      (),
        functionSolver_Type    ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5024 ) << "BCInterface1DOperatorFunctionFile::BCInterface1DOperatorFunctionFile()" << "\n";
#endif

}

template< class PhysicalSolverType >
BCInterface1DFunctionFileSolver< PhysicalSolverType >::BCInterface1DFunctionFileSolver( const data_Type& data ) :
        function_Type          (),
        functionFile_Type      (),
        functionSolver_Type    ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5024 ) << "BCInterface1DOperatorFunctionFile::BCInterface1DOperatorFunctionFile( data )" << "\n";
#endif

    this->setData( data );
}

// ===================================================
// Set Methods
// ===================================================
template< class PhysicalSolverType >
void
BCInterface1DFunctionFileSolver< PhysicalSolverType >::setData( const data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5024 ) << "BCInterface1DOperatorFunctionFile::setData" << "\n";
#endif
    functionFile_Type::setData( data );

    //functionSolver_Type::setData( data ); Cannot call directly, because it call again BCInterface1D1D_Function::setup( data )
    functionSolver_Type::M_side = data.side();

    functionSolver_Type::createAccessList( data );
}

} // Namespace LifeV

#endif /* BCInterface1DFunctionFileSolver_H */
