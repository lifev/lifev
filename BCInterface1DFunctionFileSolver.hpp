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
 *  @brief File containing the BCInterface1D_OperatorFunctionFile class
 *
 *  @date 10-05-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterface1D_OperatorFunctionFile_H
#define BCInterface1D_OperatorFunctionFile_H 1

#include <lifemc/lifesolver/BCInterface1DFunctionFile.hpp>
#include <lifemc/lifesolver/BCInterface1DFunctionSolver.hpp>

namespace LifeV
{

//! BCInterface1D_OperatorFunctionFile - LifeV bcFunction wrapper for BCInterface1D (with Operators)
/*!
 *  @author Cristiano Malossi
 *
 *  This class is an interface between BCInterface1D, SpiritParser and and a general
 *  LifeV operator (such as Oseen or FSIOperator). It allows to construct LifeV
 *  functions type for boundary conditions, using a GetPot file containing a function string and a
 *  table of discrete data (for example a discrete Flux or Pressure depending on time).
 *  The function string can contain Operator parameters.
 */
template< class PhysicalSolverType >
class BCInterface1D_OperatorFunctionFile: public virtual BCInterface1D_FunctionFile< PhysicalSolverType > ,
        public virtual BCInterface1D_OperatorFunction< PhysicalSolverType >
{
public:

    //! @name Type definitions
    //@{

    typedef PhysicalSolverType                                                    physicalSolver_Type;
    typedef BCInterface1D_Data                                                    data_Type;
    typedef BCInterface1D_Function< physicalSolver_Type >                         function_Type;
    typedef BCInterface1D_FunctionFile< physicalSolver_Type >                     functionFile_Type;
    typedef BCInterface1D_OperatorFunction< physicalSolver_Type >                 functionSolver_Type;
    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterface1D_OperatorFunctionFile();

    //! Constructor
    /*!
     * @param data BC data loaded from GetPot file
     */
    explicit BCInterface1D_OperatorFunctionFile( const data_Type& data );

    //! Destructor
    virtual ~BCInterface1D_OperatorFunctionFile() {}

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

    BCInterface1D_OperatorFunctionFile( const BCInterface1D_OperatorFunctionFile& function );

    BCInterface1D_OperatorFunctionFile& operator=( const BCInterface1D_OperatorFunctionFile& function );

    //@}

};

// ===================================================
// Factory
// ===================================================
//! Factory create function
template< typename PhysicalSolverType >
inline BCInterface1D_Function< PhysicalSolverType >* createBCInterface1D_OperatorFunctionFile()
{
    return new BCInterface1D_OperatorFunctionFile< PhysicalSolverType > ();
}

// ===================================================
// Constructors
// ===================================================
template< class PhysicalSolverType >
BCInterface1D_OperatorFunctionFile< PhysicalSolverType >::BCInterface1D_OperatorFunctionFile() :
        function_Type          (),
        functionFile_Type      (),
        functionSolver_Type    ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5024 ) << "BCInterface1D_OperatorFunctionFile::BCInterface1D_OperatorFunctionFile()" << "\n";
#endif

}

template< class PhysicalSolverType >
BCInterface1D_OperatorFunctionFile< PhysicalSolverType >::BCInterface1D_OperatorFunctionFile( const data_Type& data ) :
        function_Type          (),
        functionFile_Type      (),
        functionSolver_Type    ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5024 ) << "BCInterface1D_OperatorFunctionFile::BCInterface1D_OperatorFunctionFile( data )" << "\n";
#endif

    this->setData( data );
}

// ===================================================
// Set Methods
// ===================================================
template< class PhysicalSolverType >
void
BCInterface1D_OperatorFunctionFile< PhysicalSolverType >::setData( const data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5024 ) << "BCInterface1D_OperatorFunctionFile::setData" << "\n";
#endif
    functionFile_Type::setData( data );

    //functionSolver_Type::setData( data ); Cannot call directly, because it call again BCInterface1D1D_Function::setup( data )
    functionSolver_Type::M_side = data.side();

    functionSolver_Type::createAccessList( data );
}

} // Namespace LifeV

#endif /* BCInterface1D_OperatorFunctionFile_H */
