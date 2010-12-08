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
 *  @brief File containing the BCInterface1D_Function class
 *
 *  @date 10-05-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterface1D_Function_H
#define BCInterface1D_Function_H 1

#include <lifemc/lifesolver/BCInterface1D_Definitions.hpp>
#include <lifemc/lifesolver/BCInterface1D_Data.hpp>

#include <lifemc/lifecore/Parser.hpp>

namespace LifeV
{

//! BCInterface1D_Function - LifeV bcFunction wrapper for BCInterface
/*!
 *  @author Cristiano Malossi
 *
 *  This class is an interface between BCInterface1D and the grammar parser. It allows to construct LifeV
 *  functions type for boundary conditions, using a functions string loaded from a GetPot file.
 *
 *  <b>DETAILS:</b>
 *
 *  By default the stringSeparator is set to semicolon ";".
 *
 *  The function string has to be in this form:
 *
 *  function = q
 *
 *  where q(t).
 *  Here there is an example:
 *
 *  function = '2*sin(2*pi*t)'
 *
 *  To set a constant for complicate expression it is possible to add them before the expression
 *  using a semicolon ";":
 *
 *  function = 'a=5.67436; (a*sin(2*pi*t))'
 */
template< typename PhysicalSolverType >
class BCInterface1D_Function
{
public:

    //! @name Type definitions
    //@{

    typedef PhysicalSolverType                                                    physicalSolver_Type;
    typedef BCInterface1D_Data                                                    data_Type;
    typedef OneDimensionalModel_BCFunction                                        bcFunction_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    BCInterface1D_Function();

    //! Constructor
    /*!
     * @param data BC data loaded from GetPot file
     */
    BCInterface1D_Function( const data_Type& data );

    //! Destructor
    virtual ~BCInterface1D_Function() {}

    //@}


    //! @name Set Methods
    //@{

    //! Set data
    /*!
     * @param data BC data loaded from GetPot file
     */
    virtual void setData( const data_Type& data );

    //@}


    //! @name Get Methods
    //@{

    //! Get the base of the boundary condition
    /*!
     * @return boundary condition base
     */
    bcFunction_Type& base() { return M_base; }

    //@}

protected:

    //! @name Protected Methods
    //@{

    //! dataInterpolation
    virtual void dataInterpolation() {}

    //@}

    boost::shared_ptr< Parser >       M_parser;

private:

    //! @name Unimplemented Methods
    //@{

    BCInterface1D_Function( const BCInterface1D_Function& function );

    BCInterface1D_Function& operator=( const BCInterface1D_Function& function );

    //@}


    //! @name Private Methods
    //@{

    //! setFunction
    void setFunction() { M_base.setFunction( boost::bind( &BCInterface1D_Function::function, this, _1 ) ); }

    //! function
    Real function( const Real& t );

    //@}

    bcFunction_Type                  M_base;

};

// ===================================================
// Factory
// ===================================================
//! Factory create function
template< typename PhysicalSolverType >
inline BCInterface1D_Function< PhysicalSolverType >* createBCInterface1D_Function()
{
    return new BCInterface1D_Function< PhysicalSolverType > ();
}

// ===================================================
// Constructor
// ===================================================
template< typename PhysicalSolverType >
BCInterface1D_Function< PhysicalSolverType >::BCInterface1D_Function() :
        M_parser    (),
        M_base      ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "BCInterface1D_Function::BCInterface1D_Function()" << "\n";
#endif

}

template< typename PhysicalSolverType >
BCInterface1D_Function< PhysicalSolverType >::BCInterface1D_Function( const data_Type& data ) :
        M_parser    (),
        M_base      ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "BCInterface1D_Function::BCInterface1D_Function( data )" << "\n";
#endif

    this->setData( data );
}

// ===================================================
// Set Methods
// ===================================================
template< typename PhysicalSolverType >
void
BCInterface1D_Function< PhysicalSolverType >::setData( const data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5022 ) << "BCInterface1D_Function::setData" << "\n";
#endif

    if ( M_parser )
        M_parser->SetString( data.baseString() );
    else
        M_parser.reset( new Parser( data.baseString() ) );

    setFunction();
}

// ===================================================
// Private Methods
// ===================================================
template< typename PhysicalSolverType >
Real
BCInterface1D_Function< PhysicalSolverType >::function( const Real& t )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "BCInterface1D_Function::Function: " << "\n";
    Debug( 5021 ) << "                                                           t: " << t << "\n";
#endif

    M_parser->SetVariable( "t", t );

    this->dataInterpolation();

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5021 ) << "                                                evaluate(" << 1 << ") : " << M_parser->Evaluate( 1 ) << "\n";
#endif

    return M_parser->Evaluate( 1 );
}

} // Namespace LifeV

#endif /* BCInterface1D_Function_H */
